// swarm_mc.cpp -- 0D electron swarm Monte Carlo in a uniform field, using the
// chemistry v3 tables and the shared collision kinematics (event_kinematics).
// usage: swarm_mc BLOB OUTPUT.json EN_LIST_TD N_m3 T_gas_K electrons seed [t_relax_ns t_measure_ns]
//
// Method. Electrons in H2O gas (background target; Maxwellian at T_gas) with a
// uniform field E = (E/N) N along -z (force +e E along +z). Null-collision free
// flights with one constant majorant nu_max = N * max_{E' <= Ecap} sigma_tot v
// (Ecap enclosing every energy reached; a run whose electrons exceed Ecap is
// refused). Flights are integrated exactly for the proper momentum
// (dp/dt = eE) and with the nonrelativistic position update (energies here are
// < 1 keV). At a candidate the target velocity is drawn, the target-frame energy
// formed and real/null + channel chosen exactly as in the PIC kernel; real
// events use event_kinematics. Ionization adds the ejected electrons,
// attachment and recombination remove the projectile. Population control:
// uniform keyed resampling to the initial count whenever the population leaves
// [N0/2, 2 N0] (unbiased for the mean drift and energy).
// Observables over the measurement window: flux drift W = <v_z>, bulk drift
// W_b = d<z>/dt (least squares over samples), mean energy, net electron
// production and loss frequencies nu_i, nu_a, alpha/N = nu_i/(N W_b),
// eta/N = nu_a/(N W_b), transverse D_T = d<x^2+y^2>/dt / 4 and D_T/mu.
#include "chem_host.hpp"

#include <fstream>
#include <iostream>
#include <sstream>

using namespace pfc;

namespace {
struct El {
    pfc_u64 id;
    double u[3];   // proper velocity (m/s)
    double x[3];   // position (m)
};
double sq(double a) { return a * a; }
}  // namespace

int main(int argc, char** argv) {
    if (argc < 8) {
        std::cerr << "usage: swarm_mc BLOB OUTPUT.json EN_LIST_TD N_m3 T_gas_K electrons seed [t_relax_ns t_measure_ns]\n";
        return 2;
    }
    try {
        const Blob b = load_blob(argv[1]);
        std::vector<double> ens;
        {
            std::stringstream ss(argv[3]);
            std::string tok;
            while (std::getline(ss, tok, ',')) ens.push_back(std::stod(tok));
        }
        const double N = std::stod(argv[4]), Tg = std::stod(argv[5]);
        const std::size_t N0 = std::stoull(argv[6]);
        const pfc_u64 seed = std::stoull(argv[7]);
        const double t_relax_user = argc > 8 ? std::stod(argv[8]) * 1e-9 : 0.0;
        const double t_meas_user = argc > 9 ? std::stod(argv[9]) * 1e-9 : 0.0;
        const pfc_u32 e_sp = species_index(b, "e"), w_sp = species_index(b, "H2O");
        const pfc_u32 g = group_index(b.view, e_sp, w_sp);
        require(g != 0xFFFFFFFFu, "no e+H2O group");
        const GroupRec& gr = b.view.groups[g];
        const double me = b.view.species[e_sp].mass_eV, mt = b.view.species[w_sp].mass_eV;
        const double c = PFC_C_M_S, qe = 1.602176634e-19, me_kg = 9.1093837139e-31;
        const double kB = 8.617333262e-5;
        const double s_t = std::sqrt(kB * Tg / mt) * c;   // target 1D thermal std (m/s)
        const double trunc = 7.0;
        std::ofstream out(argv[2]);
        out << "{\"blob_sha256\": \"" << b.sha256 << "\", \"N_m3\": " << N << ", \"T_gas_K\": " << Tg << ", \"electrons\": " << N0 << ", \"seed\": " << seed
            << ", \"points\": [";
        bool first_point = true;
        for (double en : ens) {
            const double E = en * 1e-21 * N;   // V/m
            const double acc = qe * E / me_kg;   // m/s^2 along +z
            // majorant cap: energies up to Ecap
            double Ecap = std::max(200.0, 6.0 * en);   // eV
            const double vcap_t = trunc * s_t;
            auto make_majorant = [&](double Ec) {
                // bound over target-frame energies reachable from projectile energies [0, Ec]
                const double wc = vcap_t / c, gc = std::sqrt(1 + wc * wc);
                const double p = std::sqrt(Ec * (Ec + 2 * me));
                const double Ehi = (gc * Ec + gamma_minus_one(wc * wc) * me + wc * p) * (1 + 1e-12);
                return N * majorant(b.view, gr, 0.0, Ehi);
            };
            double numax = make_majorant(Ecap);
            // time scales
            const double t_relax = t_relax_user > 0 ? t_relax_user : 4000.0 / numax;
            const double t_meas = t_meas_user > 0 ? t_meas_user : 16000.0 / numax;
            const int n_samples = 200;
            const double dts = (t_relax + t_meas) / n_samples;
            std::vector<El> el(N0);
            pfc_u64 next_id = 1;
            for (auto& e : el) {
                e.id = next_id++;
                // 1 eV isotropic start
                ChemKey k{seed, 0, e.id, PFC_PURPOSE_FIXTURE_SWARM};
                const V3 d = iso_dir(chem_u01(k, 0, 0), chem_u01(k, 0, 1));
                const V3 u = u_from_p(me, scale(d, std::sqrt(1.0 * (1.0 + 2 * me))));
                e.u[0] = u.x;
                e.u[1] = u.y;
                e.u[2] = u.z;
                e.x[0] = e.x[1] = e.x[2] = 0.0;
            }
            std::uint64_t candidates = 0, reals = 0, above_cap = 0, below_grid_time = 0;
            std::vector<std::uint64_t> per_proc(b.n_proc, 0);
            double sum_vz_t = 0, sum_e_t = 0, sum_n_t = 0;
            std::vector<double> ts, zc, r2c;
            double z_offset = 0;   // population control keeps positions; no offset needed
            (void)z_offset;
            std::uint64_t meas_created = 0, meas_removed = 0;
            double meas_electron_time = 0;
            for (int smp = 0; smp < n_samples; ++smp) {
                const bool measuring = (smp + 1) * dts > t_relax;
                std::vector<El> next;
                next.reserve(el.size() + 16);
                double int_vz = 0, int_e = 0, int_n = 0;   // time integrals over this interval
                for (auto& e0 : el) {
                    El e = e0;
                    ChemKey key{seed, pfc_u64(smp + 1), e.id, PFC_PURPOSE_FIXTURE_SWARM};
                    double t = 0;
                    pfc_u32 cand = 0;
                    bool alive = true;
                    while (true) {
                        const double tau = -pfc::log(chem_u01(key, cand, PFC_W_TIME)) / numax;
                        const double rem = dts - t;
                        const bool hit = tau < rem;
                        const double fl = hit ? tau : rem;
                        // flight integrals (nonrelativistic): vz(t) = vz0 + a t
                        const double vz0 = e.u[2];
                        int_vz += vz0 * fl + 0.5 * acc * fl * fl;
                        // energy integral: 0.5 m |v|^2 integrated (eV s)
                        const double v2_0 = sq(e.u[0]) + sq(e.u[1]) + sq(e.u[2]);
                        const double ke_int = 0.5 * me_kg / qe * (v2_0 * fl + vz0 * acc * fl * fl + acc * acc * fl * fl * fl / 3.0);
                        int_e += ke_int;
                        int_n += fl;
                        e.x[0] += e.u[0] * fl;
                        e.x[1] += e.u[1] * fl;
                        e.x[2] += vz0 * fl + 0.5 * acc * fl * fl;
                        e.u[2] = vz0 + acc * fl;
                        t += fl;
                        if (!hit) break;   // no candidate within the interval
                        ++candidates;
                        // candidate
                        const V3 u = v3(e.u[0], e.u[1], e.u[2]);
                        const V3 p = p_from_u(me, u);
                        const double K = ke_from_u(me, u);
                        if (K > Ecap) ++above_cap;
                        V3 ut = v3(0, 0, 0);
                        for (pfc_u32 a = 0; a < 4u; ++a) {
                            const pfc_u32 w0 = PFC_W_TVEL + 4u * a;
                            const double r1 = std::sqrt(-2.0 * pfc::log(chem_u01(key, cand, w0)));
                            const double r2 = std::sqrt(-2.0 * pfc::log(chem_u01(key, cand, w0 + 2u)));
                            double c1, s1, c2, s2;
                            sincos_2pi(chem_u01(key, cand, w0 + 1u), &c1, &s1);
                            sincos_2pi(chem_u01(key, cand, w0 + 3u), &c2, &s2);
                            const double z0 = r1 * c1, z1 = r1 * s1, z2 = r2 * c2;
                            if (z0 * z0 + z1 * z1 + z2 * z2 <= trunc * trunc) {
                                ut = scale(v3(z0, z1, z2), s_t);
                                break;
                            }
                        }
                        const V3 wt = w_from_u(ut);
                        const double gt = std::sqrt(1.0 + norm2(wt));
                        const V3 pt = boost_in(wt, gt, p, me + K);
                        const double T = kinetic_from_p2(me, norm2(pt));
                        const double vt = std::sqrt(norm2(pt)) / (me + T) * c;
                        const int gi = grid_locate(b.view.grid + gr.grid_off, gr.grid_n, T);
                        if (gi < 0) ++below_grid_time;
                        const double y = chem_u01(key, cand, PFC_W_CHANNEL) * (numax / N);
                        pfc_u32 sel = 0xFFFFFFFFu;
                        double F = 0;
                        for (pfc_u32 j = 0; j < gr.proc_n; ++j) {
                            F += sigma_eval(b.view, gr, j, gi, T) * vt;
                            if (sel == 0xFFFFFFFFu && y < F) sel = j;
                        }
                        require(F <= numax / N, "swarm majorant violated");
                        if (sel != 0xFFFFFFFFu) {
                            const pfc_u32 pidx = gr.proc_off + sel;
                            const ProcessRec& pr = b.view.procs[pidx];
                            V3 plab[PFC_MAX_OUT];
                            double mk[PFC_MAX_OUT], eint;
                            if (event_kinematics(b.view, pr, gr, gi, key, cand, me, mt, pt, T, wt, gt, plab, mk, &eint)) {
                                ++reals;
                                per_proc[pidx]++;
                                int survivor = -1;
                                for (pfc_u32 i = 0; i < pr.n_out; ++i) {
                                    if (pr.out_species[i] != e_sp) continue;
                                    const V3 uo = u_from_p(mk[i], plab[i]);
                                    if (pr.out_dest[i] == PFC_DEST_SURVIVOR) {
                                        survivor = int(i);
                                        e.u[0] = uo.x;
                                        e.u[1] = uo.y;
                                        e.u[2] = uo.z;
                                    } else {
                                        El ne;
                                        ne.id = 0;   // assigned below (deterministic order)
                                        ne.u[0] = uo.x;
                                        ne.u[1] = uo.y;
                                        ne.u[2] = uo.z;
                                        ne.x[0] = e.x[0];
                                        ne.x[1] = e.x[1];
                                        ne.x[2] = e.x[2];
                                        next.push_back(ne);
                                        if (measuring) ++meas_created;
                                    }
                                }
                                if (survivor < 0) {
                                    alive = false;
                                    if (measuring) ++meas_removed;
                                    break;
                                }
                            }
                        }
                        ++cand;
                    }
                    if (alive) next.push_back(e);
                }
                if (measuring) {
                    sum_vz_t += int_vz;
                    sum_e_t += int_e;
                    sum_n_t += int_n;
                    meas_electron_time += int_n;
                }
                // ids for newborns (deterministic order of creation)
                for (auto& ne : next)
                    if (ne.id == 0) ne.id = next_id++;
                // population control
                if (next.size() > 2 * N0 || next.size() < N0 / 2) {
                    std::vector<El> res;
                    res.reserve(N0);
                    ChemKey pk{seed, pfc_u64(smp + 1), 0xFFFFFFFFFFFFull, PFC_PURPOSE_FIXTURE_SWARM};
                    for (std::size_t k = 0; k < N0; ++k) {
                        const std::size_t j = std::size_t(chem_u01(pk, pfc_u32(k >> 8), pfc_u32(k & 255) * 2) * double(next.size()));
                        El x = next[std::min(j, next.size() - 1)];
                        x.id = next_id++;
                        res.push_back(x);
                    }
                    next.swap(res);
                }
                require(!next.empty(), "swarm extinct");
                el.swap(next);
                // samples for bulk drift and diffusion
                double zm = 0, r2 = 0;
                for (auto& e : el) zm += e.x[2];
                zm /= double(el.size());
                for (auto& e : el) r2 += sq(e.x[0]) + sq(e.x[1]);
                r2 /= double(el.size());
                if (measuring) {
                    ts.push_back((smp + 1) * dts);
                    zc.push_back(zm);
                    r2c.push_back(r2);
                }
                // raise the cap if needed (restart not required: the bound only grows)
                (void)Ecap;
            }
            auto slope = [&](const std::vector<double>& y) {
                const std::size_t n = ts.size();
                double st = 0, sy = 0, stt = 0, sty = 0;
                for (std::size_t i = 0; i < n; ++i) {
                    st += ts[i];
                    sy += y[i];
                    stt += ts[i] * ts[i];
                    sty += ts[i] * y[i];
                }
                return (n * sty - st * sy) / (n * stt - st * st);
            };
            const double W_flux = sum_vz_t / sum_n_t;
            const double W_bulk = slope(zc);
            const double mean_e = sum_e_t / sum_n_t;
            const double nu_i = double(meas_created) / meas_electron_time, nu_a = double(meas_removed) / meas_electron_time;
            const double DT = slope(r2c) / 4.0;
            const double mu = W_bulk / E;
            out << (first_point ? "" : ",") << "\n  {\"E_over_N_Td\": " << en << ", \"W_flux_m_s\": " << W_flux << ", \"W_bulk_m_s\": " << W_bulk
                << ", \"mean_energy_eV\": " << mean_e << ", \"nu_ion_s\": " << nu_i << ", \"nu_att_s\": " << nu_a << ", \"alpha_over_N_m2\": " << nu_i / (N * W_bulk)
                << ", \"eta_over_N_m2\": " << nu_a / (N * W_bulk) << ", \"alpha_minus_eta_over_N_m2\": " << (nu_i - nu_a) / (N * W_bulk)
                << ", \"k_ion_m3_s\": " << nu_i / N << ", \"k_att_m3_s\": " << nu_a / N << ", \"mobility_N_1_m_V_s\": " << mu * N
                << ", \"DT_over_mu_eV\": " << DT / mu << ", \"candidates\": " << candidates << ", \"real_events\": " << reals
                << ", \"candidates_above_Ecap\": " << above_cap << ", \"candidates_below_grid\": " << below_grid_time << ", \"Ecap_eV\": " << Ecap
                << ", \"nu_max_s\": " << numax << ", \"t_relax_s\": " << t_relax << ", \"t_measure_s\": " << t_meas << ", \"process_events\": {";
            for (pfc_u32 p = gr.proc_off; p < gr.proc_off + gr.proc_n; ++p)
                out << (p == gr.proc_off ? "" : ", ") << "\"" << b.process_names[p] << "\": " << per_proc[p];
            out << "}}";
            out.flush();
            first_point = false;
            std::cout << "E/N " << en << " Td: W " << W_bulk << " m/s, <e> " << mean_e << " eV, alpha/N " << nu_i / (N * W_bulk) << ", eta/N "
                      << nu_a / (N * W_bulk) << " m2, above cap " << above_cap << std::endl;
        }
        out << "\n]}\n";
    } catch (const std::exception& e) {
        std::cerr << "ERROR " << e.what() << "\n";
        return 1;
    }
    return 0;
}
