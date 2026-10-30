// gpu_parity.cpp -- runs the GPU push engine on one exported fixture and writes the
// same canonical encoding as cpu_reference, plus a JSON receipt (device, NVRTC
// record, timings per phase, memory, teardown).
// Usage: gpu_parity <fixture dir> <output stem> <nvrtc library> <device uuid> [repeats]
// Writes <stem>.bin and <stem>.json. Exit status: 0 success, 2 in-process repeat
// mismatch, 3 teardown quarantine, 1 other failure.
#include "push/hall_gpu_push.hpp"
#include "push/push_abi.hpp"
#include "tests/encoding.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>

namespace hp = hall_gpu::push;
using Clock = std::chrono::steady_clock;

namespace {
double ms_since(Clock::time_point t) { return std::chrono::duration<double, std::milli>(Clock::now() - t).count(); }

std::string encode(const hp::MotionView& v, const std::vector<std::uint64_t>& id, const std::vector<double>& q, const std::vector<double>& m,
                   const std::vector<double>& w, double dt, bool diagnostics) {
    hgp_enc::RecordArrays r;
    r.refusal = v.refusal.refused ? v.refusal.message : "";
    r.refusal_rank = std::int32_t(v.refusal.rank);
    r.charge_refusal = v.charge_refusal.refused ? v.charge_refusal.message : "";
    r.charge_refusal_rank = std::int32_t(v.charge_refusal.rank);
    const std::size_t n = v.n;
    r.id = id; r.charge = q; r.mass = m; r.weight = w;
    r.flight.assign(n, dt);
    if (!v.refusal.refused) {
        r.end_pos.reserve(3 * n); r.end_u.reserve(3 * n); r.start.reserve(3 * n); r.disp.reserve(3 * n);
        for (std::size_t i = 0; i < n; ++i) {
            r.end_pos.insert(r.end_pos.end(), {v.end_x[i], v.end_y[i], v.end_z[i]});
            r.end_u.insert(r.end_u.end(), {v.end_ux[i], v.end_uy[i], v.end_uz[i]});
            r.start.insert(r.start.end(), {v.start_x[i], v.start_y[i], v.start_z[i]});
            r.disp.insert(r.disp.end(), {v.disp_x[i], v.disp_y[i], v.disp_z[i]});
        }
        r.current.assign(v.current_A_m.begin(), v.current_A_m.end());
        r.initial_charge.assign(v.initial_charge_C.begin(), v.initial_charge_C.end());
        r.endpoint_charge.assign(v.endpoint_charge_C.begin(), v.endpoint_charge_C.end());
        const auto& b = v.budget;
        r.budget = {b.kinetic_change_J, b.field_work_J, b.work_defect_J, b.continuity_C, b.plasma_dt, b.maximum_gyro, b.maximum_displacement_cells};
        r.net_hi.assign(v.free_net_hi.begin(), v.free_net_hi.end());
        r.net_lo.assign(v.free_net_lo.begin(), v.free_net_lo.end());
        r.gross.assign(v.free_gross_C.begin(), v.free_gross_C.end());
    }
    hgp_enc::Writer wr;
    r.write(wr, true);
    if (diagnostics) {
        wr.i32("diag.code1", std::vector<std::int32_t>(v.code1.begin(), v.code1.end()));
        wr.i32("diag.code2", std::vector<std::int32_t>(v.code2.begin(), v.code2.end()));
        wr.i32("diag.code3", std::vector<std::int32_t>(v.code3.begin(), v.code3.end()));
        for (int k = 0; k < hgp::diag_real_fields; ++k)
            wr.f64(std::string("diag.real.") + hgp_enc::diag_real_names[k], v.diag_real.subspan(std::size_t(k) * n, n));
        for (int k = 0; k < hgp::diag_int_fields; ++k) {
            const auto s = v.diag_int.subspan(std::size_t(k) * n, n);
            wr.i32(std::string("diag.int.") + hgp_enc::diag_int_names[k], std::vector<std::int32_t>(s.begin(), s.end()));
        }
    }
    return wr.bytes();
}
}  // namespace

int main(int argc, char** argv) {
    if (argc < 5) {
        std::cerr << "usage: gpu_parity <fixture dir> <output stem> <nvrtc library> <device uuid> [repeats]\n";
        return 1;
    }
    const std::string fixture = argv[1], stem = argv[2];
    const int repeats = argc > 5 ? std::max(1, std::atoi(argv[5])) : 5;
    std::ostringstream js;
    int rc = 0;
    hall_gpu::TeardownReport teardown;
    std::string error;
    try {
        const auto t_read = Clock::now();
        hgp_enc::Reader in(fixture + "/input.bin");
        hp::MeshTables mt;
        const auto mi = in.i32("mesh.int");
        const auto mr = in.f64("mesh.real");
        mt.mesh = {mi[0], mi[1], mi[2], mi[3], mi[4], mi[5], mi[6], mi[7], mi[8], mi[9], mr[0], mr[1], mr[2], mr[3], mr[4], mr[5]};
        mt.VN = in.f64("mesh.VN");
        const auto ea = in.i32("mesh.edge_a"), eb = in.i32("mesh.edge_b");
        mt.edge_a.assign(ea.begin(), ea.end());
        mt.edge_b.assign(eb.begin(), eb.end());
        mt.edge_length = in.f64("mesh.edge_length");
        hp::Geometry geo;
        const auto p = in.f64("policy");
        geo.policy = {p[0], p[1], p[2], p[3], {p[4], p[5], p[6], p[7]}, {p[8], p[9], p[10], p[11]}};
        const auto fr = in.f64("frame");
        for (int d = 0; d < 3; ++d) geo.native_to_field[d] = fr[d];
        const double dt = in.f64("dt")[0];
        const auto E = in.f64("field.E"), B = in.f64("field.B");
        const auto id = in.u64("carrier.id");
        const auto x = in.f64("carrier.x"), y = in.f64("carrier.y"), z = in.f64("carrier.z");
        const auto ux = in.f64("carrier.ux"), uy = in.f64("carrier.uy"), uz = in.f64("carrier.uz");
        const auto q = in.f64("carrier.charge"), m = in.f64("carrier.mass"), w = in.f64("carrier.weight");
        const bool diagnostics = in.text("mode") == "diagnostics";
        const double read_ms = ms_since(t_read);
        const std::size_t n = id.size();

        const auto t_open = Clock::now();
        hall_gpu::SessionOptions so;
        so.require_uuid = argv[4];
        so.nvrtc_library = argv[3];
        hall_gpu::Session session(so);
        const double open_ms = ms_since(t_open);
        const auto t_engine = Clock::now();
        hp::EngineOptions eo;
        eo.capacity = std::max<std::size_t>(n, 1);
        eo.diagnostics_capable = diagnostics;
        std::string first, info_runs;
        bool repeat_identical = true;
        {
            hp::MotionEngine engine(session, mt, geo, eo);
            const double engine_ms = ms_since(t_engine);
            const hp::CarrierView view{id, x, y, z, ux, uy, uz, q, m, w};
            std::ostringstream runs;
            for (int r = 0; r < repeats; ++r) {
                const auto t0 = Clock::now();
                engine.set_carriers(view);
                const auto t1 = Clock::now();
                engine.set_field(E, B);
                auto v = engine.run(dt, diagnostics);
                const double call_ms = ms_since(t1), total_ms = ms_since(t0);
                const auto t_enc = Clock::now();
                const auto bytes = encode(v, id, q, m, w, dt, diagnostics);
                const double encode_ms = ms_since(t_enc);
                if (r == 0) first = bytes;
                else if (bytes != first) repeat_identical = false;
                runs << (r ? "," : "") << "{\"repeat\":" << r << ",\"set_carriers_host_ms\":" << hall_gpu::json_number(std::chrono::duration<double, std::milli>(t1 - t0).count())
                     << ",\"field_and_run_wall_ms\":" << hall_gpu::json_number(call_ms) << ",\"total_wall_ms\":" << hall_gpu::json_number(total_ms)
                     << ",\"run_wall_ms\":" << hall_gpu::json_number(v.run_wall_ms) << ",\"host_fold_ms\":" << hall_gpu::json_number(v.host_fold_ms)
                     << ",\"encode_ms\":" << hall_gpu::json_number(encode_ms) << ",\"edge_records\":" << v.edge_records << ",\"gpu_phases_ms\":{";
                for (std::size_t k = 0; k < v.gpu_phases.size(); ++k)
                    runs << (k ? "," : "") << hall_gpu::json_string(v.gpu_phases[k].name) << ":" << hall_gpu::json_number(v.gpu_phases[k].ms);
                runs << "},\"refusal\":" << hall_gpu::json_string(v.refusal.message) << ",\"refusal_stage\":" << hall_gpu::json_string(v.refusal.stage)
                     << ",\"refusal_rank\":" << v.refusal.rank << ",\"refusal_id\":" << v.refusal.id << ",\"refusal_code\":" << v.refusal.code
                     << ",\"charge_refusal\":" << hall_gpu::json_string(v.charge_refusal.message) << "}";
            }
            info_runs = runs.str();
            js << "{\"fixture\":" << hall_gpu::json_string(fixture) << ",\"carriers\":" << n << ",\"mode\":" << hall_gpu::json_string(diagnostics ? "diagnostics" : "record")
               << ",\"repeats_in_process\":" << repeats << ",\"repeat_outputs_identical\":" << (repeat_identical ? "true" : "false")
               << ",\"output_sha256\":" << hall_gpu::json_string(hall_gpu::sha256_hex(first.data(), first.size()))
               << ",\"device\":" << hall_gpu::device_json(session.device()) << ",\"compile\":" << hall_gpu::compile_json(engine.compile_record())
               << ",\"engine_device_bytes\":" << engine.device_bytes() << ",\"memory\":" << hall_gpu::memory_json(session.memory())
               << ",\"setup_ms\":{\"read_fixture\":" << hall_gpu::json_number(read_ms) << ",\"session_open\":" << hall_gpu::json_number(open_ms)
               << ",\"engine_compile_and_allocate\":" << hall_gpu::json_number(engine_ms) << "},\"runs\":[" << info_runs << "]";
        }
        {
            std::ofstream out(stem + ".bin", std::ios::binary | std::ios::trunc);
            out.write(first.data(), std::streamsize(first.size()));
            if (!out) throw std::runtime_error("cannot write " + stem + ".bin");
        }
        teardown = session.close();
        js << ",\"teardown\":" << hall_gpu::teardown_json(teardown) << "}";
        if (!repeat_identical) rc = 2;
        if (teardown.quarantined || !teardown.released) rc = 3;
    } catch (const std::exception& e) {
        error = e.what();
        std::cerr << "gpu_parity failed: " << error << '\n';
        js.str("");
        js << "{\"fixture\":" << hall_gpu::json_string(fixture) << ",\"error\":" << hall_gpu::json_string(error) << "}";
        rc = 1;
    }
    std::ofstream(stem + ".json") << js.str() << '\n';
    return rc;
}
