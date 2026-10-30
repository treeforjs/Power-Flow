# Solid thermoelastic work and heat interface

`mhdlab/thermoelastic_energy.py` supplies two distinct accounting models. Both use fixed reference volume and small-strain, energy-conjugate stress/strain coordinates. The six-component convention is `[xx, yy, zz, gamma_xy, gamma_yz, gamma_xz]`, with `gamma_ij=2 epsilon_ij`; stress components are `[sigma_xx,sigma_yy,sigma_zz,tau_xy,tau_yz,tau_xz]`. Thus the ordinary vector dot product is stress work. Never insert another shear factor or current-volume multiplier.

The constitutive starting point is the Helmholtz potential. Stress and entropy follow its strain and negative temperature derivatives; internal energy is free energy plus temperature times entropy. These relations follow the continuum free-energy inequality described in [A. F. Bower, Applied Mechanics of Solids, section3.1](https://solidmechanics.org/Text/Chapter3_1/Chapter3_1.htm). The equations below are derived for this implementation, including temperature-dependent stiffness and integrated eigenstrain.

## Exact algorithmic exchange

At endpoints `i=0,1`, let `e_i = epsilon_i - theta_i`, where `theta_i` is thermal eigenstrain, `sigma_i=C_i e_i`, and `Phi_i=e_i dot C_i e_i/2`. Every `C_i` is symmetric positive definite in the selected strain basis. Define endpoint-average stress `sigma_bar=(sigma_0+sigma_1)/2` and `Wstress=V0 sigma_bar dot Delta epsilon`. Direct expansion gives

```text
Delta Phi = sigma_bar dot Delta e + 0.5 e1 dot (Delta C) e0

Q_alg / V0 = Delta Phi - sigma_bar dot Delta epsilon
           = -sigma_bar dot Delta theta + 0.5 e1 dot (Delta C) e0.
```

`Q_alg` is positive when energy is removed from a numerical thermal reservoir. With a midpoint-compatible momentum update, `Delta K=-Wstress+Wexternal`, so

```text
Q_alg = Delta(K + V0 Phi) - Wexternal.
Delta U_reservoir = Qsupplied - Q_alg.
```

This is an exact endpoint algebra identity, including variable modulus. It is not automatically a physical thermoelastic heat source: `Phi` is recoverable free energy, not generally elastic internal energy. Endpoint averaging is essential to the stated identity. Do not combine this work with a different mechanical work quadrature. Reaction work at moving supports belongs in `Wexternal`; reactions at fixed supports perform zero work. An independently checked global momentum/kinetic ledger is required.

## Conditional thermodynamic model

Use temperature in kelvin, per-reference-volume energies and the potential

```text
psi(epsilon,T) = psi0(T) + Phi(epsilon,T)
Phi = 0.5 e dot C(T) e; e=epsilon-theta(T)
U0(T) = psi0(T) - T psi0'(T); C0(T)=dU0/dT.
```

Here `U0` is the caloric internal energy along the stress-free state and `C0` its capacity. A measured zero-stress heat capacity can support that role under its measurement assumptions. A numerical Cp-as-Cv proxy remains a declared approximation. The absolute reference entropy in `psi0` is immaterial to the energy port and is not supplied.

At fixed total strain, the independently implemented derivatives are

```text
sigma       = C e
sigma_T     = C_T e - C theta_T
Phi_T       = 0.5 e dot C_T e - sigma dot theta_T
Phi_TT      = 0.5 e dot C_TT e - 2 theta_T dot C_T e
              + theta_T dot C theta_T - sigma dot theta_TT
eta_elastic = -Phi_T
U_elastic   = Phi - T Phi_T
U_total     = U0 + U_elastic
C_epsilon   = C0 - T Phi_TT.
```

`U_elastic` is a signed excess energy and may be negative. A nonnegative-free-energy guard must not be reused for it. The implementation checks positive fixed-strain capacity, finite values, positive absolute temperature, positive reference volume and positive-definite stiffness.

The first law `Udot=sigma dot epsilon_dot+Qdot`, with `Qdot` positive heat input, gives

```text
C_epsilon Tdot = T sigma_T dot epsilon_dot + Qdot.
```

Equivalently, `-T sigma_T dot epsilon_dot` belongs on the left of the heat equation. For constant isotropic stiffness and positive expansion, `sigma_T` is negative hydrostatic: expansion cools an adiabatic body. This sign and constant-coefficient limit agree with [Jeremy Bleyer's fully coupled thermoelasticity implementation](https://bleyerj.github.io/comet-fenicsx/tours/multiphysics/thermoelasticity_full/thermoelasticity_full.html), which uses fixed-strain heat capacity and an entropy-based heat equation. [COMSOL's constitutive theory](https://doc.comsol.com/6.3/doc/com.comsol.help.sme/sme_ug_theory.06.026.html) similarly distinguishes elastic entropy and caloric energy. Neither constant-coefficient example supplies our variable-modulus material data.

The conservative discrete base-energy debit is

```text
Q_conditional = V0 Delta U_elastic - Wstress
Delta(V0 U0) = Qsupplied - Q_conditional
R_energy = Delta(V0 U0) + Q_conditional - Qsupplied.
```

Given the compatible momentum identity, this also equals `Delta(K+V0 U_elastic)-Wexternal`. Sum Gauss-point debits using their fixed reference quadrature volumes. Do not use both `Q_alg` and `Q_conditional` in one heat update. The finite-step first-law identity is exact; physical entropy conservation is approached with time refinement, not asserted from energy closure alone. Independent tests compare adiabatic expansion against the constant-entropy endpoint and show second-order refinement.

## Binding the thermal owner

The reviewed thermal port is `/home/treefor/Documents/Codex/2026-10-28/we/outputs/ss304-thermal-properties-20261029/thermal-material-port-v1.json`, SHA256 `1f8d374f4073fdde7e935eedca0bd3d081345475dba08e7995269ee0d8f4e722`. Its header `include/ss304_thermal.hpp` has SHA256 `2fb6d6bd16be0ae0696aac64dbb7e9f61b74eeb7d81b0828059f5d496140fa1c` at review. Its `energy_j_m3_reference` and `capacity_j_m3_reference_k` provide candidate `U0(T)` and `C0(T)` callbacks only under explicit conditional identification with this Helmholtz model. The owner calls the model a fixed-reference-volume Cp-as-Cv approximation and supplies no compressible physical EOS. The coupling cannot upgrade that qualification.

Pass exactly the same accepted reference density/mass/volume as used by mechanics. The thermal owner already multiplies specific energy by reference density; the port multiplies by reference volume once. Do not multiply by density again or replace reference volume with deformed volume. Preserve owner source/scenario hashes, energy reference at300K and explicit scenario selection in the caller's transaction. The thermal owner's range begins at300K, while mechanics stops at773.15K; both gates apply. Cooling below300K must stop rather than clip. Phase change is outside this mechanical port. No latent energy is added here.

The port does not invert `U0` alone after an unconverged debit. `Q_conditional` depends on trial temperature and strain. A coupled corrector must solve the mechanical equations and `R_energy=0` together, or iterate them to a declared residual. Commit only the converged pair. On rejection restore accepted temperature, displacement, velocity, geometry and all energy ledgers together.

## Piecewise material laws and the admissible next step

The frozen expansion record interpolates cumulative relative length and Young's modulus linearly. Within a segment, `C_T` and `theta_T` are constant (for explicit constant nu), while `C_TT=theta_TT=0`. At interior table knots, those first derivatives jump. Although `Phi` stays continuous, `Phi_T`, entropy and `Phi-T Phi_T` generally jump under nonzero elastic strain. Those artificial entropy/internal-energy jumps are interpolation artifacts, not latent heat.

`ss304_thermodynamic_increment` therefore rejects a step crossing or landing on any expansion or modulus knot. It currently accepts one six-vector per Gauss point, with both endpoint temperatures strictly inside the same open segments. The first shared segment supports the planned conditional early-temperature example, with the thermal-owner floor300K and a margin below373.15K. The algorithmic exchange remains valid across knots, but its output must retain its numerical-reservoir label.

A physical-form extension across knots requires an explicitly versioned differentiable constitutive fit, derivative provenance and stability checks. At least continuous first derivatives are necessary to avoid artificial entropy jumps; a smooth second derivative makes the capacity unambiguous. Such a fit must preserve data uncertainty and admissible ranges, be tested for positive fixed-strain capacity over the intended stress domain, and undergo independent temperature/time refinement. Silently smoothing or extrapolating the frozen tables would change the material model.

Nominal `nu=0.29` is separately frozen in `data/materials/ss304-poisson-supplement-v1.json` from the primary [Lamineries MATTHEY 2022 precision-strip sheet](https://www.matthey.ch/fileadmin/user_upload/downloads/fichetechnique/FR/1.4301_v22F.pdf), pages1-2, which explicitly maps1.4301/S30400/304. Its temperature validity is unspecified. Constant nu in these conditional tests is an assumption, not measured nu(T). The source's different melting interval is retained as a discrepancy, without changing either owner's frozen material record.

## Executable interfaces and verification

- `algorithmic_exchange(...)`: endpoint free-energy exchange, exact modulus/eigenstrain decomposition, joule ledger.
- `thermodynamic_state(...)`: explicit generic smooth derivatives; entropy, elastic internal energy and fixed-strain capacity.
- `ss304_thermodynamic_increment(...)`: bounded single-point material adapter with caller-owned caloric callbacks; separate algorithmic and conditional debits and supplied-heat residual.

The module returns `physical_EOS_ready=False` throughout. The tests cover arbitrary symmetric positive-definite variable stiffness, reverse-step work, engineering shear factors, independent finite differences of `Phi_T`, `Phi_TT`, `sigma_T` and internal energy, compatible kinetic balance, signed heat differences, adiabatic cooling and refinement, positive-capacity rejection, domain checks and artificial knot jumps. Run with frozen userspace, CPU7, nice15, one thread and the shared bounded-test lock.
