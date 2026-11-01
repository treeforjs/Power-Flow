Atomic H/O elastic transport diagnostic, version 1 (2026-11-01)

This isolated component reads the 20 previously recovered, SHA-pinned NIST SRD 64
atomic hydrogen and oxygen differential elastic cross-section exports. It supports
only neutral ground-state H and O, and target-rest electron energies from 1 through
10 keV. This is not a complete collision group, a native chemistry-table addition,
a physical-model qualification, or production permission. All such flags are false.

Primary provenance
  Database: https://srdata.nist.gov/srd64/ ; DOI 10.18434/T4T88K.
  Guide: https://doi.org/10.6028/NIST.NSRDS.64-2023.
  Source-pins.json binds filenames, SHA256, Z, incident energy, query pages, and the
  acquisition receipt. Each file contains 181 angles from 0 to 180 degrees, units
  a0^2/sr, and a separately reported total cross section. The source's stated
  conversion a0^2 = 2.8002852e-21 m^2 is preserved.
  These are calculated relativistic free-atom scattering data. They are neither
  measured fragment spectra nor molecular water/OH/H2 scattering data.
  Raw exports and any arrays derived from them remain outside Git because NIST
  standard-reference-data copyright is distinct from general federal text rights.

Explicit model
  The reported total controls cross_section(), which returns m^2 (not a rate
  coefficient). Totals and DCS operands are linearly interpolated in incident
  energy between actual 1-keV source knots; neither energy tail is extrapolated.
  The 181-point DCS is piecewise linear in mu=cos(theta), with exact full angular
  endpoints. Its integral is used to normalize the angular probability only.
  Source inspection finds that this coarse-grid integral exceeds the independently
  reported total by 0.2072–3.1985 percent. We do not replace the reported rate with
  that coarse integral. The angular shape is therefore an explicitly normalized
  approximation. No angular-convergence or momentum-transfer accuracy is claimed.
  Finer NIST angular sampling is a future comparison, not part of this qualification.

  The sampled polar angle is interpreted in the target-rest laboratory frame,
  with uniform azimuth. A finite-mass relativistic two-body solution is imposed at
  that angle, then both products are boosted to the laboratory using the actual
  target proper velocity. This finite-recoil lift of calculated free-atom angles
  is an explicit modelling assumption. The electron and neutral retain their
  species, charge, nuclei and reference energy; no unspecified internal sink is
  created. The target is returned, not consumed. A later native adoption must use
  authenticated sampled-target momentum/energy feedback, not a second inventory
  debit or a fabricated target packet.

  With electron mass m, target mass M, target-rest incident kinetic energy T,
  p=sqrt(T*(T+2m)), E=m+T, C=M+E, the outgoing electron momentum magnitude is
    p_out = p*((M*E+m*m)*mu + C*sqrt(M*M-m*m*(1-mu*mu))) /
                  (C*C-p*p*mu*mu).
  This is the positive root of energy/momentum conservation at fixed target-rest
  angle. Its radicand avoids subtracting nearly equal rest-energy squares.
  The independent maximum recoil at backscatter is
    2*M*T*(T+2m)/(M*M+m*m+2*M*(m+T)).

  Registry masses are loaded from either explicitly admitted legacy70fcc... or
  CDF8435... blob identities; physically identical species do not alias provenance.
  Arbitrary construction or empty-model evaluation is unavailable. The component
  requires the explicit diagnostic mode, zero target internal energy, supported
  target/energy, unit incident direction, finite proper velocity, and fixture RNG
  purpose. RNG layout v2 words32/33 are disjoint from target rejection and dynamic
  H2 partition. Old fixed-seed trajectory parity is not a claim.

Qualification workflow
  INPUTS.json pins source, independent checks, shared CDF headers, arrival RNG,
  registry bytes, source inspection and frozen numerical environment. A one-core,
  bounded SimQueue job builds and runs normal and ASan/UBSan fixtures, repeats
  deterministic output, and verifies the immutable input set before and after.
  Independent Python/libm integration and 60-digit Decimal relativistic checks do
  not import the sampler's implementation. Endpoint maximum/zero recoil controls,
  target-rest angle recovery under target motion, source total/interpolation,
  normalized inverse CDF, distribution sanity, source tamper/missing/read-limit,
  energy/internal-state/key-domain refusal and cross-blob provenance are checked.
  Statistical controls establish implementation sanity, not physical DCS accuracy.
  Receipts and any preliminary failures are preserved by attempt name.

Pilot scope and adoption
  This component does not turn absent e+O or incomplete e+H groups into complete
  groups. An actual pilot encounter outside these energy/state bounds must still
  refuse; low-energy H/O elastic, excitation, ionization and other missing channel
  requirements are separate. Products become relevant only when their positive
  target census occurs. No old electron-only histogram is a qualified coupled EEDF.
  See artifacts/early-water-pilot-chemistry-map-v1.1-20261101 for the finite encounter
  requirements and separately conditional charged-target reaction continuity gate.
  The existing post-stage2 final-charge guard remains necessary: cross sections and
  global exact charge cannot supply missing local reaction relocation current/work.

  Future native use requires an explicit typed angular descriptor and source/model
  digest, authentic gas target update, conservative returned-target accounting,
  CPU/device qualification and unchanged per-group/process/energy guards. This
  standalone header does not mutate any shared native runtime or policy.
