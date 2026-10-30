SS304 3-D foil supports and thermoelastic work — 2026-10-29

Source commit: f6504b619f13ac435b6547d38af3e12f61b55c6a
301 selected tests pass, including72 new3D geometry/mechanics/work/coupling checks.
Frozen CPU7,nice15,one thread,shared testing-core7.lock,noGPU or production action.

New executable components:3-D Q1 transient solid; free/short-clamped/both-end
support scenarios; deformed areas/normals/gap/density; signed internal-energy
exchange and conservative3-D-cellJ to native thermalXY-cellJ projection; local
coupled caloric/mechanical corrector using the separately sealed thermal C ABI.

Matched256x32 source prefix ends8.156064549ns, below368K surface screening cap.
Thermodynamic mode rejects coefficient-knot crossings; mechanics773.15K is a
selected data limit,not melting. Actual supports/structural short remain unknown.

Edge-resolved free result:peak gap narrowing28.005nm; further mesh profile change
1.814nm; fixed-V parallel-plate indicator5.601ppm. No deformed Maxwell solve.
Coarse support variation is much smaller than mesh uncertainty, but the refined
clamped case fails the nominal proof screen at8.156ns:ratio1.00763,VM166.77MPa.
Its geometry/work handoff is rejected. This is not proof of actual foil yielding.
Thin-foil bending/full physical convergence remains unqualified.

Thermodynamic internal energy is Phi−T*Phi_T, distinct from recoverablePhi.
Both debits are saved and never mixed. The coupled-cell test iterates a nominal
310K fixed-geometry heat input to about310.297K with internal+kinetic closure.
This is a conditional Cp-as-Cv caloric model, not a physicalEOS certification.

Native heat integration requires a dedicated signed internal mechanical-exchange
ledger. Do not put it in authored_heat_J, Joule heat or PIC electric work. Parent
must iterate native thermal/mechanical stages and atomically preflight/commit;
conduction, moving native thermal/EM operators and restart are not integrated here.

See docs/solid-foil-supports-and-work.md and docs/interfaces/solid-thermoelastic-port-v2.json.
Source PDFs remain outside Git; hashes/grade/nu/melting discrepancies are frozen.
