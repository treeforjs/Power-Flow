Read-only native spatial and energy diagnostics

Scope
This isolated source consumes an immutable, synchronized accepted particle census.
It contains no mover, event selection, mutation, Poisson construction, GPU call,
source injection or heat feedback. It does not establish a breakdown mechanism,
an evolving plasma result, wall heat, or physical/domain convergence.

Native port
source/native_diagnostic_view.hpp exports PicArena::rate_population(maximum_entries)
from the immutable v5 public API. It uses actual positive native carriers rather
than the rounded Arena::view() packet weights. It exports packet Count hi/lo,
native IDs, source identity, prefix/chain/revision and accepted sink receipts.
The caller must hold/freeze the owner against concurrent changes for the read;
the before/after pointer and prefix checks are additional checks, not a mutex.
The caller must supply the EXACT runtime state registry and physical geometry
used by that arena. A source hash alone cannot bind runtime registry operands.
The output input hash, registry hash and geometry hash bind those supplied values;
parent integration must compare them to the owner's independent expected values.
This port was compiled against v5, but never called on a live native arena here.
Later v6.1 reportedly retains rate_population; it has not been compiled/qualified
by this handoff. No native arena sources have been edited.

Spatial and energy accounting
source/diagnostics.py validates unique native IDs, state mass/charge, positive
weights and finite positions/proper velocities. Canonical packets additionally
require identical carrier kinematics and an exact Fraction sum matching hi+lo.
Fraction sums retain every binary-double weight and weighted native kinetic
energy product, including products below double range. Unrepresentable rounded
values are null with exact receipts, never false zero. This is exact accounting
of the supplied native floating operands, not arbitrary precision dynamics.
Native kinetic energy is m*u^2/(sqrt(1+u^2/c^2)+1); supplied native energy is
checked against those operands with a 32-epsilon arithmetic check. Histogram
count and weighted-energy partitions are checked exactly. The EEDF is the
number-weighted kinetic-energy PDF per eV, not a Maxwellian fit or velocity PDF.
Under/overflow and unresolved native-energy weights/fractions remain explicit;
all species weight stays in the denominator. No in-range renormalization.
The final edge is included in the last bin. Zero-energy particles remain in
the first bin. A nonzero velocity with rounded zero native energy is unresolved.
Internal/reference energy are separate catalog quantities; missing values yield
null plus unsupported positive counts/fractions. Field potential energy and a
complete electromagnetic energy balance are not inferred.

Declared Cartesian 3D reaction cells subtract explicit nonoverlapping solid
boxes. Species number density is count/physical-vacuum-volume/1e6 in cm^-3.
Signed net charge is separate from all positive species densities. The optional
z observation slab selects actual carriers and uses its own clipped vacuum
volume; it NEVER changes the physical reaction volume or rate denominator.
Outside-grid and zero-volume positive weights/fractions are reported, not moved
to a boundary cell. Only the declared upper physical grid face is included in
the last cell; other cell ownership uses right-sided axis lookup.

Surface distances refer to explicitly supplied finite rectangular face patches.
Normals point in the caller-declared direction; the fixture normals point from
solid into vacuum. Face-band sum(W*v_normal)/volume is an instantaneous flux
estimator with units m^-2 s^-1, and q times it is A/m^2. It is not an observed
crossing current or an emission-origin/return classifier. Distances use the
whole supplied population and never identify source origin. Actual current is
signed charge entering accepted sinks divided by a supplied [start,end) time
window. No window means no actual current. Returns require an explicit verified
origin-to-electrode binding and provenance hash; unknown origins remain unknown.
Native incident sink energy is reported as kinetic energy, not wall heat.

H2O rate moment port
source/rate_moments.hpp is a pure bounded actual-carrier quadrature. The caller
must assign EVERY actual census entry to its independently bound physical cell
and provide the matching physical 3D vacuum volume and geometry hash. This port
does not infer a cell, target velocity, neutral density or cross section. The
same accepted census/prefix, exact registry, coordinate map and volume must be
used by spatial and rate ports. Unassigned entries must fail upstream; dropping
them before this function is not permitted by this contract.

For each explicit electron/actual-target pair, the pinned existing collision
relative() computes target-rest-frame energy and Moeller relative flux. Per cell:
  raw = sum_e,sum_n We * Wn * sigma(E_target) * relative_flux
  physical number rate [s^-1] = raw / V_physical_3D
  sum_e We * <sigma*flux>_actual_target = raw / sum_n Wn
The target average divides by ALL target weight, never only supported targets.
The pair cap is checked before invoking any query; no partial/truncated report.
Long-double Neumaier sums are compensated quadrature, not exact cross-section
arithmetic. Decimal long-double values plus double/null/underflow status remain.
Unsupported pair weight/fraction, affected projectile/target weight/fractions
and reasons are explicit. Affected fractions mean any unsupported pair for that
carrier, not the same thing as pair-weight fractions. No target gives an explicit
status and null target average; no target velocity is invented. Identical-species
queries are rejected because their self-pair convention is outside this port.

RateQuery requires dataset and kernel hashes plus model/channel/observable.
Its read-only evaluate(E) callback MUST bind the physical data owner’s exact
support/threshold semantics and return a reason for every unsupported interval.
It must distinguish inclusive/yield cross sections from supported exclusive
event kernels. Supported sigma with missing exclusive kernel has its own
positive pair fraction. A supported zero cross section is only valid where
the data owner explicitly supports that zero; gaps cannot be converted to zero.
This handoff tests analytic queries only. It does not evaluate a physical H2O
table, select reactions, generate products or certify any full mechanism.

Bounds and evidence
CLI input/output cap: 128 MiB; <=2M native entries; <=1.01M diagnostic cells;
<=128 states; <=64 disjoint solid boxes; <=256 face patches; <=100k sinks;
<=10k energy bins. Rates: <=64 queries and <=1M actual pairs. Use an outer
memory/time cap for arbitrary inputs; the CPU contract harness uses 4 GiB and
90 s per test. All work here used CPU 7, nice 15, one thread and the shared
testing-core7.lock. Frozen dependencies only. The compile-qualified port and
analytic tests are separate from any owner/live-arena integration qualification.

The preserved fixture is the original 26.478443519271536 ns post-half-drift,
post-first-solve/pre-kick electron-only midpoint, not a restart checkpoint or
common-arena census. Its 46,726 native carriers have actual fractional weights.
Its source hashes are in the output. Synthetic packet labels only identify old
native IDs; no arena lineage, gas/ion population, sink history, emission origin
or returns are invented. Eight in-plane foil faces are sampled; this is not an
axial-face census. Its physical z extent is 25 mm and the observation slab is
1 mm. No new trajectory or field solve occurs while parsing the fixture.

Build and reproduction
See build-and-check.sh for the frozen single-job command. Source archive and
dependency/input hashes are recorded in the output bindings.json and handoff.
Original exploratory outputs remain in place; handoff names final receipts.
