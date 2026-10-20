# Manufactured single-patch thermal feedback reference

The five new public files are byte-identical to `e38b819c77ccf9439c82732987e4551c722fbd2d`. The patch starts at
`b27b04fb38238c831e12bf2f3cdfe945e5834d58` and requires that complete checkout. Two wrappers provide a README and
the unchanged baseline requirements. The source subset is not a standalone
package. No initializer replacement or portable execution is supplied.

The reference evolves wall temperature and one finite parent inventory under
zero external power with explicitly manufactured inputs. Supplied areal heat
capacity, fixed kinetic reference coverage, coverage-dependent activation barrier,
constant removal cost, fixed internal energy and expected thermal translation
remain separately identified. Patch area converts areal budgets once.

New targeted tests, saved trajectories, an independent analytic temperature–
inventory/time-map oracle, temporal refinement and per-channel signed residuals
are retained with code, input, command, exit and hash associations. Estimates and
finite changes are not certified error bounds or physical validation. The archive
does not install feedback into any source history, material model or owner code.
No actual-water, foil, phase, power-feed or spectral result is inferred.

The initial 34 tests ran on the preserved initial core, not the corrected bytes.
Five new negative-value guards and one affected empty control checked the final
core after a two-line input-domain correction. The historical test source and
both core versions are retained. The initial nine CLI parser checks and five
final numeric-lexeme checks likewise retain separate source pins. The initial
CLI snapshot is an exact hash-checked reconstruction, not an original pre-edit
capture. The one saved batch and independent oracle bind the corrected final code.

Assembly does not rerun the new reference, tests, oracle, old suites or source
producers. Machine-directory projections, when needed, are reversible and leave
local raw records unchanged. Failed commands and focused corrections are retained.

Next selected task: Implement and independently verify a separate manufactured single-patch prescribed-power reference with an explicit external-input energy ledger, first qualified at zero activation barrier.
