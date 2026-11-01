Primary numerical OH supplement, isolated inspection/reader v1

The final published 2019 Chakrabarti, Laporta and Tennyson paper provides a
numerical workbook, unlike the earlier arXiv conclusion. The workbook was
retrieved by following the actual XLSX link on the publisher's supplementary
data page. Its filename contains ab264c; the article DOI is ab364c. The exact
published bytes are pinned, and no filename correction or data substitution
is performed. DOI: 10.1088/1361-6595/ab364c, PSST 28 (2019) 085013.
https://iopscience.iop.org/article/10.1088/1361-6595/ab364c/data
https://discovery.ucl.ac.uk/id/eprint/10081520/

The data page explicitly applies CC BY 3.0 to original content, with author,
title, journal and DOI attribution. Raw workbook, publisher HTML, and PDF
remain outside git. Repository artifacts contain only reader, tests, source
hashes and derived metadata/compatibility reports, never complete arrays.

This reader uses bounded ZIP/XML reads; it does not execute workbook formulas,
macros, or links. Exact source identity, sheet names, column state labels,
literal units, numeric rows and monotonicity are checked. Report output lists
actual sheet ranges and counts. Bare zero tabulated values are reported as
source values, never promoted to a proof that a missing process is absent.

The elastic sheet is an integral fixed-geometry calculation with Born
correction, not a qualified angular sampler. The excitation sheet distinguishes
three final electronic states. The single product-state-unlabelled dissociation column cannot be split
into O(3P)+H versus O(1D)+H branches: the paper's Figure 5 discusses both and
the current registry lacks a resolved excited O(1D) state. Section 4.2.1 calls
the dissociation cross section an estimate and assumes excitation into
repulsive states dissociates above D0=3.88 eV. This is not adopted as the
registry ground-product reference cost. The actual registry cost is computed
from its pinned binary64 reference energies in the queued report.

Section 4.2.2 describes DEA as an approximate resonance model, with assumed
only O-minus production and chosen Morse/exponential potentials. The report
counts source nodes below the current O-minus + H reference cost; no shift,
clamp or threshold tolerance is introduced. Its resonance columns must not
be summed again with the Total column. Excitation and dissociation overlap is
also unqualified: adding both as independent events risks double counting.

The thermal Rate sheet literally labels a rate coefficient bohr^2. This is a
dimensional conflict with paper Figure 4 (cm^3/s); it is inventoried but not
consumed as a rate pending source-column/normalization verification. No OH
ionization data appear in this workbook. Below/above-grid tails, finite-target
energy conventions, state mapping, scattering/fragment energy laws and
channel completeness remain unresolved. Existing native tables, production
scope gates, unknown channels, and energy/thermal support guards are unchanged.

Short controls and repeatable metadata extraction use frozen Python and one
SimQueue core. No simulation is launched and no native chemistry is adopted.
production_qualified=false
physical_model_validated=false
network_complete=false
native_adopted=false

Revision v1.1: the actual workbook has missing cell B546 in the A-doublet-Sigma-plus excitation column, at abscissa A546=3.8298999999999999 eV. The original strict reader refused this in job 441. inspect() now independently reports complete columns and a RefusedChannel for this incomplete column. load_channel() explicitly refuses it; no value is filled, row skipped, or missing interval interpolated. Blank formatting-only rows at ends of Rate/DEA are not data. This does not qualify the complete columns as collision event laws.

Revision v1.2: numeric_column_complete names structural completeness only. No claim is made that the unlabelled dissociation column is a sum; product assignment is unresolved. Orphan ordinate controls were added; actual source bytes and all refusal policies are unchanged.
