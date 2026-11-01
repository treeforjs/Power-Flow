Source v7.6: authentic prepared exposure introspection

Qualified job 680: 152 checks in normal and 152 in UBSan builds. Legacy raw
emission/packet/layer traces are byte-identical to the original baseline.

PreparedPatch::exposure() returns the immutable actual exposure interval;
patch_index() returns its actual owner-local patch index, distinct from caller
field patch IDs. Both require a live, unmoved token; borrowed exposure views
expire on token move/destruction. Owner::preflight(prepared, field, actual) checks
exact owner lifetime/generation, bound thermal pointer/time/revision/config and
exact immutable field/view. A thermal-only publication invalidates old exposure.

No model, precision, guard, count law, cache or existing state encoding changed.
The new source-attempt/layer owner is separate; this component alone does not
qualify retained-layer native handoff, Hall/source charge/current/energy coupling,
magnetic-barrier applicability, a populated plasma pilot, or 60 ns evolution.

source/README.txt and earlier API documents are preserved historical lineage.
The mass-aware cache is inherited from v7.4 and explicit old/CDF blob admission
from v7.5; this final suite uses the old-blob actual mass fixture. Artifact evidence
retains exact source, final/prior receipts, logs and default parity byte streams;
build binaries and large generated component checkpoints are not duplicated.
