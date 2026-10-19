# Read-only admission of exported source ledgers

`mhdlab.window_source_ledger_validation.validate_window_source_ledger(document,
expected_source=...)` checks a complete serialized source ledger using JSON
types, exact rational arithmetic and the declared original source association.
The validator uses the standard library. It does not invoke the source law,
canonical-history constructors, moment quadrature, packet sampling, transport
or a host transaction.

The expected source is the original conditioned-source manifest retained by a
preserved study or packet export. Matching it binds the ledger's original
components, geometry, intervals, pools, state identifiers, temperatures,
qualifications and provenance to that declared source. Omitting it allows an
algebra check and reports `source_association_bound: false`. Algebra alone does
not establish the provenance or physical truth of supplied source amounts.

Admission requires full component records. Separate width/current/state/time
marginals cannot reconstruct the joint emission distribution. Exact zero source
skips and declared zero-release geometry/state groups remain explicit. The
validator checks quantity completeness independently for parent number,
translation, internal excitation, removal cost and wall debit. An incomplete
channel retains an unavailable complete total and a separately named qualified
subtotal; it never becomes a structural zero.

Exact group partitions, area denominators, factored width redistribution,
unqualified-catalog references, thermal/wall residuals and SI projections are
checked together. Source availability and scalar representation remain separate.
An available rational total can still have an unsupported floating-point
projection. Signed rounding errors and positive underflow/overflow diagnostics
are retained. JSON booleans cannot substitute for numeric fields, signed zero
is checked, and rational strings must have their canonical form. Explicit work
limits bound document structure and fraction sizes.

The validation summary's `ledger_sha256` hashes the canonical JSON encoding
(`sort_keys=True`, compact separators, finite numbers). It identifies the
validated declaration. Export file hashes identify the original bytes and may
differ because whitespace and key order are preserved in those files.

The reduced-window coupling cache can round slightly above its original
coupling after a log/exp round trip. Admission permits that case only when both
cached scalars reproduce the producer's declared scalar transform exactly; it
does not use a general closeness tolerance or reevaluate the release law.

Admit a previously exported study or packet-source ledger:

```sh
python tools/admit_source_ledger_export.py --export saved_export --expected-report-sha256 APPROVED_REPORT_SHA256 --output-dir new_admission
```

The expected SHA256 is mandatory. Obtain it from the previously reviewed export
record or handoff, then compare the actual report bytes. Calculating a digest
from an unreviewed input and treating that same value as approval supplies no
independent provenance check. A report pin is a data association, rather than
an authenticity signature or proof of physical accuracy.

The command checks the pinned outer report and the exact nested input/code/
output associations. It reads JSON and CSV and hashes source and packet files;
it does not load HDF arrays, decode particle realizations or rerun the numerical
producer. Each ledger is checked against the original conditioned-source
manifest, with exact case/window identity and scalar CSV projections retained.
Producer code lists are inspected as data, rather than executed to repeat the
calculation. Changed contracts require their own qualified admission adapter.

The output contains `admission.json`, `owner_handoff.json` and exact original
export files under `export/`. The handoff identifies each ledger and original
joint-source location with paths, file hashes and JSON pointers. It retains
source labels, geometry, pool/state/temperature/timing associations, quantity
availability, projection diagnostics and the width-only comparison definition.
The original source remains authoritative for downstream launch distributions.

A packet-source export also retains hash-bound references to its original
packet/configuration files, supplied state masses and launch geometry. Those
references preserve the existing realization and its declarations. A source
study that contains no packet, mass or launch-normal context leaves those inputs
pending. The adapter does not invent values or derive a molecular level
distribution from an internal-energy cost.

Both reports keep `host_budget_accepted` and `owner_implementation_verified`
false. Source algebra and association checks do not select a host precision
policy, apply wall cooling or accept an inventory transaction. Actual state/
charge assignments, molecular energy pathways, fields, Hall response, losses,
collection geometry and observation mapping remain owner inputs. Multiple
inventory pools or state labels do not establish fragment yields or resolved
rotational/vibrational populations.

Existing destinations and unsafe or indirect input paths are refused. Staging
and final input/code rehashes assume a coordinated single writer. The portable
subset explicitly labels its empty package initializer; this isolates admission
from eager legacy imports without changing the public numerical files. A
standard-library-only execution can verify that admission does not require the
source producer's numerical runtime.

The report's `runtime.third_party_imports_required: false` describes the
admission algorithm. In the full checkout, the existing package initializer
imports legacy numerical modules and requires their installed dependencies.
The recorded `python -S` qualification applies to the explicitly isolated
portable subset with its empty initializer.

Manufactured examples remain manufactured. Admission does not recover physical
H2O coverage, absolute VUV normalization, missing experimental diagnostics or
spectral discrimination. Quark and separately owned transport/chemistry/
radiative implementations consume this local handoff through their owners.
