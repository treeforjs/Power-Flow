Proposed adoption into chemistry v3.2 -- NOT DEPLOYED
===================================================
v3.2-audit-fixes.patch was generated against the exact source files listed in
base-sha256.json. The active v3.2 runtime was not edited. Run verify_base.py on
the intended staging checkout first; it is read-only and requires every base
hash to match. If a hash changed, rebase/review rather than force the patch.

The patch preserves v3.2 reaction channels, Rudd/ion kinematics and model-estimate
metadata. It transfers conservative thermal coverage, conservative actual-parcel
release, GAS-return nondebit and strict default unknown-channel policy. Partial
placeholders below placeholder_below_eV also make the whole region unbounded.
The host ledger also records region_support counts from CarrierOut region
presence, including cold co-moving states with X=B=0; strict refusal must use
support, not a fabricated exposure floor or X>0. Kernel record ABI is unchanged.

Owner integration still required after applying to a new staging version:
  1. Regenerate the v3.2 blob, manifest, checks and strict policy with patched
     builder/policy tools. Old blob metadata is not compatible with the fix.
  2. Adopt pfc-neutral-parcels-v1 checkpoint/digest schema in native; delegate
     release to pfc::release_reservoirs and refuse legacy moments-only stocks.
     Persist/hash/merge region_support in CPU/device totals and checkpoints;
     require reachability_schema=pfc-region-reachable-support-v1.
  3. Copy/adapt audit regression tests from this v3.1.1 artifact and rerun v3.2
     CPU, sanitizer, GPU parity, partition and swarm qualification. The included
     syntax/dry-run checks are NOT that qualification.
  4. Bind new code/blob/policy identities. Production presets must require
     strict-production mode; diagnostic admission must remain explicit.

This patch does not validate v3.2 equal-speed cross-section estimates, sampled
model extrapolations, or its assumed attachment tails. Those remain physics
assumptions needing separate assessment. No incompatible patch hunk was observed
against the captured base, but checkpoint/API callers need the migration above.
