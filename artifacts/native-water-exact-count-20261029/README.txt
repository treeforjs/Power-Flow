Actual finite-rate parent integration with exact native gas inventory

All 442 existing parent-subflow controls pass against immutable native count
source b8a9a27e6a0fcd8c54e49dd341503fed165da460180b608aa0ae024b0436fb0e.
Only three Count2 member-access sites were adapted. Rates, collision grid,
product kernel, scheduler, seeds and numerical inputs are unchanged. Both saved
seven-event fixtures reproduce the earlier event IDs, volumes, weights, rates,
energies and local waiting times exactly; the first-event sample remains 73/128.

The exact count backend retains very small represented debits, exposes positive
native carriers and refuses lossy two-limb export. Its checkpoint schema is new;
legacy state is not relabeled. Native independent review/sanitizer qualification
is tracked separately from this integration regression. All logs and exact input
archive are preserved, including compiler source identity and binary receipt.

These are frozen-position operator fixtures. No physical trajectory time was
advanced, no surface source or field/heating coupling was executed, and no full
plasma breakdown or 60 ns result is claimed. The original run stays paused.
Thermal and strong-field emission remain required in the integrated model.
