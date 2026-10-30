# SS304 melt-ablation/layer proposal v1

The stateless component is implemented and verified for manufactured proposals.
Melting alone gives no removal allowance. Every requested weight must fit the
present local thermal budget including its actual product reference, internal
and kinetic energies. A refusal leaves the common-owner snapshot and the request
unchanged. Water release remains independent below melt, retained oxide survives
water depletion, and explicit oxide depletion changes structural exposure.

The component has no accepted state, native birth or publication API. It is a
proposal/evaluator for the one common owner. It does not choose foil layer
thicknesses, inventories, triggers, removal energies or release rates. All such
foil inputs remain unresolved in publication v2. The numeric test entities and
oxide properties are manufactured conservation fixtures only.

Final verification: 74 checks per normal/UBSan run, 3 compile-time contract
checks, 699 independent Fraction/Decimal checks across matching receipts. No
coupled test, physics validation, rate calibration or full kinetics claim.
See `COMPONENT-CONTRACT.md`, `REVIEW.md`, `independent-audit.json`,
`normal-002/receipt.json` and `ubsan-001/receipt.json` for scope and evidence.

The original readiness audit and appendix remain sealed. The successor
publication corrects the requested wording while preserving the distinctions
between melt eligibility, finite removal energy, below-melt desorption and
retained oxide exposure. Their binding verification is recorded separately.

Production remains stopped and frozen. Existing source/material/mechanics guards
remain untouched. Full integration stays held until the actual common owner,
native products/source law, field/geometry/thermal transactions and restart are
bound and validated under their existing constraints.
