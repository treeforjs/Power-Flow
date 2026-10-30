Heat/source owner v5

See API-v5.txt for the normative contract. v5 is owner v4 (c65ac51) with the A08
gaps G1/G2 closed by refusal: weighted packet values (exchange, ledger, paths,
deposited current/charge arrays, continuity, handoff remainder, neutral leaves)
must be exact binary rescalings of their normal-range evaluation, and the exposure
expected count/diagnostic must lie within T = 2^-40 of the exact finite-log value.
Stock drain and sealed energy-closure failures at the actual weight are typed,
scoped refusals. Everything else is v4 behavior; the v4 focused suite (97 checks)
is ported unchanged in test_source_owner.cpp, followed by the v5 checks.
precision_oracle.py and count_oracle.py are independent test oracles only.
