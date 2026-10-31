Heat/source owner v6.1 workspace -- 2026-10-31

Contract: API-v6.1.txt (the owner v6 contract of pic-electrode-heat-source-v6-20261030/API-v6.txt with
identical digests and decisions; cached mesh identity, exact matches() by comparison, exact fast emission
spectrum, state digest cache). Sealed inputs (owner v6, owner v5, source kernel v3, heat-v2 headers and
source kernel, arrival v1/v1.1 headers, native candidate003) are read in place and hashed by the scripts.

Build and run (receipted; every compile and run inside arch-base-20261028/run-frozen, core leases):
  core-lease-20261030/any-core.sh --wait /usr/bin/python -I run_source_v6_1.py build LABEL
  core-lease-20261030/any-core.sh --wait /usr/bin/python -I run_source_v6_1.py run LABEL
Outputs: /home/treefor/projects/.workspace/outputs/pic-water-common-commit-v6.1-20261031/source-v6.1/LABEL/

Results of record: s4/run-004620 (build s4), all passed
  test_source_v6_1
    emission    13191 admitted inputs (6 configurations: work function 2.0/3.1/4.3/4.5/6.0 eV, panels 16..128,
                action nodes 16..64; T 0..3000 K; signed F up to 5e9 V/m; 120 random pairs): every Distribution
                field and node bit-identical to the sealed evaluate; 489 refused inputs with identical messages;
                35.8M of 94.8M extracting-field nodes dropped as exact zeros; no configuration fell back
    a02         changed particle weight, substituted spectrum, reused owner address, mutated mesh with old
                metadata, mutation between preflight and publish (source v5 audit port rows, v6 API): v6.1
                verdicts and every recorded observation equal v6. Mutated mesh: matches/preflight/import/handoff/
                prepare/commit refused for dx, zmin, VE[0], edges[0].length (and VN[last], short_node[last]);
                owner unchanged; restored mesh matches and commits
    identity    configuration, state (before and after publication), field content and mesh digests, prepared
                patch exposure key and expected count, birth and continued candidate content, trial content
                (emission and quantized neutral release), release quotes, component checkpoint bytes: equal to v6
    timing      45x4 mesh: Field build 10.3 -> 0.94 ms, matches 9.4 -> 0.19 ms, mesh digest 8.7 -> 3.8 ms,
                prepare_patch (300 K, 2e6 V/m) 2.36 -> 0.58 ms
  test_source_v6_on_v6_1 (the sealed source v6 qualification, unchanged source, compiled against v6.1 through
                compat-v6/): 9/9 (binding, temperature, exactly-once, double-joule, typed-refusals,
                checkpoint, protocol, quantized-release, unbiasedness)
Coordinator use: pic-water-common-commit-v6.1-20261031 (README.txt there).
