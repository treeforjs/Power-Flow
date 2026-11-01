# Primary fragment table candidate with explicit missing-process refusals

This isolated table successor preserves the sealed v3.1.1 species, process prefix and PFCHEMV3 version-1 ABI, and adds one diagnostic OH dissociation process. It does not replace any live table or native binding.

The added e+OH channel uses 389 actual primary calculated knots from Chakrabarti, Laporta and Tennyson, Plasma Sources Science and Technology 28 085013 (2019), DOI [10.1088/1361-6595/ab364c](https://doi.org/10.1088/1361-6595/ab364c). The explicit diagnostic model uses a 7.36 eV vertical loss, ground O+H products, full excess energy in fragment recoil, and isotropic electron/fragment angles. Supported incident energies are 7.3605–10 eV; no source knot is shifted. The underlying workbook is CC BY 3.0; no raw workbook or numeric array is redistributed here. The isolated event model was qualified separately in oh-a4sigma-diagnostic-v1.1-20261101.

Six missing OH process categories are encoded as unknown placeholders, not zero cross sections: elastic/momentum transfer; rotational/vibrational/other electronic excitation; other neutral dissociation; ionization; dissociative attachment; radiative attachment. All 16 new group regions remain unbounded under the strict policy; no unbounded region is admitted. A supported individual process does not qualify a populated encounter group.

The H2 MCCC integral has 53 retained primary knots from 4.5–150 eV in a separate local rate bank. It has NO native sampling row or encounter group in this candidate. Thus neither a direct evaluator nor a permissive table traversal can mistake an unavailable H2 differential closure for a known zero. The authentic minimum-water scope guard must still refuse positive e+H2 encounters. Subsequent isolated H2 differential work, including the source-CDF diagnostic, requires a new candidate and identity; this table is unchanged.

SimQueue job 497 used one allocated core, frozen userspace and 2 GiB. It passed 27 scope/rate controls, 256 actual OH blob events with cold and moving targets, deterministic table regeneration and event replay, and a wrong-blob-SHA refusal. Maximum relative chemical-energy defect was 1.93e-16. Native transaction integration was not tested here. The final blob has 21 groups, 87 processes and 151 coverage regions; 48 regions are bounded and 103 refused as unbounded, with zero unbounded admissions.

Exact pins are in evidence/job-497/contract.json and receipt.json. Candidate blob SHA is 50a57ead9e71fd30d60924a492fee140197e1f0d016af40afe77b777ec4143f4; policy SHA is dacdc3a7005ee47e32da768e1bfb8f053a9cc14d875b784ca5d97d13c9f860db. This candidate uses the historical event kernel in its component qualification; the later disjoint RNG correction must be adopted under a new native kernel identity before new coupling. No old fixed-seed parity is implied.

Reproduce through SimQueue with the exact INPUTS.json dependencies:

```sh
simq submit -p powerflow -c 1 --mem 2G -t 00:05:00 -n fragment-table-candidate-v1.2 -e CUDA_VISIBLE_DEVICES= -- /home/treefor/projects/runtimes/arch-base-20261028/run-frozen /usr/bin/python3 -I /home/treefor/projects/runtimes/water-fragment-table-candidate-v1.2-20261101/tools/check.py
```

Generated blob/specification, primary rate bank, arrays and binaries remain outside git. MCCC source redistribution restrictions remain in force; this artifact includes code, metadata/hashes and aggregate evidence only. Original v1 (H2 unknown zero rows) and v1.1 (compile harness failure) are preserved on disk and are not promoted.

`production_qualified=false`, `physical_model_validated=false`, `network_complete=false`, `native_adopted=false`. Full H2O fragmentation and the broader supported-energy network remain incomplete.
