"""Validate/hash an observational capture; never run PIC, a solver, or a GPU API."""
import argparse
import array
import hashlib
import json
import math
import pathlib

parser = argparse.ArgumentParser()
parser.add_argument("directory", type=pathlib.Path)
parser.add_argument("--allow-contract-fixture", action="store_true")
parser.add_argument("--audit-prefix", type=pathlib.Path)
parser.add_argument("--write-receipt", action="store_true")
args = parser.parse_args()
root = args.directory.resolve()
meta = json.loads((root / "capture.json").read_text())
assert meta["schema"] == "pic_rhs_observation_v1"
assert meta["status"] == "complete", meta["status"]
assert not meta["synthetic_contract_fixture"] or args.allow_contract_fixture
assert meta["requested_populated_steps"] == meta["captured_populated_steps"] == 8
assert meta["captured_solves"] == len(meta["calls"]) == 16
assert meta["startup_and_unit_bias_excluded"]
assert meta["solver_contract"]["matrix_and_vector_precision"] == "FP64"
assert meta["solver_contract"]["ir_n_steps"] == 2
assert meta["no_extra_solver_or_residual_evaluation"]
assert meta["no_physics_or_checkpoint_data_changed"]
assert 0 < meta["output_limit_bytes"] <= 128 * 1024 * 1024
matrix = meta["matrix"]
n, nnz = matrix["rows"], matrix["stored_entries_including_zeros"]
assert matrix["present"] and n > 0
if not meta["synthetic_contract_fixture"]:
    assert n == 130321, "Unexpected matrix dimension for the frozen 360x360 diagnostic"
expected_files = {"capture.json"}
payload_bytes = 0


def read(kind, filename, count):
    global payload_bytes
    assert pathlib.Path(filename).name == filename
    data = (root / filename).read_bytes()
    result = array.array(kind)
    assert len(data) == count * result.itemsize, filename
    result.frombytes(data)
    payload_bytes += len(data)
    expected_files.add(filename)
    return result


rows = read("i", matrix["row_file"], n + 1)
cols = read("i", matrix["column_file"], nnz)
values = read("d", matrix["value_file"], nnz)
essential = read("i", matrix["essential_file"], matrix["essential_list_entries"])
mask = read("B", matrix["essential_mask_file"], n)
assert rows[0] == 0 and rows[-1] == nnz
assert all(a <= b for a, b in zip(rows, rows[1:]))
assert all(0 <= col < n for col in cols)
assert all(math.isfinite(value) for value in values)
essential_set = set(essential)
assert all(0 <= dof < n for dof in essential_set)
assert all(value == (index in essential_set) for index, value in enumerate(mask))
names = {"projected_previous_phi", "full_charge_rhs", "constrained_rhs",
         "incoming_zeroed_X", "accepted_solution", "guard_residual"}
last_solution = None
first_system = meta["calls"][0]["linear_system_id"]
first_field = meta["calls"][0]["field_state_address"]
previous_step = None
projected_nonzero_interiors = []
projected_vs_previous_solution_differences = []
for index, call in enumerate(meta["calls"]):
    assert call["sequence"] == index + 1
    assert call["captured_step_ordinal"] == index // 2 + 1
    assert call["half"] == index % 2 + 1
    assert call["phase"] == ("after_first_half_drift" if index % 2 == 0 else "after_second_half_drift")
    assert call["system_role"] == "main_system_within_populated_advance"
    assert call["linear_system_id"] == first_system
    assert call["field_state_address"] == first_field
    assert call["matrix_address"] == matrix["matrix_address"]
    assert call["particles_at_advance_entry"] > 0
    assert math.isfinite(call["dt_s"]) and call["dt_s"] > 0
    for field in ["driver_start_time_s", "ledger_time_s_at_advance_entry",
                  "lower_potential_V", "upper_potential_V", "residual_norm", "rhs_norm"]:
        assert math.isfinite(call[field]), field
    assert call["unchanged_cpu_guard_passed"]
    assert call["residual_norm"] <= max(1e-30, 2e-13 * call["rhs_norm"])
    if index % 2 == 0:
        if previous_step is not None:
            assert call["driver_step_index"] > previous_step
        previous_step = call["driver_step_index"]
    else:
        previous = meta["calls"][index - 1]
        for field in ["driver_step_index", "driver_start_time_s", "dt_s",
                      "ledger_time_s_at_advance_entry", "lower_potential_V",
                      "upper_potential_V", "particles_at_advance_entry"]:
            assert call[field] == previous[field], field
    assert set(call["vectors"]) == names
    vectors = {}
    for name, info in call["vectors"].items():
        assert info["entries"] == n and info["bytes"] == n * 8
        vectors[name] = read("d", info["file"], n)
        assert all(math.isfinite(value) for value in vectors[name]), name
    projected, incoming = vectors["projected_previous_phi"], vectors["incoming_zeroed_X"]
    assert all(incoming[i] == 0 for i in range(n) if not mask[i])
    assert all(projected[i] == incoming[i] for i in essential_set)
    projected_nonzero_interiors.append(sum(projected[i] != 0 for i in range(n) if not mask[i]))
    if index % 2:
        projected_vs_previous_solution_differences.append(
            sum(projected[i] != last_solution[i] for i in range(n) if not mask[i]))
    last_solution = vectors["accepted_solution"]
assert payload_bytes == meta["payload_bytes"]
entries = list(root.iterdir())
assert all(path.is_file() and not path.is_symlink() for path in entries), "Capture directory must contain only its flat regular files"
actual_files = {path.name for path in entries}
assert actual_files <= expected_files | {"validation.json"}, actual_files - expected_files
audit_parity = None
if args.audit_prefix:
    upper_rows, upper_cols, upper_values = array.array("i", [0]), array.array("i"), array.array("d")
    for row in range(n):
        entries = [(cols[k], values[k]) for k in range(rows[row], rows[row + 1])
                   if cols[k] >= row and values[k] != 0]
        entries.sort()
        upper_cols.extend(col for col, value in entries)
        upper_values.extend(value for col, value in entries)
        upper_rows.append(len(upper_cols))
    for suffix, data in [(".rows.i32", upper_rows), (".cols.i32", upper_cols),
                         (".values.f64", upper_values)]:
        assert data.tobytes() == pathlib.Path(str(args.audit_prefix) + suffix).read_bytes(), suffix
    audit_parity = True
hashes = {}
for name in sorted(expected_files):
    data = (root / name).read_bytes()
    hashes[name] = {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}
total = sum(path.stat().st_size for path in root.iterdir() if path.is_file())
assert total <= meta["output_limit_bytes"]
report = {
    "status": "validated", "synthetic_contract_fixture": meta["synthetic_contract_fixture"],
    "physical_steps": 8, "main_system_solves": 16, "all_six_vectors_finite": True,
    "incoming_X_zero_interiors_and_essential_values_checked": True,
    "second_half_projected_vs_first_half_accepted_nonessential_difference_counts":
        projected_vs_previous_solution_differences,
    "projection_note": "Original boundary projection can change DOFs outside the final essential list; those differences are observed, not treated as an invalid or reconstructed prior.",
    "source_pre_form_capture_provenance": "See source-removal proof and CPU deep-copy contracts; the validator does not synthesize a prior.",
    "projected_nonzero_interiors_per_call": projected_nonzero_interiors,
    "original_runtime_guard_metadata_checked_without_extra_MatVec": True,
    "raw_native_CSR_including_zeros": True, "audit_upper_CSR_bitwise_parity": audit_parity,
    "bytes_before_receipt": total, "limit_bytes": meta["output_limit_bytes"],
    "files": hashes, "validator_executed_GPU_or_solver": False}
text = json.dumps(report, indent=2) + "\n"
if args.write_receipt:
    assert not (root / "validation.json").exists()
    assert total + len(text.encode()) <= meta["output_limit_bytes"]
    (root / "validation.json").write_text(text)
print(json.dumps({key: value for key, value in report.items() if key != "files"}, indent=2))
