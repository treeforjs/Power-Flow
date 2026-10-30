"""Validate phase accounting without running a replay or changing its outputs."""
import json
import math
import pathlib
import sys

path = pathlib.Path(sys.argv[1])
record = json.loads(path.read_text())
assert record["schema"] == "pic_additive_step_phases_v1"
assert record["step_count"] == len(record["steps"])
names = set(record["definitions"])
assert set(record["totals"]) == names
summed = {name: {"calls": 0, "seconds_inclusive": 0.0, "seconds_exclusive": 0.0}
          for name in names}
wall = 0.0
previous = None
for step in record["steps"]:
    assert set(step["phases"]) == names
    assert math.isfinite(step["wall_seconds"]) and step["wall_seconds"] >= 0
    assert math.isfinite(step["dt_s"]) and step["dt_s"] > 0
    if previous is not None:
        assert step["step"] == previous["step"] + 1
        assert step["physical_start_s"] >= previous["physical_start_s"]
    previous = step
    wall += step["wall_seconds"]
    exclusive = 0.0
    for name, phase in step["phases"].items():
        assert isinstance(phase["calls"], int) and phase["calls"] >= 0
        inc, exc = phase["seconds_inclusive"], phase["seconds_exclusive"]
        assert math.isfinite(inc) and math.isfinite(exc)
        assert inc >= -1e-12 and exc >= -1e-12 and inc + 1e-12 >= exc
        exclusive += exc
        for key in summed[name]:
            summed[name][key] += phase[key]
    assert math.isclose(exclusive, step["wall_seconds"], rel_tol=1e-12, abs_tol=1e-12)
    assert step["phases"]["driver_step"]["calls"] == 1
    assert math.isclose(step["additive_error_seconds"],
                        exclusive - step["wall_seconds"], abs_tol=1e-12)
assert math.isclose(wall, record["step_wall_seconds"], rel_tol=1e-12, abs_tol=1e-12)
for name in names:
    assert summed[name]["calls"] == record["totals"][name]["calls"]
    for key in ["seconds_inclusive", "seconds_exclusive"]:
        assert math.isclose(summed[name][key], record["totals"][name][key],
                            rel_tol=1e-12, abs_tol=1e-12)
exclusive = sum(phase["seconds_exclusive"] for phase in record["totals"].values())
assert math.isclose(exclusive, wall, rel_tol=1e-12, abs_tol=1e-12)
terminal = record["terminal_loop_condition"]
if terminal is not None:
    assert terminal["phases"]["driver_step"]["calls"] == 1
    assert math.isclose(sum(p["seconds_exclusive"] for p in terminal["phases"].values()),
                        terminal["wall_seconds"], rel_tol=1e-12, abs_tol=1e-12)
print(json.dumps({"valid": True, "steps": record["step_count"],
                  "step_wall_seconds": wall, "additive_error_seconds": exclusive - wall,
                  "terminal_condition_separate": terminal is not None}, indent=2))
