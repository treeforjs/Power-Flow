"""Independent test oracle only; never used as production count arithmetic."""
from fractions import Fraction
from pathlib import Path
import json, sys
rows = {}
for line in (Path(sys.argv[1]) / "count-oracle.csv").read_text().splitlines():
    fields = line.split(",")
    rows[fields[0]] = sum((Fraction(float.fromhex(x)) for x in fields[1:]), Fraction())
initial = Fraction(1e15) + Fraction(-.01)
if rows["initial"] != initial:
    raise RuntimeError("initial exact Count2 import differs")
checks = 1
previous_source, previous_gas = initial, Fraction()
for step in range(1, 13):
    source, gas = rows[f"source-{step}"], rows[f"gas-{step}"]
    if not (source + gas == initial and source < previous_source and gas > previous_gas):
        raise RuntimeError(f"source debit/gas credit oracle failed at {step}")
    if previous_source - source != gas - previous_gas:
        raise RuntimeError("source and emitted leaf deltas differ exactly")
    previous_source, previous_gas = source, gas
    checks += 4
result = {"status": "passed", "checks": checks, "source_plus_gas_equals_initial_exactly": True, "oracle": "Python fractions.Fraction; production backend not imported", "final_gas_exact": str(previous_gas)}
print(json.dumps(result, indent=2))
