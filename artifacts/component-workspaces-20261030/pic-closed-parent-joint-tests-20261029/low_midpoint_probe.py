# Primitive binary64 arithmetic only; no native, Hall, field solve or coupled step.
import json
a=1.0
lo=2.0**-54
increment=2.0**-52
# Exact values for Fields::add_field with these operands.
new=a+increment
new_lo=lo
joint_mid=.5*(a+new)
receipt_mid=a+(lo+.5*increment)
assert joint_mid!=receipt_mid
print(json.dumps({"schema":"compensated-midpoint-arithmetic-v1",
"old_high_hex":a.hex(),"old_low_hex":lo.hex(),"increment_hex":increment.hex(),
"new_high_hex":new.hex(),"new_low_hex":new_lo.hex(),
"joint_high_only_midpoint_hex":joint_mid.hex(),"field_receipt_midpoint_hex":receipt_mid.hex(),
"difference":receipt_mid-joint_mid,"status":"arithmetic_counterexample",
"actual_field_or_particle_run":False},indent=2))
