#!/usr/bin/env python3
"""CPU-only numerical scale/provenance checks; no acceptance threshold fitting."""
import sys,pathlib,math,array,json
O=pathlib.Path(__file__).resolve().parents[1];sys.path.insert(0,str(O/'tools'))
from static_contract import F,sha,record,read_json,write_json,doubles,ints,check_result,verify_pins
from compare_fields import compare

def main():
 R=O/'uniform9-additive-review-v2';old=O.parent/'cathode-emission-audit-20261029/static-domain/9mm';new=O/'parent-static-runs/uniform9';control=O/'uniform-control-sealed'
 checked=check_result(new,True)
 comparison=compare(F/'particles.bin',old/'particle_fields.bin',new/'particle_fields.bin',46726)
 effect=[]
 for width in [13,17,25]:
  domain=O.parent/f'cathode-emission-audit-20261029/static-domain-extended/{width}mm'
  actual=compare(F/'particles.bin',old/'particle_fields.bin',domain/'particle_fields.bin',46726)
  ratios={}
  for name,region in comparison['regions'].items():
   numerator=region['weighted_RMS_delta_E_V_m'];denominator=actual['regions'][name]['weighted_RMS_delta_E_V_m']
   ratios[name]=dict(uniform_old_new_RMS_V_m=numerator,actual_9mm_to_domain_RMS_V_m=denominator,ratio=numerator/denominator if denominator else None)
  effect.append(dict(domain_full_width_mm=width,ratios=ratios,comparison=actual))
 csr={}
 for p in (R/'node-delta-cpu').iterdir():
  if p.name.startswith(('full_operator','constrained_operator')) or p.name in ['charge_rhs','phi_rhs','psi_rhs']:
   csr[p.name]=dict(counterfactual_sha256=sha(p),qualified_sha256=sha(control/p.name));assert csr[p.name]['counterfactual_sha256']==csr[p.name]['qualified_sha256']
 assert len(csr)==9
 I=ints(control/'constrained_operator_I.bin');J=ints(control/'constrained_operator_J.bin');V=doubles(control/'constrained_operator_V.bin');q=doubles(control/'charge_rhs');manifest=read_json(F/'manifest.json');voltage=manifest['potential_upper_boundary_V']
 essential=[all(J[k]==i or V[k]==0 for k in range(I[i],I[i+1])) for i in range(len(q))];assert sum(essential)==3618
 source=[0. if essential[i] else v for i,v in enumerate(q)];source_norm=math.sqrt(math.fsum(v*v for v in source));charge_scale=[]
 for label,folder in [('old_fresh',old),('new_fresh',new)]:
  phi=doubles(folder/'phi_nodes.bin');psi=doubles(folder/'psi_nodes.bin');charge_phi=[a-voltage*b for a,b in zip(phi,psi)]
  residual=math.sqrt(math.fsum((math.fsum(V[k]*charge_phi[J[k]] for k in range(I[i],I[i+1]))-source[i])**2 for i in range(len(q))))
  charge_scale.append(dict(name=label,free_true_charge_rhs_norm_C_per_m=source_norm,subtractive_charge_potential_max_abs_V=max(abs(v) for v in charge_phi),subtractive_charge_solution_residual_C_per_m=residual,ratio_to_free_charge_rhs=residual/source_norm,scope='Diagnostic of rounded phi-Vpsi and independent fsum CSR action; no additional solve and no newly imposed physical guard'))
 result=dict(status='CPU_causal_scale_audit_complete_no_GPU_or_followup',pins=verify_pins(),independent_scalar_guard_checks=checked,exact_reconstructed_CSR_RHS_parity=csr,old_new_uniform_fields=comparison,actual_domain_effect_fractions=effect,true_charge_scale_diagnostic=charge_scale,limits=['An actual GPU CSR/RHS dump was not produced; exact CPU reconstruction plus sealed shared-source/binary/input bindings establish the input path','No single pair identifies the exact internal atomic/reordering schedule','Ratios to actual old static domain effects are descriptive, not fitted acceptance tolerances','Driven-bias residual alone does not prove tiny charge-field accuracy; subtractive charge diagnostic is reported separately','Original exact fresh-output recipe remains failed; no GPU continuation authorized here'])
 write_json(R/'numerical-scale-and-parity.json',result)
 print(json.dumps({'free_true_charge_rhs_norm_C_per_m':source_norm,'true_charge':charge_scale,'domain_effect_fractions':[{k:v for k,v in row.items() if k!='comparison'} for row in effect]},indent=2))

if __name__=='__main__':main()
