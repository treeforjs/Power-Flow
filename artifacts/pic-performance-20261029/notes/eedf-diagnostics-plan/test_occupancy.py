import copy,json,math,pathlib,types
import snapshot_occupancy as diag
c=diag.state_io
HERE=pathlib.Path(__file__).parent

def fixture():
 o={"nx":2,"ny":2,"width_m":.001,"gap_m":.001,"depth_m":.025,"finite_foils":1,"finite_axial_extent":1,"short_z_m":-.0125,"feed_z_m":.0125,"lower_foil.left_m":.0002,"lower_foil.right_m":.0008,"lower_foil.bottom_m":.0001,"lower_foil.top_m":.0002,"upper_foil.left_m":.0002,"upper_foil.right_m":.0008,"upper_foil.bottom_m":.0008,"upper_foil.top_m":.0009}
 particles=[]
 for i,(w,z,u)in enumerate(((2.,0.,(3e7,0.,0.)),(3.,.005,(0.,0.,4e7)),(1.,0.,(0.,0.,0.)))):
  p={"id":i+1,"charge_C":-c.E_CHARGE,"mass_kg":diag.ME,"weight":w}
  for k,v in enumerate((.0003+.0001*i,.0005,z)):p["position_m["+str(k)+"]"]=v
  for k,v in enumerate(u):p["proper_velocity_m_s["+str(k)+"]"]=v
  particles.append(p)
 return types.SimpleNamespace(options=o,particles=particles)

def run_tests():
 cases=[]
 def passed(name):cases.append({"name":name,"passed":True})
 s=fixture();r=diag.analyze(s,grid=1)
 assert r["global"]["full_height"]["weighted_electrons"]==6. and r["global"]["midpoint_slab"]["weighted_electrons"]==3.
 assert math.isclose(r["global"]["full_height"]["nominal_volume_cm3"],.025,rel_tol=1e-14)
 assert math.isclose(r["global"]["full_height"]["vacuum_volume_cm3"],.022,rel_tol=1e-14)
 assert math.isclose(r["global"]["midpoint_slab"]["mean_nominal_density_cm3"],3000.,rel_tol=1e-14)
 passed("weighted count, cm3 conversion, conductor subtraction and strict midpoint slab")
 for sel in ("full_height","midpoint_slab"):
  rows=[x for x in r["energy_rows"]if x["selection"]==sel]
  integral=sum(x["number_pdf_per_eV"]*(x["energy_hi_eV"]-x["energy_lo_eV"])for x in rows)
  assert abs(integral-1.)<1e-14
  assert r["global"][sel]["angular_plus_stationary_weight"]==r["global"][sel]["weighted_electrons"]
 passed("energy PDF and angular plus stationary normalization")
 assert r["global"]["full_height"]["stationary_weight"]==1.
 assert any(row["mu_hi"]==1. and row["weighted_electrons"]==3. for row in r["angle_rows"]if row["selection"]=="full_height")
 passed("zero-speed angular exclusion and positive-z endpoint")
 overflow=diag.analyze(s,grid=1,edges=[0.,.1])
 assert overflow["global"]["full_height"]["overflow"]["weighted_electrons"]==5.
 assert not overflow["global"]["full_height"]["finite_support_complete"]
 assert overflow["global"]["full_height"]["energy_integrated_weight"]==6.
 passed("overflow retained with exact count/energy closure")
 edge=c.kinetic(s.particles[0])/c.E_CHARGE
 bounded=diag.analyze(types.SimpleNamespace(options=s.options,particles=[s.particles[0]]),grid=1,edges=[0.,edge])
 assert bounded["global"]["full_height"]["finite_support_complete"]
 passed("last finite upper energy edge is included")
 empty=diag.analyze(types.SimpleNamespace(options=s.options,particles=[]),grid=1)
 assert empty["global"]["full_height"]["weighted_electrons"]==0 and empty["global"]["full_height"]["mean_energy_eV"]is None
 passed("empty occupancy has defined zero density and undefined probability/mean")
 for name,mutation in (("nonfinite_weight",lambda x:x.particles[0].update(weight=math.nan)),("nonfinite_velocity",lambda x:x.particles[0].update({"proper_velocity_m_s[0]":math.inf})),("nonfinite_position",lambda x:x.particles[0].update({"position_m[0]":math.nan})),("negative_weight",lambda x:x.particles[0].update(weight=-1.)),("outside_z",lambda x:x.particles[0].update({"position_m[2]":1.})),("ion_rejected",lambda x:x.particles[0].update(charge_C=c.E_CHARGE)),("wrong_mass",lambda x:x.particles[0].update(mass_kg=0.))):
  bad=copy.deepcopy(s);mutation(bad)
  try:diag.analyze(bad,grid=1)
  except (ValueError,OverflowError,ZeroDivisionError):passed(name+" rejected")
  else:raise AssertionError(name+" unexpectedly accepted")
 for name,kwargs in (("nan_edge",{"edges":[0.,math.nan]}),("duplicate_edge",{"edges":[0.,0.,1.]}),("invalid_grid",{"grid":0}),("invalid_slab",{"slab_m":0.})):
  try:diag.analyze(s,**kwargs)
  except ValueError:passed(name+" rejected")
  else:raise AssertionError(name+" unexpectedly accepted")
 assert not r["crossing_flux_computed"]and not r["residence_time_weighted"]and not r["chemistry_rates_computed"]
 passed("occupancy-only labels; no inferred crossings, time average or rates")
 report={"status":"passed","cases":len(cases),"results":cases}
 (HERE/"fixture-tests.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
 return report
if __name__=="__main__":
 with c.cpu_lease():print(json.dumps(run_tests()))
