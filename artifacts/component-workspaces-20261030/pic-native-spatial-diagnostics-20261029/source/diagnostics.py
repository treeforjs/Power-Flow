"""Bounded read-only native spatial/energy diagnostics. No solver or mover."""
from fractions import Fraction
import bisect,hashlib,json,math,struct

QE=1.602176634e-19;C=299792458.;EPS=2.220446049250313e-16
def need(ok,message):
 if not ok:raise ValueError(message)
def finite(x):
 try:return isinstance(x,(int,float)) and not isinstance(x,bool) and math.isfinite(x)
 except OverflowError:return False
def f(x):
 need(finite(x),'nonfinite numeric operand');return Fraction(float(x))
def number(x):
 try:y=float(x)
 except OverflowError:return None
 return y if math.isfinite(y) and (x==0 or y!=0) else None
def count(x):
 value=number(x);return dict(numerator=str(x.numerator),denominator=str(x.denominator),rounded=value,positive_below_double=(x>0 and value is None and x<Fraction(float.fromhex('0x0.0000000000001p-1022'))))
def canonical_hash(value):return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def digest(value):return isinstance(value,str) and len(value)==64 and all(c in '0123456789abcdef' for c in value)
def canonical_count(p):
 hi,lo=p['hi'],p['lo'];need(finite(hi) and finite(lo) and hi>=0 and (hi>0 and abs(lo)<=.5*math.ulp(hi) or hi==0 and lo==0),'noncanonical Count')
 return f(hi)+f(lo)
def kinetic(m,u):
 n=(u[0]*u[0]+u[1]*u[1])+u[2]*u[2];need(math.isfinite(n),'proper velocity norm overflow')
 return m*n/(math.sqrt(1+n/(C*C))+1)
def overlap(a,b):return math.prod(max(0.,min(a[2*d+1],b[2*d+1])-max(a[2*d],b[2*d])) for d in range(3))
class Geometry:
 def __init__(self,j):
  self.axes=[j[k] for k in ['xs_m','ys_m','zs_m']]
  for a in self.axes:need(2<=len(a)<=100000 and all(finite(v) for v in a) and all(x<y for x,y in zip(a,a[1:])),'invalid diagnostic axis')
  self.dims=[len(a)-1 for a in self.axes];need(math.prod(self.dims)<=1010000,'diagnostic cell cap')
  self.bounds=[v for a in self.axes for v in [a[0],a[-1]]];self.solids=j.get('solid_boxes_m',[])
  need(len(self.solids)<=64,'solid box cap')
  for b in self.solids:need(len(b)==6 and all(finite(v) for v in b) and all(b[2*d]<b[2*d+1] for d in range(3)),'invalid solid box')
  for i,a in enumerate(self.solids):
   for b in self.solids[i+1:]:need(overlap(a,b)==0,'overlapping solid volumes require an explicit volume owner')
  self.sha=canonical_hash(j)
 def cell(self,r):
  indices=[]
  for a,x in zip(self.axes,r):
   if x<a[0] or x>a[-1]:return None
   indices.append(min(len(a)-2,bisect.bisect_right(a,x)-1))
  i,j,k=indices;return (k*self.dims[1]+j)*self.dims[0]+i
 def box(self,cell):
  nx,ny,_=self.dims;i=cell%nx;j=(cell//nx)%ny;k=cell//(nx*ny)
  return [self.axes[0][i],self.axes[0][i+1],self.axes[1][j],self.axes[1][j+1],self.axes[2][k],self.axes[2][k+1]]
 def volume(self,box,slab=None):
  b=box.copy()
  if slab:b[4]=max(b[4],slab[0]);b[5]=min(b[5],slab[1])
  if any(b[2*d]>=b[2*d+1] for d in range(3)):return 0.
  v=math.prod(b[2*d+1]-b[2*d] for d in range(3))-math.fsum(overlap(b,s) for s in self.solids)
  need(finite(v) and v>=0,'invalid physical vacuum volume');return v
def validate(data):
 need(data['schema']=='native_spatial_diagnostics_input/v1','input schema')
 need(digest(data['state_prefix']) and finite(data['time_s']) and data['time_s']>=0,'unbound state prefix/clock')
 entries=data['entries'];need(isinstance(entries,list) and len(entries)<=data['maximum_entries']<=2000000,'bounded immutable census required')
 states={s['id']:s for s in data['states']};need(len(states)==len(data['states'])<=128,'state identity cap/duplicates')
 for s in states.values():
  need(isinstance(s['id'],str) and s['id'] and finite(s['mass_kg']) and s['mass_kg']>0 and isinstance(s['charge_e'],int) and not isinstance(s['charge_e'],bool) and abs(s['charge_e'])<=128,'state mass/charge')
  for k in ['internal_J','reference_J']:need(k not in s or finite(s[k]),'nonfinite state energy')
 seen=set();by_packet={};native_ids={};packet_states={};packet_operands={};energy=[]
 for row in entries:
  need(row['state'] in states and isinstance(row['packet_id'],str) and row['packet_id'],'canonical state/packet identity')
  s=states[row['state']];n=row['native'];nid=n['id'];need(isinstance(nid,int) and nid>0 and nid not in seen,'duplicate/invalid native identity');seen.add(nid)
  need(len(n['position_m'])==len(n['proper_velocity_m_s'])==3 and all(finite(v) for v in n['position_m']+n['proper_velocity_m_s']),'invalid native kinematics')
  need(finite(n['weight']) and n['weight']>0 and n['mass_kg']==s['mass_kg'] and n['charge_C']==s['charge_e']*QE,'carrier weight/catalog operands')
  K=kinetic(n['mass_kg'],n['proper_velocity_m_s']);supplied=row.get('kinetic_J_per_particle',K);need(finite(supplied) and supplied>=0,'invalid native energy')
  need(abs(supplied-K)<=32*EPS*max(K,supplied,1e-300),'native energy/kinematics inconsistent');energy.append(supplied)
  packet=row['packet_id'];by_packet[packet]=by_packet.get(packet,Fraction())+f(n['weight']);native_ids.setdefault(packet,[]).append(nid)
  need(packet not in packet_states or packet_states[packet]==row['state'],'mixed states in one packet');packet_states[packet]=row['state']
  operands=(n['position_m'],n['proper_velocity_m_s'],n['mass_kg'],n['charge_C']);need(packet not in packet_operands or packet_operands[packet]==operands,'inconsistent canonical packet kinematics');packet_operands[packet]=operands
 packets=data.get('canonical_packets')
 if packets is not None:
  need(set(packets)==set(by_packet),'canonical packet census mismatch')
  for packet,p in packets.items():
   need(canonical_count(p['count'])==by_packet[packet],'positive carriers do not exactly represent canonical Count')
   need(sorted(p['native_ids'])==sorted(native_ids[packet]),'canonical native ids mismatch')
   need(packet_states[packet]==p['state'],'canonical packet state mismatch')
 return states,energy

def diagnose(data):
 states,energies=validate(data);geo=Geometry(data['geometry']);edges=data['energy_edges_eV'];need(2<=len(edges)<=10001 and all(finite(v) and v>=0 for v in edges) and all(a<b for a,b in zip(edges,edges[1:])),'energy edges')
 slab=data.get('observation_slab_m');need(slab is None or len(slab)==2 and all(finite(v) for v in slab) and slab[0]<slab[1],'observation slab')
 surfaces=data.get('surfaces',[]);need(len(surfaces)<=256 and len({s['id'] for s in surfaces})==len(surfaces),'surface cap/duplicate identity')
 for face in surfaces:
  axis=face['axis'];box=face['bounds_m'];need(isinstance(face['id'],str) and face['id'] and axis in [0,1,2] and face['normal_sign'] in [-1,1] and finite(face['plane_m']) and finite(face['vacuum_band_m']) and face['vacuum_band_m']>0,'surface definition')
  need(len(box)==6 and all(finite(v) for v in box) and box[2*axis]==box[2*axis+1]==face['plane_m'] and all(box[2*d]<box[2*d+1] for d in range(3) if d!=axis),'invalid finite face patch')
 totals={};cells={};obs={};hist={};outside={};unresolved={};total_energy=Fraction();global_charge=Fraction();global_charge_C=Fraction();surface={}
 def add(mapping,key,w):mapping[key]=mapping.get(key,Fraction())+w
 for row,K in zip(data['entries'],energies):
  s=row['state'];n=row['native'];w=f(n['weight']);r=n['position_m'];add(totals,s,w);global_charge+=w*states[s]['charge_e'];global_charge_C+=w*f(n['charge_C']);total_energy+=w*f(K)
  h=hist.setdefault(s,{'counts':[Fraction() for _ in edges[:-1]],'energy_J':[Fraction() for _ in edges[:-1]],'under':Fraction(),'over':Fraction(),'unresolved':Fraction(),'under_energy_J':Fraction(),'over_energy_J':Fraction(),'unresolved_energy_J':Fraction()})
  if K==0 and any(v!=0 for v in n['proper_velocity_m_s']):h['unresolved']+=w;h['unresolved_energy_J']+=w*f(K)
  else:
   E=K/QE
   if E<edges[0]:h['under']+=w;h['under_energy_J']+=w*f(K)
   elif E>edges[-1]:h['over']+=w;h['over_energy_J']+=w*f(K)
   else:
    i=min(len(edges)-2,bisect.bisect_right(edges,E)-1);h['counts'][i]+=w;h['energy_J'][i]+=w*f(K)
  cell=geo.cell(r)
  if cell is None:add(outside,s,w)
  elif geo.volume(geo.box(cell))==0:add(unresolved,s,w)
  else:
   add(cells,(cell,s),w)
   if slab and slab[0]<=r[2]<slab[1]:add(obs,(cell,s),w)
  # Geometric face-band estimators are separate from accepted wall crossings.
  u=n['proper_velocity_m_s'];gamma=math.sqrt(1+math.fsum(v*v for v in u)/(C*C));v=[q/gamma for q in u]
  for face in surfaces:
   axis=face['axis'];normal=face['normal_sign'];lo=face['bounds_m'];band=face['vacuum_band_m'];need(axis in [0,1,2] and normal in [-1,1] and finite(band) and band>0,'surface definition')
   signed=(r[axis]-face['plane_m'])*normal;tangential=[d for d in range(3) if d!=axis]
   distance=math.sqrt(signed*signed+math.fsum(max(lo[2*d]-r[d],0.,r[d]-lo[2*d+1])**2 for d in tangential))
   value=surface.setdefault((face['id'],s),{'weight':Fraction(),'normal_velocity_weighted':Fraction(),'out':Fraction(),'in':Fraction(),'distance_weighted':Fraction(),'min_distance_m':None,'distance_population_weight':Fraction()})
   value['distance_weighted']+=w*f(distance);value['distance_population_weight']+=w;value['min_distance_m']=distance if value['min_distance_m'] is None else min(distance,value['min_distance_m'])
   if signed>=0 and signed<band and all(lo[2*d]<=r[d]<=lo[2*d+1] for d in tangential):
    value['weight']+=w;vn=v[axis]*normal;value['normal_velocity_weighted']+=w*f(vn);value['out']+=w*f(max(vn,0));value['in']+=w*f(max(-vn,0))
 hist_out={}
 for s,h in hist.items():
  need(sum(h['counts'],h['under']+h['over']+h['unresolved'])==totals[s],'histogram exact count inconsistency')
  hs=sum(h['energy_J'],h['under_energy_J']+h['over_energy_J']+h['unresolved_energy_J']);expected=sum((f(r['native']['weight'])*f(K) for r,K in zip(data['entries'],energies) if r['state']==s),Fraction());need(hs==expected,'histogram exact weighted-energy inconsistency')
  hist_out[s]=dict(physical_count=count(totals[s]),bin_counts=[count(v) for v in h['counts']],weighted_energy_J=[number(v) for v in h['energy_J']],dN_dE_per_eV=[number(v/f(edges[i+1]-edges[i])) for i,v in enumerate(h['counts'])],ensemble_EEDF_per_eV=[number(v/totals[s]/f(edges[i+1]-edges[i])) for i,v in enumerate(h['counts'])],underflow_count=count(h['under']),overflow_count=count(h['over']),unresolved_energy_count=count(h['unresolved']),normalization='All physical species weight, including out-of-range/unresolved weight; no in-range renormalization',weighted_kinetic_J=number(expected),exact_weighted_kinetic_receipt=count(expected))
  hist_out[s]['complete_energy_representation']=h['unresolved']==0;hist_out[s]['native_number_weighted_mean_energy_eV']=number(expected/totals[s]/f(QE)) if h['unresolved']==0 else None
  hist_out[s]['population_fractions']=dict(below_declared_bins=number(h['under']/totals[s]),above_declared_bins=number(h['over']/totals[s]),unresolved_native_energy=number(h['unresolved']/totals[s]))
 density=[];net_cells={};net_cells_C={}
 for (cell,s),w in sorted(cells.items()):
  volume=geo.volume(geo.box(cell));entry=dict(cell_id=cell,state=s,physical_count=count(w),physical_reaction_volume_m3=volume,density_cm3=number(w/f(volume)/1000000),charge_density_C_cm3=number(w*f(states[s]['charge_e']*QE)/f(volume)/1000000))
  if slab:
   sv=geo.volume(geo.box(cell),slab);ow=obs.get((cell,s),Fraction());entry['observation']=dict(slab_m=slab,selected_count=count(ow),observation_volume_m3=sv,density_cm3=number(ow/f(sv)/1000000) if sv else None)
  density.append(entry)
  net_cells[cell]=net_cells.get(cell,Fraction())+w*states[s]['charge_e']
  net_cells_C[cell]=net_cells_C.get(cell,Fraction())+w*f(states[s]['charge_e']*QE)
 for s,w in totals.items():need(sum((n for (cell,st),n in cells.items() if st==s),outside.get(s,Fraction())+unresolved.get(s,Fraction()))==w,'spatial exact count inconsistency')
 faces=[]
 for (fid,s),v in sorted(surface.items()):
  face=next(a for a in data['surfaces'] if a['id']==fid);box=face['bounds_m'].copy();axis=face['axis'];box[2*axis]=min(face['plane_m'],face['plane_m']+face['normal_sign']*face['vacuum_band_m']);box[2*axis+1]=max(face['plane_m'],face['plane_m']+face['normal_sign']*face['vacuum_band_m']);box=[max(box[k],geo.bounds[k]) if k%2==0 else min(box[k],geo.bounds[k]) for k in range(6)];vol=geo.volume(box);need(vol>0,'undefined physical face-band volume')
  faces.append(dict(surface=fid,state=s,band_population=count(v['weight']),minimum_geometric_distance_m=v['min_distance_m'],mean_geometric_distance_m=number(v['distance_weighted']/v['distance_population_weight']),band_volume_m3=vol,snapshot_signed_normal_flux_m2_s=number(v['normal_velocity_weighted']/f(vol)),snapshot_outward_flux_m2_s=number(v['out']/f(vol)),snapshot_inward_flux_m2_s=number(v['in']/f(vol)),snapshot_charge_current_density_A_m2=number(v['normal_velocity_weighted']*f(states[s]['charge_e']*QE)/f(vol)),scope='Instantaneous face-band estimator; not actual crossing current, source origin or return classification'))
 sink_report=diagnose_sinks(data,states)
 result=dict(schema='native_spatial_diagnostics/v1',state_prefix=data['state_prefix'],input_sha256=canonical_hash(data),time_s=data['time_s'],arena_v5_source_identity=data.get('arena_v5_source_identity'),species_counts={s:count(w) for s,w in totals.items()},net_charge_equivalent_e=count(global_charge),net_charge_C=number(global_charge_C),exact_native_charge_C=count(global_charge_C),weighted_kinetic_energy_J=number(total_energy),exact_energy_receipt=count(total_energy),energy_edges_eV=edges,energy_distributions=hist_out,species_cell_densities=density,net_charge_cells=[dict(cell_id=cell,equivalent_e=count(w),net_charge_density_C_cm3=number(net_cells_C[cell]/f(geo.volume(geo.box(cell)))/1000000)) for cell,w in sorted(net_cells.items())],outside_diagnostic_grid_counts={s:count(w) for s,w in outside.items()},zero_volume_cell_counts={s:count(w) for s,w in unresolved.items()},surface_geometry_and_flux=faces,sink_diagnostics=sink_report,physical_volume_m3=geo.volume(geo.bounds),observation_slab_m=slab,geometry_sha256=geo.sha,canonical_count_receipts_checked=data.get('canonical_packets') is not None,rate_moments={'status':'separate bound rate_moments.hpp output required; no target or cross-section inferred'},qualification=dict(full_mechanism=False,physical_breakdown=False,evolving_plasma=False,new_particle_evolution=False,GPU=False),limits=['Geometric face proximity never identifies emission origin','Actual returns require explicit verified origin-to-electrode binding','Reaction volumes are physical 3D; optional observation slab only changes observation density','Missing energy/rate support is not clamped or renormalized'])
 for field in ['internal_J','reference_J']:
  unknown=sum((w for s,w in totals.items() if field not in states[s]),Fraction());value=sum((w*f(states[s][field]) for s,w in totals.items() if field in states[s]),Fraction())
  result['weighted_'+field]=number(value) if unknown==0 else None;result['unsupported_'+field+'_positive_count']=count(unknown)
  result['unsupported_'+field+'_positive_fraction']=number(unknown/sum(totals.values(),Fraction())) if totals else None
 result['state_registry_sha256']=canonical_hash(data['states']);result['spatial_unassigned_fractions']={s:dict(outside_grid=number(outside.get(s,Fraction())/w),zero_physical_volume=number(unresolved.get(s,Fraction())/w)) for s,w in totals.items()}
 result['diagnostic_geometry']=data['geometry'];result['spatial_volume_convention']='Declared physical Cartesian3D vacuum cells minus explicit nonoverlapping solid boxes; observation slab has its own clipped volume'
 return result

def diagnose_sinks(data,states):
 sinks=data.get('sinks');window=data.get('sink_time_window_s');origins=data.get('origin_electrode_bindings',{})
 if sinks is None:return {'status':'unavailable; snapshot has no accepted sink history','actual_current_A':None,'return_count':None}
 need(len(sinks)<=100000,'sink cap');need(window is None or len(window)==2 and all(finite(v) for v in window) and window[0]<window[1],'sink time window')
 grouped={};seen=set();excluded=Fraction()
 for row in sinks:
  need(row['receipt_id'] not in seen,'duplicate sink receipt');seen.add(row['receipt_id']);state=row['state'];need(state in states,'sink state registry')
  w=canonical_count(row['count']);need(w>0,'sink count');time=row['time_s'];need(finite(time),'sink time')
  need(finite(row['incident_K_J']) and row['incident_K_J']>=0,'sink kinetic energy')
  carriers=row.get('native_carrier_receipts')
  if carriers is not None:
   need(len(carriers)<=42 and len({a['native_id'] for a in carriers})==len(carriers),'sink carrier identity/cap')
   need(sum((f(a['weight']) for a in carriers),Fraction())==w and all(a['weight']>0 and finite(a['incident_K_J']) and a['incident_K_J']>=0 for a in carriers),'sink exact count/energy operands')
   energy=math.fsum(a['incident_K_J'] for a in carriers);bound=(2*len(carriers)+32)*EPS*max(energy,row['incident_K_J'],1e-300);need(abs(energy-row['incident_K_J'])<=bound,'sink energy receipt inconsistent')
  if window is not None and not(window[0]<=time<window[1]):excluded+=w;continue
  key=(row['sink'],state);g=grouped.setdefault(key,{'w':Fraction(),'q':Fraction(),'K':Fraction(),'returned':Fraction(),'return_unknown':Fraction()});g['w']+=w;g['q']+=w*f(states[state]['charge_e']*QE);g['K']+=f(row['incident_K_J'])
  origin=origins.get(row['origin'])
  if origin is None:g['return_unknown']+=w
  else:
   need(digest(origin['provenance_sha256']),'unbound origin classification')
   if origin['electrode']==row['sink']:g['returned']+=w
 return dict(status='accepted_sink_receipt_diagnostics',time_window_s=window,excluded_time_window_count=count(excluded),actual_current_convention='Signed charge entering the declared sink per elapsed time; distinct from outward-normal snapshot flux',rows=[dict(sink=k[0],state=k[1],physical_count=count(v['w']),charge_C=number(v['q']),incident_kinetic_J=number(v['K']),actual_current_A=number(v['q']/f(window[1]-window[0])) if window else None,returned_to_verified_origin_electrode_count=count(v['returned']),unknown_return_origin_count=count(v['return_unknown']),return_classification_complete=v['return_unknown']==0) for k,v in sorted(grouped.items())],scope='Signed charge/time from accepted sinks only; no heat or inferred source-face ownership')
