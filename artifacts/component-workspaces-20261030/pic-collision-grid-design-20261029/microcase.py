#!/usr/bin/env python3
"""Manufactured geometric binning prototype only. Integer-nm test geometry;
never quantizes a native particle, constructs PIC or evaluates a reaction."""
from fractions import Fraction
import argparse,bisect,hashlib,json,math,os,pathlib,random,resource,signal,time
BOX=(0,9000000,0,9000000);Z=(0,25000000)
SOLIDS=[(2000000,7000000,1800000,2000000),(2000000,7000000,7000000,7200000)]
ROOT=pathlib.Path('/home/treefor/projects/runtimes/pic-water-native-subflow-20261029')
def need(ok,text):
 if not ok:raise ValueError(text)
def inside(box,p):return box[0]<=p[0]<=box[1] and box[2]<=p[1]<=box[3]
def material(p):
 for b in SOLIDS:
  if inside(b,p):return 'material_boundary' if p[0] in b[:2] or p[1] in b[2:] else 'solid'
 return 'vacuum'
def area(b):return (b[1]-b[0])*(b[3]-b[2])
def canonical(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
class Grid:
 def __init__(self,normal_nm=1000,normal_growth=.3,tangent_nm=125000,side_tangent_nm=25000,bulk_nm=500000,nz=50,cap=1000000,near_nz=None,near_band_nm=None):
  need(normal_nm>0 and normal_growth>0 and tangent_nm>0 and side_tangent_nm>0 and bulk_nm>0 and nz>0,'positive grid parameters')
  need((near_nz is None and near_band_nm is None) or (isinstance(near_nz,int) and near_nz>0 and near_band_nm>0),'complete local-z family required')
  self.options=dict(normal_nm=normal_nm,normal_growth=normal_growth,tangent_nm=tangent_nm,side_tangent_nm=side_tangent_nm,bulk_nm=bulk_nm,nz=nz,cap=cap,near_nz=near_nz,near_band_nm=near_band_nm);self.leaves=[];self.trees=[];self.nodes=0
  def target(b):
   tx=ty=bulk_nm
   for s in SOLIDS:
    xover=min(b[1],s[1])-max(b[0],s[0]);yover=min(b[3],s[3])-max(b[2],s[2]);dx=max(s[0]-b[1],b[0]-s[1],0);dy=max(s[2]-b[3],b[2]-s[3],0)
    if xover>0 and (b[3]<=s[2] or b[2]>=s[3]):ty=min(ty,normal_nm+normal_growth*dy);tx=min(tx,tangent_nm)
    elif yover>0 and (b[1]<=s[0] or b[0]>=s[1]):tx=min(tx,normal_nm+normal_growth*dx);ty=min(ty,side_tangent_nm)
    else:
     corner=normal_nm+normal_growth*math.hypot(dx,dy);tx=min(tx,corner);ty=min(ty,corner)
   return tx,ty
  def split(b):
   self.nodes+=1;need(self.nodes<=200000,'2D geometry node cap before further splitting');tx,ty=target(b);ratios=((b[1]-b[0])/tx,(b[3]-b[2])/ty)
   if max(ratios)<=1:
    need(len(self.leaves)<100000,'2D leaf cap');idx=len(self.leaves);self.leaves.append(b);return idx
   axis=0 if ratios[0]>=ratios[1] else 1;lo,hi=b[2*axis:2*axis+2];mid=(lo+hi)//2;need(lo<mid<hi,'unresolved split');left=list(b);right=list(b);left[2*axis+1]=mid;right[2*axis]=mid;return (axis,mid,split(tuple(left)),split(tuple(right)))
  xs=[0,2000000,7000000,9000000];ys=[0,1800000,2000000,7000000,7200000,9000000]
  for j in range(len(ys)-1):
   for i in range(len(xs)-1):
    b=(xs[i],xs[i+1],ys[j],ys[j+1]);cx=(b[0]+b[1])//2;cy=(b[2]+b[3])//2
    if material((cx,cy))=='solid':continue
    self.trees.append((b,split(b)))
  self.xy_area=sum(map(area,self.leaves));need(self.xy_area==79000000000000,'exact vacuum partition area')
  def distance(b):return min(math.hypot(max(s[0]-b[1],b[0]-s[1],0),max(s[2]-b[3],b[2]-s[3],0)) for s in SOLIDS)
  self.z_counts=[near_nz if near_nz is not None and distance(b)<=near_band_nm else nz for b in self.leaves];self.offsets=[0]
  for n in self.z_counts:self.offsets.append(self.offsets[-1]+n)
  self.total=self.offsets[-1];need(self.total<=cap,'3D bin cap; no partial grid or automatic coarsening')
  self.z_families={n:[Z[0]+(Z[1]-Z[0])*k//n for k in range(n+1)] for n in set(self.z_counts)};need(all(all(a<b for a,b in zip(z,z[1:])) for z in self.z_families.values()),'z resolution');self.identity=canonical(dict(schema='manufactured_aligned_extruded_bins/v2',coordinate_unit='integer_nm_manufactured_only',box=BOX,z_families=self.z_families,z_counts_by_xy_leaf=self.z_counts,solids=SOLIDS,leaves=self.leaves,options=self.options,ownership='geometry only; Active microscopic layer collision eligibility is common-stage-owned'))
 def locate(self,p):
  need(len(p)==3 and all(isinstance(x,(int,Fraction)) for x in p),'manufactured exact coordinates required; no rounding of native positions')
  need(BOX[0]<=p[0]<=BOX[1] and BOX[2]<=p[1]<=BOX[3] and Z[0]<=p[2]<=Z[1],'outside physical3D collision box');need(material(p)=='vacuum','solid/material surface requires explicit side/stage')
  tree=None
  for box,node in self.trees:
   if box[0]<=p[0] and (p[0]<box[1] or box[1]==BOX[1]) and box[2]<=p[1] and (p[1]<box[3] or box[3]==BOX[3]):tree=node;break
  need(tree is not None,'vacuum partition hole')
  while isinstance(tree,tuple):axis,mid,left,right=tree;tree=left if p[axis]<mid else right
  z=self.z_families[self.z_counts[tree]];k=min(len(z)-2,bisect.bisect_right(z,p[2])-1);return self.offsets[tree]+k
 def volume(self,cell):
  need(0<=cell<self.total,'invalid cell identity');i=bisect.bisect_right(self.offsets,cell)-1;k=cell-self.offsets[i];z=self.z_families[self.z_counts[i]];return Fraction(area(self.leaves[i])*(z[k+1]-z[k]),10**27)
 def side_membership(self,p,foil,axis,normal,valid_stage):
  need(valid_stage=='side_resolved_vacuum_geometry_only','missing side/stage');s=SOLIDS[foil];plane=s[2*axis+1] if normal==1 else s[2*axis];tangent=1-axis
  need(normal in [-1,1] and p[axis]==plane and s[2*tangent]<p[tangent]<s[2*tangent+1] and Z[0]<p[2]<Z[1],'not an open finite material face patch')
  candidates=[]
  for i,b in enumerate(self.leaves):
   boundary=b[2*axis] if normal==1 else b[2*axis+1]
   if boundary==plane and b[2*tangent]<=p[tangent]<b[2*tangent+1]:candidates.append(i)
  need(len(candidates)==1,'ambiguous side-resolved membership');i=candidates[0];z=self.z_families[self.z_counts[i]];k=min(len(z)-2,bisect.bisect_right(z,p[2])-1)
  return dict(cell_id=self.offsets[i]+k,physical_point_unchanged=True,operator_collision_admission='not supplied; geometry membership only')
def rejected(fn):
 try:fn()
 except ValueError:return True
 return False
def main():
 p=argparse.ArgumentParser();p.add_argument('--out',type=pathlib.Path,required=True);a=p.parse_args();need(not a.out.exists(),'fresh output required');need(os.sched_getaffinity(0)=={7} and os.getpriority(os.PRIO_PROCESS,0)>=15,'CPU7 nice15 required');resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));signal.alarm(90);start=time.monotonic();a.out.mkdir();checks=[]
 source_files=[ROOT/'source/native_parent_subflow.hpp',ROOT/'inputs/collision_scheduler.hpp',ROOT/'tests/test_native_parent_subflow.cpp'];pins={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files}
 cases=[('normal2um',dict(normal_nm=2000)),('normal1um',{}),('normal0p5um',dict(normal_nm=500)),('normal0p25um',dict(normal_nm=250)),('tangent-refine',dict(tangent_nm=62500,side_tangent_nm=12500)),('axial250um',dict(nz=100)),('axial125um',dict(nz=200)),('growth-refine',dict(normal_growth=.15)),('local-z500um',dict(nz=25,near_nz=50,near_band_nm=25000)),('local-z250um',dict(nz=25,near_nz=100,near_band_nm=25000)),('local-z125um',dict(nz=25,near_nz=200,near_band_nm=25000)),('local-z-band50um',dict(nz=25,near_nz=100,near_band_nm=50000))];results=[];reference=None;local_reference=None
 for name,opts in cases:
  t=time.monotonic()
  try:
   g=Grid(**opts);broad=[b[3]-b[2] for b in g.leaves if b[0]<7000000 and b[1]>2000000 and (b[3] in [1800000,7000000] or b[2] in [2000000,7200000])];result=dict(name=name,status='CPU_GEOMETRY_PASS',options=g.options,xy_leaves=len(g.leaves),tree_nodes=g.nodes,axial_bin_counts=sorted(g.z_families),xy_rectangles_per_z_family={n:g.z_counts.count(n) for n in g.z_families},physical_3D_vacuum_bins=g.total,physical_vacuum_volume_m3=float(Fraction(g.xy_area*(Z[1]-Z[0]),10**27)),broad_face_first_normal_width_nm_min=min(broad),broad_face_first_normal_width_nm_max=max(broad),estimated_cpp_Cell_vector_bytes=16*g.total,estimated_2D_bounds_plus_tree_bytes=32*len(g.leaves)+32*g.nodes,identity=g.identity,wall_s=time.monotonic()-t)
   if name in ['normal1um','local-z250um']:
    if name=='normal1um':reference=g
    else:local_reference=g
    (a.out/(name+'-layout.json')).write_text(json.dumps(dict(identity=g.identity,coordinate_unit='integer_nm_manufactured_only',options=g.options,xy_leaves=g.leaves,z_families=g.z_families,z_counts_by_xy_leaf=g.z_counts,solid_boxes_xy=SOLIDS),indent=2)+'\n')
  except ValueError as error:
   if 'cap' not in str(error):raise
   result=dict(name=name,status='REFUSED_RESOURCE_CAP',reason=str(error),requested_options=opts,wall_s=time.monotonic()-t)
  results.append(result)
 need(reference is not None,'reference microcase refused');g=reference
 # Exact volume, source boundaries and no artificial vacuum strip.
 need(sum((g.volume(i) for i in range(g.total)),Fraction())==Fraction(1975,10**9),'full physical volume including25mm once');checks.append('Exact3D volume closure1.975e-6m3; no observation slab or cut cells')
 for point in [(4000000,1900000,12500000),(4000000,2000000,12500000),(4000000,7000000,12500000)]:need(rejected(lambda:g.locate(point)),'solid or exact material point admitted')
 pface=(4000000,2000000,12500000);side=g.side_membership(pface,0,1,1,'side_resolved_vacuum_geometry_only');need(side['physical_point_unchanged'] and side['operator_collision_admission'].startswith('not supplied'),'boundary point incorrectly evolved/admitted');need(g.locate((4000000,2000001,12500000))==side['cell_id'],'one-nm vacuum point lost to an excluded strip');checks.append('Solid/exact surface refused; explicit side identifies adjacent vacuum bin without moving point;1nm vacuum retained')
 # Geometry does not infer or alter microscopic Active operator ownership.
 checks.append('Contract: geometry returns membership only; Active microscopic-layer gas-reaction admission is not implemented or decided here')
 rng=random.Random(20261029);located=0;occupied={}
 for _ in range(20000):
  point=tuple(rng.randrange(1,high) for high in [9000000,9000000,25000000])
  if material(point)!='vacuum':continue
  cid=g.locate(point);need(g.volume(cid)>0,'nonpositive vacuum bin');located+=1;occupied[cid]=occupied.get(cid,0)+1
 need(g.locate((4000000,4000000,1))!=g.locate((4000000,4000000,24999999)),'axial slices collapsed');checks.append(f'{located}manufactured vacuum points uniquely located; full-height axial groups distinct')
 need(sum(occupied.values())==located and len(occupied)<=located and all(g.volume(i)>0 for i in occupied),'sparse occupied census mismatch');checks.append('Occupied-cell view retains every manufactured carrier and its original physical cell volume; only exact empty bins are omitted')
 need(local_reference is not None,'local-z reference refused');local=local_reference;need(sum((local.volume(i) for i in range(local.total)),Fraction())==Fraction(1975,10**9),'local-z volume closure');checks.append('Per-rectangle local z families preserve exact full-height vacuum volume; refinement band changes resolution only, never excludes gas')
 # Existing uniform grid's manufactured45x45x5 geometry, not native FP parity.
 need((45*45-2*25)*5==9875 and Fraction(9000000,45)**2*Fraction(25000000,5)/10**27==Fraction(2,10**10),'uniform geometry control');checks.append('Original45x45x5 manufactured volume/count control; live C++ uniform parity remains unrun')
 # Uniform densities: every pair moment is proportional to volume, exactly.
 ne=Fraction(7,3);nn=Fraction(11,5);k=Fraction(13,7);rate=sum((ne*nn*k*g.volume(i) for i in range(g.total)),Fraction());need(rate==ne*nn*k*Fraction(1975,10**9),'uniform rate volume invariant');checks.append('Exact uniform density integral invariant across variable-volume bins; no event selected')
 # A localized overlapping layer is diluted; disjoint layers can falsely overlap.
 A=Fraction(1,10**6);delta=Fraction(1,10**6);coarse=Fraction(25,10**6);Ne=Fraction(3);Nn=Fraction(5);fine=Ne*Nn/(A*delta);mixed=Ne*Nn/(A*coarse);need(mixed/fine==Fraction(1,25),'layer overlap microcase');need(Ne*Nn/(A*coarse)>0,'disjoint coarse false-overlap microcase');checks.append('1um co-located layer coarse product is25x too small; separated sublayers can acquire false positive overlap')
 need(rejected(lambda:Grid(cap=1)),'cell cap did not refuse');need(all(hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()==h for p,h in pins.items()),'owner source changed during read')
 report=dict(schema='collision_grid_design_cpu_microcase/v2',status='PASS_GEOMETRY_DESIGN_ONLY',cases=results,checks=checks,inputs=pins,manufactured_lookup_points=located,sparse_occupied_example=dict(complete_manufactured_census_count=located,occupied_cells=len(occupied),full_physical_geometry_bins=g.total,all_positive_entries_retained=True,physical_volume_never_renormalized_to_occupied_volume=True),reference_identity=g.identity,reference_geometry_frame='Pilot canonical [0,9]mmx[0,9]mmx[0,25]mm; source/field adapters must bind their actual frame once',physical_dimensions=dict(foil_width_m=.005,foil_thickness_m=.0002,gap_m=.005,axial_height_m=.025,vacuum_xy_width_m=.009),no_native_particle_quantization=True,native_driver_port_compiled=False,live_driver_edit=False,GPU=False,Poisson=False,evolution=False,Debye_or_sheath_convergence=False,execution=dict(affinity=sorted(os.sched_getaffinity(0)),nice=os.getpriority(os.PRIO_PROCESS,0),wall_s=time.monotonic()-start,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,address_space_cap_bytes=4*1024**3,time_cap_s=90))
 (a.out/'receipt.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(status=report['status'],cases=results,execution=report['execution']),indent=2))
if __name__=='__main__':main()
