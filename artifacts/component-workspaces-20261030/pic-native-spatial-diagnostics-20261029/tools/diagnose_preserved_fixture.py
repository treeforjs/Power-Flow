#!/usr/bin/env python3
"""Read an immutable legacy native midpoint; do not construct PIC or solve."""
import argparse,pathlib,sys,json,hashlib,struct,math
P=pathlib.Path(__file__).resolve().parents[1];sys.path.insert(0,str(P/'source'))
from diagnostics import diagnose,kinetic,QE,canonical_hash
PARTICLE=struct.Struct('<Q9d');OPTIONS=struct.Struct('<ii8d??6x8d?7xdd');CELL=struct.Struct('<4i8d')
def main():
 p=argparse.ArgumentParser();p.add_argument('--fixture',type=pathlib.Path,required=True);p.add_argument('--out-dir',type=pathlib.Path,required=True);a=p.parse_args();assert not a.out_dir.exists()
 files={name:(a.fixture/name).read_bytes() for name in ['options.bin','particles.bin','self_cells.bin','manifest.json']};assert len(files['options.bin'])==OPTIONS.size==168
 o=OPTIONS.unpack(files['options.bin']);nx,ny=o[:2];width,height,depth=o[2:5];assert nx==ny==360 and width==height==.009 and depth==.025 and o[11] and o[20]
 lower=list(o[12:16]);upper=list(o[16:20]);short,feed=o[21:23];manifest=json.loads(files['manifest.json']);assert len(files['particles.bin'])%PARTICLE.size==0 and len(files['particles.bin'])//PARTICLE.size==46726
 entries=[];mass=None
 for values in PARTICLE.iter_unpack(files['particles.bin']):
  nid,*v=values;r=v[:3];u=v[3:6];q,m,w=v[6:9];assert q==-QE and (mass is None or mass==m);mass=m;entries.append(dict(packet_id=f'legacy-native-id/{nid}',state='electron',native=dict(id=nid,position_m=r,proper_velocity_m_s=u,charge_C=q,mass_kg=m,weight=w),kinetic_J_per_particle=kinetic(m,u)))
 rawcells=files['self_cells.bin'];assert len(rawcells)==nx*ny*CELL.size
 xs=[CELL.unpack_from(rawcells,i*CELL.size)[4] for i in range(nx)];xs.append(CELL.unpack_from(rawcells,(nx-1)*CELL.size)[5]);ys=[CELL.unpack_from(rawcells,j*nx*CELL.size)[8] for j in range(ny)];ys.append(CELL.unpack_from(rawcells,(ny-1)*nx*CELL.size)[11])
 solids=[[f[0],f[1],f[2],f[3],short,feed] for f in [lower,upper]];surfaces=[]
 for electrode,foil in [('lower',lower),('upper',upper)]:
  l,r,b,t=foil
  for name,axis,plane,sign in [('back',1,b,-1),('gap',1,t,1)] if electrode=='lower' else [('gap',1,b,-1),('back',1,t,1)]:surfaces.append(dict(id=electrode+'/'+name,axis=axis,normal_sign=sign,plane_m=plane,bounds_m=[l,r,plane,plane,short,feed],vacuum_band_m=25e-6))
  for name,plane,sign in [('left',l,-1),('right',r,1)]:surfaces.append(dict(id=electrode+'/'+name,axis=0,normal_sign=sign,plane_m=plane,bounds_m=[plane,plane,b,t,short,feed],vacuum_band_m=25e-6))
 pins={name:hashlib.sha256(raw).hexdigest() for name,raw in files.items()}
 data=dict(schema='native_spatial_diagnostics_input/v1',state_prefix=canonical_hash(pins),time_s=manifest['absolute_native_time_s'],maximum_entries=46726,states=[dict(id='electron',mass_kg=mass,charge_e=-1)],entries=entries,geometry=dict(xs_m=xs,ys_m=ys,zs_m=[short,feed],solid_boxes_m=solids),surfaces=surfaces,energy_edges_eV=[0.]+[10**(-6+i/10) for i in range(121)],observation_slab_m=[-.0005,.0005])
 result=diagnose(data);result['legacy_fixture']={'input_hashes':pins,'scope':'Actual old electron-only midpoint, no arena lineage/sink history or gas/ion population inferred','native_ABI':'little endian Particle80 Options168 Cell80','actual_fixture_evolution':False};a.out_dir.mkdir();(a.out_dir/'diagnostics.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
 with (a.out_dir/'species-density-cm3.csv').open('w') as out:
  out.write('cell_id,state,physical_3D_density_cm3,observation_slab_density_cm3\n')
  for row in result['species_cell_densities']:out.write(f"{row['cell_id']},{row['state']},{row['density_cm3']},{row.get('observation',{}).get('density_cm3','')}\n")
 print('Preserved midpoint diagnostics only; no PIC construction/solve/evolution')
if __name__=='__main__':main()
