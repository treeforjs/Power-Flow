"""Recover 50 published mass-channel readouts; retain encoded PDF primitives separately."""
from pathlib import Path
import argparse,csv,hashlib,json
import numpy as np
import pymupdf

PDF_SHA='1d7b1f9c191a37e94f1cc3634587750f1af35b07eb13e0b94c60f3da47495ed7'
ASSIGNMENTS={1:'H-bearing fragment channel',2:'H2; possible D ambiguity',12:'C-bearing fragment; marker not listed in caption',14:'N or hydrocarbon fragment',16:'O-bearing fragment / CH4 ambiguity',17:'OH-bearing fragment / NH3 ambiguity',18:'H2O parent channel',28:'CO / N2 / hydrocarbon ambiguity',32:'O2 / other parent ambiguity',44:'CO2 / hydrocarbon ambiguity'}
def channel_readouts(rows):
 channels=[]
 for mass in range(1,51):
  covering=[r for r in rows if r['bar_left_amu']<=mass<=r['bar_right_amu']]
  nearby=[r for r in rows if abs(r['mass_to_charge_amu']-mass)<=.5]
  channels.append(dict(nominal_mass_channel_amu=mass,published_plot_envelope_pressure_torr=max((r['plotted_channel_pressure_torr'] for r in covering),default=None),encoded_center_nearest_pressure_torr=min(nearby,key=lambda r:abs(r['mass_to_charge_amu']-mass))['plotted_channel_pressure_torr'] if nearby else None,readout_at_plot_boundary=mass==50,neutral_species_partial_pressure_established=False))
 return channels

def recover(pdf):
 if hashlib.sha256(pdf.read_bytes()).hexdigest()!=PDF_SHA:raise ValueError('PDF differs from the source-bound dissertation')
 page=pymupdf.open(pdf)[43];drawings=page.get_drawings();xticks=[];yticks=[];bars=[];markers=[]
 for d in drawings:
  if d['fill'] and d['fill'][2]>d['fill'][0]+.2:
   if len(d['items'])!=1 or d['items'][0][0]!='re':raise ValueError('Unexpected blue source geometry')
   bars.append(d['rect'])
  if len(d['items'])!=1 or d['items'][0][0]!='l':continue
  _,a,b=d['items'][0]
  if d['color']==(0.,0.,0.) and .67<d['width']<.69:
   if abs(a.x-b.x)<1e-5 and 2.9<abs(a.y-b.y)<3.1:xticks.append(a.x)
   if abs(a.y-b.y)<1e-5 and 2.9<abs(a.x-b.x)<3.1:yticks.append(a.y)
  if d['color']==(1.,1.,1.) and abs(a.x-b.x)<1e-5 and abs(a.y-b.y)>200:markers.append(a.x)
 xticks=np.array(sorted(xticks));yticks=np.array(sorted(yticks))
 if len(xticks)!=6 or len(yticks)!=8 or len(bars)!=450:raise ValueError('Source geometry count differs')
 xf=np.polyfit(xticks,np.arange(0,51,10),1);yf=np.polyfit(yticks,-np.arange(4,12),1)
 rows=[]
 for i,r in enumerate(bars):
  m=float(np.polyval(xf,(r.x0+r.x1)/2));logp=float(np.polyval(yf,r.y0));p=10**logp
  rows.append(dict(encoded_bar_index=i,mass_to_charge_amu=m,plotted_channel_pressure_torr=p,log10_plotted_channel_pressure_torr=logp,bar_left_amu=float(np.polyval(xf,r.x0)),bar_right_amu=float(np.polyval(xf,r.x1)),pdf_bar_top_y_pt=r.y0,center_inside_displayed_mass_range=0<=m<=50,raw_instrument_data=False))
 marker_masses=[float(np.polyval(xf,x)) for x in markers]
 peaks=[]
 for mass,label in ASSIGNMENTS.items():
  near=[r for r in rows if abs(r['mass_to_charge_amu']-mass)<=.5]
  peak=max(near,key=lambda r:r['plotted_channel_pressure_torr']);center=min(near,key=lambda r:abs(r['mass_to_charge_amu']-mass))
  peaks.append(dict(nominal_mass_to_charge_amu=mass,interpretation=label,local_peak_mass_amu=peak['mass_to_charge_amu'],local_peak_channel_pressure_torr=peak['plotted_channel_pressure_torr'],nearest_center_channel_pressure_torr=center['plotted_channel_pressure_torr'],neutral_species_partial_pressure_established=False))
 channels=channel_readouts(rows)
 meta=dict(schema='mhdlab-dissertation-rga-recovery-v1',status='published_vector_geometry_recovered_not_raw_RGA',figure='2.3',printed_page=19,pdf_page_one_based=44,source_pdf_sha256=PDF_SHA,source_pdf=str(pdf.resolve()),mass_channel_count=50,encoded_graphics_rectangle_count=len(rows),independent_observation_count_inferred_from_rectangles=False,marker_masses_amu=marker_masses,axis_fit=dict(mass_polynomial=xf.tolist(),log10_pressure_polynomial=yf.tolist(),max_mass_tick_residual_amu=float(np.max(abs(np.polyval(xf,xticks)-np.arange(0,51,10)))),max_log10_pressure_tick_residual=float(np.max(abs(np.polyval(yf,yticks)+np.arange(4,12))))),user_pressure_context=dict(value_torr_approximate=5e-5,exact=False,normalization_applied=False),limitations=['Encoded bar heights reproduce the published plot, not raw acquisition records.','Encoded PDF rectangles overlap and are graphics provenance, not450 independent measured bars. User identifies50 plotted mass channels. Channel readouts use the visible encoded envelope at integer masses; these do not recover raw instrument records or calibrated neutral partial pressures.','Mass channels can include ionizer fragmentation and overlapping species. Peaks are not neutral molar fractions.','RGA model, electron energy, calibration, relative sensitivities, cracking matrix, measurement date, gas temperature and independent total-gauge reading were not recovered from this figure.','Geometric fit residuals are not experimental uncertainty.'],physical_neutral_composition_ready=False)
 return rows,channels,peaks,meta

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--pdf',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args()
 rows,channels,peaks,meta=recover(a.pdf);a.output_dir.mkdir(parents=True,exist_ok=False)
 for name,data in [('figure_2_3_rga_mass_channels.csv',channels),('figure_2_3_rga_pdf_rectangles.csv',rows),('figure_2_3_rga_peak_channels.csv',peaks)]:
  with (a.output_dir/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=data[0]);w.writeheader();w.writerows(data)
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 fig,ax=plt.subplots(figsize=(11,5.5));x=np.array([r['nominal_mass_channel_amu'] for r in channels]);y=np.array([r['published_plot_envelope_pressure_torr'] if r['published_plot_envelope_pressure_torr'] is not None else np.nan for r in channels]);ax.bar(x,y,width=.8,color='#1f77b4');ax.set_yscale('log');ax.set(xlim=(0,50),ylim=(2e-12,2e-4),xlabel='Mass-to-charge channel (amu/e; dissertation axis: AMU)',ylabel='Published channel pressure (Torr)',title='Dissertation Figure2.3 | 50 nominal mass-channel readouts')
 for row in peaks:
  m=row['nominal_mass_to_charge_amu'];ax.axvline(m,color='white',ls='--',lw=1);ax.text(m,1.3e-4,str(m),ha='center',fontsize=8)
 fig.text(.12,.01,'50 channel readouts from the published envelope; missing readouts stay unknown. Neutral-species fractions require RGA response correction.',fontsize=9);fig.tight_layout(rect=(0,.035,1,1));fig.savefig(a.output_dir/'figure_2_3_rga_recovery.png',dpi=160)
 meta['outputs_sha256']={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in a.output_dir.iterdir() if f.is_file()};(a.output_dir/'figure_2_3_rga_recovery.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(peaks,indent=2))
if __name__=='__main__':main()
