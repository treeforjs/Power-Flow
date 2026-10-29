"""Render actual MFEM cathode cell averages and material-face temperature traces."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageDraw


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_fields(directory, configuration):
    path = directory/"foil_snapshots.csv"
    sidecar = json.loads((directory/"foil_snapshots.json").read_bytes())
    if sidecar["sampling"] != "cell_average" or sha(path) != sidecar["csv_sha256"]:
        raise ValueError("snapshot sampling/hash mismatch")
    nx, ny = configuration["width_cells"], configuration["thickness_cells"]
    width, thickness = configuration["foil_width_m"], configuration["foil_thickness_m"]
    frames = {}
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            stamp = float(row["time_s"])
            ix = round((float(row["x_left_m"])+width/2)*nx/width)
            iy = round(float(row["depth_left_m"])*ny/thickness)
            value = float(row["temperature_k"])
            if not (0<=ix<nx and 0<=iy<ny and math.isfinite(value)):
                raise ValueError("unsupported snapshot cell")
            expected = {"x_left_m":-width/2+width*ix/nx,"x_right_m":-width/2+width*(ix+1)/nx,
                        "depth_left_m":thickness*iy/ny,"depth_right_m":thickness*(iy+1)/ny}
            if any(not math.isclose(float(row[key]),target,abs_tol=1e-14,rel_tol=1e-9) for key,target in expected.items()):
                raise ValueError("snapshot geometry does not match rectangular material mesh")
            field = frames.setdefault(stamp,np.full((ny,nx),np.nan))
            if not np.isnan(field[iy,ix]):
                raise ValueError("duplicate material cell")
            field[iy,ix] = value
    if not frames or any(not np.all(np.isfinite(field)) for field in frames.values()):
        raise ValueError("incomplete actual material snapshots")
    return dict(sorted(frames.items()))


def read_surfaces(directory):
    surfaces = {}
    for face in ("top","bottom","left","right"):
        config_path = directory/(face+"_surface.json")
        if not config_path.exists():
            continue
        config = json.loads(config_path.read_bytes())
        path = directory/config["table"]["path"]
        if sha(path) != config["table"]["sha256"]:
            raise ValueError("material-face table hash mismatch")
        times, edges = config["time_s"],np.asarray(config["x_edges_m"])
        values = np.full((len(times),len(edges)-1),np.nan)
        lookup = {stamp:i for i,stamp in enumerate(times)}
        baseline = config["temperature"]["initial_temperature_k"]
        if config["temperature"]["value_basis"] != "temperature_increment_kelvin":
            raise ValueError("expected native temperature increments")
        with path.open(newline="") as stream:
            for row in csv.DictReader(stream):
                i,j = lookup[float(row["time_s"])],int(row["x_index"])
                if int(row["z_index"]) != 0 or not np.isnan(values[i,j]):
                    raise ValueError("unsupported or repeated boundary patch")
                values[i,j] = baseline+float(row["temperature_value"])
        if not np.all(np.isfinite(values)):
            raise ValueError("incomplete material-face history")
        surfaces[face] = {"times":times,"edges":edges,"temperature_k":values,"config_sha256":sha(config_path)}
    return surfaces


def selected_indices(times,stride):
    selected = set(range(0,len(times),stride)) | {0,len(times)-1}
    selected.update(i for i,stamp in enumerate(times) if abs(stamp-40e-9)<1e-18 or abs(stamp-60e-9)<1e-18)
    return sorted(selected)


def capture(figure):
    figure.canvas.draw()
    return Image.fromarray(np.asarray(figure.canvas.buffer_rgba())[:,:,:3].copy())


def save_gif(frames,path,fps):
    palette = frames[-1].quantize(colors=256)
    indexed = [frame.quantize(palette=palette,dither=Image.Dither.NONE) for frame in frames]
    indexed[0].save(path,save_all=True,append_images=indexed[1:],duration=max(20,round(1000/fps)),loop=0,
                    optimize=False,disposal=2)
    with Image.open(path) as gif:
        count = gif.n_frames
        if count != len(frames):
            raise ValueError("GIF writer changed the number of native frames")
        positions = [0,count//2,count-1]
        previews = []
        for label,index in zip(("first","middle","last"),positions):
            gif.seek(index)
            frame = gif.convert("RGB")
            frame.save(path.with_name(path.stem+"_decoded_"+label+".png"))
            preview = frame.copy();preview.thumbnail((660,440));previews.append((label,index,preview))
        contact = Image.new("RGB",(660*3,470),"white")
        draw = ImageDraw.Draw(contact)
        for i,(label,index,preview) in enumerate(previews):
            contact.paste(preview,(660*i,25));draw.text((660*i+10,5),f"Decoded {label} frame {index}",fill="black")
        contact.save(path.with_name(path.stem+"_decoded_contact.png"))
    return count


def render_fields(fields,configuration,output,vmax,stride,fps):
    times = list(fields)
    indices = selected_indices(times,stride)
    width,thickness = configuration["foil_width_m"]*1e3,configuration["foil_thickness_m"]*1e3
    fig,ax = plt.subplots(figsize=(11.6,5.8),dpi=110)
    fig.subplots_adjust(left=.10,right=.87,bottom=.22,top=.80)
    image = ax.imshow(fields[times[0]],origin="lower",extent=(-width/2,width/2,0,thickness),
                      interpolation="nearest",cmap="inferno",vmin=300,vmax=vmax,aspect=10)
    ax.set_xlabel("Foil width x (mm)")
    ax.set_ylabel("Depth inward from gap face (mm)")
    ax.set_yticks([0,.05,.10,.15,.20])
    ax.text(0,thickness+.011,"Back face: physical y = +2.7 mm",ha="center",fontsize=10)
    ax.text(0,-.033,"Gap-facing face: physical y = +2.5 mm",ha="center",fontsize=10)
    fig.colorbar(image,ax=ax,pad=.035,fraction=.04,label="Foil temperature (K), fixed linear scale")
    heading = fig.suptitle("",fontsize=15,y=.96)
    fig.text(.5,.885,"Upper cathode | actual MFEM cell averages | thickness expanded 10x",ha="center",fontsize=11)
    note = ("Only two recorded interior snapshots; no intervening field is reconstructed."
            if len(times)==2 else "Selected actual native timesteps; no temporal interpolation.")
    fig.text(.5,.10,note,ha="center",fontsize=10)
    fig.text(.5,.047,"Provisional constant-property stationary solid; hot regions exceed model validity. Cathode: dissertation Figs. 4.27-4.28.",
             ha="center",fontsize=9,color="#8c240c")
    frames = []
    for i in indices:
        stamp = times[i];image.set_data(fields[stamp])
        heading.set_text(f"C.19 shot 17133 | cathode temperature | t = {stamp*1e9:.3f} ns")
        frames.append(capture(fig))
    plt.close(fig)
    path = output/"foil_temperature.gif"
    count = save_gif(frames,path,min(fps,1) if len(times)==2 else fps)
    return {"path":str(path),"frame_count":count,"actual_native_times_s":[times[i] for i in indices],
            "available_native_snapshots":len(times),"temperature_min_k":min(float(np.min(v)) for v in fields.values()),
            "temperature_max_k":max(float(np.max(v)) for v in fields.values()),"sha256":sha(path)}


def render_surfaces(surfaces,output,vmax,stride,fps):
    times = surfaces["top"]["times"]
    if any(data["times"] != times for data in surfaces.values()):
        raise ValueError("material faces have different native clocks")
    indices = selected_indices(times,stride)
    fig,axes = plt.subplots(1,2,figsize=(11.6,5.8),dpi=110)
    fig.subplots_adjust(left=.08,right=.98,bottom=.22,top=.77,wspace=.25)
    heading = fig.suptitle("",fontsize=15,y=.96)
    lines = {}
    labels = {"top":"Back face","bottom":"Gap-facing face","left":"Left width edge","right":"Right width edge"}
    for face,data in surfaces.items():
        ax = axes[0 if face in ("top","bottom") else 1]
        centers = (data["edges"][:-1]+data["edges"][1:])/2*1e3
        if face in ("left","right"):
            centers -= 2.5
        lines[face], = ax.plot(centers,data["temperature_k"][0],label=labels[face],lw=2)
    for i,ax in enumerate(axes):
        ax.set_ylim(300,vmax);ax.set_ylabel("Material-face temperature (K)");ax.grid(alpha=.2)
        ax.set_xlabel("Foil width x (mm)" if i==0 else "Depth inward from gap face (mm)")
        if ax.lines:
            ax.legend(loc="upper center",fontsize=9)
    if not {"left","right"}.issubset(surfaces):
        axes[1].set_axis_off()
        axes[1].text(.5,.5,"Side traces were not saved\nin this earlier run.",transform=axes[1].transAxes,ha="center")
    fig.text(.5,.85,"Upper cathode | actual one-sided material boundary patch averages",ha="center",fontsize=11)
    fig.text(.5,.10,"Exact recorded native frames; fixed linear Kelvin axes. Interior field is not reconstructed from these traces.",ha="center",fontsize=10)
    fig.text(.5,.047,"Provisional constant-property stationary solid; hot regions exceed model validity.",ha="center",fontsize=9,color="#8c240c")
    frames = []
    for i in indices:
        for face,line in lines.items():
            line.set_ydata(surfaces[face]["temperature_k"][i])
        heading.set_text(f"C.19 shot 17133 | cathode surface temperature | t = {times[i]*1e9:.3f} ns")
        frames.append(capture(fig))
    plt.close(fig)
    path = output/"foil_surface_temperature.gif"
    count = save_gif(frames,path,fps)
    return {"path":str(path),"frame_count":count,"actual_native_times_s":[times[i] for i in indices],"sha256":sha(path),
            "faces":{face:{"temperature_min_k":float(np.min(data["temperature_k"])),"temperature_max_k":float(np.max(data["temperature_k"]))}
                     for face,data in surfaces.items()}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-run",type=Path,required=True)
    parser.add_argument("--output-dir",type=Path,required=True)
    parser.add_argument("--frame-stride",type=int,default=1)
    parser.add_argument("--surface-frame-stride",type=int,default=3)
    parser.add_argument("--fps",type=float,default=10)
    parser.add_argument("--vmax-k",type=float)
    args = parser.parse_args()
    if args.output_dir.exists() or args.frame_stride<1 or args.surface_frame_stride<1 or args.fps<=0:
        parser.error("use a fresh output directory, positive strides and positive fps")
    run_path = args.source_run/"source_run.json"
    run = json.loads(run_path.read_bytes())
    if not run["passed"]:
        parser.error("source native solve did not pass")
    fields,surfaces = read_fields(args.source_run,run["configuration"]),read_surfaces(args.source_run)
    maximum = max(max(float(np.max(v)) for v in fields.values()),max(float(np.max(d["temperature_k"])) for d in surfaces.values()))
    vmax = args.vmax_k or math.ceil(maximum/1000)*1000
    if not math.isfinite(vmax) or vmax<maximum:
        parser.error("fixed color maximum must include all computed temperatures")
    args.output_dir.mkdir(parents=True)
    record = {"schema":"mhdlab-mfem-foil-temperature-render-v1","status":"provisional",
        "physical_experiment_prediction":False,"source_run":str(args.source_run),"source_run_sha256":sha(run_path),
        "snapshot_csv_sha256":sha(args.source_run/"foil_snapshots.csv"),"temperature_scale":{"basis":"absolute_kelvin","vmin_k":300,"vmax_k":vmax,"normalization":"fixed_linear_all_frames"},
        "electrode":{"modeled_foil":"cathode","geometry":"upper foil, physical y=+2.5..+2.7 mm",
            "evidence":"Dissertation printed pp.93-94, Figs.4.27-4.28: explicit upper Cathode label and cathode profiles at positive depth",
            "pdf_sha256":"1d7b1f9c191a37e94f1cc3634587750f1af35b07eb13e0b94c60f3da47495ed7"},
        "temporal_interpolation":False,"thickness_display_expansion":10,
        "foil":render_fields(fields,run["configuration"],args.output_dir,vmax,args.frame_stride,args.fps),
        "surfaces":render_surfaces(surfaces,args.output_dir,vmax,args.surface_frame_stride,args.fps)}
    (args.output_dir/"render.json").write_text(json.dumps(record,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(args.output_dir/"render.json")


if __name__ == "__main__":
    main()
