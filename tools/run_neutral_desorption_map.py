"""Compose native MFEM surface inputs, conservative launches and vacuum maps."""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.meta_path[:] = [f for f in sys.meta_path if type(f).__module__ != "_editable_skbc_mhdlab"]
from mhdlab.constants import SPECIES_MASS_KG
from mhdlab.mfem_surface import load_mfem_surface_samples
from mhdlab.neutral_vacuum_map import SPECIES, observe_vacuum_packets
from mhdlab.neutral_source_sampling import prepare_sampling_interval
from mhdlab.surface_coordinates import to_physical_coordinates, validate_mfem_face_axes
from mhdlab.temkin_planar import evolve_temkin_planar_release
from mhdlab.temkin_planar_boundary import sample_temkin_planar_emission_packets
from mhdlab.temkin_planar_io import write_temkin_planar_release_history
from mhdlab.temkin_release import TemkinReservoir


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def open_draws(rng, number):
    values = rng.random((number, 7))
    while np.any(values == 0):
        selected = values == 0
        values[selected] = rng.random(np.count_nonzero(selected))
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--top-surface", type=Path, required=True)
    parser.add_argument("--bottom-surface", type=Path, required=True)
    parser.add_argument("--left-surface", type=Path)
    parser.add_argument("--right-surface", type=Path)
    parser.add_argument("--case", type=Path, default=ROOT / "configs/neutral_desorption_map_provisional.json")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--resolved-only-source-history", action="store_true",
                        help="archive validated local release without optional spatially uniform comparison controls")
    args = parser.parse_args()
    output = args.output_dir
    if output.exists():
        parser.error("use a new output directory")
    case = json.loads(args.case.read_bytes())
    if case.get("schema") != "mhdlab-neutral-desorption-map-case-v1" or case.get("status") != "provisional":
        parser.error("this composed demonstration requires explicit provisional status")
    if tuple(s["name"] for s in case["species"]) != SPECIES:
        parser.error("four independent source reservoirs H,H2,H2O,OH are required in that order")
    if (args.left_surface is None) != (args.right_surface is None):
        parser.error("left and right material faces must be supplied together")
    input_paths = [args.top_surface, args.bottom_surface]
    labels = ["top", "bottom"]
    if args.left_surface is not None:
        input_paths += [args.left_surface, args.right_surface]
        labels += ["left", "right"]
    loaded = [load_mfem_surface_samples(path) for path in input_paths]
    surfaces = [item[0] for item in loaded]
    top, bottom = surfaces[:2]
    if any(not np.array_equal(getattr(top, key), getattr(bottom, key)) for key in ("time_s", "x_edges_m", "z_edges_m")):
        parser.error("both material faces must share native times and x/z partitions")
    faces = [surface.attributes["mfem_face"] for surface in surfaces]
    if any(validate_mfem_face_axes(face) != (0, 1, 2) for face in faces[:2]):
        parser.error("top/bottom must use the physical y normal")
    if faces[0]["outward_normal_sign"] != 1 or faces[1]["outward_normal_sign"] != -1:
        parser.error("top/bottom configurations must retain their respective outward signs")
    xmin, xmax = map(float, top.x_edges_m[[0, -1]])
    ymin, ymax = float(faces[1]["normal_position_m"]), float(faces[0]["normal_position_m"])
    if ymin >= ymax:
        parser.error("top material face must lie above bottom face")
    for label, surface, face, coordinate, sign in zip(labels[2:], surfaces[2:], faces[2:], (xmin, xmax), (-1, 1)):
        if (validate_mfem_face_axes(face) != (1, 0, 2)
                or face["normal_position_m"] != coordinate or face["outward_normal_sign"] != sign
                or not np.array_equal(surface.x_edges_m[[0, -1]], [ymin, ymax])
                or any(not np.array_equal(getattr(top, key), getattr(surface, key)) for key in ("time_s", "z_edges_m"))):
            parser.error(f"{label} face geometry/time must match the broad faces and retain its physical x normal")
    if case.get("emitting_faces", labels) != labels:
        parser.error("supplied material faces must match the case's emitting_faces")
    length = float(top.z_edges_m[-1] - top.z_edges_m[0])
    pools = []
    law = case["kinetic_hypothesis"]
    for species in SPECIES:
        pools.append(TemkinReservoir(name=f"assumed_{species}_pool", state_id=f"{species}:assumed_source_state",
            **{key: law[key] for key in ("attempt_frequency_s", "activation_energy_ev", "desorption_energy_ev",
                "internal_energy_ev", "coverage_coefficient", "reference_coverage_m2", "max_coverage_m2")},
            kinetic_provenance=law["provenance"], energy_provenance="Illustrative independent removal/internal energy accounts; no reaction energy or fragment recoil model"))
    state_masses = {p.state_id: SPECIES_MASS_KG[name] for p, name in zip(pools, SPECIES)}
    state_indices = {p.state_id: index for index, p in enumerate(pools)}
    rng = np.random.default_rng(case["seed"])
    output.mkdir(parents=True)
    write(output / "case.json", case)
    blocks, source_records = [], []
    initial_particles = np.zeros(4)
    for face_index, (label, surface, face, (_, extraction)) in enumerate(zip(labels, surfaces, faces, loaded)):
        shape = (4, len(surface.z_edges_m)-1, len(surface.x_edges_m)-1)
        coverage = np.stack([np.full(shape[1:], row["initial_coverage_m2"]) for row in case["species"]])
        source = evolve_temkin_planar_release(surface, coverage, pools, parameter_provenance=law["provenance"])
        area = source.patch_areas_m2
        initial_particles += np.sum(source.coverage_m2[0] * area, axis=(1, 2))
        write_temkin_planar_release_history(output / f"{label}_release.h5", source,
            include_width_uniform=not args.resolved_only_source_history)
        expected, sampled = 0., 0.
        numerical_omissions = []
        face_block_start = len(blocks)
        for interval in range(len(source.time_edges_s)-1):
            sampling_source, sampling_index, sampling_report = prepare_sampling_interval(source, interval,
                minimum_component_particles=case.get("minimum_packet_component_particles", 1e-250))
            if sampling_report["omitted_component_count"]:
                numerical_omissions.append(sampling_report)
            if sampling_source is None:
                continue
            positive = int(np.count_nonzero(sampling_source.released_m2[sampling_index] > 0))
            if not positive:
                continue
            count = max(int(case.get("side_packets_per_interval_per_face", case["packets_per_interval_per_face"])
                            if label in {"left", "right"} else case["packets_per_interval_per_face"]), positive)
            packets = sample_temkin_planar_emission_packets(sampling_source, state_masses, sampling_index,
                open_draws(rng, count), face_normal_coordinate_m=face["normal_position_m"],
                outward_normal_sign=face["outward_normal_sign"])
            blocks.append({"birth_time_s": packets.time_s, "position_m": to_physical_coordinates(packets.resolved_position_m, face),
                "velocity_m_s": to_physical_coordinates(packets.velocity_m_s, face), "particle_weight": packets.particle_weight,
                "mass_kg": packets.mass_kg,
                "face_index": np.full(len(packets.time_s), face_index, dtype=np.int64),
                "temperature_k": packets.temperature_k,
                "species_index": np.asarray([state_indices[state] for state in packets.state_ids], dtype=np.int64)})
            expected += packets.expected_translational_energy_j
            sampled += packets.sampled_translational_energy_j
            if interval % 10 == 0:
                print(f"{label}: sampled interval {interval+1}/{len(source.time_edges_s)-1}", flush=True)
        source_records.append({"face": label, "source_attributes": surface.attributes, "extraction": extraction,
            "emitting_face_area_m2": float(np.sum(area)),
            "local_to_physical_component_permutation": list(validate_mfem_face_axes(face)),
            "packet_numerical_omissions": numerical_omissions,
            "initial_particles": np.sum(source.coverage_m2[0]*area, axis=(1, 2)).tolist(),
            "remaining_particles": np.sum(source.coverage_m2[-1]*area, axis=(1, 2)).tolist(),
            "released_particles": np.sum(source.released_m2*area, axis=(0, 2, 3)).tolist(),
            "expected_launch_translational_energy_j": expected, "actual_sampled_launch_translational_energy_j": sampled,
            "thermal_sampling_energy_residual_j": sampled-expected, "release_file_sha256": sha(output/f"{label}_release.h5")})
        # Preserve completed faces before starting another source or observer.
        face_blocks = blocks[face_block_start:]
        if face_blocks:
            checkpoint = output / f"{label}_source_packets.npz"
            np.savez_compressed(checkpoint, **{key: np.concatenate([block[key] for block in face_blocks])
                                              for key in face_blocks[0]})
            write(output / f"{label}_source_record.json", {**source_records[-1],
                "packet_checkpoint_sha256": sha(checkpoint), "rng_state_after_face": rng.bit_generator.state})
    if not blocks:
        raise ValueError("no source packets were emitted")
    packets = {key: np.concatenate([block[key] for block in blocks], axis=0) for key in blocks[0]}
    np.savez_compressed(output / "source_packets.npz", **packets)
    mx, my = case["outside_width_margin_m"], case["outside_normal_margin_m"]
    x_edges = np.concatenate((np.linspace(xmin-mx, xmin, case["width_bins_per_margin"]+1)[:-1],
        np.linspace(xmin, xmax, case["width_bins_inside"]+1)[:-1],
        np.linspace(xmax, xmax+mx, case["width_bins_per_margin"]+1)))
    y_edges = np.concatenate((np.linspace(ymin-my, ymin, case["normal_bins_per_margin"]+1)[:-1],
        np.linspace(ymin, ymax, case["normal_bins_inside_solid"]+1)[:-1],
        np.linspace(ymax, ymax+my, case["normal_bins_per_margin"]+1)))
    times = np.linspace(float(top.time_s[0]), float(top.time_s[-1]), case["frames"])
    print(f"Observing {len(packets['birth_time_s'])} packets at {len(times)} frames", flush=True)
    flight_packets = {key: value for key, value in packets.items() if key not in {"face_index", "temperature_k"}}
    arrays, ledgers = observe_vacuum_packets(**flight_packets, time_s=times, x_edges_m=x_edges, y_edges_m=y_edges,
        foil_bounds_m=[xmin, xmax, ymin, ymax], invariant_current_length_m=length,
        initial_surface_particles=initial_particles)
    np.savez_compressed(output / "neutral_maps.npz", **arrays)
    write(output / "particle_ledgers.json", ledgers)
    metadata = {"schema": "mhdlab-neutral-desorption-vacuum-map-v1", "status": "provisional",
        "case_title": case.get("case_title", "Provisional MFEM-heated foil | independent neutral surface reservoirs"),
        "time_registration": case.get("time_registration", "Native MFEM source clock; no experimental trigger registration"),
        "invariant_current_length_m": length,
        "normal_axis_magnification": case.get("normal_axis_magnification", 1.0),
        "model_description": "Native constant-property MFEM Joule heating, finite-inventory first-order desorption, flux-weighted 3V thermal launches and exact neutral vacuum flight",
        "assumptions": [law["provenance"], case["transport_hypothesis"], case["chemistry_hypothesis"],
            ("All four perimeter material faces emit with the same per-area kinetic and initial-coverage assumptions; side areas are not inflated."
             if len(labels) == 4 else "Only the two broad material faces emit; foil edges do not emit in this case."),
            "Current direction is invariant/periodic: z velocity is retained for energy, while maps integrate over the declared current length.",
            case.get("heating_assumption", "The native MFEM geometry, material parameters and drive require their separately declared source qualification; no shot-specific SS304 material fit is established."),
            "Surface removal energy is a passive diagnostic and is not fed back into MFEM temperature.",
            "Conditional source-event times use the existing exact frozen-temperature Temkin/first-order law; no frame-midpoint injection.",
            "Histograms use gas-side ownership for points exactly on material faces; trajectories themselves are not shifted.",
            "Packet representation omits only per-component tails below the explicitly recorded numerical particle cutoff; original law histories are preserved, and omitted count bounds are reported separately.",
            "No water chemistry, fragment recoil, collisions, radiation, plasma forces or calibrated density/temperature inference is claimed."],
        "source_records": source_records, "emitting_faces": labels,
        "packet_coordinate_system": "x=foil_width;y=gap_normal;z=current",
        "packet_face_index_labels": labels, "reservoirs": [asdict(pool) for pool in pools],
        "species_mass_kg": {name: SPECIES_MASS_KG[name] for name in SPECIES},
        "mass_provenance": "Existing project approximate species-mass table, recorded explicitly; not isotope-resolved",
        "initial_surface_particles": initial_particles.tolist(), "seed": case["seed"],
        "packet_count": len(packets["birth_time_s"]), "frame_count": len(times),
        "width_uniform_comparison_included": not args.resolved_only_source_history,
        "inputs": {str(path): sha(path) for path in (*input_paths, args.case)},
        "output_sha256": {p.name: sha(p) for p in output.iterdir() if p.is_file()},
        "source_sha256": {name: sha(ROOT/name) for name in ("tools/run_neutral_desorption_map.py",
             "mhdlab/neutral_vacuum_map.py", "mhdlab/temkin_planar_boundary.py", "mhdlab/planar_boundary.py",
             "mhdlab/temkin_planar.py", "mhdlab/temkin_boundary.py", "mhdlab/temkin_emission_time.py",
             "mhdlab/surface_coordinates.py", "mhdlab/neutral_source_sampling.py",
             "mhdlab/mfem_surface.py", "mhdlab/planar_release.py")}}
    write(output / "metadata.json", metadata)
    print(json.dumps({"output_dir": str(output), "packets": len(packets["birth_time_s"]), "frames": len(times),
                      "status": metadata["status"]}), flush=True)


if __name__ == "__main__":
    main()
