"""Run a manufactured prescribed-field transport/state-exchange reference.

The molecule M and every collision coefficient in this case are manufactured
verification inputs. This is not a measured shot, a water mechanism, a field
solve or a prediction of Hall feedback. It exercises executable coupling.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mhdlab.kinetic_gap import GapState, Particle, advance_gap, cell_moments, total_energy_j
from mhdlab.state_kinetics import ELEMENTARY_CHARGE_C as Q, InternalState, ReactionChannel, Species, StateCatalog


def reference_case():
    me, mass = 9.1093837139e-31, 3e-26
    catalog = StateCatalog(
        [Species("e", me, -1, {}), Species("M", mass, 0, {"X": 2, "Y": 1}),
         Species("M+", mass - me, 1, {"X": 2, "Y": 1})],
        [InternalState("e", "e", 0, {}),
         InternalState("M0", "M", 0, {"modes": [0, 0, 0], "J": 0}),
         InternalState("M1", "M", 1, {"modes": [1, 0, 0], "J": 0}),
         InternalState("Mi", "M+", 2, {"modes": [0, 0, 0], "N": 0})],
    )
    provenance = "Manufactured analytic verification channel; no physical molecular rate asserted"
    metadata = {"kinematics_provenance": "Explicit isotropic COM reference assumption, not measured molecular scattering"}
    channels = (
        ReactionChannel("excitation", ("e", "M0"), ("e", "M1"),
                        (0, 1, 2, 100), (0, 0, 1e-18, 1e-18),
                        "jacobi_energy_partition", provenance, metadata),
        ReactionChannel("deexcitation", ("e", "M1"), ("e", "M0"),
                        (0, 100), (1e-18, 1e-18),
                        "jacobi_energy_partition", provenance, metadata),
    )
    state = GapState((
        Particle("e", (0, 0, 0), (0, me * math.sqrt(2 * 5 * Q / me), 0)),
        Particle("M0", (0, 0, 0), (0, 0, 0)),
        Particle("Mi", (0, 0, 0), (0, 0, 0)),
    ))
    return catalog, state, channels


def explicit_reference_kernel(channel, catalog, incoming, rng):
    direction = rng.normal(size=3)
    direction /= np.linalg.norm(direction)
    return [1.0], [direction]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("runs/kinetic_gap_reference"))
    parser.add_argument("--steps", type=int, default=200)
    parser.add_argument("--dt-s", type=float, default=1e-12)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    if not 1 <= args.steps <= 10_000 or not np.isfinite(args.dt_s) or args.dt_s <= 0:
        parser.error("steps must be 1..10000 and dt-s must be finite and positive")
    catalog, state, channels = reference_case()
    rng = np.random.default_rng(args.seed)
    volume = 5e-24
    electric, magnetic = [0, 1000.0, 0], [0.5, 0, 0]
    rows, counts = [], {channel.name: 0 for channel in channels}
    field_gain, quadrature_work = 0.0, 0.0
    max_step_residual = 0.0
    initial_energy = total_energy_j(state, catalog)

    def snapshot():
        moments = cell_moments(state, catalog, volume)
        electron = next(p for p in state.particles if p.state_id == "e")
        ion = next(p for p in state.particles if p.state_id == "Mi")
        molecular_internal_ev = math.fsum(catalog.energy_ev(p.state_id) for p in state.particles)
        rows.append(dict(time_s=state.time_s, total_energy_ev=total_energy_j(state, catalog) / Q,
                         internal_energy_ev=molecular_internal_ev, field_gain_ev=field_gain / Q,
                         cumulative_balance_ev=(total_energy_j(state, catalog) - initial_energy - field_gain) / Q,
                         electron_vz_m_s=electron.momentum_kg_m_s[2] / catalog.mass_kg("e"),
                         ion_vz_m_s=ion.momentum_kg_m_s[2] / catalog.mass_kg("Mi"),
                         excited_molecule_count=sum(p.state_id == "M1" for p in state.particles),
                         current_z_a_m2=moments["current_density_a_m2"][2],
                         excitation_events=counts["excitation"], deexcitation_events=counts["deexcitation"]))

    snapshot()
    for _ in range(args.steps):
        result = advance_gap(state, catalog, args.dt_s, collision_volume_m3=volume,
                             electric_field_v_m=electric, magnetic_field_t=magnetic,
                             channels=channels, differential_kernels={c.name: explicit_reference_kernel for c in channels},
                             rng=rng)
        state = result.state
        field_gain = math.fsum((field_gain, result.field_kinetic_gain_j))
        quadrature_work = math.fsum((quadrature_work, result.electric_displacement_work_j))
        max_step_residual = max(max_step_residual, abs(result.total_energy_residual_j))
        for name, count in result.reaction_counts.items():
            counts[name] += count
        snapshot()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_dir / "history.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    root = Path(__file__).resolve().parents[1]
    source_paths = ["mhdlab/kinetic_gap.py", "mhdlab/kinetic_orbits.py", "mhdlab/state_kinetics.py",
                    "mhdlab/reaction_kinematics.py", "tools/run_gap_reference.py"]
    summary = {
        "scope": "manufactured prescribed-field homogeneous-cell transport and state-exchange reference",
        "is_measured_shot": False, "is_physical_water_mechanism": False,
        "steps": args.steps, "dt_s": args.dt_s, "seed": args.seed,
        "electric_field_v_m": electric, "magnetic_field_t": magnetic,
        "collision_volume_m3": volume, "catalog_fingerprint": catalog.fingerprint,
        "reaction_counts": counts, "initial_total_energy_ev": initial_energy / Q,
        "final_total_energy_ev": total_energy_j(state, catalog) / Q,
        "field_kinetic_gain_ev": field_gain / Q,
        "electric_displacement_work_ev": quadrature_work / Q,
        "maximum_step_balance_residual_ev": max_step_residual / Q,
        "final_cumulative_balance_residual_ev": rows[-1]["cumulative_balance_ev"],
        "source_sha256": {p: hashlib.sha256((root / p).read_bytes()).hexdigest() for p in source_paths},
        "not_implemented": ["self-consistent field/Hall feedback", "spatial collision mesh and walls",
                            "physical water state-to-state and differential reaction data",
                            "radiative population evolution", "FLASH/Quark runtime coupling"],
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
        time_ns = np.array([r["time_s"] for r in rows]) * 1e9
        axes[0].plot(time_ns, [r["total_energy_ev"] for r in rows], label="Total particle energy")
        axes[0].plot(time_ns, [initial_energy / Q + r["field_gain_ev"] for r in rows], "--", label="Initial + field work")
        axes[0].plot(time_ns, [r["internal_energy_ev"] for r in rows], label="Stored internal energy")
        axes[0].set_ylabel("Energy [eV]")
        axes[0].legend()
        axes[1].plot(time_ns, [r["electron_vz_m_s"] for r in rows], label="Electron")
        axes[1].plot(time_ns, [r["ion_vz_m_s"] for r in rows], label="Ion")
        axes[1].set_ylabel("Out-of-plane velocity [m/s]")
        axes[1].set_xlabel("Time [ns]")
        axes[1].legend()
        fig.suptitle("Manufactured reference: prescribed fields and state exchange")
        fig.tight_layout()
        fig.savefig(args.output_dir / "overview.png", dpi=150)
        plt.close(fig)
    print(json.dumps({"output_dir": str(args.output_dir), "reaction_counts": counts,
                      "maximum_step_balance_residual_ev": max_step_residual / Q}))


if __name__ == "__main__":
    main()
