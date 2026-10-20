#!/usr/bin/env python3
"""Export conservative source ledgers associated with preserved launch packets."""
from __future__ import annotations

import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))

from mhdlab.window_source_ledger import window_source_ledger
from tools.export_window_source_ledger import regular_file, same_json, write_ledger_csv
from tools.inspect_windowed_source_moments import associated_source, decode_packet_metadata
from tools.sample_temkin_source_windows import CODE_PATHS as PRODUCER_CODE_PATHS, sha

CODE_PATHS=(
    "tools/export_packet_source_ledger.py", "mhdlab/window_source_ledger.py",
    "tools/export_window_source_ledger.py", "tools/inspect_windowed_source_moments.py",
    "tools/compare_windowed_source_resolution.py", "mhdlab/source_window_plan.py",
)


def pin_implementation():
    # Imported producer and writer dependencies have their own literal lists;
    # include them in the association rather than hashing this wrapper alone.
    from tools.export_window_source_ledger import CODE_PATHS as EXPORT_CODE_PATHS
    from tools.inspect_windowed_source_moments import CODE_PATHS as ADMISSION_CODE_PATHS
    return {name:sha(ROOT/name) for name in sorted(set(CODE_PATHS)|set(PRODUCER_CODE_PATHS)
                                                |set(EXPORT_CODE_PATHS)|set(ADMISSION_CODE_PATHS))}


def recheck_inputs(directory, hashes):
    for name,digest in hashes.items():
        if Path(name).name != name:
            raise ValueError("packet input names must be plain local filenames")
        path=directory/name
        if sha(regular_file(path))!=digest:
            raise ValueError(f"packet input is indirect, missing or changed: {name}")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--packet-export",type=Path,required=True)
    parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args()
    if args.output_dir.exists() or args.output_dir.is_symlink():
        parser.error("--output-dir must be a new path")
    try:
        if args.packet_export.is_symlink() or (args.packet_export/"report.json").is_symlink():
            raise ValueError("packet export and report must be direct local inputs")
        regular_file(args.packet_export/"report.json")
        code_hashes=pin_implementation()
        source,packet_report,input_hashes=associated_source(args.packet_export)
        recheck_inputs(args.packet_export,input_hashes)
        manifest=source.manifest()
        if (not same_json(packet_report["conditioned_source"],manifest)
                or not same_json(packet_report["packets"]["source"],manifest)):
            raise ValueError("saved packet source JSON types differ from canonical replay")
        with np.load(args.packet_export/"boundary_packets.npz",allow_pickle=False) as packets:
            if not same_json(decode_packet_metadata(packets["metadata_json"]),packet_report["packets"]):
                raise ValueError("saved NPZ metadata JSON types differ from the packet report")
        ledger=window_source_ledger(source,include_components=True)
        report={
            "schema":"mhdlab-packet-source-ledger-export-v1",
            "parameter_status":source.status,"source_parameter_provenance":source.parameter_provenance,
            "conditioning_provenance":source.provenance,
            "original_mass_status":packet_report["mass_status"],
            "original_mass_provenance":packet_report["mass_provenance"],
            "input_export_sha256":input_hashes,"source_sha256":code_hashes,
            "host_budget_accepted":False,"source_evolution_repeated":False,
            "packet_sampling_repeated":False,"source_moments_recomputed":False,
            "runtime":{"python":sys.version,**{name:importlib.metadata.version(name) for name in ("numpy","scipy","h5py")}},
            "limits":[
                "The original packet export is read and associated with its canonical conditioned source. No packets, source histories or source moments are regenerated.",
                "Packet schema, dimensions and source association are checked; individual saved random realizations are not independently recomputed.",
                "Ledger values are deterministic expected source budgets. They are separate from each saved packet's sampled translation, statistical uncertainty or a host transaction.",
                "Original interval origins, local selections, state identifiers and pool costs remain distinct. Composite source selections are not photon recording gates.",
                "Each extensive channel retains independent source qualification and separate projection flags. Unqualified positive contributors never become structural zeros.",
                "Original state mass declarations remain provenance; no mass-dependent velocity or new molecular internal-state distribution is computed.",
                "Source labels do not validate physical water parameters, monolayers, absolute spectral calibration or photon response.",
                "Selected numerical code pins exclude eager legacy package imports. Staging and rehashes assume a coordinated single writer.",
            ],
        }
        args.output_dir.parent.mkdir(parents=True,exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=".packet-source-ledger-",dir=args.output_dir.parent) as temporary:
            output=Path(temporary)/"result"
            output.mkdir()
            original=output/"packet_export"
            original.mkdir()
            for name in input_hashes:
                shutil.copyfile(args.packet_export/name,original/name)
            (output/"ledger.json").write_text(json.dumps(ledger,indent=2,allow_nan=False)+"\n")
            write_ledger_csv(output/"ledger.csv",[{"case_id":"packet_source","window_index":None,
                                                "window":None,"ledger":ledger}])
            recheck_inputs(args.packet_export,input_hashes)
            recheck_inputs(original,input_hashes)
            if pin_implementation()!=code_hashes:
                raise ValueError("implementation changed during packet-source ledger export")
            report["output_sha256"]={str(path.relative_to(output)):sha(path) for path in output.rglob("*") if path.is_file()}
            (output/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
            if args.output_dir.exists() or args.output_dir.is_symlink():
                raise ValueError("output path appeared during export; refusing replacement")
            os.rename(output,args.output_dir)
    except (OSError,ValueError,TypeError,KeyError,OverflowError,RuntimeError) as error:
        parser.error(str(error))
    print(json.dumps({"output_dir":str(args.output_dir),"parameter_status":source.status,
                      "host_budget_accepted":False,"source_evolution_repeated":False,
                      "packet_sampling_repeated":False}))


if __name__=="__main__":
    main()
