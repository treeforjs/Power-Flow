#!/usr/bin/env python3
"""Fit explicit monolayer coverage against supplied binned predictions."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from mhdlab.coverage_regression import (
    BinnedObservation, CalibrationScale, CoveragePredictionGrid, fit_surface_coverage,
)

CONFIG_SCHEMA = 'mhdlab-surface-coverage-fit-config-v1'
SCHEMA = 'mhdlab-surface-coverage-fit-run-v1'
CODE_PATHS = ('mhdlab/coverage_regression.py', 'tools/fit_surface_coverage.py')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f'duplicate configuration key {key!r}')
        result[key] = value
    return result


def load_config(path):
    def invalid_constant(value):
        raise ValueError(f'nonfinite JSON constant {value!r}')
    config = json.loads(Path(path).read_bytes(), object_pairs_hook=_object, parse_constant=invalid_constant)
    required = {'schema', 'grid', 'observation', 'calibration'}
    if not isinstance(config, dict) or not required <= set(config) or set(config) - required - {'manufactured_truth'}:
        raise ValueError('configuration must contain schema, grid, observation and calibration, with optional manufactured_truth')
    if config['schema'] != CONFIG_SCHEMA:
        raise ValueError(f'configuration schema must be {CONFIG_SCHEMA}')
    if not all(isinstance(config[name], dict) for name in ('grid', 'observation', 'calibration')):
        raise ValueError('grid, observation and calibration must be JSON objects')
    grid = CoveragePredictionGrid(**config['grid'])
    observation = BinnedObservation(**config['observation'])
    calibration = CalibrationScale(**config['calibration'])
    if 'manufactured_truth' in config and (grid.status != 'manufactured' or observation.status != 'manufactured'):
        raise ValueError('manufactured_truth is permitted only for explicitly manufactured grid and observations')
    return config, grid, observation, calibration


def write_outputs(output, config, grid, observation, result, config_hash, code_hashes, *, plot):
    (output / 'input_config.json').write_text(json.dumps(config, indent=2, allow_nan=False) + '\n')
    profile = result['grid_profile']
    with (output / 'coverage_profile.csv').open('w', newline='') as stream:
        columns = tuple(profile)
        writer = csv.writer(stream)
        writer.writerow(columns)
        writer.writerows(zip(*(profile[name] for name in columns)))
    best = result['best_fit']
    with (output / 'best_fit_bins.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['bin_index', 'bin_left', 'bin_right', 'observed_signal', 'known_background',
                         'scaled_prediction_signal', 'prediction_with_background', 'observed_minus_prediction'])
        for index in range(len(observation.values)):
            writer.writerow([index, grid.bin_edges[index], grid.bin_edges[index + 1], observation.values[index],
                             observation.known_background[index], best['prediction_signal'][index],
                             best['prediction_with_background'][index],
                             observation.values[index] - best['prediction_with_background'][index]])
    if plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        figure, axes = plt.subplots(1, 2, figsize=(11, 4.3), constrained_layout=True)
        centers = .5 * (grid.bin_edges[:-1] + grid.bin_edges[1:])
        axes[0].errorbar(centers, observation.values, yerr=np.sqrt(np.diag(observation.covariance)),
                         fmt='o', label='Supplied observations and marginal errors')
        axes[0].stairs(best['prediction_with_background'], grid.bin_edges, label='Fitted prediction plus known background')
        axes[0].stairs(observation.known_background, grid.bin_edges, label='Known background', color='gray')
        axes[0].set(xlabel=f'Bin coordinate [{grid.bin_axis_units}]', ylabel=grid.signal_units,
                    title=f'Coverage {best["coverage_ml"]:.6g} ML; scale {best["calibration_scale"]:.6g}')
        axes[0].legend(fontsize=7)
        axes[1].scatter(profile['coverage_ml'], profile['delta_chi2'], label='Supplied coverage knots')
        candidates = result['search_candidates']
        axes[1].scatter([item['coverage_ml'] for item in candidates],
                         [item['delta_chi2'] for item in candidates],
                         marker='x', label='Search candidates')
        axes[1].plot([best['coverage_ml']], [0.], '*', markersize=12, label='Selected minimum')
        axes[1].set(xlabel='Coverage under supplied monolayer definition [ML]', ylabel='Profile objective minus minimum',
                    title='Finite grid support; no automatic confidence interval')
        axes[1].legend(fontsize=7)
        figure.suptitle(f'Coverage regression: {grid.status} predictions / {observation.status} observations')
        figure.savefig(output / 'coverage_fit.png', dpi=160)
        plt.close(figure)
    report = {
        'schema': SCHEMA, 'input_config_sha256': config_hash, 'source_sha256': code_hashes,
        'fit': result,
        'input_labels': {'prediction_status': grid.status, 'observation_status': observation.status,
                         'monolayer_provenance': grid.monolayer_provenance, 'calibration_provenance': config['calibration']['provenance']},
        'output_sha256': {path.name: sha(path) for path in sorted(output.iterdir()) if path.is_file()},
        'limits': ['Predictions and diagnostic response are supplied inputs; no chemistry, radiation or gap evolution is calculated.',
                   'The monolayer definition and calibration constraint are explicit; the example is not a physical water-coverage estimate.',
                   'Finite interpolation, supplied Gaussian covariance and scale-penalty assumptions require separate physical qualification.',
                   'An observation fit does not repair unknown calibration units or recover missing original detector data.'],
    }
    if 'manufactured_truth' in config:
        report['manufactured_truth_not_used_by_fit'] = config['manufactured_truth']
    (output / 'fit_report.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True, type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--plot', action='store_true')
    args = parser.parse_args()
    if args.output_dir.exists() or args.output_dir.is_symlink():
        parser.error('--output-dir must be a new path')
    try:
        config_hash = sha(args.config)
        code_hashes = {name: sha(ROOT / name) for name in CODE_PATHS}
        config, grid, observation, calibration = load_config(args.config)
        result = fit_surface_coverage(grid, observation, calibration)
        if sha(args.config) != config_hash or any(sha(ROOT / name) != digest for name, digest in code_hashes.items()):
            raise ValueError('configuration or implementation changed during fitting')
        args.output_dir.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix='.coverage-fit-', dir=args.output_dir.parent) as directory:
            output = Path(directory) / 'result'
            output.mkdir()
            report = write_outputs(output, config, grid, observation, result, config_hash, code_hashes, plot=args.plot)
            if args.output_dir.exists() or args.output_dir.is_symlink():
                raise ValueError('output path appeared during fitting')
            os.rename(output, args.output_dir)
        print(json.dumps({'output_dir': str(args.output_dir), 'best_fit': report['fit']['best_fit'],
                          'diagnostics': report['fit']['diagnostics']}, indent=2, allow_nan=False))
    except (OSError, ValueError, TypeError, KeyError, OverflowError) as error:
        parser.exit(2, f'coverage fit: {error}\n')


if __name__ == '__main__':
    main()
