"""CCD calibration, LOS integration, and synthetic spectra."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class CCDCalibration:
    wavelength_poly_nm: list[float]
    gate_s: float
    radiometric_response: float = 1.0
    background: float = 0.0
    instrument_fwhm_nm: float = 0.2

    def wavelength_axis(self, n_pixels: int) -> np.ndarray:
        pixels = np.arange(n_pixels, dtype=float)
        coeff = list(reversed(self.wavelength_poly_nm))
        out = np.zeros_like(pixels)
        for c in coeff:
            out = out * pixels + c
        return out


def load_ccd_image(path: str | Path) -> np.ndarray:
    import imageio.v3 as iio

    return np.asarray(iio.imread(path), dtype=float)


def synthesize_spectrum(
    wavelength_nm: np.ndarray,
    lines: list,
    line_strengths: dict[str, float],
    calibration: CCDCalibration,
    scale: float = 1.0,
    wavelength_offset_nm: float = 0.0,
) -> np.ndarray:
    sigma = calibration.instrument_fwhm_nm / max(2.0 * np.sqrt(2.0 * np.log(2.0)), 1e-30)
    spectrum = np.zeros_like(wavelength_nm, dtype=float)
    for line in lines:
        strength = line_strengths.get(line.name, 0.0)
        center = line.wavelength_nm + wavelength_offset_nm
        spectrum += strength * np.exp(-0.5 * ((wavelength_nm - center) / sigma) ** 2)
    return scale * calibration.radiometric_response * calibration.gate_s * spectrum + calibration.background


def integrate_los(emissivity: np.ndarray, raster, los: dict) -> float:
    """Integrate a 2D emissivity map along a straight chord using nearest cells."""
    start = np.asarray(los["start"], dtype=float)
    end = np.asarray(los["end"], dtype=float)
    samples = int(los.get("samples", 200))
    pts = start[None, :] + (end - start)[None, :] * np.linspace(0.0, 1.0, samples)[:, None]
    ix = np.clip(np.searchsorted(raster.x, pts[:, 0]), 0, raster.x.size - 1)
    iy = np.clip(np.searchsorted(raster.y, pts[:, 1]), 0, raster.y.size - 1)
    length = float(np.linalg.norm(end - start))
    return float(emissivity[iy, ix].sum() * length / max(samples, 1))


def fit_absolute_spectrum(measured: np.ndarray, synthetic: np.ndarray) -> dict[str, float]:
    a = np.column_stack([synthetic.ravel(), np.ones(measured.size)])
    coeff, *_ = np.linalg.lstsq(a, measured.ravel(), rcond=None)
    pred = (a @ coeff).reshape(measured.shape)
    rms = float(np.sqrt(np.mean((pred - measured) ** 2)))
    return {"scale": float(coeff[0]), "offset": float(coeff[1]), "rms": rms}
