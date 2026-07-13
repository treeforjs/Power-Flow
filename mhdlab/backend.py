"""Array backend selection for optional GPU acceleration."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class ArrayBackend:
    name: str
    xp: object

    @classmethod
    def from_preference(cls, preference: str = "numpy") -> "ArrayBackend":
        pref = (preference or "numpy").lower()
        if pref in {"cuda", "cupy", "gpu", "auto"}:
            try:
                import cupy as cp

                # Touch the current device and compile a tiny kernel so broken
                # header/toolkit installs fail before a long simulation starts.
                _ = cp.cuda.runtime.getDeviceCount()
                _ = int(cp.arange(2).sum().get())
                return cls(name="cupy", xp=cp)
            except Exception:
                if pref != "auto":
                    raise
        return cls(name="numpy", xp=np)

    def asnumpy(self, value):
        if self.name == "cupy":
            return self.xp.asnumpy(value)
        return np.asarray(value)

    @property
    def is_gpu(self) -> bool:
        return self.name == "cupy"
