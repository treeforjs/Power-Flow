"""Reduced MHD / kinetic-neutral / CR research prototype."""

from .config import load_config
from .runner import run_from_config

__all__ = ["load_config", "run_from_config"]
