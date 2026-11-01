"""Audit qualification entry point: cpu|gpu LABEL, in a matching SimQueue allocation.
CPU and GPU receipts are separate; neither alone claims complete qualification.
"""
from pathlib import Path
import runpy
runpy.run_path(str(Path(__file__).with_name("audit-tests.py")), run_name="__main__")
