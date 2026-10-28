"""Require canonical research interchange arrays to live in their hashed file."""

from __future__ import annotations

import h5py


def local_dataset(handle: h5py.Group, name: str) -> h5py.Dataset:
    """Resolve only hard-linked paths with self-contained dataset storage.

    External/soft links and virtual/external storage make a hash of the parent
    file insufficient to identify the actual arrays. Canonical interchange
    readers reject those dependencies instead of silently following them.
    """
    if not isinstance(name, str) or not name or any(part in {"", ".", ".."} for part in name.split("/")):
        raise ValueError("dataset paths must be explicit nonempty relative HDF5 paths")
    current = handle
    parts = name.split("/")
    for index, part in enumerate(parts):
        if not isinstance(current, h5py.Group):
            raise ValueError(f"non-group component in dataset path: {name}")
        if not isinstance(current.get(part, getlink=True), h5py.HardLink):
            raise ValueError(f"dataset path must use local hard links: {name}")
        current = current[part]
        if index < len(parts) - 1 and not isinstance(current, h5py.Group):
            raise ValueError(f"non-group component in dataset path: {name}")
    if not isinstance(current, h5py.Dataset):
        raise ValueError(f"required dataset is missing: {name}")
    if current.is_virtual or current.external:
        raise ValueError(f"dataset storage must be self-contained and nonvirtual: {name}")
    return current
