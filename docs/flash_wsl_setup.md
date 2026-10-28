# FLASH workstation preparation

Historical setup record. MFEM Joule now replaces FLASH in the active workflow;
see the [current backend](mfem_joule_heating.md). The private HDF5 installation
and passing MPI evidence remain available. FLASH acquisition is no longer an
active requirement.

The selected host is Ubuntu 24.04 WSL2 on this PC. The prerequisite check found
GNU C/Fortran 13.3 and OpenMPI 4.1.6. No licensed FLASH source root has been
identified. The source has not been downloaded, as confirmed during setup.
No FLASH build, heating benchmark or physical foil calculation has been run here.
The official [download page](https://flash.rochester.edu/site/flashcode/download.html)
uses the credentials supplied through its code-request process. Source acquisition
remains a separate user step; no credentials are stored in this repository.

Parallel HDF5 1.14.6 was built privately under
`/home/tsmit/.local/opt/power-flow/hdf5-1.14.6-openmpi`. Its official release
archive was checked against the upstream SHA256 list before extraction. The
build uses the installed OpenMPI C/Fortran wrappers and the v16 default API.
Compression filters were disabled. System packages and shell startup files were
not changed. This prerequisite build has not established compatibility with an
uninspected FLASH version.

From the Windows repository directory, the reproducible build command is:

```powershell
wsl.exe -d Ubuntu-24.04 --exec bash tools/build_parallel_hdf5.sh /home/tsmit/.local/opt/power-flow/hdf5-1.14.6-openmpi /home/tsmit/.local/share/power-flow/hdf5-1.14.6-openmpi-build-20261027
```

The script requires fresh work and install paths. Preserve the completed build;
do not repeat that command against the installed prefix. Configure/build/install
logs, verified archive digest, source and CMake cache remain in the external
work directory. Public sources: [HDF5 release](https://github.com/HDFGroup/hdf5/releases/tag/hdf5_1.14.6)
and [HDF5 CMake instructions](https://github.com/HDFGroup/hdf5/blob/hdf5_1.14.6/release_docs/INSTALL_CMake.txt).

The following check compiled Fortran with the MPI/HDF5 API and collectively
created an HDF5 file with two MPI ranks:

```powershell
python tools/probe_flash_environment.py --wsl-distro Ubuntu-24.04 --hdf5-prefix /home/tsmit/.local/opt/power-flow/hdf5-1.14.6-openmpi --compile-smoke --output runs/flash_preflight_hdf5_runtime/report.json
```

The private shared libraries require their runtime path. The probe scopes
`LD_LIBRARY_PATH` to its child process. An initial successful compile followed by
a missing-library execution failure is preserved in
`runs/flash_preflight_parallel_hdf5/report.json`; the runtime-path correction and
passing two-rank test are in `runs/flash_preflight_hdf5_runtime/report.json`.
The script records exact arguments, exits, tool versions and its own source hash.

Once the licensed source path is supplied, use a new report path and add
`--source-root /absolute/linux/source/path`. Inspection hashes the setup script
and checks the expected source layout without copying licensed files or invoking
the setup script. Record the actual version independently. Keep source, build,
site configuration and binaries outside the public repository.

Next, configure a version-specific site file with these MPI wrappers, HDF5
include/library directories and runtime library path. Check the installed
version's I/O API and resistivity implementation. The public
[FLASH MHD guide](https://flash.rochester.edu/site/flashcode/user_support/flash_ug_devel/node107.html)
documents magnetic diffusivity rather than electrical resistivity as its
transport coefficient and a general uniform-grid restriction for the implicit
solver. These are setup leads, not verification of an uninspected installation.

Run the [constant-property benchmark](flash_heating_verification.md) before
SS304, finite-width edges or desorption. Require magnetic/Joule energy accounting,
current closure, correct lattice temperature, unit verification, independent
mesh/time refinement and a nonsingular physical drive. The shot 17133 heating
benchmark remains distinct from the shot 15750/15765 spectroscopy cases.
