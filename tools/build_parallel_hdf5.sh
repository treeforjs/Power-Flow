#!/usr/bin/env bash
# Build a private prerequisite; does not build or qualify FLASH.
set -euo pipefail
if [[ $# -ne 2 ]]; then
  echo "usage: bash tools/build_parallel_hdf5.sh ABS_INSTALL_PREFIX ABS_NEW_WORK_DIR" >&2
  exit 2
fi
prefix=$1
work=$2
[[ $prefix = /* && $work = /* && $prefix != / && $work != / ]] || exit 2
[[ ! -e $prefix && ! -e $work ]] || {
  echo "Use new prefix and work directories; existing installations are preserved." >&2
  exit 2
}
mkdir -p "$work"
cd "$work"
release=https://github.com/HDFGroup/hdf5/releases/download/hdf5_1.14.6
curl --fail --location --retry 2 "$release/hdf5-1.14.6.tar.gz" -o hdf5-1.14.6.tar.gz
curl --fail --location --retry 2 "$release/hdf5-1.14.6.sha256sums.txt" -o upstream.sha256sums.txt
python3 - <<'PY'
from pathlib import Path
import hashlib
import tarfile
name = 'hdf5-1.14.6.tar.gz'
matches = [line.split()[0] for line in Path('upstream.sha256sums.txt').read_text().splitlines()
           if len(line.split()) == 2 and line.split()[1].lstrip('*').removeprefix('./') == name]
actual = hashlib.sha256(Path(name).read_bytes()).hexdigest()
if len(matches) != 1 or matches[0].lower() != actual:
    raise SystemExit('upstream SHA256 missing, ambiguous or mismatched')
Path('verified_source.sha256').write_text(actual + '  ' + name + '\n')
with tarfile.open(name) as archive:
    archive.extractall('source', filter='data')
PY
src=$(find "$work/source" -mindepth 1 -maxdepth 1 -type d)
[[ -f $src/CMakeLists.txt ]] || exit 3
cmake -S "$src" -B "$work/build" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$prefix" \
  -DCMAKE_C_COMPILER=mpicc -DCMAKE_Fortran_COMPILER=mpif90 \
  -DHDF5_ENABLE_PARALLEL=ON -DHDF5_BUILD_FORTRAN=ON \
  -DHDF5_BUILD_CPP_LIB=OFF -DHDF5_BUILD_JAVA=OFF \
  -DHDF5_BUILD_HL_LIB=ON -DBUILD_TESTING=OFF \
  -DHDF5_BUILD_TOOLS=OFF -DHDF5_BUILD_EXAMPLES=OFF \
  -DHDF5_ENABLE_Z_LIB_SUPPORT=OFF -DHDF5_ENABLE_SZIP_SUPPORT=OFF \
  -DDEFAULT_API_VERSION=v16 2>&1 | tee configure.log
cmake --build "$work/build" --parallel 4 2>&1 | tee build.log
cmake --install "$work/build" 2>&1 | tee install.log
"$prefix/bin/h5pfc" -showconfig > installed_hdf5_config.txt
echo "Private HDF5 installed: $prefix"
echo "Run tools/probe_flash_environment.py --hdf5-prefix $prefix --compile-smoke next."
