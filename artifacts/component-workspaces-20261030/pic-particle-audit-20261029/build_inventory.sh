#!/usr/bin/env bash
set -eu
exec taskset -c 0 nice -n 15 /usr/sbin/c++ -O3 -DNDEBUG \
  -I/home/treefor/projects/runtimes/mfem-gpu-pic-20261029/source/include \
  -isystem /home/treefor/projects/runtimes/mfem-4.10-pic-20261029/install-mpi/include \
  -isystem /home/treefor/projects/runtimes/mfem-4.10-pic-20261029/tpl/hypre/include \
  -isystem /home/treefor/projects/runtimes/mfem-4.10-pic-20261029/tpl/metis/include \
  /home/treefor/projects/runtimes/pic-particle-audit-20261029/checkpoint_inventory.cpp \
  -o /home/treefor/projects/runtimes/pic-particle-audit-20261029/checkpoint_inventory \
  /home/treefor/projects/runtimes/mfem-4.10-pic-20261029/install-mpi/lib/libmfem.a \
  /home/treefor/projects/runtimes/mfem-4.10-pic-20261029/tpl/hypre/lib/libHYPRE.a \
  /home/treefor/projects/runtimes/mfem-4.10-pic-20261029/tpl/metis/lib/libmetis.a \
  /home/treefor/projects/runtimes/mfem-4.10-pic-20261029/tpl/GKlib/lib/libGKlib.a \
  /usr/lib/libmpi.so -Wl,-rpath,/usr/lib
