#!/usr/bin/env bash
# Build Open MPI 5.0.11 and PETSc 3.26.0 (optimized and debug arches) for SCycle on a Mac without
# Homebrew (the Intel Mac mini from 2026-10-06, the iMac Pro from 2026-10-28). C and C++ come from
# the Command Line Tools (Apple clang), Fortran from conda-forge: gfortran in its own environment
# ~/miniforge3/envs/scycle-fortran, whose bin must be on PATH (its assembler runs clang) and whose
# lib ends up on SCycle's rpath (the Fortran runtime); CMake from ~/miniforge3/envs/scycle-build.
# PETSc's configure options are those of docs/SERVER.md. Logs in ~/opt: build-openmpi.log,
# build-petsc-<arch>.log, check-petsc-<arch>.log. On a new Mac, first the Command Line Tools
# (xcode-select --install, the one step that needs an administrator), then Miniforge and the two
# environments (to match another machine exactly, create each from its `conda list --explicit`):
#   curl -fsSLO https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh
#   bash Miniforge3-MacOSX-x86_64.sh -b -p ~/miniforge3
#   ~/miniforge3/bin/mamba create -y -p ~/miniforge3/envs/scycle-fortran -c conda-forge gfortran
#   ~/miniforge3/bin/mamba create -y -p ~/miniforge3/envs/scycle-build -c conda-forge cmake hdf5
# (Apple Silicon: the MacOSX-arm64 installer.) Then, from anywhere, detached:
#   nohup bash tools/build-mpi-petsc.sh > ~/opt/build.log 2>&1 &
set -euo pipefail
OPT=$HOME/opt
F=$HOME/miniforge3/envs/scycle-fortran
B=$HOME/miniforge3/envs/scycle-build
MPIV=5.0.11
PETSCV=3.26.0
NJ=12
export PATH=$PATH:$F/bin
mkdir -p "$OPT/downloads"
cd "$OPT/downloads"
[ -f openmpi-$MPIV.tar.bz2 ] || curl -fsSLO https://download.open-mpi.org/release/open-mpi/v5.0/openmpi-$MPIV.tar.bz2
[ -f petsc-$PETSCV.tar.gz ] || curl -fsSLO https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-$PETSCV.tar.gz

echo "$(date '+%F %T') Open MPI $MPIV"
if [ ! -x "$OPT/openmpi-$MPIV/bin/mpicc" ]; then
  rm -rf "openmpi-$MPIV"
  tar xjf openmpi-$MPIV.tar.bz2
  cd "openmpi-$MPIV"
  ./configure --prefix="$OPT/openmpi-$MPIV" CC=/usr/bin/clang CXX=/usr/bin/clang++ FC="$F/bin/gfortran" \
    --with-libevent=internal --with-hwloc=internal --with-pmix=internal --with-prrte=internal \
    > "$OPT/build-openmpi.log" 2>&1
  make -j$NJ >> "$OPT/build-openmpi.log" 2>&1
  make install >> "$OPT/build-openmpi.log" 2>&1
  cd "$OPT/downloads"
fi
export PATH=$OPT/openmpi-$MPIV/bin:$PATH

# PETSc copies gfortran's library directories into its link lines, and the environment's lib also
# holds conda's libc++, which C++ links would then take instead of the system's: hide its link-time
# names (the environment's clang, which gfortran assembles with, loads libc++.1.dylib)
for f in libc++.dylib libc++.a; do
  if [ -e "$F/lib/$f" ]; then mv "$F/lib/$f" "$F/lib/$f.hidden"; fi
done

[ -d "$OPT/petsc-$PETSCV" ] || tar xzf petsc-$PETSCV.tar.gz -C "$OPT"
cd "$OPT/petsc-$PETSCV"
COMMON=(--with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 --with-mpiexec="$OPT/openmpi-$MPIV/bin/mpiexec"
        --with-scalar-type=real --with-64-bit-indices=0 --with-x=0
        --download-mumps --download-scalapack --download-metis --download-parmetis --download-hypre --download-hdf5
        --with-cmake-exec="$B/bin/cmake")
for arch in opt debug; do
  echo "$(date '+%F %T') PETSc $PETSCV, arch-darwin-c-$arch"
  if [ $arch = opt ]; then extra=(--with-debugging=0 COPTFLAGS=-O3 CXXOPTFLAGS=-O3 FOPTFLAGS=-O3); else extra=(--with-debugging=1); fi
  ./configure PETSC_ARCH=arch-darwin-c-$arch "${COMMON[@]}" "${extra[@]}" > "$OPT/build-petsc-$arch.log" 2>&1
  make PETSC_DIR="$OPT/petsc-$PETSCV" PETSC_ARCH=arch-darwin-c-$arch all >> "$OPT/build-petsc-$arch.log" 2>&1
  make PETSC_DIR="$OPT/petsc-$PETSCV" PETSC_ARCH=arch-darwin-c-$arch check > "$OPT/check-petsc-$arch.log" 2>&1
  if otool -L "$OPT/petsc-$PETSCV/arch-darwin-c-$arch/lib/libpetsc.dylib" | grep '@rpath/libc++' > /dev/null; then
    echo "libpetsc links conda's libc++" >&2; exit 1
  fi
done
echo "$(date '+%F %T') done"
