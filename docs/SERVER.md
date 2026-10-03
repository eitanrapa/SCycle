# Running SCycle on a Linux server

How to build SCycle on a Linux machine and run long batches there. Use `master`, the fork's only
branch: it holds all the two-fault work (stages 1 to 5) and the fixes that let GCC build the code
(`933d1bd`, `89da3ae`); the code at the tag `stage4` does not compile with GCC.
The steps assume a machine you log into with ssh; for a cluster with a scheduler, see section 6.

## 1. Packages

A C, C++ and Fortran compiler, MPI, BLAS/LAPACK, git and Python 3. On Debian or Ubuntu:

```bash
sudo apt install build-essential gfortran openmpi-bin libopenmpi-dev libopenblas-dev \
                 git python3-venv hdf5-tools
```

`hdf5-tools` provides `h5diff`, which `tools/regress.sh` needs; the HDF5 that PETSc builds does
not include it. On a cluster, load a compiler and an MPI module instead (`module avail`).

## 2. PETSc

SCycle needs PETSc 3.14 or newer with HDF5, MUMPS and HYPRE, real scalars and 32-bit indices. The
Mac uses 3.26.0. Building it takes about an hour, most of it in the packages PETSc downloads:

```bash
git clone -b v3.26.0 --depth 1 https://gitlab.com/petsc/petsc.git ~/opt/petsc-3.26.0
cd ~/opt/petsc-3.26.0
./configure PETSC_ARCH=arch-linux-c-opt \
  --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 \
  --with-debugging=0 COPTFLAGS=-O3 CXXOPTFLAGS=-O3 FOPTFLAGS=-O3 \
  --with-scalar-type=real --with-64-bit-indices=0 --with-x=0 \
  --download-mumps --download-scalapack --download-metis --download-parmetis \
  --download-hypre --download-hdf5 --download-cmake
make PETSC_DIR=$PWD PETSC_ARCH=arch-linux-c-opt all
make PETSC_DIR=$PWD PETSC_ARCH=arch-linux-c-opt check
```

Configure uses the system BLAS/LAPACK; without one, add `--download-fblaslapack`. For working on
the code, build a debug PETSc next to it: the same configure line with
`PETSC_ARCH=arch-linux-c-debug --with-debugging=1` and without the three `OPTFLAGS`.

A PETSc that is already installed (a cluster module) will do if

```bash
grep -E "^PETSC_(SCALAR|PRECISION|INDEX_SIZE) " $PETSC_DIR/$PETSC_ARCH/lib/petsc/conf/petscvariables
grep -E "define PETSC_HAVE_(MUMPS|HYPRE|HDF5) " $PETSC_DIR/$PETSC_ARCH/include/petscconf.h
```

print `real`, `double` and `32`, and all three `HAVE` lines (`PETSC_ARCH` is empty for an install
made with `--prefix`).

## 3. SCycle

```bash
git clone https://github.com/eitanrapa/SCycle.git ~/SCycle
cd ~/SCycle
export PETSC_DIR=$HOME/opt/petsc-3.26.0 PETSC_ARCH=arch-linux-c-opt   # also put this in ~/.bashrc
make -C source -j8 WERROR=1
```

GCC 16 builds the tree without warnings. Python, for the input generators and the analysis:

```bash
python3 -m venv ~/venvs/scycle
~/venvs/scycle/bin/pip install numpy h5py
source ~/venvs/scycle/bin/activate
```

## 4. Check the build

```bash
./source/main tools/mms.in                        # convergence rates: order 4 gives about 3.5 for u, 2.5 for sxy
tools/regress.sh baseline data/regress-baseline   # about a minute with the optimized build
```

The second command stores the baseline that later changes are compared with
(`tools/regress.sh compare data/regress-baseline`). Make it on the server, with the `PETSC_ARCH`
you will use: another compiler or CPU changes the last digits, the adaptive steps then diverge, and
a bit-for-bit comparison with the Mac's baseline fails.

To check that the server reproduces the Mac, compare the physics instead. Copy the Mac's baseline
(85 MB) over, from the Mac:

```bash
scp -r data/regress-baseline <server>:SCycle/data/regress-baseline-mac
```

then, on the server:

```bash
python tools/compare_runs.py data/regress-baseline-mac data/regress-baseline
```

It lists the earthquakes and the end state of each case. On the Mac, the optimized clang and GCC
builds match the debug baseline with the same number of events, onsets within 0.021 yr over
4753 yr (ex1), 1e-6 yr (ex2) and 2.2e-4 yr (ex4s, ex4g); a server build should agree as closely.

## 5. Batches

Generate a batch on the machine that runs it: the inputs hold absolute paths.

```bash
python examples/two_faults/stage5_batch1.py data/stage5_batch1   # or the driver of a newer batch
tools/batch.sh start  data/stage5_batch1 -j 6                    # 6 runs at a time, 1 MPI rank each
tools/batch.sh status data/stage5_batch1
tools/batch.sh stop   data/stage5_batch1                         # each run checkpoints and stops
python tools/two_fault.py data/stage5_batch1/alt/ --zref 10 --window 1000
```

The runs are detached, so they outlive the ssh session. After a stop or a reboot, `start` resumes
each run from its last checkpoint; `tools/batch.sh` with no arguments prints its usage, and the top
of the script explains the rest. SCycle stops cleanly on SIGTERM or SIGINT: it finishes the step,
writes a checkpoint, closes its files and exits with status 143 (130 for SIGINT). Never `kill -9` a
run: HDF5 output killed while it is being written can become unreadable, and with it the whole
history of the run.

Sizing, from the first batch (581 x 172 nodes, 12,000 yr, about 110,000 steps per run): each run
took 16 to 17 hours on one core of an Apple M4 with six running at once, needed about 1.1 GB of
memory, and wrote about 0.9 GB. A run matrix finishes soonest as single-rank runs side by side, as
many as the cores and the memory allow. `-n R` gives each run R MPI ranks.

To fetch results, everything but the HDF5 files (the CSV files of `tools/two_fault.py`,
`faultSeries.txt`, the logs) is small. From the Mac:

```bash
rsync -av --exclude '*.h5' <server>:SCycle/data/stage5_batch1/ data/stage5_batch1_server/
```

## 6. Clusters

On a cluster with a scheduler, do not run on the login node. With Slurm, `--wait` keeps the queue
inside the job:

```bash
sbatch -N 1 -n 1 -c 6 --mem=8G -t 2-00:00:00 \
       --wrap "tools/batch.sh start data/stage5_batch1 -j 6 --wait"
```

At the time limit Slurm sends SIGTERM: every run checkpoints and stops, and submitting the same
job again resumes them. That holds for single-rank runs. With `-n 2` or more, `mpirun` receives the
signal too and kills its ranks at once, with the risk to the output described in section 5, so
keep scheduled runs at one rank or stop them with `tools/batch.sh stop` before the limit.

## 7. Troubleshooting

- Runs stop when you log out: the server ends a session's processes at logout (systemd's
  `KillUserProcesses=yes`). Ask the administrator, or start the batch outside the session with
  `loginctl enable-linger` followed by `systemd-run --user --scope tools/batch.sh start ...`.
- `h5diff: command not found`: install `hdf5-tools` (section 1).
- Several multi-rank runs at once are slow: Open MPI binds each run's ranks to the same first
  cores. Start the batch with `MPIEXEC_FLAGS="--bind-to none"`.
- A run is `failed` in `tools/batch.sh status`: its `run.log` ends with the error. Once it is
  fixed, delete `run.exit` in the run's output directory and `start` again.
