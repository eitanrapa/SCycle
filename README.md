# README #

Scyle is a parallelized finite difference code that simulates earthquake cycles on a strike-slip fault in 2D. The fault is governed by rate-and-state friction using either the aging or slip law. The coseismic phase can be fully inertial (meaning that wave-propagation is resolved) or quasi-dynamic, while the interseismic phase is assumed to be quasi-dynamic. A number of different rheological choices are supported: elasticity, linear Maxwell viscoelasticity, and power-law viscoelasticity with dislocation creep, diffusion creep, or both operating in parallel. Two strain-weakening mechanisms are supported: grain size reduction, and shear heating.

Several representative examples, including visualization files, can be found in the examples/ folder.

Please see scycle-manual.pdf for more information, including example files. Parts of the manual are out of date; `docs/AUDIT.md` lists the differences.

## Quickstart

Requirements: PETSc 3.14 or newer configured with HDF5 and MUMPS (and HYPRE for the AMG solver options), real scalars and 32-bit indices, plus an MPI implementation.

```bash
export PETSC_DIR=/path/to/petsc PETSC_ARCH=your-arch   # PETSC_ARCH empty for a --prefix install
make -C source -j8                                      # builds source/main
./source/main examples/ex1.in                           # spring slider, about 15 s
mpirun -n 4 ./source/main examples/ex2.in               # 2D elastic, quasi-dynamic
```

Output goes to the prefix given by `outputDir` in the input file (`data/ex2_data_1D.h5`, ...). If `<outputDir>checkpoint.h5` exists, a run restarts from it unless the input sets `restartFromChkpt = 0`. Stop a run with Ctrl-C or `kill` (SIGTERM): it writes a checkpoint and closes its files first, and running it again continues it. Do not use `kill -9`, which can leave the HDF5 output unreadable. Load results with `examples/loadFuncs.py` (Python, h5py) or `matlab/visualizePetsc` (MATLAB); `examples/visualize_ex*.ipynb` and `visualize_ex*.m` show how.

Unknown keys in an input file are ignored without a warning; `tools/checkkeys.py file.in` lists them.

## Documentation in this repository

- `CLAUDE.md`: build, run, input rules, code map and conventions.
- `docs/AUDIT.md`: defects found in an October 2026 audit, how each was fixed and verified, and the known limitations that remain.
- `docs/TWO_FAULT_DESIGN.md`: design for extending the code to two interacting faults.
- `docs/SERVER.md`: building on a Linux server and running batches of long runs there with `tools/batch.sh`.
- `tools/regress.sh`, `tools/mms.in`: regression and convergence checks; `tools/compare_runs.py` compares the output of different builds by its physics.
- `SEAS_benchmarks/BP1`: SEAS benchmark problem 1 (`createICs.py` writes its grid and initial conditions).
