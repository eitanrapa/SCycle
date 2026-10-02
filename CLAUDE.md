# CLAUDE.md

Guidance for working in this repository. Read `docs/AUDIT.md` for what was found and fixed, and
`docs/TWO_FAULT_DESIGN.md` for the planned extension.

## What this is

SCycle simulates earthquake sequences and aseismic slip on a vertical strike-slip fault in 2D
antiplane shear (C++, PETSc), by Kali Allison (Stanford). It solves the quasi-static or dynamic
momentum balance on a rectangular grid with summation-by-parts (SBP-SAT) finite differences and
couples it to rate-and-state friction on the fault. Optional physics: power-law viscoelasticity
(dislocation, diffusion, dissolution-precipitation creep, plasticity), grain-size evolution, heat
equation with shear heating, pore pressure with slip-dependent permeability (fault valving).

**Research goal of this fork.** Extend the single boundary fault to two parallel faults that share
the far-field load (SAF–SJF, CF–AF), the second as an interior frictional interface, and ask whether
and when they trade slip. Stages: (i) elastic, verifying stress transfer, slip partitioning and the
single-fault limit; (ii) power-law viscoelastic; (iii) grain-size evolution, fault valving and shear
heating. Output: a library of two-fault models with partitioning histories and surface-velocity
profiles. Design and verification gates: `docs/TWO_FAULT_DESIGN.md`.

Upstream: `origin` = github.com/kali-allison/SCycle (`master`). Audit fixes live on branch
`audit/fixes-2026-10`. Nothing has been pushed.

## Build

Requires PETSc >= 3.14 (HDF5 timestepping API) with HDF5 and MUMPS, real scalars, 32-bit indices;
HYPRE for the `AMG` and `CG_PCAMG` solver options. On this machine PETSc 3.26.0 is built from source
(Homebrew's `petsc` was unusable: its `hdf5-mpi` conflicts with `hdf5`, and it lacks MUMPS/HYPRE):

```bash
export PETSC_DIR=$HOME/opt/petsc-3.26.0
export PETSC_ARCH=arch-darwin-c-debug     # or arch-darwin-c-opt (-O3) for production runs
make -C source -j8                        # -> source/main
make -C source -j8 WERROR=1               # warnings as errors (the tree is warning-free)
make -C source DEBUG_MODULES=-DVERBOSE=2  # more tracing; VERBOSE=1 (default) prints one line per step
make -C source clean
```

MPI is Homebrew Open MPI (`mpirun`). `h5diff`/`h5dump`/`h5ls` come from Homebrew `hdf5`.
A scratch Python with h5py for the loaders: `python3 -m venv <dir> && <dir>/bin/pip install h5py numpy`.
MATLAB is installed but its license has expired (the MATLAB loaders are untested).

## Run

```bash
./source/main examples/ex2.in            # from the repository root
mpirun -n 4 ./source/main examples/ex2.in
```

- `outputDir` is a file-name **prefix** (`data/ex2_` writes `data/ex2_data_1D.h5`, ...); its
  directory must exist. `data/` is shipped and ignored by git.
- **Restart is on by default** (`restartFromChkpt = 1`): if `<outputDir>checkpoint.h5` exists, the
  run resumes from it and appends to the existing output, printing `RESTARTING`. Use a fresh prefix
  or `restartFromChkpt = 0` for a clean run. After a restart, `maxStepCount` counts steps taken in
  the new run, not the total.
- Checkpoints (every `strideChkpt` steps) are written to `checkpoint.h5.tmp` and renamed when
  complete; the 1D/2D output is flushed at the same time, so a killed job keeps a usable restart.
- The exit status is nonzero after a PETSc error. Input errors stop with a message and `assert(0)`
  (exit 134); asserts are always on (no `NDEBUG`).
- Modes (`main.cpp`): `bulkDeformationType` = `linearElastic` | `powerLaw`; `momentumBalanceType`
  = `quasidynamic` | `dynamic` (elastic only) | `quasidynamic_and_dynamic`;
  `systemEvolutionType = steadyStateIts` (power law) iterates to a steady state, and each outer
  iteration runs a full transient to `maxTime`. Special drivers: `isMMS = 1`,
  `computeGreensFunction_fault = 1`, `computeGreensFunction_offFault = 1`.

## Input files

- One `key = value` per line; the delimiter is exactly ` = `. A scalar value ends at the first
  space, so anything after it (usually `# comment`) is ignored. Lists keep the whole remainder:
  `aVals = [0.0135 0.0300]`. Lines without ` = ` are skipped.
- **Unknown or misspelled keys are silently ignored.** Run `tools/checkkeys.py file.in` to list keys
  that no component reads. It cannot see a key read by the wrong component: `linSolver` is read only
  by the pressure equation; the momentum balance reads `linSolverSS`/`linSolverTrans` and the heat
  equation `linSolver_heateq`.
- Depth profiles: `<name>Vals = [...]` with `<name>Depths = [...]` (km), interpolated linearly in z;
  equal lengths, non-decreasing depths; a repeated depth makes a step (`lockedDepths = [0 40 40 500]`).
- Every component re-reads the whole file in its own `loadSettings`, so key names are per
  component: momentum-balance boundaries are `momBal_bcR_qd`, `momBal_bcT_qd`, `momBal_bcL_qd`,
  `momBal_bcB_qd` (and `_fd` variants in the combined mode); creep laws use prefixes `disl_`,
  `disl2_`, `diff_`, `dp_`.
- Units: lengths km (grid, depths), slip and Dc m, vL m/s, time s, mu GPa, rho g/cm^3, stresses MPa.
- With hydraulic coupling, `sNVals` is the **total** normal stress (sNEff = sN - p, p including the
  hydrostatic part when g > 0). The examples' comments describe the no-pressure case.
- Rejected at startup: `sbpCompatibilityType = compatible` (fails MMS), `timeControlType = PI`
  (not implemented), implicit pressure with an explicit-only integrator.
- Initial conditions can be loaded as PETSc binary Vecs named `y`, `z`, `psi`, `slip`, `prestress`,
  `tauQS`, `tau`, `fault_a`, ... from `inputDir` (see `SEAS_benchmarks/BP1/createICs.py`).
- Faults have names. The default fault, `fault`, reads the plain keys. A fault named `<name>` reads
  them too, then `<name>_<key>` overrides any of them (`fault2_aVals`, `fault2_vCreep`, ...); it uses
  the files `<name>_psi`, ..., the HDF5 groups `/<name>`, `/<name>_qd` and integrand keys
  `<name>_slip`, `<name>_psi`. `vCreep` sets the slip velocity of creeping nodes (default `vL`).
  No input key declares an additional fault yet (Stage 2 of `docs/TWO_FAULT_DESIGN.md`).

## Architecture

- `main.cpp`: builds `Domain` (grid, coordinates, scatters, run flags) and dispatches to a mediator.
- **Mediators** (`strikeSlip_*.cpp`, all `ProblemContext`): `StrikeSlip_LinearElastic_qd`,
  `_LinearElastic_fd`, `_LinearElastic_qd_fd`, `_PowerLaw_qd`, `_PowerLaw_qd_fd`. Each owns the
  material (`LinearElastic` or `PowerLaw`), the fault(s) (`Fault_qd`, `Fault_fd`), optional
  `HeatEquation`, `PressureEq`, `GrainSizeEvolution`, the time integrator, and all output.
  They are largely copy-pasted: **a fix in one usually has to be mirrored in the others.**
  `StrikeSlip_LinearElastic_qd` holds a list of faults (`_faults`, boundary fault `_fault` first) and
  loops over it for per-fault work; the other mediators still hold a single fault.
- **Integrand**: `map<string,Vec>`. Explicit (`varEx`): `slip`, `psi`, `gVxy`/`gVxz` (power law),
  `grainSize`, and `pressure`/`permeability` when `hydraulicTimeIntType = explicit`. Implicit
  (`varIm`, backward Euler inside IMEX): `Temp`, and `pressure`/`permeability` when implicit.
- **`d_dt` sequence (quasi-dynamic)**: bcL = slip/faultTypeScale (symmFault: slip = 2u at y = 0),
  bcR = vL t/2 + bcRShift (remoteLoading) -> solve A u = rhs (KSP; MUMPS Cholesky by default,
  factored once) -> stresses -> sxy on the fault + prestress = tauQS -> rate-and-state root solve
  for V (`ComputeVel_qd`, bracketed Newton) -> state-law rate -> optional heat/pressure/grain size.
  IMEX steps then solve the implicit fields and recompute explicit rates with them.
- **SBP operators**: `sbpOps_m_varGrid` (curvilinear, used whenever `bCoordTrans > 0` or y/z come
  from files) and `sbpOps_m_constGrid`; `spmat.cpp` builds the 1D blocks; only `fullyCompatible`
  operators are allowed. Order 2 or 4 (`order`).
- **Integrators**: `odeSolver.cpp` (RK32, RK43, FEuler), `odeSolverImex.cpp` (RK32_WBE, RK43_WBE),
  `odeSolver_WaveEq.cpp` (fully dynamic leapfrog). Step control uses `timeIntInds` (empty = all
  explicit variables) with `scale`, norm `normType`, controller P or PID. Quasi-dynamic steps never
  go below `minDeltaT` (default: one grid wave-transit time) and are accepted there even when the
  error exceeds `timeStepTol`.
- **Combined mode** switches to fully dynamic when max V exceeds `trigger_qd2fd` and back below
  `trigger_fd2qd` (`checkSwitchRegime`, `prepare_qd2fd`, `prepare_fd2qd`).
- **Output** (HDF5, PETSc timestepping, one index per output step): `data_context.h5` (coordinates,
  parameters), `data_1D.h5` (fault and boundary series every `stride1D` steps: `/time`, `/fault`,
  `/momBal`, `/heatEquation`, `/pressureEq`), `data_2D.h5` (body fields every `stride2D`),
  `data_steadyState.h5`, `checkpoint.h5`, and text context files (`domain.txt`, `mediator.txt`, ...).
  Each mediator's `writeStep1D/2D` must run first in `timeMonitor`: it opens the files and advances
  the HDF5 time index that the components write into.
- **Loaders**: `examples/loadFuncs.py` (all datasets by default; reshapes to (Ny, Nz[, Nt]) and
  (N[, Nt]) with time last) and `matlab/visualizePetsc/load*_hdf5.m`.

## Conventions and traps

- y is distance across the fault, z is depth (positive down). Body index `Ii = iy*Nz + iz`
  (z fastest). The fault is the **left boundary** y = 0; there is no interior fault yet.
- `Domain::_y0` is the size-Nz left-boundary template and holds **z** there; `_z0` is the size-Ny
  top-boundary template and holds **y** there.
- `Nz = 1` is a spring slider (ex1). The combined mode refuses it.
- `lockedVals`: > 0.5 locked, < -0.5 creeps at vL, otherwise frictional.
- The aging law freezes psi where b <= 1e-3 (and where exp((f0-psi)/b) overflows).
- Explicit slip-dependent permeability is unstable when |V| dt > ~2.8 kL_p; with quasi-dynamic
  minimum steps this happens during events. Use `hydraulicTimeIntType = implicit` (exact relaxation).
- `momBal_bcT_qd`/`momBal_bcB_qd = remoteLoading` stay at their initial displacement.

## Testing

```bash
tools/regress.sh baseline <dir>      # run ex1, ex2 (fresh, outputDir redirected) and store output
tools/regress.sh compare  <dir>      # rerun and h5diff every dataset; exit 1 on any difference
SCYCLE_BIN=... REGRESS_NP=2 REGRESS_CASES="ex1 ex2" REGRESS_DELTA=1e-12 REGRESS_WORK=... tools/regress.sh ...
./source/main tools/mms.in           # MMS convergence, Ny = Nz = 21, 41, 81 (order 4: u ~3.5, sxy ~2.5)
tools/checkkeys.py examples/*.in     # keys no component reads
python3 SEAS_benchmarks/BP1/createICs.py [--dz 0.1]   # BP1 grid and initial conditions
./source/main <input> -objects_dump  # PETSc objects never freed (none expected)
```

- Runs are bit-reproducible for a given build and rank count; 1 vs 2 ranks agree to round-off
  (~1e-10). Debug and optimized builds differ slightly once adaptive steps diverge (ex1: 2201 vs
  2202 steps, event onsets within 0.02 yr over 4500 yr), so compare a build only with a baseline
  made by the same PETSC_ARCH.
- A baseline from the end of the audit (debug build, 1 rank) is in `data/regress-baseline/`
  (ignored by git; `README.txt` there records its commit). Regenerate it after any intended change
  of results and say so in the commit message.
- ex2 takes about 70 s with the debug build; ex1 about 13 s.

## Working agreements

- One branch per task; one commit per fix, with the defect, the failure it caused and the test
  evidence in the message. Never push without asking.
- Run `tools/regress.sh compare` after any change under `source/`. ex1/ex2 must stay bit-identical
  unless the change is meant to alter results; then say how much and why.
- Mirror fixes across the mediator classes; build with `WERROR=1`.
- Physics changes that alter published behaviour get their own commit, marked in the message, so
  they can be reverted alone.
