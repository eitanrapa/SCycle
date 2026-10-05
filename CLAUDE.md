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

Upstream: `origin` = github.com/kali-allison/SCycle (`master`). All work is on `master` of the
user's fork (remote `fork`, github.com/eitanrapa/SCycle), ahead of upstream's: the audit fixes,
two-fault stages 1 to 4, then stage 5, which continues there. Tags mark where each part ended:
`audit-2026-10` and `stage1` to `stage4` (they were branches until 2026-10-03, when the user asked
for a single branch). `docs/REVERSIBLE_STRENGTH_PLAN.md` plans and tracks stage 5.

## Build

Requires PETSc >= 3.14 (HDF5 timestepping API) with HDF5 and MUMPS, real scalars, 32-bit indices;
HYPRE for the `AMG` and `CG_PCAMG` solver options. Since 2026-10-04 the work runs on a Mac Studio
(M2 Ultra, 24 cores, 192 GB; before that an M4 Mac). PETSc 3.26.0 is built from source, both arches,
by `~/opt/build-mpi-petsc.sh` (logs in `~/opt`), with the configure options of `docs/SERVER.md`:

```bash
export PETSC_DIR=$HOME/opt/petsc-3.26.0
export PETSC_ARCH=arch-darwin-c-debug     # or arch-darwin-c-opt (-O3) for production runs
export PATH=$HOME/opt/openmpi-5.0.11/bin:$PATH   # mpirun
make -C source -j8                        # -> source/main
make -C source -j8 WERROR=1               # warnings as errors (the tree is warning-free)
make -C source DEBUG_MODULES=-DVERBOSE=2  # more tracing; VERBOSE=1 (default) prints one line per step
make -C source clean                      # also after switching PETSC_ARCH (objects are per arch)
```

Two traps on this machine. Homebrew (`/opt/homebrew`) belongs to another account, so `brew install`
fails; Open MPI 5.0.11 is therefore built into `~/opt/openmpi-5.0.11`, and PETSc uses Homebrew's
CMake 4.4.3 (`--with-cmake-exec`) instead of `--download-cmake`. And the Command Line Tools keep a
stale `/Library/Developer/CommandLineTools/usr/include/c++/v1` (four files from 2020-22) that hides
the SDK's libc++: plain `clang++` cannot compile `#include <vector>`. `mpicxx` passes the SDK's
directory itself (`-stdlib++-isystem`, set in `~/opt/openmpi-5.0.11/share/openmpi/mpic++-wrapper-data.txt`),
so builds through PETSc work; deleting that stale directory (sudo) would fix `clang++` for everyone.
`h5diff`/`h5dump`/`h5ls` come from Homebrew `hdf5`, `gfortran` from Homebrew `gcc` (16).
Linux servers: `docs/SERVER.md` (PETSc configure line, build, checks, batches). GCC builds the tree
without warnings too (checked with g++ 16 against libstdc++).
Python with numpy, h5py, scipy and matplotlib for the tools and loaders: `~/venvs/scycle/bin/python`.
MATLAB R2023a is installed (its license not checked); the MATLAB loaders are untested.

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
  complete; the 1D/2D output is flushed at the same time.
- **Stop a run with SIGTERM or SIGINT** (`kill <pid>`, Ctrl-C): it finishes the step, writes a
  checkpoint, closes its files and exits with 128 + the signal (143, 130); a second signal ends it
  at once. With MPI, signal the ranks, not `mpirun` (Open MPI's kills them at once). A run killed
  outright (`kill -9`, power loss) can leave its HDF5 output unreadable between checkpoints.
- `tools/batch.sh start|status|stop <dir>` runs the inputs `<dir>/*.in` detached (outliving the
  shell), resumes them after a stop or reboot, and records each run's end in `<outputDir>run.exit`.
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
  by the pressure equation (`AMG`, the default, or `MUMPSLU`); the momentum balance reads
  `linSolverSS`/`linSolverTrans` and the heat equation `linSolver_heateq`.
- Depth profiles: `<name>Vals = [...]` with `<name>Depths = [...]` (km), interpolated linearly in z;
  equal lengths, non-decreasing depths; a repeated depth makes a step (`lockedDepths = [0 40 40 500]`).
- Every component re-reads the whole file in its own `loadSettings`, so key names are per
  component: momentum-balance boundaries are `momBal_bcR_qd`, `momBal_bcT_qd`, `momBal_bcL_qd`,
  `momBal_bcB_qd` (and `_fd` variants in the combined mode); creep laws use prefixes `disl_`,
  `disl2_`, `diff_`, `dp_`.
- Units: lengths km (grid, depths), slip and Dc m, vL m/s, time s, mu GPa, rho g/cm^3, stresses MPa.
  Pore pressure: k m^2, eta_p Pa s, beta_p 1/MPa, rho_f g/cm^3, g m/s^2, p MPa (consistent with z
  in km, since m^2/Pa = km^2/MPa); the recharge flux `bcB_q0` m/s.
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
- Interior faults (both quasi-dynamic mediators, elastic and power law): `interiorFaults = [name ...]`,
  `<name>_y` (km) for each, and `momBal_bcL_qd = remoteLoading` for the full domain. Each is placed
  midway between two grid rows, at least 6 rows from the y-boundaries; `mediator.txt` records where.
  Not with a boundary fault (`symmFault` would mirror it) or the steady-state initial guess (or
  `steadyStateIts`, `computeSSHeatEq`). `interiorFaultKinkLift = 1` (default) gives second-order fault
  traction. Heat needs `wVals > 0` (one Gaussian per fault, normalized on the grid; set
  `bcLType_trans = Dirichlet`). With pore pressure each fault has its own `PressureEq`: a fault
  `<name>` reads `<name>_` overrides of the pressure keys and writes `/<name>_pressureEq`.
  `examples/interior_fault/make_inputs.py` writes one interior fault and its half-space twin;
  `examples/two_faults/make_inputs.py` two faults (seismic, creeping or locked second fault) on a
  grid refined around both, each fault exactly midway between two rows; `--rheology powerlaw` gives
  the stage 4 baseline (ex4's creep and geotherm to 60 km, coupled heat, optional grain size).
- With interior faults, two outputs are on by default (off otherwise): `computeSurfVel = 1` writes
  the instantaneous surface velocity `/momBal/surfVel` (m/s) with each 1D output, and
  `strideSeries = N` appends per-fault max slip rate, its depth, potency rate and potency to
  `faultSeries.txt` every N steps (restarts append; keep the last line per step).
  `data_context.h5` gives each interior fault's `y`, `iRow` and `dy` as attributes of its group.
  `tools/two_fault.py <outputDir>` writes `events.csv`, `partition.csv` and `surfvel.csv`.
- The base (`momBal_bcB_qd`) is traction-free by default, so with loading only at the sides the
  depth-averaged velocity is a straight line to the loaded side and the surface velocity meets vL/2
  there with a slope (docs/TWO_FAULT_DESIGN.md, section 7). `movingBase` (quasi-dynamic mediators)
  moves the base with the plates: (vL t/faultTypeScale) tanh((y - yc)/w) plus its initial
  displacement, `momBal_bcB_center` (km, default mid-domain; a boundary fault: 0) and
  `momBal_bcB_width` (km, holds 90% of the change, 0 a step; about 60 km for the two-fault models).
  Faults must be locked where they meet it (a boundary fault under a step must not be).
- A `timeIntInds` that lists `slip` or `psi` means those of every fault: `<name>_slip`, `<name>_psi`
  are added with the same scale (a note prints the final list).
- Stage 5 mechanisms (`docs/REVERSIBLE_STRENGTH_PLAN.md`, section 6 lists what is done). Bulk state
  fields, `StrikeSlip_PowerLaw_qd` only (the other mediators refuse them): `<prefix>type = off |
  transient | constant` with strain hardening `hard_` (needs `wDislCreep = yes`) and water content
  `water_` (through `wDislWetDry = yes`, which mixes `disl_` as the dry and `disl2_` as the wet
  law, or pressure solution); each writes `/<name>` to data_2D.h5, joins a non-empty
  `timeIntInds`, and `<prefix>eTest` drives it with a prescribed strain rate for tests. Grain size:
  the cataclastic sink `grainSizeEv_fCatVals` (a fraction of the faults' work spread over the
  frictional-heat kernel, so it needs `wVals > 0`; the heat equation gets the rest) and the Zener
  cap `grainSizeEv_dZVals`. Pore pressure (every fault's `PressureEq`): `bcB_q0` prescribes the
  basal flux, `bcB_pulse*` pulse it, `bcB_sourceDepth`/`bcB_sourceWidth` deliver it at depth.
  Fault cohesion (both quasi-dynamic mediators, per fault): `cohesionEvolution = explicit |
  implicit` with `cohesionMaxVals`, `cohesionTauHeal`, `cohesionDheal` (and `cohesionReseal`);
  implicit needs an IMEX integrator and bounds the step (`cohesionStepFrac`), explicit puts
  `<name>_cohesion` under error control and is slow in events when `cohesionDheal` is millimetres.

## Architecture

- `main.cpp`: builds `Domain` (grid, coordinates, scatters, run flags) and dispatches to a mediator.
- **Mediators** (`strikeSlip_*.cpp`, all `ProblemContext`): `StrikeSlip_LinearElastic_qd`,
  `_LinearElastic_fd`, `_LinearElastic_qd_fd`, `_PowerLaw_qd`, `_PowerLaw_qd_fd`. Each owns the
  material (`LinearElastic` or `PowerLaw`), the fault(s) (`Fault_qd`, `Fault_fd`), optional
  `HeatEquation`, `PressureEq`, `GrainSizeEvolution`, the time integrator, and all output.
  They are largely copy-pasted: **a fix in one usually has to be mirrored in the others.**
  `StrikeSlip_LinearElastic_qd` and `StrikeSlip_PowerLaw_qd` hold a list of faults (`_faults`,
  boundary fault `_fault` first) and loop over it for per-fault work; the multi-fault pieces they
  share are in `source/multiFault.hpp`. The other mediators still hold a single fault.
- **Bulk state fields** (`bulkStateField.hpp`; `hardeningState`, `waterState`): the power-law
  mediator holds them in `_bulkStates`, pushes each into `PowerLaw` before the momentum balance
  and evaluates its rate after the faults' (so is the grain size, whose sink takes their work).
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
  explicit variables) with `scale`, norm `normType`, controller P or PID. Steps never go below
  `minDeltaT` and are accepted there even when the error exceeds `timeStepTol` (a warning says so
  and the run summary counts such steps). The quasi-dynamic default is 1e-3 of the shear-wave time
  of the smallest cell; the old default, that time itself, overrode error control in events (A-95).
- **Combined mode** switches to fully dynamic when max V exceeds `trigger_qd2fd` and back below
  `trigger_fd2qd` (`checkSwitchRegime`, `prepare_qd2fd`, `prepare_fd2qd`).
- **Output** (HDF5, PETSc timestepping, one index per output step): `data_context.h5` (coordinates,
  parameters), `data_1D.h5` (fault and boundary series every `stride1D` steps: `/time`, `/fault`,
  `/momBal`, `/heatEquation`, `/pressureEq`; `/<name>` per fault), `data_2D.h5` (body fields every
  `stride2D`), `faultSeries.txt` (per step, multi-fault runs),
  `data_steadyState.h5`, `checkpoint.h5`, and text context files (`domain.txt`, `mediator.txt`, ...).
  Each mediator's `writeStep1D/2D` must run first in `timeMonitor`: it opens the files and advances
  the HDF5 time index that the components write into.
- **Loaders**: `examples/loadFuncs.py` (all datasets by default; reshapes to (Ny, Nz[, Nt]) and
  (N[, Nt]) with time last) and `matlab/visualizePetsc/load*_hdf5.m`.

## Conventions and traps

- y is distance across the fault, z is depth (positive down). Body index `Ii = iy*Nz + iz`
  (z fastest). A fault is either the **left boundary** y = 0 (symmFault: a half-space, mirrored) or
  an interior fault of the full domain (`InteriorFaultLift`, a jump-corrected right-hand side).
- `Domain::_y0` is the size-Nz left-boundary template and holds **z** there; `_z0` is the size-Ny
  top-boundary template and holds **y** there.
- `Nz = 1` is a spring slider (ex1). The combined mode refuses it.
- `lockedVals`: > 0.5 locked, < -0.5 creeps at vL, otherwise frictional.
- The aging law freezes psi where b <= 1e-3 (and where exp((f0-psi)/b) overflows).
- Explicit slip-dependent permeability needs |V| dt < ~2.8 kL_p, so an event with kL_p = 1 mm takes
  ~1e5 steps. Use `hydraulicTimeIntType = implicit` (exact relaxation).
- Near an interior fault (rows iRow-2 .. iRow+3) the y-strain comes from u without the fault's jump
  (`InteriorFaultLift::correctStress` in `LinearElastic::computeStresses`, `correctStrain` in
  `PowerLaw::computeTotalStrains`); anything new that differentiates u in y must do the same.
- `momBal_bcT_qd`/`momBal_bcB_qd = remoteLoading` stay at their initial displacement (`movingBase`
  is the bottom boundary that moves).

## Testing

```bash
tools/regress.sh baseline <dir>      # run ex1, ex2, ex4s, ex4g (fresh, outputDir redirected), store output
tools/regress.sh compare  <dir>      # rerun and h5diff every dataset; exit 1 on any difference
SCYCLE_BIN=... REGRESS_NP=2 REGRESS_CASES="ex1 ex2" REGRESS_DELTA=1e-12 REGRESS_WORK=... tools/regress.sh ...
./source/main tools/mms.in           # MMS convergence, Ny = Nz = 21, 41, 81 (order 4: u ~3.5, sxy ~2.5)
tools/checkkeys.py examples/*.in     # keys no component reads
tools/compare_runs.py REF OTHER      # regression outputs of different builds: events and end state
python3 SEAS_benchmarks/BP1/createICs.py [--dz 0.1]   # BP1 grid and initial conditions
./source/main <input> -objects_dump all  # every PETSc object never freed (none expected; without
                                         # "all" only objects made by a Create call are listed)
```

- Runs on 1 rank are bit-reproducible for a given build, also across a stop and restart. Runs on 2
  ranks are not: two identical 2-rank runs of ex2 differ by 2.6e-9 (relative, slip), so compare them
  with a tolerance. Debug and optimized builds differ slightly once adaptive steps diverge (ex1:
  2201 vs 2202 steps, event onsets within 0.02 yr over 4500 yr), so compare a build only with a
  baseline made by the same PETSC_ARCH, and different builds with `tools/compare_runs.py`.
  Another machine is another build too: on the Mac Studio ex1 reproduces the M4 Mac's baseline bit
  for bit, but ex2, ex4s and ex4g differ in the last digits (probably MUMPS through Accelerate's
  BLAS, which is tuned per chip) while agreeing in the physics (`tools/compare_runs.py`: the
  same events, onsets within 1.5e-6 yr in ex2, 6.5e-5 yr in ex4s, 1.5e-4 yr in ex4g).
- The baseline (debug build, 1 rank) is in `data/regress-baseline/` (ignored by git), made on the
  Mac Studio at `64a5211` (`README.txt` there); the M4 Mac's is kept in `data/regress-baseline-m4/`
  and the one from the end of the audit in `data/regress-baseline-35ba52f/`. Regenerate it after any
  intended change of results and say so in the commit message.
- With the debug build ex1 takes about 30 s, ex2 140 s, ex4s and ex4g 50 s each on the Mac Studio
  while a batch runs (13, 70 and 30 s on the M4); the optimized build runs ex2 in about 65 s and
  ex4s, ex4g in 5 s. ex4s and ex4g are short power-law runs (coupled heat with the
  implicit-explicit integrator; grain size coupled to diffusion creep with the explicit one), the
  only cases that exercise `StrikeSlip_PowerLaw_qd`.

## Working agreements

- Work on `master`, without a branch per task or stage; one commit per fix, with the defect, the
  failure it caused and the test evidence in the message. Never push without asking.
- Run `tools/regress.sh compare` after any change under `source/`. ex1/ex2 must stay bit-identical
  unless the change is meant to alter results; then say how much and why.
- Mirror fixes across the mediator classes; build with `WERROR=1`.
- Physics changes that alter published behaviour get their own commit, marked in the message, so
  they can be reverted alone.
