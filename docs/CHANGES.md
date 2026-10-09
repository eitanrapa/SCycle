# Changes from upstream SCycle

This fork (github.com/eitanrapa/SCycle, `master`) starts from upstream master
(github.com/kali-allison/SCycle, `74a132f`, 19 December 2024), which has had no commits since. By
`df6e3c4` (6 October 2026) the fork is 147 commits ahead: in `source/`, 60 files with 6,551 lines
added and 1,588 removed; elsewhere, 77 files with 5,216 added and 10,025 removed (most of the
removals are the stale `tests/` directory).

This file is an index. The details are elsewhere:

- `docs/AUDIT.md`: every fix, with the defect, the failure it caused and the test evidence
  (findings A-01 to A-107), how much each changes results (its section 4) and open items (O-01 to
  O-15).
- `docs/TWO_FAULT_DESIGN.md`: the two-fault extension, stages 1 to 4 (tags `stage1` to `stage4`).
- `docs/REVERSIBLE_STRENGTH_PLAN.md`: stage 5, the reversible-strength mechanisms and the runs.
- Each commit message (`git show <hash>`).

## 1. Fixes that change results

Most fixes turn a crash, a hang or lost output into a result. These change the numbers of runs that
already worked. Anyone comparing with results made with the original code, such as the published
papers, should check which apply. Sizes are in `docs/AUDIT.md`, section 4.

- **Every run.** The Jacobian of the velocity solve (A-23; the converged velocity is unchanged to
  `rootTol`), roots accepted once the bracket collapses (A-24; a few ulps), the 4th-order operators'
  R coefficient (A-33, A-34; 4e-8 relative in ex2), steps that reach `minDeltaT` (A-50), and the
  default `minDeltaT`, now 1e-3 of the smallest cell's shear-wave time (A-95; ex2's peak slip rate
  +3.4%, onset and slip unchanged).
- **Heat equation.** The heat source dropped on constant grids (A-27); the MUMPS factorization reused
  with a stale time step (A-28); the fault heat flux's sign in the steady-state and explicit solvers
  (A-29); two-point geotherms lost (A-30); **viscous dissipation 41% too large in the heat source and
  half its true value in the grain-size sink (A-31, marked PHYSICS CHANGE)**; implicit-explicit
  coupling one step late (A-55); the frictional-heat kernel not normalized when the grid does not
  resolve the shear-zone width (A-97).
- **Power law.** The second dislocation mechanism (A-41) and the units of the initial viscosity
  guess (A-42); dissolution-precipitation creep wrong by a factor of the stress, and off by default
  (A-43, A-44); the steady-state viscosity iteration's damping (A-45); the grain-size step cap
  (A-48); the atan top-boundary loading (A-46); the written fault stress, output only (A-96); the
  time-step bound across restarts (A-102).
- **Pore pressure.** The initial effective normal stress (A-39), the order of the rate evaluations
  (A-40), the Picard stop (A-47), slip-dependent permeability (A-57), `pVals` (A-77), the bottom
  flux (A-100), `linSolver` (A-101).
- **Fully dynamic phases.** Absorbing boundaries (A-26), cohesion (A-37, also quasi-dynamic),
  creeping nodes (A-38), the hand-off of the particle velocity (A-56).
- **Inputs.** The shear modulus profile, now in depth (A-32), list parsing (A-70), and the example
  inputs (A-88 to A-91; ex2 now has an earthquake at 54.5 yr).

## 2. All fixes

107 findings: 34 Critical, 27 High, 30 Medium, 16 Low (severity scale in `docs/AUDIT.md`,
section 2). "Latent" means found by reading the code; the failure was not reproduced.

### Crashes, undefined behaviour and memory

| | What was wrong | Severity | Commit |
|---|---|---|---|
| A-01 | Unsupported problem types ran through an uninitialized pointer | Critical | `e9cdc49` |
| A-02 | BracketedNewton evaluated its first residual at an unset member | Critical (latent) | `b7aa7b7` |
| A-03 | Flash heating gave NaN at locked fault nodes | Critical | `5c2606b` |
| A-04 | Fault_fd::setPhi re-acquired arrays instead of restoring them | Critical | `0802469` |
| A-05 | The pressure state handed to the integrator was freed | Critical | `5c87a81` |
| A-06 | Grain-size flags and output strides unset in the power-law qd_fd class | Critical | `4931d1c` |
| A-07 | Fully dynamic runs wrote to time Vecs that were never created | Critical | `a54f046` |
| A-08 | Members read before being set in six classes | Critical | `bcf2fed` |
| A-09 | Radiogenic decay length unset; ex4 produced NaN temperatures | Critical | `0e844ac` |
| A-10 | A zero shear-zone width divided by zero | Critical | `25e29c9` |
| A-11 | Matrix preallocation arrays on the stack | Critical (latent) | `009d188` |
| A-12 | The heat equation was used without being created | Critical | `aece9f9` |
| A-13 | Power-law qd_fd built the ice-stream forcing before the material existed | Critical | `e77dcda` |
| A-14 | The ice-stream context was written to a destroyed viewer | Critical | `2febba7` |
| A-15 | Pressure fields missing from the integrand at phase changes and in steady-state iterations | Critical | `04298a8` |
| A-16 | Implicit pore pressure aborted on constant grids | Critical | `ca11678` |
| A-17 | updateBCMats destroyed the wrong boundary matrices | Medium (latent) | `ca11678` |
| A-18 | LinearElastic::setupKSP created one solver and configured another | Medium | `ca11678` |
| A-19 | The slip-rate ceiling was applied while the array was checked out | Medium | `8954864` |
| A-20 | Power-law qd_fd copied a temperature entry that does not exist | Medium | `28b40f9` |
| A-21 | The fault Green's-function driver failed in parallel | Critical | `89a6b0b` |
| A-22 | PETSc objects leaked | Medium | `2febba7`, `9866e9b`, `f8fa0e4`, `4942bd8`, `aece9f9`, `e77dcda`, `89a6b0b`, `8916df7` |
| A-98 | The implicit-explicit integrators never freed their implicit-variable Vecs | Low | `8df49c8` |

### Physics and numerics (wrong results)

| | What was wrong | Severity | Commit |
|---|---|---|---|
| A-23 | Wrong Jacobian in the quasi-dynamic velocity solve | Medium | `da28ff5` |
| A-24 | The root finder never accepted a bracket collapsed to round-off | Critical | `74fdbb9` |
| A-25 | Bisect::setBounds lost the bracket when the bounds were reversed | Low | `9be7f3d` |
| A-26 | Absorbing boundaries were disabled in dynamic phases | High | `4777e9a` |
| A-27 | The transient heat solve dropped the heat source on constant grids | High | `e665cde` |
| A-28 | The transient heat solve reused the first step's factorization | High | `b8f197e` |
| A-29 | The fault heat flux had the wrong sign in two heat solvers | High | `24f6732` |
| A-30 | Two-point temperature profiles lost the geotherm | High | `06bada2` |
| A-31 | Viscous dissipation used wrong factors | High | `574a54b` |
| A-32 | The shear modulus profile was interpolated across the fault-normal direction | High | `9608cec` |
| A-33 | Fourth-order operators depended on the number of MPI ranks | High | `b65b5e0` |
| A-34 | R_y averaged the wrong neighbours | High | `b65b5e0` |
| A-35 | Compatible 4th-order boundary rows had coefficient errors | High, partly fixed | `e892279` |
| A-36 | getMus returned mu for the metric-weighted moduli | Low (latent) | `dd9de2f` |
| A-37 | The documented cohesion had no effect | High | `6b3dad7`, `aa53679` |
| A-38 | Creeping nodes slipped freely in dynamic phases | High | `8954864` |
| A-39 | The initial effective normal stress ignored the pore pressure | High | `61bbac0` |
| A-40 | Power-law runs did not evolve slip-dependent permeability | High | `325945f` |
| A-41 | The second dislocation mechanism was only partly wired in | Medium | `9745bb2` |
| A-42 | The initial viscosity guess was off by a units factor | Low | `9745bb2` |
| A-43 | Dissolution-precipitation creep returned a strain rate as an inverse viscosity | High | `b42c53d` |
| A-44 | Dissolution-precipitation creep did nothing by default | Medium | `b42c53d` |
| A-45 | The steady-state viscosity iteration diverged without warning | High, partly fixed | `c366a1b` |
| A-46 | Top-boundary loading used depths as horizontal coordinates | High | `e77dcda` |
| A-47 | The implicit pressure Picard loop ignored its tolerance | Medium | `4942bd8` |
| A-97 | The boundary fault's frictional-heat kernel was not normalized on the grid | High | `3f402f3` |
| A-100 | The pressure equation's bottom flux used the ends of the depth lists, not the bottom node | High | `e0ae91a` |

### Time integration, checkpoints and I/O

| | What was wrong | Severity | Commit |
|---|---|---|---|
| A-48 | The grain-size limit could override the Maxwell time-step limit | Medium | `00b3b33` |
| A-49 | The next step size could be unset or stale | Critical | `39f8663` |
| A-50 | Rejected steps were accepted with the wrong step size | High | `fbe6ee0` |
| A-51 | Checkpoints stored a never-set error estimate | Low | `fbe6ee0` |
| A-52 | The relative error norms were broken | Critical | `d823c7b` |
| A-53 | The PID controller produced NaN in integrators created mid-run | Critical | `d823c7b` |
| A-54 | FEuler ignored stop requests and could not checkpoint in parallel | Medium | `d823c7b` |
| A-55 | IMEX coupling lagged one step | High | `185a2ac` |
| A-56 | The quasi-dynamic to dynamic hand-off lost the particle velocity | High | `b78b1cb` |
| A-57 | Slip-dependent permeability went unstable during events | Critical, partly fixed | `04298a8`, `35ba52f` |
| A-58 | RK32 restarts read the wrong HDF5 group | Critical | `ab5973c` |
| A-59 | Checkpoint and steady-state files stayed open after loading | Medium | `ab5973c`, `9866e9b` |
| A-60 | qd_fd regime flags were checkpointed through the wrong type | Critical | `da39e0e` |
| A-61 | Restart flags without a checkpoint file aborted the run | Critical | `14a1331` |
| A-62 | Pore pressure was never written or checkpointed | Critical | `4942bd8` |
| A-63 | Steady-state restarts overwrote the first dislocation mechanism | Medium | `9866e9b` |
| A-64 | sxz was never written; sdev was misnamed | Critical | `f8fa0e4` |
| A-65 | A job killed during a checkpoint write could lose the only restart file | Critical (latent) | `3949db3` |
| A-66 | The off-fault Green's function stopped after 10 columns | Critical | `89a6b0b` |
| A-67 | The 2D output stride was taken from the 1D stride | Medium | `28b40f9`, `2febba7` |
| A-68 | The power-law qd_fd class ignored initTime | Medium | `28b40f9` |
| A-95 | The default minimum time step overrode error control during earthquakes | Critical | `5bb2fee` |
| A-96 | The power-law explicit rate function overwrote the written fault stress | Medium | `ae6c6a8` |
| A-99 | Restarts changed the boundary fault's frictional-heat kernel by roundoff | Low | `fc8938f` |
| A-102 | Restarts lost the time-step bound set by the monitor | Medium | `b03b76e` |
| A-105 | A run killed between checkpoints could lose its output | Critical | `4df9668` |
| A-106 | A restart from a checkpoint taken before the first output shifted the output by one entry | Medium | `8986bdb` |
| A-107 | Restarted power-law runs wrote a steady-state index that fresh runs do not | Low | `1a81344` |

### Input handling and error reporting

| | What was wrong | Severity | Commit |
|---|---|---|---|
| A-69 | Failed runs exited with status 0 | Medium | `e9cdc49` |
| A-70 | The list parser dropped or invented entries | High | `7232307` |
| A-71 | Depth profiles were not validated | Critical | `aece9f9` |
| A-72 | A missing input file was not reported | Medium | `aece9f9` |
| A-73 | The NaN and Inf checks in dissolution-precipitation creep could not fire | Medium | `2005eed` |
| A-74 | A rigidFault left boundary was always rejected in fully dynamic runs | Critical | `28b40f9` |
| A-75 | The regime-switch fallback compared limit_qd with the wrong trigger | Medium | `b78b1cb` |
| A-76 | Root-finder iteration counts were never recorded | Low | `9be7f3d` |
| A-77 | The initial pore pressure pVals was ignored | High | `4942bd8` |
| A-78 | The mediators did not parse hydraulicTimeIntType; implicit pressure could freeze | High | `04298a8` |
| A-79 | The mediators ignored errors from the pressure and grain-size rates | Critical | `04298a8` |
| A-101 | The pore-pressure solver ignored `linSolver` | Medium | `6e39b4d` |

### Build, repository and tooling

| | What was wrong | Severity | Commit |
|---|---|---|---|
| A-80 | The Makefile named a nonexistent linker and hid link failures | Low | `44a81b1` |
| A-81 | Stale header dependencies and mandatory -Werror | Low | `44a81b1` |
| A-82 | Binaries, OS files and stale test sources were tracked | Low | `8c44ddc`, `f82d430` |
| A-83 | There was no automated regression check | Low | `26139e7`, `9c0cc61` |
| A-84 | The manufactured-solution test restarted from the previous resolution | Medium | `14a1331`, `e0b9041` |
| A-85 | Misspelled or retired input keys were ignored without a warning | Low | `df64c96` |
| A-86 | Unreachable debug drivers with hard-coded paths | Low | `89a6b0b` |
| A-87 | Misleading comment on Domain::_y0 and _z0 | Low | `c844c59` |
| A-103 | GCC could not compile the multi-fault helpers | Medium | `933d1bd` |
| A-104 | GCC warned at ten unchecked PetscFree calls | Low | `89da3ae` |

### Examples, benchmarks and post-processing

| | What was wrong | Severity | Commit |
|---|---|---|---|
| A-88 | ex1 and ex2 started from zero shear stress | High | `6b2315e` |
| A-89 | ex3 ran the wrong problem class | High | `6b2315e` |
| A-90 | ex4 and ex5 failed input validation | Medium | `6b2315e` |
| A-91 | The examples' momentum-balance solver keys were ignored | Medium | `6b2315e` |
| A-92 | The Python loader did not read current output | Medium | `f3b5829` |
| A-93 | The MATLAB loaders called functions that were never committed | Medium | `f3b5829` |
| A-94 | SEAS benchmark BP1 could not be run | Medium | `df64c96` |

## 3. Behaviour that changed on purpose

- **Stopping and restarting.** Checkpoints are written to a temporary file and renamed (A-65).
  SIGTERM and SIGINT finish the step, write a checkpoint, close the files and exit with 143 or 130
  (A-105). Restarts continue bit for bit, also with the implicit-explicit integrators and the
  power-law time-step bound (A-55, A-58, A-99, A-102, A-106).
- **Errors.** The exit status is nonzero after an error (A-69). Unsupported problem types (A-01),
  compatible SBP operators (O-02), the PI controller (O-04), implicit pore pressure with an
  explicit-only integrator (A-78) and the combined mode with Nz = 1 (O-13) are refused at startup.
  Depth profiles are validated (A-71), and a missing input file is named (A-72).
- **Defaults.** `minDeltaT` (A-95); dissolution-precipitation creep is wet everywhere without a
  `pl_wetDist` file (A-44); `he_Lrad` is 10 km (A-09); the pore-pressure solver takes `linSolver =
  AMG` or `MUMPSLU` (A-101).
- **Startup notes** explain behaviour left as it was: top and bottom `remoteLoading` stay at their
  initial displacement (O-03), `sNVals` is the total normal stress with hydraulic coupling (O-07),
  and explicit slip-dependent permeability is slow in events (O-01).
- **Removed.** The stale `tests/` directory (A-82) and unreachable debug drivers with hard-coded
  paths (A-86).
- **Build.** The Makefile links with PETSc's variables and generates dependencies; `WERROR=1` is
  opt-in, and the tree builds without warnings with clang and GCC (A-80, A-81, A-103, A-104).

## 4. New capabilities

**Two faults** (stages 1 to 4, `docs/TWO_FAULT_DESIGN.md`):

- Named faults: a fault `<name>` reads `<name>_<key>` overrides and has its own files and HDF5
  groups; `vCreep` per fault (`8d78596`).
- Interior faults in a full domain: `interiorFaults`, `<name>_y`, a jump-corrected right-hand side,
  stresses near the fault from u without the jump, second-order fault traction (`c602b99`,
  `44312d7`). The linear-elastic and power-law quasi-dynamic mediators hold a list of faults
  (`64c6ea0`, `e1a40ae`, `18bb957`, shared helpers in `source/multiFault.hpp`, `e6e260b`).
- Two-fault outputs: the surface velocity `/momBal/surfVel`, `faultSeries.txt` and the fault
  positions (`cd2fac5`).
- Frictional heat of interior faults, one Gaussian per fault normalized on the grid (`ae2591c`).
- A pore-pressure equation per fault (`b25167d`).
- `movingBase`: a bottom boundary that moves with the plates (`38a2437`).

**Reversible-strength mechanisms** (stage 5, `docs/REVERSIBLE_STRENGTH_PLAN.md`; power-law
quasi-dynamic mediator):

- Bulk state fields and strain hardening with recovery, `hard_` (`e884fae`).
- Water content of the shear zone, `water_`, with wet-dry dislocation creep (`a2b91d2`).
- A cataclastic grain-size sink and a Zener pinning cap (`672b999`).
- Prescribed and pulsed basal recharge of pore fluid, at the base or at depth (`d6a9556`).
- Evolving fault cohesion with reseal hardening (`c2f556e`).
- Fabric and cement states (`986f39c`); directional fabric in Hill's form (`7b43bd0`, `1d04ecd`).
- Phase segregation, a memory of pressure solution, `seg_` (`48984ca`).
- Pseudotachylite products on grain-size-sensitive creep, `pt_` (`61a5019`).
- Root-strength probes in `faultSeries.txt` (`6e48e56`).

**Other:** the SEAS benchmark BP1 runs, with initial conditions from `createICs.py` (`df64c96`).

## 5. Tools, examples and documentation

- `tools/regress.sh` (bit-for-bit regression on ex1, ex2, ex4s, ex4g), `tools/mms.in` (convergence
  test), `tools/checkkeys.py` (input keys no component reads), `tools/compare_runs.py` (builds
  compared by their physics), `tools/batch.sh` (detached, resumable batches), `tools/two_fault.py`
  (events, slip partitioning and surface velocity of two-fault runs), `tools/zero_d_trading.py`
  (the two-root screen of stage 5), `tools/build-mpi-petsc.sh` (Open MPI and PETSc on a Mac).
- Examples: ex1 to ex5 read the keys the code uses (A-88 to A-91); ex4s and ex4g are short
  power-law regression cases; generators in `examples/interior_fault/` and `examples/two_faults/`
  write the interior-fault, two-fault and stage 5 inputs.
- Loaders: `examples/loadFuncs.py` reads current output (A-92); the MATLAB loaders were rewritten
  but not run (A-93).
- Documentation: `CLAUDE.md`, a quickstart in `README.md`, `docs/AUDIT.md`,
  `docs/TWO_FAULT_DESIGN.md`, `docs/REVERSIBLE_STRENGTH_PLAN.md`, `docs/SERVER.md` and this file.

## 6. Known limitations

From `docs/AUDIT.md`, section 5: explicit slip-dependent permeability is slow in quasi-dynamic events
(O-01); compatible SBP operators fail the convergence test and are refused (O-02); top and bottom
remote loading stay fixed (O-03); the PI controller is not implemented (O-04); `steadyStateIts` is
impractically slow for ex4 (O-05); the aging law freezes the state variable where b <= 1e-3 (O-06);
`sNVals` is the total normal stress with hydraulic coupling (O-07); `maxStepCount` counts steps after
a restart (O-08); the MATLAB tools were not run (O-09); about half of `genFuncs.cpp` is unused and the
mediators are copy-pasted (O-10); some Vecs are allocated every step (O-11); constant pore-pressure
fields are written at every output (O-12); the combined mode refuses Nz = 1 (O-13); the off-fault
Green's function cannot resume a crashed run (O-14); `checkkeys.py` cannot see a key read by the
wrong component (O-15).
