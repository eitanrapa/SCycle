# Two-fault extension: design note

Status: Stages 0 to 4 are done, in that order on `master`; the tags `audit-2026-10` and `stage1` to
`stage4` mark where each ended.
`docs/REVERSIBLE_STRENGTH_PLAN.md` plans the next stage, reversible strength mechanisms.
Line numbers refer to upstream commit `74a132f` and will drift; function names are the stable
reference.

## 1. Goal

Determine with forward models whether two neighbouring sub-parallel strike-slip faults that share
the far-field load trade slip, and under what conditions, and produce the matching surface-velocity
predictions (SAF–SJF, CF–AF). Three stages of increasing physics:

1. **Elastic.** Two rate-and-state faults in an elastic half-space. Verify stress transfer and slip
   partitioning, and recover the single-fault limit when one fault is suppressed.
2. **Power-law viscoelastic.** Off-fault power-law creep only. Test whether trading emerges from
   viscoelastic relaxation alone.
3. **Reversible strength change.** Grain-size evolution, fault valving (pore pressure) and shear heating,
   each in turn.

Deliverable: a library of two-fault models spanning the trading regimes, each with a long-term
partitioning history and present-day surface-velocity profiles.

## 2. Geometry: full domain, both faults interior

Today the single fault is the left boundary `y = 0` with the `symmFault` condition `u(0+) = slip/2`,
which is an antisymmetric half-space: the solution satisfies `u(-y) = -u(y)`.

Adding a second fault at `y = +y_f` to that model silently adds its mirror image at `-y_f`, with the
same slip. The model becomes a three-fault system whose long-term rates satisfy `V1 + 2 V2 = vL`,
not `V1 + V2 = vL`. The `rigidFault` option avoids the image but doubles the stiffness of fault 1 and
flattens the surface-velocity profile on one side. Neither represents SAF–SJF.

**Decision:** model the full domain `y ∈ [-L1, L2]` with **both** faults interior and treated by the
same code. Load both far boundaries with remote displacement:

    bcL = -vL t / 2 + bcLShift        bcR = +vL t / 2 + bcRShift

Consequences:
- Unknowns double relative to the half-space. The momentum-balance matrix is still factored once
  (MUMPS Cholesky with a reused preconditioner), so cost per step grows roughly linearly.
- The built-in `sinh` grid stretching (`Domain::setFields`) refines only near `y = 0`. Supply the
  y-coordinates from a file instead (the existing `inputDir + "y"` PETSc-binary hook in
  `Domain::setFields`), refined near both faults. Metrics are computed numerically from `y(q)`
  (`SbpOps_m_varGrid::constructJacobian`), so any smooth monotone `y(q)` works.
- Each fault must sit at the midpoint between two grid lines `i_f` and `i_f + 1`, at least 6 nodes
  from the left and right boundaries (the Neumann closure of `D1` spans 6 columns).
- The left-boundary assert in `StrikeSlip_LinearElastic_qd::checkInput` currently rejects
  `momBal_bcL_qd = remoteLoading`; `d_dt` must update `bcL` like `bcR`.

## 3. Mechanics of an interior fault: jump-corrected right-hand side ("lift")

There is no interface, multi-block or jump support anywhere in the code today. Two routes exist.

### Route B (chosen for stages 1–3): single grid, jump-corrected right-hand side

The momentum-balance operator is

    A = A_y + A_z + A_L + A_R + A_T + A_B
    A_y = H z_r (D_q μ q_y D_q − H_y⁻¹ R_y)       (y part, built in constructDyymu)
    A_z = H y_q (D_r μ r_z D_r − H_z⁻¹ R_z)       (z part, built in constructDzzmu)

plus SAT boundary terms. For fault k with slip `δ_k(z)` define the step field

    U_k = δ_k(z) · step_k(y),   step_k = 1 for grid rows iy > i_f(k), 0 otherwise
    Ũ_k = δ_k(z) for every row (its smooth extension)

and solve

    A u = rhs_BC + Σ_k J_k,       J_k = A U_k − diag(step_k) A Ũ_k  (= A_y U_k)

The identity `J_k = A_y U_k` is exact because `D_q` and the undivided `R_y` stencils annihilate
constants (so `A_y Ũ_k = 0`), the z-operators and the top/bottom SATs act within one y-row (so they
commute with `diag(step_k)`), and the left/right SATs see a constant field when the fault is at least
6 nodes away. `J_k` costs two `MatMult`s with the existing `A`; **`A` and its factorization are
unchanged**, and with `δ_k ≡ 0` the code path is bit-identical to today's.

Then:
- The physical displacement is `u`. The continuous part is `w = u − Σ U_k`.
- Every **y-strain** must be computed from `w`, not `u`: `σ_xy = μ D_y w`
  (`LinearElastic::computeStresses`, `PowerLaw::computeTotalStrains`). The z-strain keeps `u`
  (the jump `∂z δ_k` is physical).
- Fault traction: average `μ D_y w` on rows `i_f` and `i_f + 1` (two row scatters).
- Accuracy: `w` is C¹ but not C² across the fault when slip varies with depth
  (`[w_yy] = −(μ δ')'/μ`). Displacement converges at O(h²); traction at O(h), with error
  ≈ 0.24·h·|[w_yy]|. With `δ` linear in depth the traction is O(h²) (a useful test).
- If O(h) traction is not good enough, lift the kink too ("B+"): add `½ c(z) (y − y_f)²` with
  `c = −(μ δ')'/μ` to `U_k` for rows past the fault; traction becomes O(h²). Implement only if
  the convergence test in stage 2 requires it.
- Rejected alternative: representing fault slip as an eigenstrain in one column of cells through the
  power-law viscous source term. It omits the `R_y` contribution and leaves an O(δ/h) error in the
  fault traction.

Limits of Route B:
- **μ must be continuous across each fault** (no bimaterial contrast).
- **Not for fully dynamic phases** (`_fd` and `_qd_fd` mediators); those need characteristic
  interface conditions.

### Route A (deferred): split-node multi-block SBP-SAT

Blocks concatenated in y keep the global index layout `Ii = iy·Nz + iz`. Needs per-block 1D
operators, interface selectors, symmetric energy-stable interface SATs for `[u] = δ` and
`[μ u_y] = 0`, a fifth data vector in `setRhs`, and a per-block `R` term. References:
Kozdon, Dunham & Nordström (2012); Erickson & Day (2016); Kozdon, Erickson & Wilcox (2021,
hybridized SBP-SAT for multiple faults). Needed only for a material contrast across a fault, fully
dynamic ruptures, or if Route B traction accuracy proves insufficient. Encapsulate Route B behind a
small interface (`InteriorFaultLift`: `setSlip`, `addToRhs`, `strainField`, `traction`) so Route A
can replace it without touching the mediators.

## 4. Code changes by stage

### Stage 0: regression baseline (done in the audit, tag `audit-2026-10`)
- `tools/regress.sh baseline|compare <dir>` runs `examples/ex1.in` and `examples/ex2.in` and
  compares every HDF5 dataset with `h5diff`. Every later stage must keep these bit-identical in
  single-fault mode. A baseline from the end of the audit is in `data/regress-baseline/`.
- `tools/mms.in`: manufactured-solution convergence of the elastic quasi-dynamic solver on
  Ny = Nz = 21, 41, 81 (order 4: u converges at about 3.5, sxy at about 2.5; order 2: about 2).
- `SEAS_benchmarks/BP1`: runnable (`createICs.py` writes the grid and initial conditions; 25 m
  spacing gives 401 x 2401 nodes). This is the reference for gate 2d.
- The partition-dependent coefficient of the 4th-order `R` operator is fixed (`b65b5e0`): 1 and
  2 ranks agree to about 1e-10, so convergence tests may run in parallel.
- `tools/checkkeys.py` lists input keys that no component reads; run it on every new input.
- 94 findings were fixed first (`docs/AUDIT.md`), including ones on the paths stages 2 to 4
  rely on: the IMEX coupling of temperature and pressure, permeability stability through events,
  creeping nodes in fully dynamic phases, viscous dissipation in the heat and grain-size equations,
  and atomic checkpoints for long runs.

### Stage 1: generalize the fault class, still one fault, bit-identical (done)
Done (tag `stage1`; commits `8d78596`, `3967df3`, `64c6ea0`); ex1, ex2 and
a spot suite (ex5 with heat and flash heating, implicit pore pressure through an event, BP1, the
combined quasi-dynamic and dynamic mode) stay bit-identical, and so does a restart.
- **Named faults.** `Fault`, `Fault_qd` and `Fault_fd` take a name, default `"fault"`, which keeps
  everything as before. A fault named `<name>`:
  - reads every unprefixed key, then lets `<name>_<key>` override it; lists are replaced, not
    appended. For example `fault2_aVals`, `fault2_lockedVals`, `fault2_stateLaw`;
  - loads its initial conditions from `<name>_psi`, `<name>_slip`, `<name>_prestress`, ...,
    `<name>_a`, `<name>_b`, `<name>_Dc` in `inputDir`;
  - writes the HDF5 groups `/<name>`, `/<name>_qd` and the file `<name>.txt`;
  - integrates `<name>_slip` and `<name>_psi` (`Fault::_slipKey`, `_psiKey`).
- **Per-fault settings.** `vCreep` (new key, default `vL`) is the slip velocity of creeping nodes,
  in `ComputeVel_qd` and `Fault_fd::d_dt`. The radiation damping is `_etaScale*sqrt(mu*rho)`, with
  default `1/faultTypeScale`; `Fault_qd::setEtaScale(0.5)` gives the interior-fault value.
- **Mediator.** `StrikeSlip_LinearElastic_qd` owns `std::vector<Fault_qd*> _faults`, with `_fault`,
  the boundary fault at y = 0, first. Stress scatter, rates, output and checkpoints loop over the
  faults, each through its own scatter. The boundary fault alone sets `bcL` and carries the
  steady-state solve, the heat source and the pore pressure (`PressureEq::setSlipKey`).
- **Rows.** `Domain::makeRowScatter(iy, scatter)` returns a scatter from a body field to grid row
  `iy`, created once and kept in `_scatters` as `"body2row<iy>"`; row 0 equals `body2L`, row Ny-1
  `body2R`, also in parallel.
- Not done, by design: the power-law and combined-mode classes still hold one fault, and no input
  key yet declares an additional fault. Stage 2 adds the declaration (and teaches
  `tools/checkkeys.py` the prefixed keys), creates interior faults from it, and treats a model
  without a boundary fault (`_fault` then is NULL and the boundary-only couplings must be skipped).

### Stage 2: one interior fault (done)
Done (tag `stage2`; commit `c602b99` and the example `cc9625d`).
- **Input.** `interiorFaults = [name ...]`, `<name>_y = <km>` for each, and
  `momBal_bcL_qd = remoteLoading` (the full domain: bcL = -vL t/2, bcR = +vL t/2). The fault is
  placed midway between the two grid rows around `<name>_y`, at least 6 rows from the y-boundaries.
  A fault named `fault` reads the plain keys and files, so a single interior fault needs no prefixes.
- **Method.** `InteriorFaultLift` implements Route B with the B+ kink lift
  (`interiorFaultKinkLift = 1`, the default; 0 gives the first-order traction).
  `LinearElastic::computeStresses` forms mu D_y u and lets each lift replace it, on the rows
  i_f - 2 .. i_f + 3 whose D_y stencil crosses the fault, with mu D_y w + mu c (y - y_f) past it
  (`correctStress`); the fault traction subtracts that term's own contribution on row i_f + 1.
  c is taken from the operator, c = -(A delta)/(A q) on row i_f with q = (y - y_f)^2/2, so that
  A Ut vanishes on the fault rows. (Both refinements came in stage 3, `44312d7`: the first version
  used w on every row past the fault, which put errors up to 1 MPa into sxy near the remote
  boundaries, and computed c with D1 twice, which made the order-2 traction first order at the free
  surface.)
- **Refused.** Interior faults with a `symmFault` boundary (mirror images) or `rigidFault`, with
  `guessSteadyStateICs = 1`, with heat or pore pressure (stage 4), with `isMMS` or a body force.
  Interior faults take their initial state from `<name>_stateVals`, `<name>_prestressScalar` or
  the files `<name>_psi`, `<name>_prestress` (see `examples/interior_fault/make_inputs.py`).
- **Results.** See gates 2a, 2b and 2d below; 2c is covered by 2a, which compares with the
  half-space solution of the existing, separately verified code.

### Stage 3: two interior faults, elastic (done)
Done (tag `stage3`; `44312d7`, `cd2fac5`, `5bb2fee` and the example and tool).
- **Input.** `interiorFaults = [fault f2]` with `fault_y`, `f2_y`; the second fault inherits the
  plain keys and overrides them with `f2_` keys. `examples/two_faults/make_inputs.py` writes a
  case with a seismic (`--f2 rs`), creeping (`vs`) or locked (`locked`) second fault, or none
  (`--single`) on the same grid: uniform spacing `--h0` within `--band` of each fault, each fault
  exactly midway between two rows, geometric growth elsewhere.
- **Step size.** A `timeIntInds` listing `slip` or `psi` now covers every fault (`<name>_slip`,
  `<name>_psi` added with the same scale). The default `minDeltaT` was the shear-wave time of the
  smallest cell, and steps at the floor are accepted whatever their error; with 100 m cells it
  stopped the first two-fault runs with an overflow of the state variable. It is now 1e-3 of that
  time (audit finding A-95, `5bb2fee`; ex2's peak slip rate changes by 3.4%).
- **Outputs** (section 6): `/momBal/surfVel` (`computeSurfVel`), `faultSeries.txt`
  (`strideSeries`), fault positions in `data_context.h5`; `tools/two_fault.py` turns a run into
  `events.csv`, `partition.csv` and `surfvel.csv`.
- **Results.** Gates 3a, 3b and 3c below all pass.

### Stage 4: viscoelastic and reversible-strength physics (done)
Done (tag `stage4`). Interior faults now work in both quasi-dynamic mediators
with every physics module; the gates are in section 5.
- **4a, power law** (`e1a40ae`, `e6e260b`, `18bb957`). `StrikeSlip_PowerLaw_qd` has the fault
  list, interior faults, remote loading on the left and the stage 3 outputs, through the helpers
  of `source/multiFault.hpp` that the elastic mediator now uses too. `PowerLaw::computeTotalStrains`
  applies each lift's `correctStrain` near its fault, so sxy, sdev, the viscosity and the viscous
  strain rates see u without the jump. The surface velocity adds the viscous source of the viscous
  strain rates. `interiorFaultKinkSource = 1` (off by default) adds the jump of the viscous source
  to the B+ curvature; it did not improve the match with the half-space beyond that model's own
  discretization error, and with it a fault that never slips still perturbs the solution.
- **4b, heat** (`ae2591c`). `FaultWorkKernel` spreads each interior fault's work tau V over a
  Gaussian of width `wVals` centred on it, normalized on the grid so that the body integral equals
  the work exactly (gate 4b-1); `HeatEquation::setFaultHeatSource` takes it as Qfric. Interior
  faults need `wVals > 0`; set `bcLType_trans = Dirichlet` so the far left side matches the right.
  The same normalization fixed the half-space kernel (`3f402f3`, audit A-97).
- **4c, grain size**: no fault-specific code; the gates pass without changes.
- **4d, pore pressure** (`b25167d`). One `PressureEq` per fault, named after it: the default
  fault's keeps the plain keys, `/pressureEq`, `pressure` and `permeability`; that of a fault
  `<name>` reads `<name>_` overrides and uses `/<name>_pressureEq`, `<name>_pressure`,
  `<name>_permeability`. Each fault's effective normal stress follows its own pore pressure. Use
  `hydraulicTimeIntType = implicit`: the explicit form needs |V| dt < 2.8 kL_p, about 1e5 steps
  per event with kL_p = 1 mm.
- Still refused with interior faults: the steady-state initial guess, `steadyStateIts` and the
  steady-state heat solve (they need a boundary fault), `isMMS`, a body force, the `atan_u` top
  boundary of the power law.
- **Regression**: `examples/ex4s.in` (power law, coupled heat, implicit-explicit) and
  `examples/ex4g.in` (power law, grain size, explicit) join ex1 and ex2 in `tools/regress.sh`.
- **Baseline** for the next stage: `examples/two_faults/make_inputs.py --rheology powerlaw`
  (ex4's dislocation creep and geotherm, Lz = 60 km, coupled heat with a 10 m kernel, optional
  grain size). Cost, measured as the reversible-strength plan asks (581 x 172 nodes, `--grainsize`,
  optimized build, 1 process, 200 yr): 3873 steps in 812 s of integration, 0.21 s per step. Of
  these, about 1200 are coseismic per earthquake (a fault faster than 1e-3 m/s) and about 500
  postseismic; interseismic steps run from 8e5 to 3e7 s, about 500 per century, the Maxwell-time
  cap binding only in the decade after an event (it fell to 2e5 s, then recovered to 2.7e7 s).
  One cycle of about a century is thus about 2000 steps and 7 minutes; 100 cycles about 12 hours.
  The run itself: fault ruptured first (the trigger), and fault f2, with the same friction,
  ruptured next, at 116 yr.

## 5. Verification gates

| Stage | Test | Pass criterion |
|---|---|---|
| 1 | ex1, ex2 via `tools/regress.sh` | bit-identical |
| 2a | Prescribed smooth `δ(z)`, μ uniform, fault at `y = 0` of a symmetric full domain, vs the existing `symmFault` half-space run with `bcL = δ/2` | `u` O(h²), traction O(h) over three grids. **Passed**: Gaussian slip, Ny = 32..256 against Ny = 513: `u` rates 3.7, 3.6; traction rates 1.95, 1.89, 1.70 with B+ (0.93-0.99 without) |
| 2b | Same with `δ` linear in depth | traction O(h²). **Passed** in the stronger form of uniform slip: exact to round-off (4e-13 in `u`) |
| 2c | Analytic antiplane screw dislocation with free-surface image | Covered by 2a |
| 2d | Rate-and-state on the interior fault, other boundary far | recurrence interval and peak V within 1% of the half-space BP1-style run. **Passed** with a triggered first event (steady-state start nucleates from round-off) and `minDeltaT = 1e-4`: at Ny_half = 151 onset +0.68%, peak V +0.54%, slip 0.08%; at 301 onset +0.18%, peak V +0.15%, slip 0.02%. Restarts bit-identical; 2 ranks match onset and slip. (With the old default step floor, A-95, peak V differed by 1.8% and 0.29%.) |
| 3a | Fault 2 locked (`lockedVals`) | reproduces 2d. **Passed** bit for bit: with f2 locked and `constantState`, every dataset of `fault` and `/momBal` (surfDisp, surfVel) and the series equal the run without f2 on the same grid, over 1533 steps through an event |
| 3b | Fault 2 velocity-strengthening | long-term `V1 + V2 = vL` from cumulative slip. **Passed**: f2 with a - b = 0.004, 20 km from fault (ex2 friction), 2980 yr, 14 events; the cycle becomes periodic (recurrence 224.3 yr, peak V 7.22 m/s); over the last cycle, mid-interseismic to mid-interseismic, (slip_1 + slip_2)/(vL dt) = 0.9994 to 1.0045 over depth, and within 0.9% over the last 2, 3 and 5 cycles. f2 takes 0.6 to 0.8% (its friction barely weakens at slow rates, while fault's deep part, with a - b up to 0.2, does). Earlier cycles show the initial stress relaxing (sum up to 1.03) |
| 3c | Coseismic stress change on fault 2 from a fault-1 event | matches the 2D screw-dislocation formula at distance `y_f`. **Passed**: against the antiplane solution for the strip (traction-free top and bottom, held sides; cosine series in depth, exact in y) from the computed slip of fault, the change on the locked f2 20 km away (max 1.29 MPa) agrees to 2.6e-5 MPa (0.002%) |
| 4 | Power-law single-fault limit | matches the existing `StrikeSlip_PowerLaw_qd` run; heat budget `∫Q = Σ τ_k V_k` |
| 4a-1 | Interior fault in the power-law full domain against the half-space power-law run (ex4s material, Nz = 161 to 50 km, cold start, to 1e9 s) | **Passed**: at Ny_half = 101 / 201 onset +1.37% / +0.37%, peak V +1.0% / +0.31%, slip 0.17% / 0.055%, slip at 30 km (viscous afterslip) 0.07% / 0.03%, surface velocity 0.039% / 0.008%; fault traction below 15 km within 2e-3 MPa at 201; near-fault viscous strain within 8% only at the fault's bottom end on the traction-free base, < 0.3% beyond 2 km |
| 4a-2 | Second fault locked (power law) | **Passed** bit for bit: all 1D and 2D datasets (u, sxy, gVxy, effVisc) |
| 4b-1 | Heat budget, two slipping interior faults, uniform 0.5 km cells | **Passed**: body integral of Qfric equals the faults' work to 1.1e-14 (w = 10 m) and 2.7e-15 (w = 1 km) |
| 4b-2 | Single-fault limit with coupled heat (power law, w = 1 km) | **Passed**: temperature within 4.3e-5 K (0.16% of the rise) |
| 4b-3 | Second fault locked, coupled heat | **Passed** bit for bit, T and Qfric included |
| 4c | Single-fault limit with grain-size evolution (coupled to diffusion creep); second fault locked | **Passed**: grain size within 0.03% beyond 2 km (0.7% at the bottom corner), mechanics as 4a-1; locked case bit for bit |
| 4d-1 | Single-fault limit with valving pore pressure (implicit, slip-dependent permeability, an initial pressure diffusing out) | **Passed**: onset +0.61%, peak V +0.49%, pore pressure within 0.13%, log10 k within 3.6e-3 during the event |
| 4d-2 | Second fault locked with its own pressure equation | **Passed** bit for bit, /pressureEq included |

Keep μ uniform across both faults (Route B does not handle a modulus jump at a fault). The
partition-dependent `R` coefficient that once made results depend on the rank count is fixed.

## 6. Output for the research goal

Implemented in stage 3 (`cd2fac5`), on by default when interior faults are declared:
- Per fault, per 1D stride: slip, V, tau, tau_QS, psi (`Fault::writeStep`, group `/<name>`).
- Per step (`strideSeries`): `faultSeries.txt` with, per fault, max V, its depth, potency rate and
  potency, so event catalogs and peak slip rates do not depend on `stride1D`.
- Surface velocity (`computeSurfVel`): `/momBal/surfVel` with each 1D output. Because the
  right-hand side is linear in (bcL, bcR, slip of each fault), one extra back-substitution with the
  factored A, A v = rhs(-vL/2, +vL/2) + sum_k J(V_k), gives the instantaneous velocity field; its
  top row is written. No differencing of displacement across adaptive steps. The far field is
  exactly +-vL/2.
- Context: each interior fault's `y`, `iRow` and `dy` as attributes of its group in
  `data_context.h5`.
- `tools/two_fault.py <outputDir>`: event catalog (`events.csv`: onset and end, duration, peak V,
  depth, potency, moment per unit length), slip at a reference depth and each fault's share over a
  trailing window (`partition.csv`), and the time-weighted mean and the latest interseismic surface
  velocity (`surfvel.csv`, outputs where every fault is slower than `--vinter`).

## 7. Bottom boundary: free or moving base (added 2026-10-04)

**The free base and the far field.** With traction-free top and bottom (`momBal_bcB_qd =
freeSurface`, the default) and the load applied only at the sides, equilibrium makes the
depth-integrated shear force `F = ∫ σ_xy dz` the same at every `y`, so the depth-averaged velocity
is a straight line from the faults' depth-averaged slip rate to `vL/2` at the loaded side. The
surface velocity follows that line beyond a few times `Lz/π` from the faults and meets `vL/2` only at
the boundary, with a slope: from below while the faults are locked, from above after an
earthquake. In ex2 (one fault, 100 km by 30 km) the line holds to a few millionths of `vL/2`
beyond 20 km; 5 yr before its earthquake the surface moves at 0.74 `vL/2` at 50 km and 0.95 at
90 km, 5 yr after at up to 1.79 `vL/2`. Making the domain wider only flattens the line toward the
faults' rate; a steady basal traction would change no velocity (the elastic response is linear
in the rates of the loads). The base has to carry load beneath the far field.

**`momBal_bcB_qd = movingBase`** (`MovingBase` in `multiFault.hpp`, both quasi-dynamic
mediators): the base displacement is imposed, `u(y, Lz, t) = (vL t/faultTypeScale) s(y) +
shift(y)` with `s(y) = tanh((y - yc)/w)`, from `-vL/2` to `+vL/2` across a transition centred at
`yc` (`momBal_bcB_center`, km; default the middle of the domain) whose width `W`
(`momBal_bcB_width`, km) holds 90% of the change, `w = W/(2 atanh 0.9)`; `W = 0` is a step. With a
boundary fault the transition is centred on it and a step moves the whole base at `vL/2`. `shift`
is the base row of the initial displacement (steady-state guess, file or zero), kept in the
checkpoint; the surface-velocity solve gets the base's rate. Where a fault meets the base, the base
pins it: interior faults must be locked at their deepest node (the two-fault runs are, below
`--lock-depth`), a boundary fault too unless `W = 0`, when it must slip there; the code checks.
`remoteLoading` keeps its old meaning (the base held at its initial displacement). The other
mediators refuse `movingBase`. Generator: `make_inputs.py --moving-base W [--base-center Y]`.
A step beneath a boundary fault in a power-law model is costly unless the fault creeps at `vL`
where it meets the base: in ex4s the mismatch at that corner shears the hot mantle so hard that
the smallest Maxwell time drops from 3e7 s to 3e4 s and bounds the steps. Give the transition a
width and lock the fault at the base instead (as in the two-fault models).

Verification: off by default, `tools/regress.sh compare` is bit-identical on all 16 datasets;
with the fault creeping at `vL` at every depth the surface moves at `vL/2` to 1e-12; ex2 with a
moving base matches an independent finite-difference solve of the same problem (its own fault
slip rate, Dirichlet base) to 0.014 `vL/2`, while a free-base solve with that slip differs by
0.16; a run stopped after 300 steps and restarted is bit-identical to the uninterrupted one. With
the moving base, ex2's surface moves at 0.944, 0.981 and 0.996 `vL/2` at 50, 70 and 90 km before
its earthquake and has no overshoot after it; an unbounded elastic half-space with the same fault
slip levels off more slowly (a `1/y` tail), since a base pinned at 30 km couples more stiffly.

**Width of the transition.** The base stands for the mantle below the model: far from the plate
boundary it moves with the plates, beneath it the mantle shear zone carries the plate motion.
Observations of that zone at the base of the lithosphere: beneath the central San Andreas fault
the lithosphere-asthenosphere boundary changes character across less than 50 km, read as shear
localized through the whole lithosphere (Ford, Fischer & Lekic 2014, Geology,
doi:10.1130/G35128.1); SKS anisotropy gives each fault of the San Andreas system a shear zone about
40 km wide at the base of the lithosphere, intense shear over 50-100 km beneath the system, and
broadening into horizontal asthenospheric shear below (Bonnin, Tommasi et al. 2012, GJI); beneath
New Zealand the mantle deforms across about 200 km (Molnar et al. 1999, Science), while at the
initiation of the Alpine fault deformation at 25-50 km depth spread across at least 60 km in
strands (Kidder et al. 2021, Geology). For two faults 20 km apart (SAF-SJF, the Marlborough
faults) two such 40-km zones overlap into one about 60 km wide, and the model's base at 60 km lies
just below its lithosphere (the ex4 geotherm reaches 1488 K at 50 km). Estimate: `W = 60 km`,
centred between the faults, with 40 to 100 km as the range to explore (200 km would need a wider
domain: the tanh tail must be negligible at the sides, 1e-4 for 60 km at 100 km from the centre).

**What the moving base changes in the trading problem.** With a free base the two fault columns
carry the same depth-integrated force, the series condition of the screen
(docs/REVERSIBLE_STRENGTH_PLAN.md, section 2). A moving base takes part of that load and fixes
where the deep motion changes from one plate's velocity to the other's: centred between the
faults it favours neither, and how much it pulls toward the symmetric partition depends on how
strongly the mantle above it couples to the lithosphere.

**Coupling through a power-law mantle.** An interseismic deficit near the faults spreads along the
lithosphere over `lambda = sqrt(mu H_e R_b t)`, with `R_b = ∫ dz/eta` through the mantle above the
base, before the base takes it up (a lithosphere coupled to a viscous channel, Elsasser's stress
diffusion). Only where `lambda` is shorter than the distance to the side boundaries does the base,
rather than the sides, carry the far field. In the two-fault power-law model (ex4 geotherm, 1488 K
at 50 km) the mantle at 45-60 km has effective viscosities of 2e18 to 9e18 Pa s 30 km from the
side in test runs at 340-410 yr (`R_b` 4e-15 to 7e-15 m/(Pa s)), so `lambda` is 170-220 km over a
250-yr cycle with `H_e` = 30 km: longer than the 90 km from the faults to the sides of the 200 km
domain. There the moving base changes the surface velocity by about 1% (free and 60 km moving
base, 50-150 yr after the first event), and the sides still set the ramp.

**A wider domain** (1000 km, Ny 651 instead of 581, faults and base transition centred) lets the
base act: 50 yr after the initial earthquake the moving base holds the surface at the plate rate
beyond about 150 km from the faults (0.97 and 1.05 `vL/2` at 150 and 270 km, against 0.47 and 0.65
with a free base). Late in the cycle it does not: 300 yr after, the lithosphere lags the plates over
hundreds of kilometres with either base (0.41 and 0.64 against 0.36 and 0.58), the faults being
locked and the mantle too weak to hold the lithosphere to the base; the profile is narrow early and
broad late, as for an elastic layer over a viscous half-space (Savage & Prescott 1978, JGR). The
wide domain also loads the faults far more slowly: after the initial earthquake `f2` first ruptured
at 1412 yr with a free base, and not by 661 yr with the moving base, against 415 and 407 yr in the
200 km domain. So with the ex4 mantle the base condition cannot give plate-rate far fields late in
the cycle. That needs a mantle that couples (an effective viscosity above about 1e20 Pa s over the
30 km above the base, which also shortens `lambda` below the distance to the sides) or shear zones
that carry the plate motion down to the base beneath the faults: choices about the model's
physics, not its boundary. Test runs and figures: `data/mb_tests/` (`plot_elastic.py`,
`plot_twofault.py`).

**After spin-up (correction, 2026-10-05).** The runs above started from zero stress and covered
300-650 yr. Run to 3000 yr and compared after 2000 yr of spin-up with the first batch's null run
(free base), the moving base alone and the moving base with a mantle shear zone
(`make_inputs.py --mantle-zone 60`: the bands' factor 1000 from 30 to 40 km, 60 km wide, tapering
into the asthenosphere by 45 km) all give late-interseismic surface velocities that match the
elastic half-space with both faults locked to 20 km and each slipping `vL/2` below (Savage &
Burford 1973) within about 75 km of the faults. 5 yr before an earthquake: 0.54, 0.62 and 0.59
`vL/2` at 30 km from the midpoint (half-space 0.60), 0.68, 0.74 and 0.72 at 50 km (0.75). They
leave it near the sides, where the boundary forces `vL/2` at 100 km and the half-space gives 0.86.
With a 20 km locking depth the half-space itself approaches the plate rate only as
`1 - 2D/(pi x)` (0.92 at 150 km, 0.95 at 250 km), so what remains of the "taper" is the domain's
width, not missing physics; velocities above the plate rate near the faults after an earthquake
(up to 1.6-1.9 `vL/2` 5 yr after) are afterslip and viscous relaxation, also physical. The mantle
zone changes nothing; the moving base costs about 40% more steps by 3000 yr (the zone 50%). The
conclusion above, that plate-rate far fields need a stiffer mantle or shear zones to the base,
came from runs that were not spun up and from the wrong criterion (the plate rate within 100 km
of the faults); it is superseded. A wider domain would remove the boundary's pull but needs a much
longer spin-up, its faults loading far more slowly. Figure: `data/mb_tests/plot_late.py`.
