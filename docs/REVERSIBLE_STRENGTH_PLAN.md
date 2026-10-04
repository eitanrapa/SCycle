# Stage 5: reversible strength mechanisms (plan)

Status: planning note, written 2026-10-02 against commit `5bb2fee`; line numbers refer to it and
will drift, function names are the stable reference. It follows `docs/TWO_FAULT_DESIGN.md`, whose
Stages 3 (two elastic interior faults) and 4 (power law, heat, grain size and pore pressure per
fault) are done (tag `stage4`). This stage is in progress on `master`; section 6 marks the items
that are done.

## 1. Goal

Cawood & Dolan (2024, Seismica, doi:10.26443/seismica.v3i2.1165) list 19 weakening (W1-W19) and
11 strengthening (S1-S11) mechanisms that could make neighbouring faults trade slip over millennia,
and keep those that are (1) reversible at constant depth, (2) strengthen a fault while it is fast
and weaken it while it is slow, (3) act over 20-25 m of slip (centuries to millennia) and (4) act at
the brittle-ductile transition (BDT), where the ductile root sets the long-term rate. Four pairs
survive: hydration against dehydration (W2, W3 vs S1, S2), fabric development against hydrothermal
cementation (W11-W13 vs S6, S9), strain hardening against recovery (S7 vs W14), and shear folding
against rotational weakening (S11 vs W17).

Stage 5 builds the mechanisms of that list that fit SCycle's physics and asks, with the two-fault
model, which of them produce alternating fast and slow periods, with what period, and what the
surface-velocity field looks like in each period. The assessment of 2026-10-01 sorted the list:

| Mechanism (paper) | In SCycle after Stage 4 | Stage 5 |
|---|---|---|
| W8 grain-size reduction, S3 static growth | yes (Austin-Evans wattmeter) | run as is; add the cataclastic sink (W6) and a pinning cap (W10) |
| W18, W19 shear heating | yes (heat equation, Arrhenius feedback) | run as is |
| W4 overpressure, hydraulic part of S6 | yes (fault valving, 1D pore pressure) | add pulsed basal recharge; test hydraulic embrittlement |
| S7 strain hardening vs W14 recovery | no | new bulk state field (section 4.1) |
| W2, W3 hydration vs S1, S2 dehydration | static wetness field only | new bulk state field replacing it (4.2) |
| W11-W13 fabric vs S6, S9 cementation | no | two proxy state fields, fabric and cement (4.5), last |
| friction-side healing (paper's introduction allows upper-crustal strength change) | cohesion exists, static | evolving cohesion (4.6), a control |
| W12 fixed viscous anisotropy | no | optional (4.7); its evolution is out of reach |
| W1, W9, S4, S5, S10 reactions; W5 bulk porosity; W7 melt; W15, W16, S8, S11, W17 geometry | no | out of reach in 2D antiplane strain (section 8) |

## 2. What alternation needs: a two-root reduction

Before any 2D work, the long-term behaviour of two faults sharing a plate rate can be reduced to
two shear-zone roots in series: they carry the same far-field shear stress and their slip rates add
to `vL`. With power-law roots of strength `S_i`, `V_i = vL S_i^-n / (S_1^-n + S_2^-n)`. Give each
root a slow strengthening state `h` that grows with strain and relaxes in time, and a fast
weakening state `D` that grows with strain and heals in time,

    S_i = (1 + h_i) exp(-beta D_i)
    dh_i/dt = (V_i/w)/gamma_h - h_i/tau_h            (hardening strain gamma_h, recovery time tau_h)
    dD_i/dt = (V_i/w)/gamma_D - D_i/tau_D            (weakening strain gamma_D, healing time tau_D)

with `w` the shear-zone width. In units of `w gamma_h / vL` the parameters are
`T_h = tau_h vL/(w gamma_h)`, `T_D = tau_D vL/(w gamma_h)`, `R = gamma_h/gamma_D` and `beta`.
A 40-line RK4 script (to become `tools/zero_d_trading.py`, section 6) integrated this system from
the symmetric state with a 0.1 % perturbation (n = 3):

| Feedbacks present | Result |
|---|---|
| strengthen-when-fast only (`beta = 0`, any `T_h`) | equal sharing, perturbation decays |
| weaken-when-fast only | equal sharing if `n beta R T_D / 2 < 1`, permanent lock-in above it |
| both | alternation in a band just above the lock-in threshold; lock-in deeper in it |

Selected cases with `R = 10` (period in units of `T_h`):

| beta | T_h | T_D | n beta R T_D/2 | regime | period / T_h |
|---|---|---|---|---|---|
| 0.2 | 3 .. 30 | 1 | 3.0 | alternation, full swings (f between 0.00 and 1.00) | 7.8, 3.1, 2.4 |
| 0.2 | 3 .. 30 | 3 | 9.0 | lock-in | |
| 0.4 | 3 .. 30 | 0.3 | 1.8 | alternation | 2.2, 1.2, 0.9 |
| 0.4 | 3 .. 30 | 1 | 6.0 | lock-in | |
| 0.8 | 3 .. 30 | 0.1 | 1.2 | alternation, partial swings (f between 0.2 and 0.8) | 0.9, 0.5, 0.3 |
| 0.8 | 10, 30 | 0.3 | 3.6 | alternation | 3.7, 2.9 |

Three consequences organise this plan:

1. **The paper's criterion 2 is necessary but not sufficient.** A mechanism that only strengthens
   the fast fault (any of the four surviving pairs on its own) drives the system to equal sharing,
   not to alternation. Alternation needs a fast weakening feedback as well, strong enough to
   localise slip on one fault, and a slow strengthening memory that later overturns the choice.
   The weakening feedbacks are the ones already in the code (grain-size reduction, shear heating)
   and the paper's W6, W7 and W11-W13; the slow memories are S7/W14, S1/W2 and S6.
2. **Every mechanism is classified by its sign and its time scale** (table in section 4), and the
   2D run matrix pairs one weakening feedback with one strengthening memory.
3. **The period is set by the recovery time of the slow memory** (1 to 8 `T_h` above), so the
   20-25 m criterion constrains `tau_h` and `gamma_h` jointly, not either alone.

The reduction ignores elastic stress transfer, depth structure, the rate-and-state fault and its own
hysteresis (velocity-weakening patches, deepening of the seismogenic zone when a fault is fast); the
2D runs test whether those add routes to alternation that the 0D model lacks. The scan over the
stress exponent and recovery time is in section 2.1.

### 2.1 Alternation band against the stress exponent

With `T_D = 0.3` and `R = 10` fixed, the localisation number `N_loc = n beta R T_D / 2` and the
recovery time `T_h` were scanned for three stress exponents (S sharing, L lock-in, A alternation
with the period in units of `T_h`):

| `N_loc` | n = 2: `T_h` = 3, 10, 30, 100 | n = 3 | n = 5 |
|---|---|---|---|
| 1.2 | S, A 0.9, A 0.6, A 0.4 | S, A 0.7, A 0.4, A 0.3 | S, A 0.6, A 0.3, A 0.2 |
| 1.5 | A 2.4, A 1.2, A 0.8, A 0.7 | A 1.8, A 0.9, A 0.6, A 0.5 | A 1.3, A 0.7, A 0.4, A 0.3 |
| 2 | A 3.6, A 2.0, A 1.5, A 1.3 | A 2.5, A 1.4, A 1.1, A 0.9 | A 1.7, A 0.9, A 0.7, A 0.6 |
| 3 | L, A 4.3, A 3.3, A 2.9 | A 5.0, A 2.7, A 2.2, A 1.9 | A 2.7, A 1.6, A 1.3, A 1.2 |
| 4 | L, L, A 5.9, A 4.9 | L, A 4.5, A 3.5, A 3.1 | A 4.4, A 2.5, A 2.0, A 1.9 |
| 6 | L, L, L, L | L, L, L, A 5.9 | L, A 4.8, A 3.7, A 3.4 |

The two edges have simple estimates. Equal sharing is linearly unstable when `N_loc > 1`
(perturbation growth rate `n beta R/2 - 1/T_D` for the weakening state alone). Lock-in is permanent
when the fully weakened fast root, strength `exp(-beta R T_D)` relative to its hardened baseline,
cannot be overturned by the most it can harden, `1 + T_h`: roughly `beta R T_D > ln(1 + T_h)`, that
is `N_loc > (n/2) ln(1 + T_h)`, which gives 2.1, 3.6, 5.1, 6.9 for n = 3 and `T_h` = 3, 10, 30,
100, in line with the table. The alternation band is therefore

    1 < N_loc < (n/2) ln(1 + T_h),

wider for longer recovery times and larger stress exponents, and the period grows from a fraction
of `T_h` at the lower edge to several `T_h` at the upper edge. Partial swings (the dominant fault
takes 80 % rather than 100 %) occur near the lower edge; full switching near the upper edge. The
script `zero_d_trading.py` reproduces every number here in under a minute.

### 2.2 The laws as implemented (added 2026-10-04)

The states as built (4.1, 4.5) saturate in [0, 1] and enter the flow law as strength factors, so the
reduction above overstates both feedbacks; the first batch of the run matrix (section 6, item 10),
mapped through it, shared slip equally in all six runs. With the implemented laws a root creeps at
`[(1 + beta_f Phi)/H]^n` times the rate of the unaltered law at the same stress,
`H = 1 + a (S - S_ref)`, and in units of the recovery time `tau_r`

    dS/dt   = Xh f (1 - S) - S,                   Xh = e_full tau_r / gamma_h
    dPhi/dt = (Xf f (1 - Phi) - Phi) / r,         Xf = e_full tau_c / gamma_f,   r = tau_c / tau_r

with `e_full = vL/w` the strain rate of a root of width `w` that takes the whole plate rate and `f`
its share (`zero_d_trading.py --laws implemented`, which also reads these keys from an input file).
Three results:

1. **The fabric's gain is bounded.** Equal sharing is unstable when
   `N_f = n beta_f u/((1 + u)(1 + u + beta_f u)) > 1`, `u = Xf/2`. Its largest value over `Xf`,
   `n beta_f/(1 + sqrt(1 + beta_f))^2` at `Xf = 2/sqrt(1 + beta_f)`, exceeds 1 only for `beta_f > 3`
   (n = 3): 0.80 at `beta_f = 2`, 1.15 at 4, 1.50 at 8, 1.61 at 10, tending to `n`. The saturation
   `(1 - Phi)` and the gentler weakening `1/(1 + beta_f Phi)` against `exp(-beta D)` remove most of
   the gain the linear model gives the same parameters.
2. **The memory's push is bounded too.** Above the threshold the fabric locks one root in, and its
   lock-in branch withstands a hardening bias `b = n ln(H_fast/H_slow)` up to a hold `b_c` (largest
   of `n ln[(1 + beta_f Phi(f))/(1 + beta_f Phi(1 - f))] - ln(f/(1 - f))` over `f > 1/2`, at `f*`).
   The memory builds at most `b_max = n ln[H(S_ss(f*))/H(S_ss(1 - f*))]`, `S_ss(f) = Xh f/(1 + Xh f)`,
   below `n ln(1 + a)` and largest for `Xh` of 2 to 3; a saturated memory (`Xh >> 1`, both roots
   near `S = 1`) pushes little.
3. **Alternation when `b_max > b_c`**, given a fast fabric. Over 135 cases (`beta_f` 4-10, `Xf`
   0.5 to 1, `a` 0.5-1, `Xh` 1.5-4) this criterion gives the regime of the two-root equations in 134
   with `r = 0.01` and 127 with `r = 0.03` (the others: `N_f` within 7 % of 1, which shares, and
   `b_max` within 1 % of `b_c`); with `r = 0.1` it fails in half of them.

Selected cases (`r = 0.03`, `Xh = 2.5`; A alternation with the period in `tau_r` and the range of
the share `f`, L lock-in with the dominant share):

| `beta_f` | `Xf` | `N_f` | `b_c` | `a = 0.5` | `a = 0.7` | `a = 1.0` |
|---|---|---|---|---|---|---|
| 4 | 1.0 | 1.14 | 0.15 | A 1.57, 0.10-0.90 | A 1.32, 0.11-0.89 | A 1.12, 0.11-0.89 |
| 5 | 0.5 | 1.20 | 0.16 | A 1.74, 0.13-0.87 | A 1.43, 0.13-0.87 | A 1.21, 0.14-0.86 |
| 6 | 0.76 | 1.35 | 0.46 | A 3.61, 0.03-0.97 | A 2.34, 0.03-0.97 | A 1.78, 0.03-0.97 |
| 8 | 0.5 | 1.48 | 0.60 | L 0.89 | A 3.33, 0.03-0.97 | A 2.22, 0.03-0.97 |
| 10 | 0.5 | 1.60 | 0.87 | L 0.95 | L 0.92 | A 3.36, 0.01-0.99 |

So with a hardened-to-recovered stress ratio of 1.5 (`a = 0.5`) alternation needs `beta_f` of 4 to
6, `beta_f = 8` needs `a` of 0.7 and `beta_f = 10` about 1; the window of `N_f` that alternates is
narrow (1.1 to 1.6) where the linear model's was wide, and the 2D model's own stabilising effects
(the substrate beneath the roots, elastic coupling) may move it up. The first batch mapped to
`N_f` = 0.80 (alt, fab), 0.65 (lockin) and 0.29 (share) with `Xh = 20`: sharing in every run,
as observed. `zero_d_trading.py --laws implemented --scan` prints the full table in about a minute.

## 3. Common infrastructure (one commit, bit-identical)

All new bulk physics follows the pattern of `GrainSizeEvolution`: a scalar body field on the
(Ny x Nz) grid, integrated explicitly as part of `varEx`, driven by the local deviatoric stress,
viscous strain rates and temperature, and fed back into `PowerLaw` as a multiplier. Everything is
off by default; `tools/regress.sh compare` must stay bit-identical.

### 3.1 A base class for bulk state fields

`GrainSizeEvolution` (`source/grainSizeEvolution.cpp`) already has the full life cycle a new field
needs: constructor (lines 9-34: `loadSettings`, `checkInput`, `allocateFields`,
`setMaterialParameters`, then `loadCheckpoint` / `loadCheckpointSS` / `loadFieldsFromFiles`),
depth-profile keys `<prefix><name>Vals/Depths` read with `loadVectorFromInputFile` (66-135),
`initiateIntegrand` adding a deep copy to `varEx` (296-314), `updateFields` copying it back
(316-330), `d_dt` (393-466), `computeMaxTimeStep` (333-391), and symmetric `writeContext`,
`writeStep`, `writeCheckpoint`, `loadCheckpoint`, `loadCheckpointSS` (647-863). Factor this into
a small base class `BulkStateField` (name, key prefix, integrand key, the Vec, I/O, checkpoint) with
virtual `rate()` and `maxTimeStep()`; derive `GrainSizeEvolution` from it without changing its
results, then the new fields. Each derived class owns its keys (`hard_`, `water_`, `fabric_`), its
HDF5 group (`/hard`, ...) and its text context file, so `tools/checkkeys.py` (regex at line 25)
picks the keys up automatically.

Units inside `PowerLaw` (from the code comments): strains are millistrain, strain rates 1e-3/s,
stress MPa, effective viscosity GPa s, temperature K, activation terms `QR` in K, and
`sdev * dgdev` is kW/m^3, the unit of every heat source. State fields are kept dimensionless and
O(1) so that the default error control (scale 1, `L2_absolute`, `odeSolver.cpp:521-535`) treats
them like `psi`; a field of magnitude 1e12 would pin the run at `minDeltaT`.

### 3.2 Multiplier hooks in `PowerLaw`

`PowerLaw::computeViscosity` (`source/powerLaw.cpp:975-1009`) calls each mechanism's
`computeInvEffVisc` (989-993) and sums the inverse viscosities with the cap (996-1002).
`computeDevViscStrainRates` (1199-1259) reads the mechanisms' member Vecs to form `dgVdev_disl`
(1231-1244), and `solveMomentumBalance` orders `computeViscosity` before `computeViscStrainRates`,
so a factor applied to a mechanism's `_invEffVisc` between lines 992 and 996 stays consistent
everywhere. Add:

- `Vec _hardFactor` with `updateHardFactor(const Vec&)` (a clone of `updateGrainSize`, 875-889):
  multiplies `_disl->_invEffVisc` and `_disl2->_invEffVisc`.
- `updateWetDist(const Vec&)`: `_wetDist` (allocated 271, loaded from `pl_wetDist` or set to 1 at
  371-374) becomes the evolving water content; it already multiplies dissolution-precipitation
  creep (`dissolutionPrecipitationCreep.cpp:322`, and 260 in the guess). Add the mix of
  `_disl->_invEffVisc` (dry) and `_disl2->_invEffVisc` (wet; the second law exists for a wet
  end-member, `wDislCreep2`, `powerLaw.cpp:160`, 206-208): by default the pointwise geometric mean
  with weights `1 - chi` and `chi` (the log-linear form of 4.2), optionally the arithmetic sum with
  the same weights; enabled by a new key `wDislWetDry = yes`.
- `Vec _fabricFactor`: scales every mechanism's `_invEffVisc` before the sum (not the total after
  `VecReciprocal`, so the mechanism partition still adds up to `dgVdev`).
- Mirror the three in `guessSteadyStateEffVisc` (1312-1356, lines 1328-1337), and write the active
  fields in `writeStep2D` (1661-1719) and the checkpoint (1721-1795) only when they are on.

### 3.3 Mediator hooks (`StrikeSlip_PowerLaw_qd`)

The grain-size module shows every touch point in `source/strikeSlip_powerLaw_qd.cpp`: keys and
asserts (170-188, 240, 289-291), construction after `PowerLaw` (78-79; the first copy into the
material happens there), `initiateIntegrand` (510-512) and `solveSStau` (1519), both `d_dt`
overloads (`updateFields` and push into the material at 1022-1023 / 1122-1125, the rate at
1030-1038 / 1135-1143), `timeMonitor` (`writeStep` 554, checkpoint 565, time-step bound 578-582),
`writeContext` (909), `writeSS` (1804), the steady-state viscosity loop (1611-1614) and the
destructor (121). A list `std::vector<BulkStateField*> _bulkStates` driven from one loop at each of
these points keeps the mediator diff small. The fully dynamic classes get no new physics: each
module refuses `_D->_momentumBalanceType != "quasidynamic"` in its `checkInput` (`domain.hpp:39`;
message and `assert(0)` as in `strikeSlip_linearElastic_qd.cpp:401-404`), which covers
`_LinearElastic_fd`, `_LinearElastic_qd_fd` and `_PowerLaw_qd_fd` with one check each, since those
classes ignore unknown keys and would otherwise run the new laws silently. Also refuse
`wLinearMaxwell = yes` (the viscosity is never recomputed, 984-986) and `isMMS`.

Order of operations: in the explicit `d_dt` the grain-size block (1030-1038) runs before the fault
block (1042-1049), so fault tractions and slip rates there are one evaluation old. The cataclastic
sink (4.3) and the cementation proxy (4.5) need the current `tau_k V_k`; move the state-field block
below the fault block in both overloads (bit-identical for the current grain-size law, which does
not use fault quantities) and update the five call sites the report lists (1031, 1136, 580,
1037/1142, 1742) plus the `_qd_fd` mirror (1633-1642).

### 3.4 Fault work spread into the body (shared by heat, grain size, cementation)

`HeatEquation::constructMapV` (`source/heatEquation.cpp:493-561`) builds `_MapV`, an
`Ny*Nz x Nz` matrix that copies a size-Nz fault vector into every grid row (511-514; it is
independent of the fault position), and `_Gw = exp(-y^2/(2 w^2)) / (sqrt(2 pi) w)` (547) from the
depth profile `wVals/wDepths` (m, converted to km, 524-528), centred at `y = 0` only.
`computeFrictionalShearHeating` (1586-1618) then forms `Qfric = MapV(tau V) .* Gw` in kW/m^3 when
`w > 0`, or a boundary flux when `w = 0`. Stage 4 must centre one kernel on each fault
(`_lifts[k]->_yFault`, `interiorFaultLift.hpp:50`) and refuse `w = 0` for interior faults. Do that
as a small class `FaultWorkKernel` (per fault: `y_f`, `w(z)`, `Gw_k`; method `spread(tau_k, V_k)`
adding `MapV(tau_k .* V_k) .* Gw_k` into a body field `Qfault`, kW/m^3), owned by the mediator and
handed to the heat equation, the grain-size module (4.3) and the cement state (4.5). Energy must be
partitioned once: the heat source already subtracts the grain-size fraction of dislocation work
(`dgV_sh = dgVdev - f dgVdev_disl`, `strikeSlip_powerLaw_qd.cpp:1176-1181`); the cataclastic
fraction `f_cat Qfault` is likewise removed from `Qfric`, so gate 4 of the design document
(`integral Q = sum tau_k V_k`) still holds.

### 3.5 Time stepping and error control

- Empty `timeIntInds` puts every explicit key under error control with scale 1
  (`odeSolver.cpp:593-597, 609-613`); a non-empty list hides unlisted keys from step control. Each
  module appends its key with a scale through the `PressureEq::addErrorControl` pattern
  (`pressureEq.cpp:1010-1048`), called from the mediator constructor as at line 73, so an input that
  lists `[psi slip]` (ex4, BP1) still controls the new field.
- Linear relaxation terms with short time scales (cohesion destroyed at `C |V|/D_heal` with
  `D_heal` of millimetres, permeability-like) are unstable in explicit RK below `minDeltaT`, which
  the integrators accept regardless (`odeSolverImex.cpp:448-455`). The exact per-step update of
  `PressureEq::relaxPermeability` (`pressureEq.cpp:964-1002`, called from `be` at 1523-1532 with
  `varImo`, `varIm`, `dt`) is the pattern: put the field in `varIm`, update it exactly in the
  mediator's IMEX `d_dt` after the fault loop, push it into the material, and require `RK32_WBE` or
  `RK43_WBE` as the implicit pressure does (257-260). Fields whose rates are slow relative to the
  step (hardening, water, fabric, grain size) stay explicit.
- `timeMonitor` tightens `maxDeltaT` from the modules' `computeMaxTimeStep` one step late
  (`odeSolver.cpp:694-704`); the relaxation times of the new fields are years, so this is harmless,
  but each module still supplies a bound.

### 3.6 Inputs, tools and regression

- Keys are literal `var.compare("hard_...")` strings in the new `.cpp` files, so `checkkeys.py`
  needs no change; keys assembled from runtime strings need a special case as `disl2_` has
  (`tools/checkkeys.py:34-37`).
- `tools/regress.sh` compares only ex1 and ex2, both elastic. Before Stage 5 touches the power-law
  path, add a power-law baseline: `REGRESS_CASES="ex1 ex2 ex4s"` with `ex4s` a short transient
  version of `examples/ex4.in` (no `steadyStateIts`, `maxTime` of a few hundred years, coarse grid)
  plus a grain-size case (`evolveGrainSize = 1`, `grainSizeEvCoupling = coupled`), and store it
  with `tools/regress.sh baseline`. New datasets are written only when a module is on, since an
  extra dataset makes `h5diff -r` differ.
- Each module gets: an analytic unit test (a 1 x 1 or 1 x Nz run with prescribed stress; section
  4), a restart test (`strideChkpt`, kill, restart, bit-identical against the uninterrupted run,
  as in Stage 1), and `./source/main <input> -objects_dump` clean.
- `CLAUDE.md` and `docs/TWO_FAULT_DESIGN.md` get a line per module when it lands.

### 3.7 The sub-grid width of the shear-zone core

Three of the literature reports arrive at the same obstacle from different sides. The paper's
20-25 m of slip per period is a shear strain of `delta/w` in a root of width `w`: 0.45 in 50 m,
0.02 in 1 km. The drying strain of 4.2 (1-3) is reached within a period only if the deforming core
is 10-30 m wide, which is the width of natural ultramylonite cores (the paper's Pofadder example
is 30 m); the hardening strains of 4.1 (0.03-0.6) allow cores of 40 m to 700 m; the comminuted
zone of 4.3 is centimetres. The two-fault grid has 50 m cells next to the faults and the
frictional-heat kernel is 10 m wide in ex4, so a state field evolved from the cell-mean strain
rate sees strains that are too small by the ratio of the cell width to the core width whenever
the core is narrower than a cell, and a cell-mean grain size cannot represent a centimetre vein at
all. The options are:

1. **Resolve the core**: 10 m cells across a 100 m band at BDT depths. The grid generator already
   refines near the faults; the cost is a few times more rows and a smaller `minDeltaT`, and it is
   the honest first check of every other option.
2. **A sub-grid core width `w_sz(z)`** (input profile, 10-100 m): the state ODEs use the core strain
   rate `e_core = e_cell dy / w_sz` where `dy` is the local cell width, and the viscosity seen by the
   cell is the series combination of core and host at equal stress (for a grain-size-sensitive law
   `d_eq^-m = d_host^-m + (w_sz/dy) d_core^-m`; for a multiplier `M` on the strain rate,
   `M_eq = 1 + (w_sz/dy)(M - 1)`). This is the two-component scheme of 4.3 generalised, and it is
   consistent with how `wVals` already treats frictional heat as a sub-grid source.
3. **Accept the cell scale** and report rates as cell-scale parameters. Cheapest, but the
   constants then depend on `dy` and cannot be compared with the laboratory values.

The plan takes option 1 for the verification runs of each module (so that the laws are tested as
written) and option 2 for the library runs, with the comparison between them as a reported
sensitivity. `w_sz` is one profile shared by every module, owned by the mediator next to `wVals`.

## 4. Mechanisms

Classification used in the run matrix (section 5): "memory" is a state that strengthens the fast
fault (the paper's criterion 2) and relaxes when it is slow; "feedback" is a state that weakens the
fast fault further. The 0D model (section 2) says alternation needs one of each.

| Module | State | Sign | Strain scale | Time scale | Section |
|---|---|---|---|---|---|
| strain hardening / recovery | `S` on dislocation creep | memory | `gamma_h` | `tau_r(T)` | 4.1 |
| water content | `chi` on wet/dry creep | memory (dry when fast) | `gamma_dry` | `tau_hyd` | 4.2 |
| grain size (existing) + cataclastic sink | `d` | feedback | wattmeter | static growth | 4.3 |
| Zener pinning cap | parameter `d_Z` | limits the feedback | | | 4.3 |
| shear heating (existing) | `T` | feedback | | thermal diffusion `w^2/kappa` | Stage 4 |
| pulsed basal recharge | forcing `q_b(t)` | neither | | pulse period | 4.4 |
| fabric | `Phi` | feedback | `gamma_f` | `tau_c` | 4.5 |
| cementation | `Psi` | memory, driven by earthquakes | `D_cem` | `tau_Psi`, `gamma_Psi` | 4.5 |
| cohesion healing | `C` on the fault | feedback (slip destroys, time heals); memory in the reseal variant | `D_heal` | `tau_heal` | 4.6 |
| fixed anisotropy | parameter | neither | | | 4.7 |

Each subsection gives the law, the inputs, the code changes, the tests and the expected role.
Stress is `sdev` (MPa), the dislocation strain rate `dgVdev_disl` (1e-3/s), temperature `_T` (K),
all available where `GrainSizeEvolution::d_dt` is called.

### 4.1 Strain hardening and recovery (S7 vs W14)

**Law.** Dolan, Bowman & Sammis (2007) proposed the mechanism in words only: during rapid creep the
ductile root strain-hardens as dislocations tangle, slip shifts to the annealed root of the
neighbouring system, and the hardened root recovers while it is slow. The paper's section 5.3
places it at the BDT because climb is suppressed there, and gives no numbers. The quantitative
templates come from metals and olivine:

- Kocks-Mecking-Estrin: `d rho/d gamma = k1 sqrt(rho) - k2 rho` with the Taylor stress
  `tau = tau_0 + alpha mu b sqrt(rho)`, saturating at `sqrt(rho_s) = k1/k2`; static recovery by
  climb adds `d rho/dt = -K_s(T) rho^2` (Mecking & Kocks 1981; Kocks & Mecking 2003; olivine:
  Karato, Rubie & Yan 1993). Strain to saturation in metals is 0.03-0.1 at tens of MPa.
- Back stress (kinematic hardening): `d sigma_b/d eps = gamma (sigma_b,max - sigma_b)` with the
  strain rate driven by `sigma - sigma_b` (Hansen et al. 2019 for olivine low-temperature
  plasticity, `gamma = 75`, hardening strain about 0.04; Hansen et al. 2021 at 1250-1300 C,
  `gamma` 5-55, steady back stress 48-61 % of the applied stress; Wallis et al. 2021 see the same
  long-range back stresses after high-temperature creep). Hansen et al. state that recovery terms
  are the missing piece at elevated temperature.
- Dislocation theory of transient creep (Breithaupt, Katz, Hansen et al. 2023): storage
  `(e/b)(m sqrt(rho) + n/d)` against grain-boundary and pipe-diffusion recovery; after a stress
  step the viscosity drops by about two orders of magnitude and re-equilibrates over months to
  years, with transient and steady-state viscosities of different stress dependence, unlike a
  Burgers body.
- Nonlinear Burgers transients used for postseismic modelling (Masuti et al. 2016; Masuti &
  Barbot 2021) saturate at strains of 0.3-0.8 % at 100 MPa with transient-to-steady viscosity
  ratios of 0.17-0.67 (Chopra 1997), too small a strain scale for a 20-25 m period.
- Quartz and feldspar at BDT conditions: Hirth & Tullis (1992) regime 1 (natural bulging
  recrystallization, about 280-400 C) is climb-limited and strain-hardens, with recovery only by
  bulging recrystallization that creates strain-free grains, exactly the paper's W14; recovery in
  quartz needs water-related defects and pressure (Tullis & Yund 1989: days at 800 C and 1.5 GPa
  wet, nothing in four days dry); kick-and-creep experiments show tangled dislocations at
  1e12-1e13 /m^2 after high-stress loading and recrystallization during relaxation (Trepmann &
  Stöckhert 2013); Wallis, Lloyd & Hansen (2018) document strain hardening by geometrically
  necessary dislocations driving marbles from creep to friction at the BDT. No transient or
  recovery law exists for quartz or feldspar at 300-450 C and 10-200 MPa; extrapolating olivine
  or hot-quartz recovery kinetics to 350-450 C gives recovery times anywhere between 10 and 1e7 yr.

So the hardening amplitude, strain scale and recovery time are free parameters to be reported as
such. Hardened-to-recovered stress ratios of 1.1-1.7 at fixed strain rate are what olivine and
Burgers data support. No published fault-system model uses hardening and recovery to switch slip;
the closest are damage rheologies with healing, which switch between localized and distributed
deformation on the healing time (Ben-Zion et al. 1999; Lyakhovsky, Ben-Zion & Agnon 2001), and
viscoelastic coupling of two segments, where a 10-30 % strength contrast gives bimodal recurrence
(Lynch et al. 2003).

**Form implemented.** One scalar `S` in [0, 1] per node (0 recovered, 1 fully hardened), as the
hardening report recommends:

    dS/dt = (dgVdev_disl / gamma_h) (1 - S) - S / tau_r(T),
    tau_r(T) = tau_r0 exp( (Q_r/R) (1/T - 1/T_ref) ),
    strain rate = A exp(-QR/T) (sdev / H(S))^n,   H(S) = 1 + a (S - S_ref),

that is a factor `H^-n` on `_disl->_invEffVisc` (and `_disl2`). In the notation of section 2 this
is the memory state with `h = a S`, the same `gamma_h`, and `tau_h = tau_r`. Steady state at strain
rate `e`: `S_ss = x/(1 + x)` with the Damköhler number `x = e tau_r / gamma_h`. The calibrated flow law is
recovered to round-off at `S = S_ref`, chosen as `S_ss` at the reference strain rate `e_ref` of the
single-fault steady state (key `hard_eRef`, or `hard_SRef` directly). Analytic solutions, the unit
tests: at constant `e`, `S(t) = S_ss + (S_0 - S_ss) exp(-t/tau_eff)` with `tau_eff = tau_r/(1 + x)`,
in strain `S(gamma) = S_ss + (S_0 - S_ss) exp(-(gamma/gamma_h)(1 + 1/x))`; at zero strain rate
`S = S_0 exp(-t/tau_r)`; the stress at constant strain rate is the power-law stress times `H(S(t))`.
Two side effects to document: at steady state the law has an apparent stress exponent
`1/n_app = 1/n + a x/((1 + x)^2 H)` (about 2.9 for `n = 4`, `a = 0.5`, `x = 1`), so calibrations
must be compared at `S_ref`; and in the bulging regime recovery is recrystallization, which the
grain-size module also represents, so runs that activate both must watch for double counting.

**Parameter window** (from the sweep in the hardening report): `a` 0.3-0.6; during the fast period
`x_f = e_f tau_r/gamma_h` of 1-3 and during the slow period `x_s` below 0.3; `gamma_h` no larger
than the strain accumulated in one fast period, `delta/w`, 0.45, 0.11 and 0.022 for 22.5 m of slip
over 50 m, 200 m and 1 km. Examples: `gamma_h = 0.05`, `tau_r` 1-3 kyr, strain rates 1e-13 to
1e-11 /s give a strength ratio of 1.37-1.39 between the periods and a hardening e-folding time of
140-150 yr; `gamma_h = 0.1`, `tau_r = 3 kyr`, 1e-14 to 1e-12 /s give 1.24 and 1.5 kyr. The rule of
thumb is `tau_r` of 1-3 times `gamma_h/e_f`, 0.2-5 kyr, the observed period itself, which is what
a memory must have. In a 50 m root the strain per period far exceeds any plausible `gamma_h`, so
the state saturates early in the period; the mechanism is most natural for roots of 200 m to 1 km.

**Keys.** `hard_type = off | transient | constant`, `hard_a`, `hard_gammaHVals/Depths`,
`hard_tauR0`, `hard_QR` (K), `hard_TRef`, `hard_SRef` or `hard_eRef`, initial `hard_SVals/Depths`
or a file `hard_S` in `inputDir`.

**Code.** New `HardeningState : BulkStateField` (rate from `dgVdev_disl` and `T`);
`PowerLaw::_hardFactor` applied at `computeViscosity` 992-996 and `guessSteadyStateEffVisc`
1328-1337; mediator loop (3.3). Requires `wDislCreep = yes` (otherwise `dgVdev_disl` is zero,
`powerLaw.cpp:1231-1244`). A back-stress form, `A exp(-QR/T) (sdev - sigma_b)^n` with
`d sigma_b/dt = H_b e (1 - sigma_b/sigma_b,max) - sigma_b/tau_r`, is the more physical alternative
(it allows the observed reverse strain after a stress drop) but needs a stress-dependent
renormalisation of `A`; keep it for a later iteration.

**Tests.** (a) Unit: a 2 x 2 grid with `effViscCap` large, prescribed constant stress through the
boundary conditions, `S(t)` against the two analytic solutions to 1e-10, and `S_ss`, `n_app` from a
long constant-rate run. (b) Single fault: `S` frozen at `S_ref` is bit-identical to Stage 4's run;
`S` free changes the recurrence interval by less than 1 % when `tau_r` is short compared with the
interseismic period and measurably when it is long. (c) Two faults: the sign test of section 2,
equal sharing when the module runs alone (`grainSizeEvCoupling = no`, `thermalCoupling = no`);
then the pairing runs of section 5.

**Role.** Memory. Within one shear zone it is a stabilising, delocalising feedback (the hardening
report's and section 2's conclusion): with `tau_r` short it is only a rate-strengthening steady
law, with `tau_r` infinite a one-shot transfer; alternation needs `tau_r` comparable to the fast
period and a weakening feedback alongside it.

### 4.2 Water content (W2, W3 vs S1, S2)

**Law.** Quartz dislocation creep in equilibrium with water scales as `A sigma^n f_H2O^r exp(-Q/RT)`
with `r` near 1 (Hirth, Teyssier & Dunlap 2001: `A = 10^-11.2 MPa^-n-r/s`, `n = 4`, `Q = 135
kJ/mol`; Tokle, Hirth & Behr 2019: `n = 4` at low stress and 2.7 at high stress and low
temperature, `Q` 125 and 115 kJ/mol, `r` 1.0-1.2; Fukuda, Holyoke & Kronenberg 2018: `r = 1.0`)
or 0.49 in the natural-plus-experimental fit of Lusk, Platt & Platt (2021: `A = 10^-9.3`,
`n = 3.5`, `Q = 118 kJ/mol`), which is 10-50 times faster than the laboratory laws at BDT conditions.
Laboratory dry quartzite is 1.5-2.5 times stronger than the same rock with water (Koch et al. 1989:
wet `A = 5.05e-6 MPa^-n/s`, `n = 2.61`, `Q = 145 kJ/mol`; dry `A = 1.16e-7`, `n = 2.72`,
`Q = 134`), while vacuum-dried quartz has `Q` near 300 kJ/mol (Kronenberg & Tullis 1984) and
extrapolates to absurd strengths at 350 C, so laboratory-dry laws are not the dry end-member.
At 50 MPa the Koch laws give a dry-to-wet viscosity ratio of 3.4 at 350 C and 4.0 at 400 C; the
fugacity route (`r = 1`, wet fugacity 30-400 MPa against a dry value of 1-10 MPa) gives 10-300,
and `r = 0.49` gives 3-17. The design range is therefore a ratio of 3-30 with 10 as default, and the
spread between laws at these conditions (20-50 between Hirth et al. and Lusk et al.) is as large as
the wet-dry effect itself. Weakening is controlled by the presence of water rather than its
concentration: natural recrystallized grains with only 20-100 H per 1e6 Si crept at low stress
because grain-boundary water was present (Kilian et al. 2016; Kronenberg et al. 2020), and Hirth
& Beeler (2015) argue that the water fugacity at the BDT is set by the pore pressure, which
supports tying the supply to the fault's `p`.

Water loss with recrystallization is quantified once: in the Moine thrust footwall the
fluid-inclusion water falls from 4080 to 1570 ppm (molar H per 1e6 Si) as the recrystallized
fraction rises from 17 to 97 %, roughly linearly, with the strongest localization where the rock
is driest (Kronenberg et al. 2020), so fully recrystallized grains keep about a third of the
protolith water; a kilometre-thick ultramylonite holds half the water of its weakly deformed
host (Finch, Weinberg & Hunter 2016); grain-boundary-migration recrystallization leaves 20-100
H per 1e6 Si (Kilian et al. 2016). Full recrystallization of sheared quartzite needs shear
strains up to 8 (Heilbronner & Tullis 2006), so a drying strain scale of 1-3 for subgrain-rotation
and migration recrystallization, plausibly 3-10 for bulging at 300-400 C, is defensible; no
source gives water against strain directly. Uptake is not diffusive: interstitial hydrogen moves a
millimetre in a thousand years at 400 C but does not weaken, molecular water needs 1e5 yr per
100 um (Kronenberg et al. 1986; Gleason & DeSisto 2008), and microcracks are the entry path,
healing within days at 400 C (Smith & Evans 1984), so re-wetting is event-limited: a cracking
episode delivers inclusions at once, and a hydration time of 1-6 kyr follows from re-wetting
5-10 % of the deficit per event every 100-300 yr. Two findings cut the other way and make the
sign of the strain-rate term a hypothesis rather than a fact: recrystallization can add water to
initially dry quartz (Palazzin et al. 2018), and creep cavitation pumps fluid into the fastest
layers (Précigout et al. 2017, 2019). For pressure solution the rate is proportional to the
lumped product of grain-boundary diffusivity, film thickness and solubility (Gratier, Dysthe &
Renard 2013; Pluymakers & Spiers 2015), which is `B D c` in SCycle's law, so a linear scaling by
a wetted-boundary fraction is the mean-field limit, with an optional percolation threshold of
0.1-0.3. The closest published models are a hydration front feeding viscosity in a
one-dimensional shear zone (Kaatz et al. 2023) and a fracture-gated reaction-progress state
(Gueydan, Leroy & Jolivet 2004); no earthquake-sequence model evolves the wetness of the viscous
root.

**Form implemented.** `chi` in [0, 1] is the shear-zone wetness, with a floor `chi_min` for the
water that recrystallization does not remove:

    dchi/dt = (chi_sup - chi) / tau_hyd  -  (dgVdev_disl / gamma_dry) (chi - chi_min),

the same structure as the aging law. Steady state `chi_ss = chi_min + (chi_sup - chi_min)/(1 + x)`
with `x = tau_hyd e / gamma_dry`; the crossover strain rate `gamma_dry/tau_hyd` must sit between the
slow-period and fast-period strain rates of the root, which gives the design rule
`tau_hyd = gamma_dry w / V_plate`. Analytic solutions, the unit tests: at constant `e`,
`chi(t) = chi_ss + (chi_0 - chi_ss) exp(-t/tau_eff)` with `1/tau_eff = 1/tau_hyd + e/gamma_dry`;
at zero strain rate `chi` relaxes to `chi_sup` on `tau_hyd`; with `tau_hyd` infinite and any
history `e(t)`, `chi - chi_min = (chi_0 - chi_min) exp(-gamma(t)/gamma_dry)` depends only on the
accumulated strain, a test with a time-varying rate. `chi_sup` is a depth profile, a function of
the fault pore pressure at that depth (`min(1, f_H2O(p, T)/f_ref)`, after Stage 4 from
`PressureEq::_p` scattered along the fault-adjacent column) or pulsed in time (4.4); an optional
dry-grain enhancement of uptake, `1/tau_hyd` multiplied by `1 + kappa (1 - chi)`, encodes the
microcracking argument and keeps the equation bounded.

Coupling: the `disl2_` law is the wet end-member and `disl_` the dry one. The default mix is
log-linear, `1/eta_disl = (1/eta_wet)^chi (1/eta_dry)^(1 - chi)`, which is the effective fugacity
`f_dry (f_wet/f_dry)^chi` when the two laws differ only in `r` and gives a viscosity that varies
smoothly between the end-members; the arithmetic mix `chi/eta_wet + (1 - chi)/eta_dry`, which
SCycle's series structure suggests, keeps the wet law dominant until `chi` falls below the inverse
viscosity ratio and is offered as an option. Pressure solution scales with `chi` as today, or with
`(chi - chi_c)/(1 - chi_c)` above a threshold.

**Values to explore.** `chi_sup = 1` or pressure-linked; `chi_min = 0.3`; `gamma_dry = 2` (1-10);
`tau_hyd = 1 kyr` (1-6 kyr); dry-to-wet viscosity ratio 10 (3-30); `chi_c = 0`. Drying time
`gamma_dry w / V`: 290-1000 yr for a 10 m core at 10-35 mm/yr, 0.9-3 kyr for 30 m, 3-10 kyr for
100 m. Slip per fast period is `gamma_dry w`, so matching 20-25 m needs a core of 10-25 m for
`gamma_dry` of 1-2: this module needs the sub-grid width of 3.7 (or 10 m cells across the band),
otherwise the drying cannot complete within a period in a 50 m to 1 km cell. Stiffness is mild
(the drying time is at least decades even at postseismic strain rates of 1e-9 /s), so `chi` is an
explicit integrand like `grainSize`.

**Keys.** `water_type`, `water_chiSupVals/Depths`, `water_chiMin`, `water_tauHyd`, `water_kappa`,
`water_gammaDry`, `water_chiC`, `water_mix = log | arithmetic`, `water_supplyFrom = profile |
pressure | pulses`, initial `water_chiVals/Depths` or the existing file `pl_wetDist`;
`wDislWetDry = yes` in `PowerLaw`.

**Code.** New `WaterState : BulkStateField`; `PowerLaw::updateWetDist`, the mix of
`_disl->_invEffVisc` and `_disl2->_invEffVisc` at `computeViscosity` 991-992 (pointwise geometric
mean for the log-linear form) and `guessSteadyStateEffVisc` 1328-1337; `_wetDist` is already in
the context, checkpoint and SS checkpoint (1602, 1756, 1825, 1898), add it to `writeStep2D` when
evolving; clamp `chi` to [`chi_min`, 1] after each `updateFields`. Note the dry end-member must be a
fugacity-scaled wet law or the Koch dry law, not a vacuum-dried one.

**Tests.** (a) Unit: the three analytic solutions above. (b) `chi = 1` everywhere reproduces the
current wet run bit-identically (weights 1 and 0); `chi = chi_min` with `chi_min = 0` the dry run.
(c) Single fault: the steady-state `chi_ss` profile against the formula; two faults: sign test
(equal sharing), then pairings.

**Role.** Memory, if expulsion is strain-driven and supply is slow; the water report's linear
analysis of two load-sharing zones gives decaying antisymmetric modes, the same conclusion as
section 2. The published wet-dry contrast bounds the strength ratio between periods to the
design range above, which is comparable to what hardening gives (4.1).

### 4.3 Cataclastic grain-size sink (W6, grain-size side of W7) and Zener pinning (W10)

**Law.** The existing wattmeter (Austin & Evans 2007) is an energy balance: the rate of change of
`1/d` equals the fraction `f` of dissipated work that becomes grain-boundary energy, divided by
`c gamma`. Coseismic comminution is the same balance with frictional work as the source, so one
event leaves `1/d_after = 1/d_before + (f_cat/(c gamma)) integral(tau V) / w_cat` over a comminuted
zone of width `w_cat`: the product depends on the deposited energy density, not on the slip rate.
Measured partitions: the Punchbowl ultracataclasite stores fracture-surface energy worth a few
tenths of a percent of the frictional work (Chester, Chester & Kronenberg 2005); Chi-Chi gouge
surface energy is 6 % of the breakdown work, itself a few percent of the frictional work (Ma et al.
2006); Gole Larghe pseudotachylytes put 97-99 % of the work into heat (Pittarello et al. 2008);
Wilson et al. (2005) find nanometre fines after a single rupture. With a surface energy of 1 J/m^2
and `c = pi`, a 1 cm zone comminuted to 1 um by 1 m of slip at 100 MPa needs `f_cat` of 3e-4, a
10 cm zone to 0.1 um 3e-2; the defensible range of `f_cat` as a fraction of total frictional work is
1e-3 to 3e-2, central 1e-2, the same order as the dynamic-recrystallization fraction `f` fitted to
the quartz piezometer (0.010-0.025, Tokle & Hirth 2021; 0.1 in Austin & Evans 2007 for olivine and
calcite). Products: pseudotachylyte veins of millimetres to centimetres, ultracataclasite of
centimetres to decimetres, damage over 0.1-1 m, BDT ultramylonite bands of centimetres to metres
(Price et al. 2012; Hawemann et al. 2019; Campbell & Menegon 2019); ruptures penetrate up to about
4 km below the seismogenic zone (Scholz 1988), with aftershock depths deepening by 3 km for four
years after Landers (Rolandone, Bürgmann & Nadeau 2004) and a transient 1-3 km deflection of the
BDT (Ellis & Stöckhert 2004).

Static growth, `d^p - d_0^p = k_g f_H2O^r exp(-Q/RT) t`, is the uncertain half. For quartz the
recommended law (Tokle & Hirth 2021, refitting Tullis & Yund 1982 and Fukuda et al. 2019) has
`p = 3`, `k_g = 0.261 um^3/s/MPa^1.38`, `r = 1.38`, `Q = 134 kJ/mol`; the Fukuda et al. (2019)
novaculite law has `p = 2.9`, `k_g = 10^-5.8`, `r = 1.9`, `Q = 60 kJ/mol`. With water fugacities of
50-140 MPa at 350-450 C the time for 1 um to grow to 10 um is 90 kyr, 6.5 kyr and 650 yr at 350,
400 and 450 C under the first law, and 1 kyr, 150 yr and 30 yr under the second; to 50 um a
hundred times longer. So at 400-450 C fine coseismic products can coarsen to 10 um within one
recurrence interval under one law and not under the other, at 350 C under neither, and dry
conditions give no growth at all: the two laws must both be run. In SCycle units (d in m, `p = 3`)
the first law is `A = 0.261e-18 f_H2O^1.38 m^3/s` (6e-17 to 2.4e-16 m^3/s over that temperature
range) with `QR = 16100 K`; the second `A = 10^-23.2 f_H2O^1.9 m^2.9/s`, `QR = 7200 K`. The
wattmeter constants for quartz are `c = pi`, `gamma = 1 J/m^2`, `f = 0.015`. Zener pinning gives
a limiting grain size `d_Z = K r_p / phi_p^m` (Smith-Zener; Herwegh et al. 2011; Linckens et al.
2011) with `K` 0.7-1.3 and `m` 0.5-1; for quartz-feldspar-mica ultramylonites with second-phase
radii of 1-10 um at 10-40 % the cap is 3-30 um, while quartz-rich mylonites have no practical cap.
Bercovici & Ricard (2012) multiply growth by a pinning factor `1 - (d/d_Z)^2`. Piezometers for the
initial grain size: Stipp & Tullis 2003 `d = 10^3.56 sigma^-1.26` (um, MPa; stresses scaled by
0.73 after Holyoke & Kronenberg 2010), Cross et al. 2017 `d = 10^3.91 sigma^-1.41`, about 12 um at
100 MPa.

**Form implemented.** `GrainSizeEvolution::d_dt` (393-466) becomes

    dd/dt = growth * max(0, 1 - (d/d_Z)^q)  -  (f/(c gamma)) d^2 (sdev dgdev_disl)  -  (f_cat/(c gamma)) d^2 Qfault,

with `Qfault` the fault work spread by the kernels of 3.4 (kW/m^3, the same unit as
`sdev * dgdev`, so no conversion), `f_cat` a depth profile, `d_Z` the pinned grain size (default
infinite, which reproduces the current law exactly) and `q` 1 or 2. `computeMaxTimeStep` (361-363)
gets the extra work term; `computeSteadyStateGrainSize` (555-569) keeps its closed form only
without the cap and with a sink still proportional to `d^2`, so with a finite `d_Z` it becomes
`min(d_ss, d_Z)` or a bracketed root (`rootFinder.hpp` is already included).

The sub-grid problem of 3.7 is sharpest here. The frictional-heat kernel width `w` is 10 m in ex4
and the cells next to the faults are 50 m, while the comminuted zone is centimetres: the cell-mean
energy density is a hundred to a thousand times too small, and a cell-mean grain size would stay
at hundreds of micrometres while the physical zone is at 1 um. For fault-parallel shear the thin
zone and its host deform in series at equal stress, so the cell behaves, for a grain-size-sensitive
law with exponent `m`, as if it had `d_eq^-m = d_host^-m + (w_cat/dy) d_f^-m`, about
`d_f (dy/w_cat)^(1/m)`, ten times `d_f` for a ratio of a thousand and `m = 3`. The recommended
form therefore keeps two grain sizes near each fault: the body field `d_host` as today, and a
size-Nz fault-zone grain size `d_f` per fault, evolved with the comminution energy density
`tau V / w_cat(z)` (input `grainSizeEv_wCatVals/Depths`, 0.01-1 m), the same static growth and
cap, and combined by the series rule of 3.7 into the `d` that the creep laws see in the cells
within the kernel width. Implement the cell-mean form (one line in `d_dt`) first as the limit
`w_cat = dy`, then the two-component form; the difference between them is itself a result.

**Keys.** `grainSizeEv_fCatVals/Depths` (default 0), `grainSizeEv_wCatVals/Depths` (default: cell
width), `grainSizeEv_dZVals/Depths` (default none), `grainSizeEv_dZExp` (1 or 2), a floor
`grainSizeEv_dMin` (0.1 um; the nanometre fines are a small mass fraction).

**Code.** `d_dt` gains a `const Vec& Qfault` argument (a zero Vec when no kernel exists); the
mediator evaluates the kernels after the fault block (3.3); the heat source subtracts
`f_cat Qfault` (3.4); the two-component form adds a per-fault `Vec _dFault` to the module with its
own checkpoint entries and a `/grainSizeEv/<name>` dataset; write `fCat`, `wCat` and `dZ` into
`/grainSizeEv` and the checkpoint only when given.

**Tests.** (a) Unit, no growth and constant `Qfault`: `1/d(t) = 1/d_0 + (f_cat/(c gamma)) Qfault t`;
growth only: `d^p = d_0^p + A exp(-QR/T) t`, checked against the annealing times above; with the
cap `d -> d_Z`. (b) `f_cat = 0`, no `d_Z`: bit-identical to the Stage 4 grain-size run. (c) Energy:
`integral (Qfric + f_cat Qfault)` equals `sum tau_k V_k` over an event to 1e-6. (d) Series rule:
two cells of prescribed grain size in series against the closed form. (e) Two faults: grain size
reduced at the deep tip of each rupture, annealing between events, with the paper's statement
that this weakens the fast fault further (lock-in) as the expected outcome.

**Role.** Feedback with a long memory (annealing), and the one module whose rates are set by
published energy partitions rather than by choice. The pinning cap limits how far static growth
can restrengthen an abandoned root and is a one-parameter sensitivity, not a mechanism.

### 4.4 Pulsed basal recharge and hydraulic embrittlement (W4, S2 at the fault)

**Law.** The module is the fault-valve model of Zhu, Allison, Dunham & Yang (2020): along-fault
flow `n beta dp/dt = d/dz[(k/eta)(dp/dz - rho g)]`, permeability `k = k_min + (k* - k_min)
exp(-(sigma - p)/sigma*)` with `dk*/dt = -(V/L)(k* - k_max) - (k* - k_min)/T`, `p = 0` at the
surface and a constant upward flux `q_0` at the base, "a crude approximation for a fluid source at
depth". Their featured parameters: `n = 0.01`, `beta = 1e-9 /Pa`, `eta = 1e-4 Pa s`,
`k_min = 1e-19 m^2`, `k_max = 1e-15 m^2`, `L = 1 m`, `sigma* = 30 MPa`, `T = 1e8 s` (about 3 yr;
0.3 to 300 yr explored), `q_0 = 3e-9 m/s` (3e-10 for the short-`T` run). Findings that matter here:

- The controlling ratio is `T` over the recurrence interval. `T` comparable to or longer than the
  recurrence gives little valving (the fault stays permeable); `T` much shorter gives reduced
  overpressure cycling but quasi-periodic slow-slip events at the base of the seismogenic zone
  (Ozawa, Yang & Dunham 2024 derive the instability of steady slip with steady flow above a critical
  flux). The cyclic strength change is 10-20 MPa.
- **No multi-event pressure cycle exists in these models**: every large rupture resets `k*` to
  `k_max` within `L` of slip and drains the column, so the pressure cycle is at most one earthquake
  cycle long. Multi-event accumulation would need incomplete coseismic enhancement (`L` much larger
  than the coseismic slip, or a smaller `k_max`) or a supply that varies on longer time scales. The
  related models (Yang & Dunham 2021; Yang, Yehya, Iwalewa & Rice 2021; Dal Zilio, Hegyi, Behr &
  Gerya 2022; Petrini et al. 2020; Perez-Silva et al. 2023; Farge, Jaupart & Shapiro 2021, where
  the input flux selects quiescent, clustered or periodic bursts of a valve chain) show no two
  parallel faults trading slip under fluid control.
- Steady state of the coupled laws at the plate rate: `k*_ss = (k_min + k_max V T/L)/(1 + V T/L)`,
  and the permeability needed to pass `q_0` at a near-lithostatic gradient fixes
  `sigma - p = sigma* ln[(k*_ss - k_min)/(k_need - k_min)]`, about 50 MPa for their numbers: a
  closed-form check of the Stage 4 single-fault valving run.
- Fluid supply from depth: the regional metamorphic flux is about 1e-12 m/s (Connolly 2010), so
  `q_0` of 3e-9 m/s is flux focused into the fault zone from a wide catchment. Episodic supply is
  documented at 1 kyr pulses every 6-8 kyr (Louis et al. 2019), 10-15 kyr cyclicity (Vignaroli et
  al. 2022), flow lasting under 15-45 kyr (Hickey et al. 2014), and porosity waves with spacing of
  a compaction length (Connolly 1997, 2010; periods of tens of kyr). Menegon & Fagereng (2021)
  describe episodic fluid escape at the base of the seismogenic zone driven by tectonic pressure
  gradients in the creeping root.
- At quasi-steady flux the fault passes the flux by adjusting `k`, so the response to a pulse of
  amplitude `A` is logarithmic, `Delta(sigma - p) = -sigma* ln(1 + A)`: with `sigma* = 30 MPa`,
  `A = 9` removes about 70 MPa and the valve must open. A pulse cannot be stored (the column holds
  about 0.02 m^3/m^2 per 100 MPa against about 100 m^3/m^2 delivered per kyr at 3e-9 m/s); it is
  transmitted.
- Hydraulic embrittlement: a strength-envelope estimate (`tau_f` proportional to `(1 - lambda) z`,
  `tau_v` Arrhenius) deepens the frictional-viscous intersection by about 3 km for `lambda` from 0.4
  to 0.8 and 5 km to 0.9 with a wet quartzite law, 1-2 km with feldspar laws. Hirth & Beeler (2015)
  caution that the effective-stress law loses efficiency near the BDT as real contact area grows.
  No cycle model that moves the BDT by changing `lambda` was found; Allison & Dunham (2018) cap the
  effective normal stress at about 50 MPa (Rice 1992) and keep it fixed.

**Form implemented.** The bottom flux of `PressureEq` is today `q_b = (1 + bcB_ratio) rho_f^2 g k/eta`,
set once in `setUpSBP` (`source/pressureEq.cpp:661-668`) and refreshed when permeability changes
(`updateBoundaryCoefficient`, 583-629); it therefore scales with the bottom-node permeability.
Make the imposed part a prescribed flux in m/s, independent of `k`, and let it vary in time,

    q_b(t) = q_0 (1 + sum_i A_i P((t - t_i)/tau_i)),

with `P` a raised cosine `(1 + cos(pi s))/2` for `|s| <= 1` (or a Gaussian; square pulses only
for unit tests), amplitude `A_i`, centre `t_i`, half-duration `tau_i`, and an optional period for
repetition. The factor is applied to the imposed part just before `setRhs` in the two explicit rate
functions (1313, 1392) and in the backward-Euler solve (1619, after the refresh at 1599 inside the
pressure-dependent loop), never in the initial steady-state solve (831), so the start is still a
steady state; a pulse that is nonzero at `initTime` is refused. `_bcB` enters only the right-hand
side (the matrix at 1629-1659 depends on `dt` and the coefficient), so no factorization changes. In
`be` the time argument is `t_{n+1}` (`odeSolverImex.cpp:462-470`), the backward-Euler-consistent
choice. The pulse is a pure function of absolute time, so restarts need no new checkpoint data.
Write `q_b(t)` as a size-1 timestepped dataset `/pressureEq/bcB` and record the pulse keys in
`p_context.txt`, which today holds only `g` and the time-integration type (1985-1986).

Where the pulse enters matters: a transient imposed at a base hundreds of kilometres down diffuses
with `D = k/(eta n beta)` of order 0.01 m^2/s at the steady permeability and is smeared over 1e5 yr
before it reaches the seismogenic zone, while Zhu et al.'s constant `q_0` is insensitive to the
base depth. Either keep the pressure column short (the Stage 4 two-fault baseline needs `Lz` of
60-100 km for the viscous roots, so the pulse boundary at 25-30 km is not available) or add a
depth-Gaussian source term `Q(z, t)` centred at 15-25 km to the pressure equation (a new term in
`dp_dt` and `be`, keyed `bcB_sourceDepth`, `bcB_sourceWidth`), which is the recommended form: the
paper's pulses are delivered to the BDT, not to the base of the lithosphere.

**Values to explore.** Background `q_0` 3e-10 to 3e-9 m/s; amplitude 3-10 (porosity waves need at
least three times the background; conserving the time-averaged flux with Louis et al.'s duty cycle
gives 5-8), up to 30; half-duration 0.5-1 kyr (Louis) and 5-7 kyr (Vignaroli); recurrence 6-8 kyr,
10-15 kyr and 30-50 kyr; ramp time at least `max(T, H^2/4D at k_max)` to keep the implicit step
well conditioned. A 1 kyr pulse spans 10-20 recurrence intervals, which is the paper's fast period
(clusters of 4-6 events, 20-25 m of slip); the hypothesis to test is whether two faults sharing one
source alternate because the pulse is captured by whichever fault is more permeable.

**Keys.** `bcB_q0` (m/s; replaces the use of `bcB_ratio` when given), `bcB_pulseShape = cosine |
gaussian | square`, `bcB_pulseAmp = [...]`, `bcB_pulseCentre = [...]`, `bcB_pulseHalfDur = [...]`,
`bcB_pulsePeriod = [...]` (0 = once), optional `bcB_sourceDepth`, `bcB_sourceWidth` (km); equal
lengths, positive durations.

**Pitfalls.** A square edge inside an adaptive step corrupts the embedded error estimate: clamp
`maxDeltaT` to the next edge in `timeMonitor`, or use `tanh` ramps. Backward Euler with steps longer
than the pulse under-resolves it: use the step-averaged flux. The background flux uses the last
entries of the depth lists (665-666) while the refresh uses the bottom node (583-629); the two agree
only when the lists end at the node value (pre-existing). Explicit pressure is conditionally stable:
use `hydraulicTimeIntType = implicit` and `RK*_WBE`. Stage 4's per-fault `PressureEq` must remove
the collisions the report lists: unprefixed keys, the `/pressureEq` group, the `pressure` and
`permeability` integrand keys, and the refusal of pore pressure with interior faults
(`strikeSlip_linearElastic_qd.cpp:401-404`), by giving `PressureEq` a name and prefix as `Fault`
has (`fault.cpp:12-13, 69-79`).

**Tests.** (a) Constant permeability, a step `Delta q` in the basal flux at `t = 0` on a column of
height `H` with `p = 0` at the top: the excess pressure is

    u(z, t) = (eta Delta q / k) { z - (8H/pi^2) sum_m (-1)^m/(2m+1)^2 sin((2m+1) pi z / 2H) exp(-(2m+1)^2 pi^2 D t / 4H^2) },

`D = k/(eta n beta)`, so the base rises to `eta Delta q H / k` with the slowest mode decaying on
`4H^2/(pi^2 D)`, half of the rise at `0.79 H^2/(4D)`, 95 % at `4.5 H^2/(4D)`, and early on
`u(H, t) = (2 eta Delta q / k) sqrt(D t / pi)`. With `n beta = 1e-11 /Pa` and `H = 15 km` the time
`H^2/(4D)` is 18 kyr at `k = 1e-19 m^2`, 180 yr at 1e-17, 1.8 yr at 1e-15; the steady basal rise
for `Delta q = 3e-10 m/s` is 4.5 MPa at `k = 1e-16 m^2` and 45 MPa at 1e-17, which is why a sealed
column must valve. The valving steady state of the Stage 4 run is checked against
`sigma - p = sigma* ln[(k*_ss - k_min)/(k_need - k_min)]` above. (b) No pulse keys: bit-identical.
(c) Single fault, valving on, `kT_p` long: the recurrence modulation over several events; then a
pulse train with a period of 1-10 kyr. (d) Two faults, a pulse on one: does the pressure rise at
12-20 km move frictional slip below the locked zone (the effective normal stress falls and the
fault's frictional strength drops under the viscous flow stress), and does the other fault slow
down. Note that with the input `a - b > 0` below about 17 km, embrittlement there converts viscous
strain into aseismic frictional slip and narrows the surface-velocity profile rather than
extending seismic rupture; a run that wants deeper seismic slip must also move the `a - b`
transition.

**Role.** A forcing, not a mechanism. It tests the paper's dismissal of W4 (cycles tied to the
seismic cycle) and whether an imposed hydration pulse on one root is enough to switch faults.

### 4.5 Fabric development against cementation (W11-W13 vs S6, S9), proxy

**Law.** The paper treats foliation, viscous anisotropy and interconnected weak layers (W11-W13)
as unidirectional products of strain that only cementation or veining can undo, and cementation
(S6, S9) as something deep ruptures inject during fast periods until, after 20-25 m of slip, the
root is strong enough to hand over. Quantitative templates:

- Strain softening with healing: damage rheology (Lyakhovsky, Ben-Zion & Agnon 1997; Lyakhovsky &
  Ben-Zion 2008) grows a damage variable with the strain-energy rate and heals it with a
  logarithmic, Dieterich-like law (healing fast at first, damage above 0.5 persisting for
  centuries); in fault-network simulations a long healing time gives regular characteristic
  events, a short one a disordered network, and an intermediate one alternating periods of
  intense activity and quiescence (Lyakhovsky, Ben-Zion & Agnon 2001; Ben-Zion et al. 1999),
  the closest established analogue of slip trading. Geodynamic codes heal accumulated strain as
  `d eps/dt = e - H(T) eps` with Arrhenius healing (Fuchs & Becker 2019; the ASPECT
  `strain_dependent` rheology), while classical tectonic strain softening (Huismans & Beaumont
  2003; Lavier, Buck & Poliakov 2000) has no healing at all.
- Fabric strain scale: aligned weak inclusions give a fivefold traction drop by a shear strain
  of 5 after hardening up to 1 (Dabrowski, Schmid & Podladchikov 2012); CPO anisotropy saturates
  near 10 (Hansen, Zimmerman & Kohlstedt 2012); interconnection of weak phases begins at low
  strain (Holyoke & Tullis 2006; Montési 2007, 2013). So `gamma_f` of 1-10 and a weakening factor
  of 1-10 (the Pofadder strain-rate contrast of a thousand between ultramylonite and mylonite is a
  strength ratio of 6-10 for `n` of 3-4, partly grain size, so an upper bound).
- Cementation rates span six orders of magnitude. Diagenetic quartz precipitation (Walderhaug
  1994, 1996; Lander, Larese & Bonnell 2008) would need 1e8 yr to seal a millimetre crack and is
  irrelevant; Williams & Fagereng (2022) find that the fast mechanisms (coseismic pressure drops,
  frictional heating, dissolution-precipitation) can cement micrometre-to-millimetre slip surfaces
  in days to years but cannot build mesoscale vein networks interseismically; hydrothermal
  experiments on faulted sandstone regain up to 35 MPa of cohesion in six hours at 927 C with an
  activation energy near 70 kJ/mol (Tenthorey & Cox 2006; Tenthorey, Cox & Todd 2003), which
  scales to 0.2 yr at 400 C and 270 yr at 150 C as an upper bound on rate; crack sealing by
  pressure solution takes years to a thousand years for sub-millimetre cracks at 100-300 C
  (Renard, Gratier & Jamtveit 2000). The field anchor is reseal hardening: dilation breccias
  resealed between earthquakes became stronger than intact rock and later rupture moved elsewhere
  (Woodcock, Dickson & Tarasewicz 2007). The cementation time is therefore a scanned parameter,
  10^2 to 10^3 yr in fast periods, with Arrhenius dependence (`Q` 50-70 kJ/mol).

**Sign analysis, which changes the form.** The obvious proxy, fabric `Phi` growing with viscous
strain and destroyed by cementation driven by the overlying fault's slip, is degenerate: its two
rates are both proportional to the slip rate, so the long-term steady state
`Phi_ss = (1/(w gamma_f)) / (1/(w gamma_f) + tau/E_c)` does not depend on how fast the fault is,
and the only term that breaks the tie, time healing `1/tau_c`, removes fabric from the slow fault
and leaves the fast one weaker. That proxy is a bounded weakening feedback (the classical
localisation feedback of Montési 2013 and Dabrowski et al. 2012), never the paper's memory. The
paper's narrative needs strength that grows with activity beyond the pristine state, as in the
reseal-hardening field case. Two forms do that: a cement state that strengthens, or a gate that
stops fabric growth above a critical strain rate (the paper's argument that fast strain rates give
recrystallization and hardening rather than planar fabrics). The continuous form is preferred.

**Form implemented.** Two states near each fault, both in [0, 1]:

    dPhi/dt = (dgVdev / gamma_f) (1 - Phi)  -  Phi / tau_c(T)                       (fabric, feedback)
    dPsi/dt = (Qfault / E_c) (1 - Psi)  -  Psi ( 1/tau_Psi + dgVdev / gamma_Psi )     (cement, memory)
    viscosity factor on every mechanism:  (1 + beta_c Psi) / (1 + beta_f Phi).

Fabric grows with total viscous strain on `gamma_f` and heals in time on `tau_c`; cement grows
with the frictional work of deep ruptures spread by the kernels of 3.4 (`Qfault`, so a slip of
`delta` at traction `tau` adds `tau delta / E_c`; with `E_c = tau D_cem` the hand-over scale
`D_cem` of 10-30 m is the input) and is removed by renewed deformation on `gamma_Psi` and by
dissolution in time on `tau_Psi`. With constant inputs both equations are linear with the
exponential solutions of 4.1 (unit tests). The gated alternative replaces the fabric source by
`(dgVdev/gamma_f)(1 - Phi) H(e_crit - dgVdev)`.

**Keys.** `fabric_type`, `fabric_betaF`, `fabric_gammaFVals/Depths`, `fabric_tauC0`, `fabric_QC`,
`fabric_TRef`, `fabric_eCrit` (0 = no gate), `cement_type`, `cement_betaC`, `cement_DCem`,
`cement_tauPsi`, `cement_gammaPsi`, initial `fabric_PhiVals/Depths`, `cement_PsiVals/Depths`.

**Values to explore.** `gamma_f = 3` (1-10), `beta_f = 4` (1-10), `tau_c` 10^2 to 10^3 yr;
`D_cem` 10-30 m, `beta_c` 0.2-1, `gamma_Psi` 1-10, `tau_Psi` 10^3 to 10^4 yr. The strain-scale
problem of 3.7 applies: at the ultramylonite rate of 3e-10 /s a fabric strain of 3 takes 300 yr,
at 1e-12 /s in a kilometre cell 1e5 yr.

**Code.** New `FabricState` and `CementState : BulkStateField` (inputs `dgVdev`, `Qfault`, `T`);
`PowerLaw::_fabricFactor` applied to each mechanism's `_invEffVisc` before the sum (3.2) and in the
steady-state guess; `Qfault` from the kernels after the fault block (3.3).

**Tests.** Unit solutions; `Phi = 0`, `Psi = 0` frozen bit-identical; two-fault sign tests: fabric
alone must reproduce the lock-in threshold of section 2 in the 2D model (a quantitative check of
the reduction itself); cement alone must give equal sharing; the pair maps onto the 0D band.

**Role.** The fabric state is the tunable weakening feedback, the cement state the only memory in
this plan that is driven by earthquakes rather than by strain, so it is the form in which the
paper's hydrothermal narrative can be tested at all. Both are the least constrained modules (the
paper itself notes that ore-deposit cementation takes far longer than the observed periods), so
they are built last.

### 4.6 Evolving cohesion on the fault (friction-side control)

**Law.** Hydrothermal healing of quartz gouge has a solution-transfer component beyond the
Dieterich log-time healing: friction gains 0.010-0.014 per decade of hold with a cutoff time that
shortens with temperature (activation energy 54 kJ/mol) and is removed over slip of about 500 um
rather than the usual 10 um (Nakatani & Scholz 2004); at 65 C healing accelerates beyond
log-linear for holds longer than a thousand seconds and follows a pressure-solution contact-growth
model (Yasuhara, Marone & Elsworth 2005); healing rates rise with temperature up to lithification
(Karner, Marone & Evans 1997); restrengthening of hydrothermal quartz gouge is largely independent
of normal stress, that is cohesive, and promotes instability (Muhuri et al. 2003); faulted
sandstone regains up to 35 MPa of cohesion in hours at 927 C with `Q` near 70 kJ/mol (Tenthorey &
Cox 2006); in halite gouge friction and cohesion grow on different clocks, cohesion of several MPa
becomes significant under near-lithostatic fluid pressure, and the slip that removes it is much
shorter than `Dc` (van den Ende & Niemeijer 2019); the Chen-Niemeijer-Spiers microphysical model
derives `a`, `b`, `Dc` and the accelerating healing from pressure-solution compaction (Chen,
Verberne & Spiers 2015; Chen, Niemeijer & Spiers 2017; van den Ende et al. 2018 compare it with
rate-and-state in cycle simulations). Plausible cohesion at 5-15 km: 0.1-1 MPa for uncemented
gouge, 1-10 MPa for hydrothermally cemented gouge (the recommended `C_max`), 20-35 MPa as the
upper bound of fully relithified fault rock, against effective normal stresses of 50-150 MPa.
The aging law already supplies Dieterich healing of order `b sigma_n ln(t/t_0)`, about 20 MPa
over 300 yr for `b = 0.01` and 100 MPa, so `C` must stand for the non-Dieterich part only, or it
double counts. Pressure-solution lithification as the stick-slip mechanism goes back to Angevine,
Turcotte & Furnish (1982); Sleep & Blanpied (1992) coupled it to fluid pressure. No rate-and-state
cycle model with evolving cohesion or healing-rate variation producing supercycles was found.

**Form implemented.** Cohesion `C` (MPa) already adds to the frictional strength
(`ComputeVel_qd::getResid`, `source/fault.cpp:1383, 1401`; the velocity solve treats
`|tauQS| <= C` as stuck, 1330) and is a static depth profile (`cohesionVals/Depths`, 105-106). Make
it a per-fault integrand `<name>_cohesion`:

    dC/dt = (C_max - C)/tau_heal - C |V| / D_heal,

healing toward `C_max` in time, destroyed by slip over `D_heal`. Steady state at slip rate `V`:
`C_ss = (C_max/tau_heal) / (1/tau_heal + V/D_heal)`, which is also the consistent initial value
under `guessSteadyStateICs` (`Fault::guessSS` adds cohesion at 784). The rate is zero on locked and
creeping nodes (`lockedVals`, 1322-1326), which never consult cohesion.

**Code.** `Fault`: `_cohesionKey = _prefix + "cohesion"` (constructor 12-13), keys in
`parseSetting` before line 142 (`cohesionEvolution`, `cohesionMaxVals/Depths`, `cohesionTauHeal`,
`cohesionDheal`; the `<name>_` override is automatic, 69-79), constraints and the
`momentumBalanceType` refusal in `checkInput` (244-292; `Fault_fd` shares it, so the leapfrog
classes are refused with one check), `Fault_qd::initiateIntegrand` (909-917 pattern),
`updateFields` (937-938), `d_dt` after 1016, `writeStep` (634-641, only when evolving),
`writeContext` (576-580, 598); checkpoints already carry `_cohesion` (676, 1263, 1157).
`StrikeSlip_LinearElastic_qd`: add `cohesion` to the `timeIntInds` expansion (90-92) and the
automatic scale (`max C_max`). `StrikeSlip_PowerLaw_qd`: mirror the zeroing at 1072-1073 and
1155-1156, and repair the block at 1065-1069, which overwrites `_strength` without cohesion after
`_fault->d_dt` (so `tauP` written and used by the heat source lacks it; pre-existing). Stiffness:
`C |V|/D_heal` with `D_heal` of millimetres is explicit-unstable in events; use the IMEX exact
update of 3.5 (`C = C_inf + (C_old - C_inf) exp(-a dt)`, `a = 1/tau_heal + |V|/D_heal`) when
`D_heal < 0.1 m`, else the explicit key.

**Values and a caveat.** `C_max` 1-10 MPa, `tau_heal` 10-1000 yr with Arrhenius dependence (`Q`
50-70 kJ/mol), `D_heal` 1e-3 to 1e-1 m. The steady state `C_ss(V) = C_max / (1 + V/V_c)` with
`V_c = D_heal/tau_heal` is an extra velocity-weakening term of slope up to `C_max/4` per e-fold of
`V` at `V = V_c`, 2.5 MPa for `C_max = 10 MPa` against a typical `(b - a) sigma_n` of 0.4 MPa: put
`V_c` far below creep rates (`tau_heal = 300 yr`, `D_heal = 1 mm` gives 1e-13 m/s), so that `C`
is zero on creeping sections and `C_max` on locked ones, and include `C` in the nucleation-size
estimate when choosing the grid.

**Tests.** (a) Spring slider (`Nz = 1`, ex1 style): the stress drop grows by
`C_max (1 - exp(-T_r/tau_heal))` with the recurrence `T_r`; compare the recurrence against the
single-degree-of-freedom prediction. (b) Keys absent: bit-identical. (c) Two faults, elastic:
the sign test.

**Role.** Feedback (slip destroys cohesion, time restores it), the friction-side analogue of fault
maturation (W15), but reversible. In the 0D classification it promotes lock-in; it is the cheapest
control for the elastic two-fault model and needs no power law. The friction-side memory would be
reseal hardening, strength that grows with activity: `dC_max/dt = (|V|/D_rh)(C_lim - C_max)
- (C_max - C_min)/tau_long`, so repeated ruptures raise the ceiling on the fast fault until it
hands over; it is a two-line addition to the same module and the only way the paper's narrative
can be tested in the elastic model, so it is included as an option (`cohesionReseal = 1`).

### 4.7 Fixed viscous anisotropy (W12), optional

**Law.** The transversely isotropic viscous rheology of Mühlhaus, Moresi, Hobbs & Dufour (2002) and
Lev & Hager (2008), with the director normal to the foliation, reduces in antiplane strain with a
fault-parallel foliation to two uncoupled scalars: `sxy = 2 eta_S exy` on the foliation plane (easy)
and `sxz = 2 eta_N exz` across it (hard), no off-diagonal terms. Measured and inferred ratios:
Lev & Hager use `eta_S/eta_N = 0.1` as conservative; olivine CPO gives a shear viscosity fifteen
times smaller than the normal one, saturating at a shear strain of 10 (Hansen, Zimmerman &
Kohlstedt 2012); Mameri et al. (2019, 2020) find a factor of 2 in strength at equal strain rate and
implement it as a Hill yield function; layered-composite theory gives a normal-to-shear ratio of
2.3-3 for a quartz-mica viscosity contrast of 10 and 5.5-8 for 30 at 20-50 % mica (Christensen
1987; Treagus 2003). For power-law flow the Hill form is `s_e^2 = sxy^2 + r sxz^2`,
`exy = A s_e^(n-1) sxy`, `exz = A s_e^(n-1) r sxz` with `0 < r <= 1`; the viscosity ratio is
`r^-(n+1)/2` at equal stress and `m = r^-(n+1)/(2n)` at equal strain rate, so a strength ratio of
2, 5 or 10 at equal strain rate is `r` of 0.34, 0.082 or 0.028 for `n = 3.5`. Recommended `m` of
2-10 for micaceous mylonite, default 5, specified at equal strain rate.

**Form implemented.** With `r = m^(-2n/(n+1))` from the input `m`, replace the deviatoric invariant
in `computeSDev` (1166-1195) by `s_e` and multiply the `xz` rate in `computeViscStrainRates`
(1028-1036) by `r`; every mechanism's `computeInvEffVisc` then receives `s_e`. `r = 1` is the
isotropic code to round-off (the regression test); pure `xy` and pure `xz` loading give the two
viscosity ratios above and mixed loading at `sxy = sxz` checks `s_e` (unit tests). The power-law
MMS needs the same change. Fixed anisotropy does not evolve, so it changes the baseline strength
of a root and its localisation, not the trading dynamics; build it only if a pairing run needs a
weaker fault-parallel direction to localise. The evolution of the fabric orientation (the paper's
W12, S8) remains out of reach.

## 5. Run matrix and diagnostics

**Baseline (Stage 4 deliverable).** The two-fault power-law model: the grid of
`examples/two_faults/make_inputs.py` (faults 20 km apart in a 200 km domain, 50 m cells within 4 km
of each fault, friction of ex2, 50 MPa normal stress) extended in depth to hold the viscous roots
(ex4 uses the LAB at 50 km, `Nz_lab`, and a thermal profile from 283 K at the surface), with the
ex4 dislocation law (`A = 1585 MPa^-n s^-1`, `n = 3`, `Q/R = 41570 K`), `wVals` of 10 m (ex4) to
1 km for the frictional heat kernel, `thermalCoupling = coupled`, `RK43_WBE`, and the
Austin-Evans grain size from a piezometer start (quartz constants of 4.3). Each root's width then
emerges from the rheology and the grid; the state-field rates scale with the local strain rate, so
the strain per metre of slip, `1/w`, is a resolution-dependent quantity that every run must report
(diagnostics below). The verification runs of each module use 10 m cells across a 100 m band at
BDT depths; the library uses 50 m cells with the sub-grid core width of 3.7, after gate 5.0b has
shown the two agree.

**Screen first.** For each module, run `tools/zero_d_trading.py` with the module's `(gamma, tau)`
pair and the feedback it is paired with (grain size, shear heating or the fabric proxy) mapped to
`(beta, R, T_D)`, and pick three to five parameter sets spanning sharing, alternation and lock-in.
Only those go to 2D.

**2D runs, in order.**

1. Stage 4 physics alone (grain size; heating; both): the null experiment. Expected lock-in or
   sharing, never alternation; record which, and the strength that a memory must overcome.
2. Each memory alone (4.1, 4.2): expected equal sharing, with a measurable transient.
3. Pairings: 4.1 with grain size, 4.1 with heating, 4.2 with grain size, 4.1 with the fabric state
   of 4.5 as a tunable feedback, and the cement state of 4.5 with grain size (the one pairing whose
   memory is driven by earthquakes). For each, the three to five screened parameter sets.
4. Controls: 4.6 on the elastic model; 4.4 pulses on one fault of the Stage 4 model.
5. The library: for the pairings that alternate, a sweep of fault separation (10, 20, 40 km), normal
   stress contrast and the second fault's frictional parameters, each run long enough for five
   switches.

**Diagnostics.** `tools/two_fault.py` already writes the event catalogue, the trailing-window slip
partition at a reference depth and the interseismic surface-velocity profiles. Add, per fault and
per `strideSeries` step, the column average of each active state field and of `effVisc` within
`w` of the fault at a reference BDT depth (new columns in `faultSeries.txt`), so the strength
history of each root is available without 2D output; and a `switches.csv` from the partition
series (times at which the dominant fault changes, period, slip per period). The 20-25 m criterion
is then read directly as slip per period.

**Cost.** Stage 4 must measure the step count per cycle and the time per step of the baseline;
until then assume several thousand steps per cycle and hours per run in the optimized build on a
few ranks, so the library is tens of runs, not hundreds, and the 0D screen is what keeps the
sweeps small. Long runs rely on the atomic checkpoints and `restartFromChkpt = 1`.

## 6. Order of work

On `master` after Stage 4 (tag `stage4`); one commit per item, each with its
tests in the message, each leaving ex1, ex2 and the power-law baseline bit-identical unless the
message says why not.

1. `tools/zero_d_trading.py` (the script of section 2 with a command line for the parameters and a
   CSV of regimes) and this plan. Small. **Done**: it reproduces every entry of the tables of
   sections 2 and 2.1 (the scan in 33 s).
2. Power-law regression baseline (`ex4s`, grain-size case), stored under `data/regress-baseline/`.
   Small. **Done in Stage 4** (`c44b308`): `examples/ex4s.in` and `examples/ex4g.in`, default cases
   of `tools/regress.sh`.
3. Infrastructure (3.1-3.5): `BulkStateField`, `GrainSizeEvolution` derived from it, the three
   `PowerLaw` hooks, `FaultWorkKernel` shared with the heat equation, the mediator loop, the
   refusals in the fully dynamic classes, state-field block moved below the fault block. Medium;
   bit-identical. **Done** with item 4 (`e884fae`): `BulkStateField`, the hardening and wet-dry
   hooks of `PowerLaw`, the mediator loop and the refusals; `FaultWorkKernel` came in Stage 4; the
   grain-size block moved below the fault block with item 6. `GrainSizeEvolution` is not derived
   from the base class: its steady-state, piezometer and checkpoint paths would all change for no
   gain yet.
4. Strain hardening (4.1). Medium. **Done** (`e884fae`): gates 5.1a (2.5e-12 to 2.9e-12) and 5.1b
   (bit-identical with S frozen at S_ref); restart bit-identical.
5. Water content (4.2). Medium. **Done** (`a2b91d2`): gate 5.2a (6.9e-13 to 1.4e-12), 5.2b (chi = 1
   and 0 bit-identical to the wet and the dry run); pressure-linked and pulsed supply not yet.
6. Cataclastic sink and pinning cap (4.3), with the energy-budget test. Small. **Done**
   (`672b999`), cell-mean form: gates 5.3a (closed forms to 1e-11, capped steady state to 2.7e-15,
   fCat = 0 bit-identical) and 5.3b (budget to 4e-16 and 1.1e-14, 30% stored). The two-component
   form (a fault-zone grain size over w_cat) is left to the sub-grid core width of 3.7. Found and
   fixed on the way: A-98 (`8df49c8`) and A-99 (`fc8938f`).
7. Pulsed recharge (4.4), with the diffusion unit test; depends on Stage 4's per-fault
   `PressureEq`. Small. **Done** (`d6a9556`): `bcB_q0`, pulses (cosine, Gaussian, box; periodic),
   the source in depth; gate 5.4a (flux step against the series to 5e-7, backward Euler first order
   and exact pulse volumes, source to 3.5e-6) and 5.4b (bit-identical without the keys). The
   dataset is `qb` (m/s), not `bcB`. Found and fixed on the way: A-100 (`e0ae91a`, the bottom flux
   from the list ends, the pitfall named above) and A-101 (`6e39b4d`, `linSolver` ignored).
8. Evolving cohesion (4.6) on the elastic mediator first (no Stage 4 dependency; can be done in
   parallel with 3-7). Small to medium. **Done** (`c2f556e`), both quasi-dynamic mediators, with reseal
   hardening: gate 5.6a (spring slider, Dheal = 0.1 m: 0.66% explicit, 0.40% implicit against the
   single-degree-of-freedom recurrence; Dheal = 1 mm: implicit within 0.40% of explicit, whose
   events cost 12 times the steps) and 5.6b (bit-identical without the keys). The implicit update
   needs a step bound (`cohesionStepFrac`): one step must not remove the cohesion at once. Found and
   fixed on the way: A-102 (`b03b76e`, the monitor's step bound lost on restart). The sign test
   with two faults (c) belongs to the run matrix.
9. Fabric and cement states (4.5). Small once 3 exists. **Done** (`986f39c`): gate 5.5a (unit
   solutions to 2e-12; frozen at zero bit-identical; the strength factor exact to 3e-16). The
   two-fault sign tests 5.5b and 5.5c belong to the run matrix.
10. Diagnostics (section 5) and the run matrix; the library. **Diagnostics done** (`6e48e56`,
   `3b2bf77`): root-strength probes in `faultSeries.txt` (`seriesDepth`), `switches.csv`, and the
   generator options `--lock-depth`, `--series-depth`, `--set`. **Run matrix: blocked on the root
   geometry.** Measured in the stage 4 baseline (two faults 20 km apart, ex4 creep and geotherm):
   the faults reach the bottom of the domain and below 15 km creep frictionally at 0.3 to 0.9 vL, so
   the viscous strain rate next to them is 1e-21 /s at 12-18 km and 4e-17 to 4e-14 /s at 25-40 km.
   Locked below 20 km (`--lock-depth 20`, 440 yr), the flow beneath is broad: 8e-15 /s below
   30 km with half-widths over 15 km and the far field at half that rate, so the two roots merge
   into one substrate and do not form the two roots in series of section 2. Bulk state fields
   driven by these rates act on 1e5-yr scales with laboratory strains. Distinct, localized roots
   need a choice: a prescribed weak band beneath each fault (`disl_A` as a body field from a file),
   a localizing feedback (grain-size-sensitive creep with the wattmeter), or a wider fault
   separation; and the state strain scales then follow from the root width through the 0D screen
   (alternation for n = 3, R = 10, beta = 2, T_D = 0.1, T_h = 10: period 24.7 w gamma/vL).
11. Optional: fixed anisotropy (4.7).

Each physics commit is marked as altering published behaviour only when it is on, so it can be
reverted alone.

## 7. Verification gates

| Gate | Test | Pass criterion |
|---|---|---|
| 5.0 | infrastructure commit, `tools/regress.sh compare` on ex1, ex2, `ex4s`, grain-size case | bit-identical |
| 5.0b | sub-grid core width (3.7) with `w_sz` equal to the cell width | bit-identical to the cell-scale form; a 10 m-cell run and a 50 m-cell run with `w_sz = 10 m` agree on the root strength history to 10 % |
| 5.1a | hardening unit test, constant and zero strain rate | `S(t)` within 1e-10 of the analytic solution |
| 5.1b | single fault, `S` frozen at `S_ref` | bit-identical to Stage 4 |
| 5.1c | two faults, hardening alone | partition returns to 0.5 after a perturbation (equal sharing) |
| 5.2a | water unit test; `chi = 1` and `chi = 0` frozen | analytic to 1e-10; bit-identical to the wet and the dry run |
| 5.3a | grain-size sink unit tests; `f_cat = 0` | analytic; bit-identical |
| 5.3b | energy budget over an event | `integral (Qfric + f_cat Qfault) = sum tau_k V_k` to 1e-6 |
| 5.4a | flux step at constant permeability | gradient `Delta q eta/k` at long times; transient against the 1D solution |
| 5.4b | no pulse keys | bit-identical |
| 5.5a | fabric and cement unit tests; `Phi = 0`, `Psi = 0` frozen | analytic; bit-identical |
| 5.5b | fabric as pure feedback, two faults | lock-in threshold within a factor 2 of `n beta R T_D/2 = 1` from the 0D model |
| 5.5c | cement alone, two faults | equal sharing |
| 5.6a | spring slider with cohesion healing | recurrence against the single-degree-of-freedom prediction within 2 % |
| 5.6b | keys absent | bit-identical |
| 5.7 | every module: restart | bit-identical to the uninterrupted run |
| 5.8 | every module: `-objects_dump`, `checkkeys.py` on all inputs | clean |
| 5.9 | first alternating 2D run | at least three switches; slip per period and period reported against the memory's strain and time scales (`gamma_h`, `tau_r` for 4.1) |

## 8. Out of scope, and why

- Mineral reactions and phase changes (W1, W9, S4, S5, S10): no mineralogy, and the paper itself
  rejects them as irreversible or needing a change of depth.
- Fluid flow through the bulk (W5 porosity as a weak phase, S2 in the bulk, porosity waves, creep
  cavitation): the pore pressure is one-dimensional along each fault; a two-phase bulk model is a
  different code. The water field of 4.2 is the phenomenological stand-in.
- Melt rheology (W7 as a persistent weak layer): flash heating gives the coseismic weakening only.
- Geometry (W15, W16, S8, S11, W17): planar vertical faults in antiplane strain carry no
  along-strike structure, fold geometry or fabric rotation. A fixed anisotropy (4.7) is the most
  that can be represented.
- A modulus contrast across a fault: Route B of the design document needs continuous `mu`.

## 9. References

Entries marked (abstract) were checked against the abstract or citing text only; the rest against
the full text or the paper's own bibliography.

Fault systems and the hypothesis
- Cawood, T. K. & Dolan, J. F. (2024). Seismica 3(2). doi:10.26443/seismica.v3i2.1165
- Dolan, J. F., Bowman, D. D. & Sammis, C. G. (2007). Geology 35, 855. doi:10.1130/G23789A.1
- Dolan, J. F. et al. (2016). EPSL 446, 123. doi:10.1016/j.epsl.2016.04.011
- Dolan, J. et al. (2024). EPSL 625, 118484. doi:10.1016/j.epsl.2023.118484
- Gauriau, J. & Dolan, J. F. (2021). G-cubed 22. doi:10.1029/2021GC009938
- Lynch, J. C., Bürgmann, R., Richards, M. A. & Ferencz, R. M. (2003). GRL 30. doi:10.1029/2002GL016765
- Kenner, S. J. & Simons, M. (2005). GJI. doi:10.1111/j.1365-246X.2005.02460.x
- Ben-Zion, Y., Dahmen, K., Lyakhovsky, V., Ertas, D. & Agnon, A. (1999). EPSL. doi:10.1016/S0012-821X(99)00187-9
- Lyakhovsky, V., Ben-Zion, Y. & Agnon, A. (2001). JGR. doi:10.1029/2000JB900218

Strain hardening and recovery (4.1)
- Mecking, H. & Kocks, U. F. (1981). Acta Metall. doi:10.1016/0001-6160(81)90112-7
- Estrin, Y. & Mecking, H. (1984). Acta Metall. doi:10.1016/0001-6160(84)90202-5
- Kocks, U. F. & Mecking, H. (2003). Prog. Mater. Sci. doi:10.1016/S0079-6425(02)00003-8
- Karato, S., Rubie, D. C. & Yan, H. (1993). JGR. doi:10.1029/93JB00472
- Hansen, L. N. et al. (2019). JGR Solid Earth. doi:10.1029/2018JB016736
- Hansen, L. N. et al. (2021). JGR Solid Earth. doi:10.1029/2020JB021325
- Wallis, D. et al. (2021). Nature Communications. doi:10.1038/s41467-021-23633-8
- Breithaupt, T., Katz, R. F., Hansen, L. N. et al. (2023). PNAS 120. doi:10.1073/pnas.2203448120
- Masuti, S. et al. (2016). Nature. doi:10.1038/nature19783
- Masuti, S. & Barbot, S. (2021). Earth, Planets and Space. doi:10.1186/s40623-021-01543-9
- Chopra, P. N. (1997). Tectonophysics. doi:10.1016/S0040-1951(97)00134-0
- Hirth, G. & Tullis, J. (1992). JSG 14, 145. doi:10.1016/0191-8141(92)90053-Y
- Tullis, J. & Yund, R. A. (1989). GRL 16, 1343. doi:10.1029/GL016i011p01343
- Trepmann, C. A. & Stöckhert, B. (2013). Solid Earth 4, 263.
- Wallis, D., Lloyd, G. E. & Hansen, L. N. (2018). JSG 107, 25. doi:10.1016/j.jsg.2017.11.008
- Montési, L. G. J. (2004). JGR. doi:10.1029/2003JB002925 (postseismic relaxation of a power-law shear zone, a benchmark for the `S`-free limit)

Grain size (4.3)
- Austin, N. J. & Evans, B. (2007). Geology 35, 343. doi:10.1130/G23244A.1
- Tokle, L. & Hirth, G. (2021). JGR Solid Earth 126. doi:10.1029/2020JB021475
- Fukuda, J. et al. (2019). Solid Earth 10, 621. doi:10.5194/se-10-621-2019
- Tullis, J. & Yund, R. A. (1982). J. Geol. 90, 301. doi:10.1086/628681 (abstract)
- Platt, J. P. & Behr, W. M. (2011). Geology 39, 127. doi:10.1130/G31561.1
- Chester, J. S., Chester, F. M. & Kronenberg, A. K. (2005). Nature 437, 133. doi:10.1038/nature03942
- Wilson, B. et al. (2005). Nature 434, 749. doi:10.1038/nature03433
- Ma, K.-F. et al. (2006). Nature 444, 473. doi:10.1038/nature05253
- Pittarello, L. et al. (2008). EPSL 269, 131. doi:10.1016/j.epsl.2008.01.052 (abstract)
- Aben, F. M., Brantut, N. & Mitchell, T. M. (2020). JGR 125. doi:10.1029/2020JB019860
- Scholz, C. H. (1988). Geol. Rundschau 77, 319. doi:10.1007/BF01848693
- Rolandone, F., Bürgmann, R. & Nadeau, R. M. (2004). GRL 31. doi:10.1029/2004GL021379
- Ellis, S. & Stöckhert, B. (2004). JGR 109. doi:10.1029/2003JB002744
- Price, N. A., Johnson, S. E., Gerbi, C. C. & West, D. P. (2012). Tectonophysics 518-521, 63. doi:10.1016/j.tecto.2011.11.011 (abstract)
- Kirkpatrick, J. D. & Rowe, C. D. (2013). JSG 52, 183. doi:10.1016/j.jsg.2013.03.003
- Campbell, L. R. & Menegon, L. (2019). JGR 124. doi:10.1029/2019JB018052
- Hawemann, F. et al. (2019). Solid Earth 10, 1635. doi:10.5194/se-10-1635-2019
- Herwegh, M., Linckens, J., Ebert, A., Berger, A. & Brodhag, S. H. (2011). JSG 33, 1728. doi:10.1016/j.jsg.2011.08.011
- Linckens, J., Herwegh, M., Müntener, O. & Mercolli, I. (2011). JGR 116. doi:10.1029/2010JB008119 (abstract)
- Bercovici, D. & Ricard, Y. (2012). PEPI 202-203, 27. doi:10.1016/j.pepi.2012.05.003
- Stipp, M. & Tullis, J. (2003). GRL 30, 2088. doi:10.1029/2003GL018444
- Holyoke, C. W. & Kronenberg, A. K. (2010). Tectonophysics 494, 17. doi:10.1016/j.tecto.2010.08.001 (abstract)
- Cross, A. J., Prior, D. J., Stipp, M. & Kidder, S. (2017). GRL 44, 6667. doi:10.1002/2017GL073836

Fault valving and fluid pulses (4.4)
- Zhu, W., Allison, K. L., Dunham, E. M. & Yang, Y. (2020). Nature Communications 11, 4833. doi:10.1038/s41467-020-18598-z
- Allison, K. L. & Dunham, E. M. (2018). Tectonophysics 733, 232. doi:10.1016/j.tecto.2017.10.021
- Ozawa, S., Yang, Y. & Dunham, E. M. (2024). JGR Solid Earth. doi:10.1029/2024JB029165
- Yang, Y. & Dunham, E. M. (2021). JGR 126. doi:10.1029/2020JB021258
- Yang, Z., Yehya, A., Iwalewa, T. M. & Rice, J. R. (2021). JGR. doi:10.1029/2021JB021787
- Dal Zilio, L., Hegyi, B., Behr, W. & Gerya, T. (2022). Tectonophysics 838, 229516. doi:10.1016/j.tecto.2022.229516
- Petrini, C. et al. (2020). Tectonophysics 791, 228504. doi:10.1016/j.tecto.2020.228504 (abstract)
- Perez-Silva, A. et al. (2023). JGR 128. doi:10.1029/2022JB026332
- Farge, G., Jaupart, C. & Shapiro, N. M. (2021). JGR. doi:10.1029/2021JB021894
- Connolly, J. A. D. (1997). JGR 102, 18149. doi:10.1029/97JB00731
- Connolly, J. A. D. (2010). Elements 6, 165. doi:10.2113/gselements.6.3.165
- Louis, S., Luijendijk, E., Dunkl, I. & Person, M. (2019). Geology 47, 938. doi:10.1130/G46254.1
- Vignaroli, G. et al. (2022). Tectonophysics 827, 229269. doi:10.1016/j.tecto.2022.229269 (abstract)
- Hickey, K. A. et al. (2014). Economic Geology 109, 1461. doi:10.2113/econgeo.109.5.1461 (abstract)
- Menegon, L. & Fagereng, A. (2021). Geology 49, 1255. doi:10.1130/G49012.1
- Sibson, R. H. (1992). Tectonophysics 211, 283. doi:10.1016/0040-1951(92)90065-E
- Sibson, R. H. (2020). Earth, Planets and Space 72, 31. doi:10.1186/s40623-020-01153-x
- Hirth, G. & Beeler, N. M. (2015). Geology 43, 223. doi:10.1130/G36361.1
- Rice, J. R. (1992). In Fault Mechanics and Transport Properties of Rocks, 475-503.
- Carslaw, H. S. & Jaeger, J. C. (1959). Conduction of Heat in Solids, sections 2.9 and 3.8 (the flux-step series solutions of 4.4).

Water content (4.2)
- Hirth, G., Teyssier, C. & Dunlap, W. J. (2001). Int. J. Earth Sci. 90, 77. doi:10.1007/s005310000152
- Tokle, L., Hirth, G. & Behr, W. M. (2019). EPSL 505, 152. doi:10.1016/j.epsl.2018.10.017
- Fukuda, J., Holyoke, C. W. & Kronenberg, A. K. (2018). JGR 123. doi:10.1029/2017JB015133
- Lusk, A. D. J., Platt, J. P. & Platt, J. A. (2021). JGR 126. doi:10.1029/2020JB021302
- Lu, L. X. & Jiang, D. (2019). JGR 124. doi:10.1029/2018JB016226
- Rutter, E. H. & Brodie, K. H. (2004). JSG 26, 259 and 2011. doi:10.1016/S0191-8141(03)00096-8 (abstract)
- Gleason, G. C. & Tullis, J. (1995). Tectonophysics 247, 1. doi:10.1016/0040-1951(95)00011-B (compilation values)
- Luan, F. C. & Paterson, M. S. (1992). JGR 97, 301. doi:10.1029/91JB01748 (abstract)
- Koch, P. S., Christie, J. M., Ord, A. & George, R. P. (1989). JGR 94, 13975. doi:10.1029/JB094iB10p13975
- Kronenberg, A. K. & Tullis, J. (1984). JGR 89, 4281. doi:10.1029/JB089iB06p04281
- Jaoul, O., Tullis, J. & Kronenberg, A. (1984). JGR 89, 4298. doi:10.1029/JB089iB06p04298
- Kronenberg, A. K., Kirby, S. H., Aines, R. D. & Rossman, G. R. (1986). JGR 91, 12723. doi:10.1029/JB091iB12p12723
- Kronenberg, A. K. et al. (2020). Geology 48, 557. doi:10.1130/G47041.1
- Finch, M. A., Weinberg, R. F. & Hunter, N. J. (2016). Geology 44, 599. doi:10.1130/G37972.1
- Kilian, R. et al. (2016). JGR 121. doi:10.1002/2015JB012771
- Palazzin, G. et al. (2018). JSG 114, 95. doi:10.1016/j.jsg.2018.05.021 (via Cawood & Dolan 2024)
- Stünitz, H. et al. (2017). JGR 122, 866. doi:10.1002/2016JB013533
- Heilbronner, R. & Tullis, J. (2006). JGR 111. doi:10.1029/2005JB004194
- Smith, D. L. & Evans, B. (1984). JGR 89, 4125. doi:10.1029/JB089iB06p04125
- Gleason, G. C. & DeSisto, S. (2008). Tectonophysics 446, 16. doi:10.1016/j.tecto.2007.09.006 (via Cawood & Dolan 2024)
- Précigout, J., Prigent, C., Palasse, L. & Pochon, A. (2017). Nature Communications 8, 15736. doi:10.1038/ncomms15736
- Précigout, J., Stünitz, H. & Villeneuve, J. (2019). Scientific Reports 9. doi:10.1038/s41598-019-40020-y
- Fukuda, J. & Shimizu, I. (2019). Earth, Planets and Space 71. doi:10.1186/s40623-019-1117-4
- Gratier, J.-P., Dysthe, D. K. & Renard, F. (2013). Advances in Geophysics 54, 47. doi:10.1016/B978-0-12-380940-7.00002-0
- Pluymakers, A. M. H. & Spiers, C. J. (2015). Geol. Soc. London Spec. Publ. 409. doi:10.1144/SP409.6
- Kaatz, L. et al. (2023). G-cubed 24. doi:10.1029/2022GC010830
- Gueydan, F., Leroy, Y. M. & Jolivet, L. (2004). JGR 109. doi:10.1029/2003JB002806

Fabric, cementation and damage rheology (4.5)
- Lyakhovsky, V., Ben-Zion, Y. & Agnon, A. (1997). JGR 102, 27635. doi:10.1029/97JB01896
- Lyakhovsky, V. & Ben-Zion, Y. (2008). GJI. doi:10.1111/j.1365-246X.2007.03652.x
- Fuchs, L. & Becker, T. W. (2019). GJI 218, 601. doi:10.1093/gji/ggz167
- Huismans, R. S. & Beaumont, C. (2003). JGR. doi:10.1029/2002JB002026
- Lavier, L. L., Buck, W. R. & Poliakov, A. N. B. (2000). JGR. doi:10.1029/2000JB900108
- Montési, L. G. J. (2007). GRL. doi:10.1029/2007GL029250
- Montési, L. G. J. (2013). JSG 50, 254. doi:10.1016/j.jsg.2012.12.011
- Holyoke, C. W. & Tullis, J. (2006). JSG. doi:10.1016/j.jsg.2006.01.008
- Dabrowski, M., Schmid, D. W. & Podladchikov, Y. Y. (2012). JGR. doi:10.1029/2012JB009183
- Williams, R. T. & Fagereng, A. (2022). Reviews of Geophysics 60. doi:10.1029/2021RG000768
- Walderhaug, O. (1994). J. Sediment. Res. 64, 324. doi:10.2110/jsr.64.324; (1996) AAPG Bull. 80, 731
- Lander, R. H., Larese, R. E. & Bonnell, L. M. (2008). AAPG Bull. doi:10.1306/07160808037
- Tenthorey, E., Cox, S. F. & Todd, H. F. (2003). EPSL. doi:10.1016/S0012-821X(02)01082-8
- Tenthorey, E. & Cox, S. F. (2006). JGR 111. doi:10.1029/2005JB004122
- Renard, F., Gratier, J.-P. & Jamtveit, B. (2000). JSG. doi:10.1016/S0191-8141(00)00064-X (abstract)
- Woodcock, N. H., Dickson, J. A. D. & Tarasewicz, J. P. T. (2007). Geol. Soc. London Spec. Publ. 270, 43. doi:10.1144/GSL.SP.2007.270.01.03

Cohesion and frictional healing (4.6)
- Nakatani, M. & Scholz, C. H. (2004). JGR 109. doi:10.1029/2001JB001522 and doi:10.1029/2003JB002938
- Yasuhara, H., Marone, C. & Elsworth, D. (2005). JGR 110. doi:10.1029/2004JB003327
- Karner, S. L., Marone, C. & Evans, B. (1997). Tectonophysics. doi:10.1016/S0040-1951(97)00077-2 (abstract)
- Muhuri, S. K., Dewers, T. A., Scott, T. E. & Reches, Z. (2003). Geology. doi:10.1130/G19601.1 (abstract)
- van den Ende, M. P. A. & Niemeijer, A. R. (2019). Scientific Reports 9. doi:10.1038/s41598-019-46241-5
- Chen, J., Verberne, B. A. & Spiers, C. J. (2015). EPSL. doi:10.1016/j.epsl.2015.03.044
- Chen, J., Niemeijer, A. R. & Spiers, C. J. (2017). JGR. doi:10.1002/2017JB014226
- van den Ende, M. P. A. et al. (2018). Tectonophysics. doi:10.1016/j.tecto.2017.11.040
- Angevine, C. L., Turcotte, D. L. & Furnish, M. D. (1982). Tectonics 1, 151. doi:10.1029/TC001i002p00151
- Sleep, N. H. & Blanpied, M. L. (1992). Nature 359, 687. (abstract)

Viscous anisotropy (4.7)
- Mühlhaus, H.-B., Moresi, L., Hobbs, B. & Dufour, F. (2002). Pure Appl. Geophys. doi:10.1007/s00024-002-8737-4
- Lev, E. & Hager, B. H. (2008). GJI. doi:10.1111/j.1365-246X.2008.03731.x
- Hansen, L. N., Zimmerman, M. E. & Kohlstedt, D. L. (2012). Nature 492, 415. doi:10.1038/nature11671
- Mameri, L. et al. (2019). PEPI. doi:10.1016/j.pepi.2018.11.002; (2020) GJI. doi:10.1093/gji/ggaa400
- Christensen, R. M. (1987). Geophys. J. R. Astron. Soc. doi:10.1111/j.1365-246X.1987.tb01666.x
- Treagus, S. H. (2003). Tectonophysics. doi:10.1016/S0040-1951(03)00239-7
- Hill, R. (1948). Proc. R. Soc. A 193, 281.

## Appendix A. The two-root screening script

The script below is the one behind section 2; it becomes `tools/zero_d_trading.py` in the first
commit of the stage. It needs only numpy.

```python
#!/usr/bin/env python3
"""Two shear-zone roots in series: do they share, lock in, or alternate?

Both roots carry the same far-field stress and their slip rates add to the plate rate, so the
partition is f_1 = S_1^-n / (S_1^-n + S_2^-n). Each root has a slow strengthening memory h (grows
with strain, relaxes in time) and a fast weakening feedback D (grows with strain, heals in time):

    S = (1 + h) exp(-beta D),   dh/dt = f - h/Th,   dD/dt = R f - D/TD

in time units of w*gamma_h/vL (w shear-zone width, gamma_h hardening strain). Th and TD are the
recovery and healing times in those units, R = gamma_h/gamma_D. Equal sharing is unstable when
N_loc = n beta R TD/2 > 1; lock-in is permanent when N_loc > (n/2) ln(1 + Th), roughly.

  python3 tools/zero_d_trading.py --n 3 --beta 0.4 --Th 10 --TD 0.3 --R 10
  python3 tools/zero_d_trading.py --scan      # the regime table of docs/REVERSIBLE_STRENGTH_PLAN.md
"""
import argparse, itertools
import numpy as np


def integrate(n, beta, Th, TD, R, T=None, dt=None, eps=1e-3):
    """RK4 from the symmetric steady state with a relative perturbation eps of D. Returns t, f_1."""
    n, beta, Th, TD, R = (np.atleast_1d(np.asarray(v, float)) for v in (n, beta, Th, TD, R))
    M = max(v.size for v in (n, beta, Th, TD, R))
    n, beta, Th, TD, R = (np.broadcast_to(v, M).copy() for v in (n, beta, Th, TD, R))
    h = np.stack([0.5*Th, 0.5*Th]); D = np.stack([0.5*R*TD*(1 + eps), 0.5*R*TD*(1 - eps)])

    def rhs(h, D):
        S = (1 + h)*np.exp(-beta*D); w = S**-n; f = w/w.sum(axis=0)
        return f - h/Th, R*f - D/TD, f[0]

    T = T or 30*Th.max(); dt = dt or 0.01*min(1.0, TD.min(), 1.0/R.max())
    N = int(T/dt); every = max(1, N//20000); ts, fs = [], []
    for k in range(N):
        k1h, k1D, f1 = rhs(h, D); k2h, k2D, _ = rhs(h + 0.5*dt*k1h, D + 0.5*dt*k1D)
        k3h, k3D, _ = rhs(h + 0.5*dt*k2h, D + 0.5*dt*k2D); k4h, k4D, _ = rhs(h + dt*k3h, D + dt*k3D)
        h = h + dt*(k1h + 2*k2h + 2*k3h + k4h)/6; D = D + dt*(k1D + 2*k2D + 2*k3D + k4D)/6
        if k % every == 0:
            ts.append(k*dt); fs.append(f1.copy())
    return np.array(ts), np.array(fs)


def classify(t, f, Th):
    """Regime from the second half of the series: sharing, lock-in or alternation (with period)."""
    half = t > 0.5*t[-1]; f = f[half]; t = t[half]
    amp = f.max() - f.min(); cross = np.flatnonzero(np.diff(np.sign(f - 0.5)) != 0)
    if amp < 0.05 and abs(f[-1] - 0.5) < 0.05: return 'sharing', amp, float('nan'), f.max(), f.min()
    if amp < 0.05: return 'lock-in', amp, float('nan'), f.max(), f.min()
    per = np.mean(np.diff(t[cross]))*2 if cross.size >= 3 else float('nan')
    return ('alternation' if cross.size >= 2 else 'transient'), amp, per, f.max(), f.min()


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--n', type=float, default=3.0); p.add_argument('--beta', type=float, default=0.4)
    p.add_argument('--Th', type=float, default=10.0); p.add_argument('--TD', type=float, default=0.3)
    p.add_argument('--R', type=float, default=10.0); p.add_argument('--T', type=float, default=None)
    p.add_argument('--scan', action='store_true', help='regime table over N_loc and Th for n = 2, 3, 5')
    p.add_argument('--csv', default=None, help='write the f_1(t) series of a single run to this file')
    a = p.parse_args()
    if a.scan:
        cases = list(itertools.product([2., 3., 5.], [1.2, 1.5, 2., 3., 4., 6.], [3., 10., 30., 100.]))
        n, Nloc, Th = (np.array(c) for c in zip(*cases)); TD, R = 0.3, 10.0; beta = 2*Nloc/(n*R*TD)
        t, f = integrate(n, beta, Th, TD, R, T=30*Th.max(), dt=0.003)
        print('N_loc = n beta R TD/2 with TD = %g, R = %g; columns Th = 3, 10, 30, 100; A = alternation (period/Th)' % (TD, R))
        for nn in [2., 3., 5.]:
            print('n = %g' % nn)
            for c in [1.2, 1.5, 2., 3., 4., 6.]:
                row = []
                for th in [3., 10., 30., 100.]:
                    j = cases.index((nn, c, th)); reg, amp, per, fmax, fmin = classify(t, f[:, j], th)
                    row.append({'sharing': '   S    ', 'lock-in': '   L    '}.get(reg, 'A(%5.2f)' % (per/th)))
                print('  N_loc %4.1f: %s' % (c, ' '.join(row)))
        return
    t, f = integrate(a.n, a.beta, a.Th, a.TD, a.R, T=a.T)
    reg, amp, per, fmax, fmin = classify(t, f[:, 0], a.Th)
    print('N_loc = %.2f, upper bound (n/2) ln(1+Th) = %.2f' % (a.n*a.beta*a.R*a.TD/2, a.n/2*np.log(1 + a.Th)))
    print('regime: %s; swing of f_1: %.2f .. %.2f; period: %.3g (%.2f Th)' % (reg, fmin, fmax, per, per/a.Th))
    if a.csv:
        np.savetxt(a.csv, np.column_stack([t, f[:, 0]]), header='t f1', comments='')


if __name__ == '__main__':
    main()
```
