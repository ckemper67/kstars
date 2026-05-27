# MPC Guider Architecture

This document describes the current design of the MPC (Model Predictive
Control) guide algorithm in `kstars/ekos/guide/internalguide/`, how it
relates to the other guiders (Linear, Hysteresis, GPG), and the known
limitations -- with possible options forward where the limitation is
worth a code change.

Audience: someone reading the MPC source for the first time and trying
to decide where to extend it. Not a tutorial on MPC theory.

## 1. Scope and intent

The internal guider exposes four algorithms behind the same `guide(offset)`
interface:

- `LinearGuider`     -- proportional + averaging window
- `HysteresisGuider` -- proportional + dead-band
- `GaussianProcessGuider` (GPG) -- non-parametric Gaussian process
                        with a periodic + Matern kernel, learns the
                        worm period online
- `MPCGuider`        -- parametric model-predictive controller with
                        Laguerre-network parametrization, online FFT
                        period detection, Luenberger observer, and an
                        Internal Model Principle (IMP) disturbance
                        model

The first three are reference baselines. MPC is the active development
target.

Design intent of MPCGuider:

- *Predict* and *cancel* periodic worm error before the camera sees it
  (one-frame anticipation under a typical 2-4s exposure).
- Learn the worm period from the open-loop residual via FFT; no user
  configuration of PE.
- Stay numerically robust under plant-model mismatch (compliance,
  backlash, stiction, transport delay) instead of issuing
  catastrophically large corrections.

## 2. Top-level data flow

    offset (arcsec, from guide star centroid)
      |
      v
    MPCGuider::guide(offset)
      |--- dt = wall-clock since last call
      |--- reconstruct delta_theta, velocity, delta_omega
      |--- push (offset, velocity) into history (50-frame ring)
      |--- passive harmonic-drive classifier (ratio of vRMS / posRMS)
      |--- FFT period estimator (estimates T1, T2 from open-loop error)
      |--- if dt changed / parameters changed / periods updated:
      |      rebuild TelescopePlant
      |      rebuild MPCSolver gains (K_r, K_x)
      |
      |--- Luenberger observer step:
      |      x_hat(k) = x_hat_pred(k) + L * (offset - C * x_hat_pred(k))
      |
      |--- delta_u = K_r * setpoint - K_x . delta_x_hat + DOB feedforward
      |--- punch-through bias on sustained direction reversals
      |--- per-step rate limit |delta_u| <= 5"
      |--- current_u += delta_u
      |
      |--- predict x_hat_pred(k+1) for next call
      v
    return -(current_u - prev_u)  -- the arcsec correction to apply

## 3. Files and responsibilities

| File                       | Role                                          |
|----------------------------|-----------------------------------------------|
| `mpcguider.{h,cpp}`        | Public entry point. Owns observer state,      |
|                            | FFT estimator, passive-detection counters,    |
|                            | timing/dt tracking. Calls into Plant/Solver.  |
| `TelescopePlant.{h,cpp}`   | Builds the discrete-time state-space matrices |
|                            | (A_aug, B_aug, C_aug) the MPC uses for both   |
|                            | prediction and gain computation. Selects      |
|                            | between three model orders (see Section 4).   |
| `MPCSolver.{h,cpp}`        | Builds the Laguerre-parametrized cost-to-go,  |
|                            | solves for feedback gains K_r, K_x via SVD.   |
|                            | Hot path: `computeDeltaU(x_aug, setpoint)`.   |
|                            | Owns rate limit, punch-through, current_u.    |
| `LaguerreNetwork.{h,cpp}`  | The Laguerre-function basis (N functions with |
|                            | pole a) used to compress the prediction       |
|                            | horizon. Reduces the QP dimension from Np to  |
|                            | N (typically Np=100, N=10).                   |
| `fftperiodestimator.{h,cpp}`| Sliding-window FFT with quadratic-interp     |
|                            | peak-picking. Returns up to two strongest     |
|                            | periods. Gated by total-time-observed.        |

## 4. Plant model selection

`TelescopePlant::discretize()` selects one of three augmented systems
based on what the caller has declared. The selection logic, in order:

    if (omega1_ > 0.0)                    -> 7-state rigid IMP
    else if (Ks_ > 1e7 || J_axis_ < 1e-6) -> 3-state rigid
    else                                  -> 5-state flexible (no IMP)

`omega1_` is set by `setDisturbanceFrequencies()` after FFT learning.
`Ks_` and `J_axis_` are set by `setSpringStiffness()` /
`setAxisInertia()` (compliance declaration -- Section 6).

### 4a. 3-state rigid (default)

State: `[delta_theta_m, delta_v_m, theta_track]`. Motor modeled as a
pure integrator with friction p = Bf / Jm. Used when no periodic
disturbance has been learned and the mount is treated as a rigid
single-mass system. Asymptotically gives `theta_track(k+1) =
theta_track(k) + B * u`.

### 4b. 7-state rigid IMP

State: `[delta_theta_m, delta_v_m, d1, d1_dot, d2, d2_dot, y_track]`.
Two damped harmonic oscillators (d1 at omega1_, d2 at omega2_) appended
to the motor block to model the periodic disturbance internally. The
controller's optimal `u(k)` then includes a feedforward cancellation
of the next-frame disturbance value. Selected once FFT has learned a
period.

The output map is `y = theta_m + d1 + d2` (rows of `Caug` are
[1, 0, 1, 0, omega2>0 ? 1 : 0, 0, 1]); the tracking row only updates
from disturbance and disturbance feedforward to avoid double-counting
motor position in prediction.

### 4c. 5-state flexible (no IMP)

State: `[delta_theta_m, delta_v_m, delta_theta_a, delta_v_a,
theta_a_track]`. Two-mass system: motor inertia Jm coupled to axis
inertia Ja through a spring of stiffness Ks and damping Bs. Selected
when compliance has been declared explicitly via
`setMechanicalParams()`.

Currently mutually exclusive with the IMP path (you cannot have both
flexible coupling and periodic-disturbance cancellation). See Section
7 for why.

### 4d. Discretization

- 3-state rigid: exact analytic discretization of the 1st-order motor.
- 7-state IMP: exact analytic motor + analytic damped-oscillator
  discretization for each harmonic block (closed-form sin/cos/exp).
- 5-state flexible: **2nd-order Taylor expansion**
  `Ad = I + Ac*dt + 0.5*Ac^2*dt^2`. This is the weak link -- see
  Section 7.

## 5. Observer, FFT, and online learning

### 5a. FFT period estimator

`mpcguider.cpp` pushes `(t, openLoopError = offset - current_u)` into
`FFTPeriodEstimator` every frame. Every 10 frames after frame 40 it
calls `estimatePeriods()` with a 100s-elapsed gate.

When the estimate is accepted (>2% change), `omega1_` and `omega2_`
are EMA-smoothed (alpha=0.1), `m_PeriodsLearned = true` is set, and
the plant is rebuilt -- which routes selection to the 7-state IMP if
no compliance has been declared.

The estimator window is 600 samples, with quadratic interpolation on
the spectrum peak for sub-bin frequency resolution.

### 5b. Luenberger observer (IMP path only)

The 7-state path runs an explicit observer:

    x_hat_pred(k) = Ad6 * x_hat(k-1) + Bd6 * delta_u(k-1)
    x_hat(k)      = x_hat_pred(k) + L * (offset - C * x_hat_pred(k))

with gains tuned for slow motor-state correction and aggressive
disturbance-state correction:

    L = [0.02, 0.002/dt, 0.60, 0.60*omega1, 0.30, 0.30*omega2]

The motor-position gain is deliberately small to prevent the observer
from "absorbing" periodic error into the motor state instead of the
disturbance states.

The 5-state and 3-state paths do **not** use an observer; they pass
`[0, 0, ..., offset]` directly as `x_aug` (incremental motor/axis
states zeroed, absolute position taken from the measurement).

### 5c. Passive harmonic detection

A second, simpler classifier compares `rms(velocity) / rms(offset)`
over the last 40 frames. When ratio > 0.08 and RMS > 0.15" sustained
for 10 frames, `m_IsHarmonicDetected` flips. It currently only logs a
diagnostic; no parameter is actually changed (a stale `raise R to 1.0`
log message remains; see Plan Section 1b).

## 6. Compliance declaration (current)

`MPCGuider::setMechanicalParams(Ks, Bs, J_axis)` is the public hook
for declaring that the mount has finite spring stiffness and a
separate axis inertia. When set such that `Ks < 1e7 && J_axis > 1e-6`,
the controller:

- forces `m_Initialized = false` (triggers plant rebuild),
- propagates Ks, Bs, J_axis to `TelescopePlant`,
- **skips** the FFT-driven `setDisturbanceFrequencies()` call even if
  periods were learned,
- builds the 5-state flexible plant.

This is configuration, not online detection. There is no system-ID
step that estimates Ks/J_axis from observed data; the caller has to
declare it (the benchmark scenarios do this for H17 and H16).

Rationale for skipping IMP: the matrices for a fused
flexible-plus-IMP system would be 9-state, and the 5-state flexible
already uses the unstable 2nd-order Taylor discretization. Adding
two harmonic-oscillator blocks on top would require a proper matrix-
exponential discretization (scaling-and-squaring) before the gains
would be trustworthy.

## 7. Limitations

### 7a. Flexible-plant discretization is fragile

`TelescopePlant::discretize()` uses 2nd-order Taylor expansion for the
5-state flexible path. The continuous Ac matrix has entries of order
`Ks / Jm`. For a typical worm-drive mount with `Jm = 0.02`, `Ks = 250`,
this is 12500. The 2nd-order term `0.5 * Ac^2 * dt^2` with `dt = 2s`
gives matrix entries of order `0.5 * 1.5e8 * 4 = 3e8`. The resulting
Ad has eigenvalues far outside the unit circle, the MPC gains diverge,
and `computeDeltaU` returns NaN.

This is why compliance declaration is currently opt-in per scenario in
the benchmark, not blanket-applied:

- H17 (Jm=0.5, dt=4s): borderline stable, 0.000" RMS with declaration
  vs. 0.908" without.
- H16 (Jm=0.02, dt=4s): also stable in practice -- 0.000" with vs.
  0.384" without -- because the slow PE (T=480s) does not excite the
  resonance enough to manifest the bad gains.
- H7  (Jm=0.02, dt=2s, fast PE T=30s): NaN under declaration; opt-out
  preserves the 7-state IMP path which gives 1.515" RMS.

Forward options:

1. **Replace the 2nd-order Taylor with matrix-exponential discretization
   via scaling-and-squaring.** Standard textbook algorithm; about 100
   lines including Pade approximation. After this, compliance can be
   declared safely for all 2-mass scenarios. *Highest leverage, ~1 day.*

2. **Pre-check `Ad` eigenvalues; fall back to rigid when unstable.**
   Cheap (one `eigenvalues().array().abs().max()` call per rebuild),
   but masks the root cause and the user does not know why their
   declared compliance was ignored. *Useful as a safety net even
   after option 1.*

3. **Adaptive sub-stepping inside `discretize()`** -- split dt into N
   sub-steps until `Ac * (dt/N)` is small enough for Taylor. Works
   but mixes two concerns (modeling and integration) in one place.

### 7b. Flexible plant has no observer

The 3-state and 5-state paths feed `[0, 0, ..., offset]` directly as
`x_aug`. The motor/axis velocity states are assumed zero each frame,
which is wrong as soon as the mount has any inertia. The 7-state IMP
path does the right thing with a Luenberger observer.

Forward options:

1. **Add a 4-state observer to the flexible path** mirroring the
   IMP one. Needs an L gain matrix and a pred/correct cycle; same
   shape as the existing IMP observer.

2. **Identify motor/axis velocity from finite differences** of the
   measured offset (already done in `mpcguider.cpp` for diagnostics
   only). Pass `[0, delta_v_m_est, 0, delta_v_a_est, offset]` as
   `x_aug` instead of zeros. Cheap, but noise-sensitive.

### 7c. Compliance and IMP are mutually exclusive

A compliant mount with periodic worm error currently has to choose:
flexible plant (no PE feedforward) or rigid IMP (no compliance model).
H7 / H8 fall in this gap -- they have both, but the controller can
only model one at a time.

Forward options:

1. **Build a 9-state flexible IMP plant** -- 2-mass motor/axis plus
   two harmonic-oscillator blocks plus tracking output. The
   continuous-time A_c is 8x8 (with the tracking integrator giving
   the 9th row); discretization requires the same matrix-exp work
   from 7a. *Largest payoff; ~2-3 days.*

2. **Layer IMP feedforward outside the plant** -- run the harmonic
   oscillator states as a separate observer (same dynamics as today's
   d1/d2 block) and add `-K_d * (d1 + d2)` as an external feedforward
   to the flexible-plant MPC's `delta_u`. Decouples the two models
   but loses the joint optimization the MPC was meant to give.

### 7d. No online identification of mechanical parameters

The user (or the benchmark) has to know Ks, Bs, J_axis to declare
them. There is no calibration step that estimates them from a chirp,
step response, or natural disturbance.

Forward options:

1. **Add an open-loop chirp during the calibration phase** and fit
   omega_n, zeta, DC gain from the frequency response. Calibration
   already exists (`calibrationprocess.cpp`); adding a chirp is a
   small extension. *Practical for users who own the mount and can
   afford an extra ~30s of calibration.*

2. **Online recursive least squares** from the closed-loop residual.
   Hard to make robust (the closed-loop transfer function masks the
   plant) but does not need a special calibration step.

### 7e. Sawtooth and other non-sinusoidal PE

The IMP model uses two damped harmonic oscillators. A sawtooth has a
`1/n` Fourier comb; two harmonics capture only the fundamental and
second harmonic, leaving the higher modes uncancelled. H10 is the
benchmark scenario that exposes this -- MPC trails the simple guiders.

Forward options:

1. **Add more harmonic blocks** (4-6 oscillators tuned to detected
   harmonics 2 omega1, 3 omega1, ...). Each additional block adds two
   states; the matrices grow but stay structured-sparse.

2. **Switch to a state-space PE model fit to the actual waveform**
   (e.g., fit a periodic spline). More general but loses the clean
   "two-line IMP" interpretability.

### 7f. Backlash punch-through is heuristic

`MPCSolver::computeDeltaU` adds `+/-backlash` to delta_u on commanded
direction reversals, gated by `frames_since_reversal_ >= 30`. The
threshold was tuned to stay dormant under typical PE (T=30s, half-
period in frames << 30 for dt >= 2s) and engage on slow PE / step
recoveries. There is no model of the gear gap state itself.

Forward options:

1. **Model the gear gap as a 6th plant state** (gap position 0..2*B).
   The state evolves under commanded motion and the spring force
   only transmits when the state is at +B or -B. Would let the MPC
   plan the gap traversal as part of the optimization rather than as
   a fixed bias.

2. **Estimate backlash size online** from the asymmetry of step
   responses. Currently the user has to call `setBacklash()`; the
   controller could learn it. The benchmark numbers show that knowing
   the value matters more than punch-through timing -- H13 hits 0.000"
   when the value is set correctly.

### 7g. Per-step rate limit is a safety net, not a feature

`MAX_DELTA_U = 5.0` per frame caps the worst-case correction under
plant-model mismatch. It was added to prevent the rigid-plant
controller from issuing 590000" corrections on stiff systems. It
bounds the damage but does not address the cause; once Section 7a
is fixed, the rate limit is mostly cosmetic.

Forward options:

1. **Replace with a slew-rate-aware MPC formulation** that includes
   `|delta_u| <= U_max` as a constraint inside the optimization
   (proper QP solve instead of unconstrained SVD).

2. **Leave as-is.** The rate limit costs nothing and is a reasonable
   sanity guard regardless of the model.

## 8. Benchmark coverage

`Tests/ekos/guide/benchmarkguidealgorithms.cpp` exercises the algorithms
under 17 scenarios (S1-S6 + H1-H17) covering:

- clean and noisy sinusoidal PE (H1, S6)
- non-sinusoidal PE: sawtooth, triangle, half-rectified (H10-H12)
- non-stationary and stochastic PE (H2, H3)
- step disturbances and mid-run plant changes (H5, H15)
- period drift (H14)
- compliance and resonance (H7, H8)
- backlash (H13, H16, H17)
- slow-response motor (H17)

H17 is the scenario the compliance routing was added for. H7 and H8
are the scenarios that block "always declare compliance" until 7a is
fixed.

## 9. Where to start when extending

- **Adding a scenario**: a new entry in `benchmarkguidealgorithms.cpp`;
  follow the H1-H17 patterns. Use `runComplianceScenario` for 2-mass
  cases, `runScenario` otherwise.
- **Tuning Q, R, alpha, N**: `MPCGuider::setParameters()`. The
  benchmark's `--tune` mode sweeps over predefined grids.
- **Adding a plant model**: extend `TelescopePlant::discretize()` and
  the matching `x_aug` block in `MPCGuider::guide()`. `MPCSolver`
  adapts via `plant.getOrder()`.
- **Adding online detection** (compliance, backlash, friction): the
  passive classifier in `MPCGuider::guide()` is the template. Keep
  the trigger conservative and the action small.
