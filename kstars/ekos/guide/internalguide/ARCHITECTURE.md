# Ekos Internal Guider Architecture

This document describes the high-level architecture of the Ekos Internal Guider system within KStars.

## Overview

The guiding system follows a classic feedback control loop: **Capture -> Register -> Transform -> Control -> Actuate**. It is modular in two independent dimensions:

- **Registration algorithm** -- how the per-frame offset (dx, dy, dtheta) is measured from the image.
- **Correction algorithm** -- how that offset is converted into mount pulse durations.

Any registration algorithm can be paired with any correction algorithm.

## Key Components

### 1. Orchestration (`Ekos::Guide`)
The `Guide` class (in `guide.cpp`) acts as the module manager. It coordinates between the UI, the camera, the mount, and the specific guider instance (Internal or External like PHD2).

### 2. Guider Logic (`Ekos::InternalGuider`)
`InternalGuider` (in `internalguider.cpp`) implements the `GuideInterface`. It manages the state machine (IDLE, CALIBRATING, GUIDING, DITHERING) and handles the timing of exposures and pulses.

### 3. Math Engine (`cgmath`)
`cgmath` (in `gmath.cpp`) is the core of the guider. It:
- Dispatches each frame to the active registration algorithm (`findLocalStarPosition`).
- Converts the resulting pixel offset to RA/Dec arcseconds via the `Calibration` object.
- Invokes the active correction algorithm to compute pulse durations.

### 4. Transformation (`Calibration`)
The `Calibration` class (in `calibration.cpp`) stores the calibration result: a rotation matrix and scaling factors that convert pixel shifts (dx, dy) into RA/Dec arcsecond drifts.

### 5. Registration Algorithms

Registration algorithms measure the offset of the current frame relative to a reference. The active algorithm is selected via `Options::guideAlgorithm()` (the "Star Detection / Registration Algorithm" combo in the UI). All algorithms return a pixel position that `cgmath` subtracts from the lock position to obtain the drift.

| Enum | Class/File | Method |
|------|-----------|--------|
| `SMART_THRESHOLD` | `GuideAlgorithms` | Adaptive threshold centroiding |
| `SEP_THRESHOLD` | `GuideAlgorithms` | SExtractor-based single-star centroid |
| `CENTROID_THRESHOLD` | `GuideAlgorithms` | Simple center-of-mass centroid |
| `AUTO_THRESHOLD` | `GuideAlgorithms` | Automatic threshold selection |
| `NO_THRESHOLD` | `GuideAlgorithms` | No thresholding |
| `SEP_MULTISTAR` | `GuideStars` | Multi-star reference geometry matching |
| `DONUTS_REGISTRATION` | `Donuts::Registrar` | 4-quadrant phase-only image correlation |

#### DONUTS Registrar (`Donuts::Registrar`)
Inspired by McCormac et al. 2013 (PASP 125, 548). Differences from the original paper:
- 4-quadrant profile split (instead of full-frame 1-D projections) for rotation sensitivity.
- Phase-only correlation (Kuglin & Hines 1975) instead of amplitude cross-correlation, giving invariance to brightness changes.
- Weighted least-squares solver for translation (dx, dy), rotation (dtheta), and optionally isotropic scale.

`Donuts::Registrar` is stateful: `setReference()` stores the reference frame and resets the accumulated field rotation to zero. Each `measure()` call de-rotates the current frame by the internally-tracked accumulated rotation before correlation, then updates the accumulated rotation with the returned `Transform::dtheta`. `cgmath` reads `dtheta` from the returned `Transform` and emits it as `newRotationDelta` for display and future de-rotator integration.

### 6. Correction Algorithms

All three operate on arcsecond drift and return an arcsecond correction; `cgmath` converts that to a pulse duration. RA and DEC can use different algorithms independently.

- **`LinearGuider`** -- Proportional-plus-trend controller. Fits a least-squares slope to the last N samples and corrects for both the current offset and the drift velocity, preventing systematic walk under constant seeing bias.

- **`HysteresisGuider`** -- Proportional controller with output smoothing. Each correction is a weighted blend of the current error and the previous output (`output = gain * ((1-h)*error + h*last_output)`), reducing pulse chatter when seeing noise dominates.

- **`GPG` (Gaussian Process Guider)** -- Predictive RA-only controller that learns the mount's periodic error. Fits a Gaussian process to the residual error history and issues corrections even between frames ("dark guiding"). Contributed by the MPI-IS group.

## Data Flow (Guiding Loop)

1. **Capture:** `InternalGuider` requests a frame from the camera.
2. **Analysis:** Upon arrival, the image (`FITSData`) is passed to `cgmath::performProcessing`.
3. **Registration:** `cgmath` calls `findLocalStarPosition`, which dispatches to the active registration algorithm and returns a pixel position (or a synthesized virtual position for DONUTS_REGISTRATION).
4. **Coordinate Mapping:** The pixel drift is rotated and scaled into RA/Dec arcseconds via the `Calibration` object.
5. **Pulse Calculation:** The arcsecond drift is processed by the active correction algorithm (`LinearGuider`, `HysteresisGuider`, or `GPG`).
6. **Execution:** The calculated pulses are emitted as signals, which `InternalGuider` forwards to the mount via INDI.

## Hardware Interaction

### 1. Pulse Guiding (Translation)
1. **`InternalGuider`**: Emits `newMultiPulse(ra_dir, ra_ms, dec_dir, dec_ms)`.
2. **`Guide`**: Receives the signal and calls `m_Guider->doPulse()`.
3. **`ISD::Guider`**: Maps directions to INDI properties (`TELESCOPE_TIMED_GUIDE_WE`, `TELESCOPE_TIMED_GUIDE_NS`).
4. **INDI Driver**: Executes the physical movement on the mount or ST4 port.

### 2. Adaptive Optics (Fast Translation)
1. **`InternalGuider`**: Processes frames at high speed (50 Hz+) and calculates (dx, dy).
2. **`Guide`**: Forwards corrections to the `ISD::AdaptiveOptics` device.
3. **`ISD::AdaptiveOptics`**: Maps corrections to INDI AO properties (`AO_NS`, `AO_WE`).
4. **AO Hardware**: Rapidly tilts the mirror to stabilize the image.

### 3. Rotator Integration (Rotation -- DONUTS_REGISTRATION)
`cgmath` already emits `newRotationDelta(double dTheta)` on each frame when DONUTS_REGISTRATION is active. Full de-rotator integration is planned:
1. **`InternalGuider`**: Will emit `newRotationCorrection(double dTheta)`.
2. **`Guide`**: Will forward to the active `ISD::Rotator`.
3. **`ISD::Rotator`**: Will update `ABS_ROTATOR_ANGLE`.
4. **INDI Rotator Driver**: Adjusts the physical sensor orientation.

## Extensibility

### Adding a new registration algorithm
1. Add an enum value after `DONUTS_REGISTRATION` in `gmath.h`.
2. Implement the measurement in a new class or within `GuideAlgorithms`.
3. Add a branch in `cgmath::findLocalStarPosition` that returns a pixel position (synthesize a virtual position relative to `targetPosition` if the algorithm measures frame-to-frame offsets rather than absolute positions).
4. Add a combo item to `kcfg_GuideAlgorithm` in `opsguide.ui`.

### Adding a new correction algorithm
1. Implement a class with a `guide(double arcsecDrift) -> double` interface (see `LinearGuider` as a template).
2. Add an enum value to `RAGuidePulseAlgorithm` / `DECGuidePulseAlgorithm` in `kstars.kcfg`.
3. Branch into the new algorithm in `cgmath::processAxis`.
