# Ekos Internal Guider Architecture

This document describes the high-level architecture of the Ekos Internal Guider system within KStars.

## Overview

The guiding system follows a classic feedback control loop: **Capture -> Measure -> Transform -> Control -> Actuate**. It is designed to be modular, allowing for different star detection algorithms and pulse calculation strategies.

## Key Components

### 1. Orchestration (`Ekos::Guide`)
The `Guide` class (in `guide.cpp`) acts as the module manager. It coordinates between the UI, the camera, the mount, and the specific guider instance (Internal or External like PHD2).

### 2. Guider Logic (`Ekos::InternalGuider`)
`InternalGuider` (in `internalguider.cpp`) implements the `GuideInterface`. It manages the state machine (IDLE, CALIBRATING, GUIDING, DITHERING) and handles the timing of exposures and pulses.

### 3. Math Engine (`ISD::cgmath`)
`cgmath` (in `gmath.cpp`) is the "brain" of the guider. It encapsulates:
- **Star Tracking:** Mapping pixels to coordinates.
- **Drift Calculation:** Comparing current position to the lock (target) position.
- **Processing:** Invoking the control algorithms to decide pulse durations.

### 4. Transformation (`Calibration`)
The `Calibration` class (in `calibration.cpp`) stores the results of the calibration process. It contains a rotation matrix and scaling factors to convert pixel shifts $(\Delta x, \Delta y)$ into RA/Dec arcsecond drifts.

### 5. Measurement Algorithms (`GuideAlgorithms`)
Static helpers that perform image analysis:
- **Centroid:** Basic center-of-mass calculation.
- **SEP:** Source Extraction and Photometry for multi-star detection.
- **Star Correspondence:** Tracking stars across frames using relative geometry.

## Data Flow (Guiding Loop)

1. **Capture:** `InternalGuider` requests a frame from the camera.
2. **Analysis:** Upon arrival, the image (`FITSData`) is passed to `cgmath::performProcessing`.
3. **Detection:** `cgmath` calls `findLocalStarPosition` which uses the selected `GuideAlgorithm`.
4. **Coordinate Mapping:** The pixel drift is rotated and scaled into RA/Dec arcseconds via the `Calibration` object.
5. **Pulse Calculation:** The arcsecond drift is processed by a pulse algorithm (e.g., `LinearGuider`, `HysteresisGuider`, or `GPG`).
6. **Execution:** The calculated pulses are emitted as signals, which `InternalGuider` forwards to the mount via INDI.

## Hardware Interaction

### 1. Pulse Guiding (Translation)
The execution of translation corrections $(\Delta x, \Delta y)$ follows this path:
1. **`InternalGuider`**: Emits `newMultiPulse(ra_dir, ra_ms, dec_dir, dec_ms)`.
2. **`Guide`**: Receives the signal and calls `m_Guider->doPulse()`.
3. **`ISD::Guider`**: An INDI-aware class that maps directions to standard INDI properties:
    - `TELESCOPE_TIMED_GUIDE_WE` (West/East)
    - `TELESCOPE_TIMED_GUIDE_NS` (North/South)
4. **INDI Driver**: Receives the `NumberVector` and executes the physical movement on the mount or ST4 port.

### 2. Adaptive Optics (Fast Translation)
For high-frequency Tip-Tilt corrections:
1. **`InternalGuider`**: Processes frames at high speed (50Hz+) and calculates $(\Delta x, \Delta y)$.
2. **`Guide`**: Forwards corrections to the `ISD::AdaptiveOptics` device.
3. **`ISD::AdaptiveOptics`**: Maps corrections to INDI AO properties (e.g., `AO_NS`, `AO_WE`).
4. **AO Hardware**: Rapidly tilts the mirror to stabilize the image before the exposure completes.

### 3. Rotator Integration (Rotation - DONUTS)
To support Alt-Az de-rotation via the DONUTS algorithm, a new control path is required:
1. **`InternalGuider`**: Will emit a new signal `newRotationCorrection(double dTheta)`.
2. **`Guide`**: Will receive this signal and coordinate with the active `ISD::Rotator`.
3. **`ISD::Rotator`**: Will handle the low-frequency "nudge" by updating the `ABS_ROTATOR_ANGLE` property.
4. **INDI Rotator Driver**: Adjusts the physical orientation of the sensor.

## Extensibility

To add a new guiding approach (like full-image correlation):
1. Add a new algorithm enum in `kstars.kcfg`.
2. Implement the measurement logic in a new class (e.g., `DonutsGuider`) or within `GuideAlgorithms`.
3. Update `cgmath::performProcessing` to branch into the new logic when selected.
4. Ensure the output is a $(\Delta x, \Delta y)$ shift and a $\Delta \theta$ rotation, enabling the tiered L1 (AO) and L2 (Rotator) control loops.
