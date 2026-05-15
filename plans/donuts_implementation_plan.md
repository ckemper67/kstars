# Plan: Implement DONUTS Autoguiding Algorithm

## Introduction
DONUTS is a high-precision autoguiding algorithm designed for wide-field surveys (originally for the Next Generation Transit Survey). Unlike traditional centroid-based guiders, DONUTS uses the entire image (or a large portion of it) to calculate a global shift $(\Delta x, \Delta y)$ between a reference frame and the current frame.

**Reference:** [McCormac et al. (2013), "Donuts: a standalone guiding algorithm for the Next-Generation Transit Survey"](https://arxiv.org/abs/1304.2405)

## Objectives
- Integrate DONUTS into the Ekos Internal Guider.
- Support high-precision guiding even with defocused stars or complex fields.
- Maintain compatibility with the existing calibration and pulse-generation pipeline.

## Implementation Steps

### 1. Configuration & UI
- **File:** `kstars/kstars.kcfg`
- **Task:** Add `GUIDE_DONUTS` (value 6) to the `GuideAlgorithm` entry.
- **Task:** Add any DONUTS-specific options (e.g., `DonutsBorder` to ignore edges).

### 2. Core Logic Implementation
- **New Files:** `kstars/ekos/guide/internalguide/donuts.h/cpp`
- **Dependency:** Integrate `pocketfft` (header-only) for Fourier-domain processing.
- **Class:** `DonutsGuider`
- **Functionality:**
    - **Quadrant Projections:** 
        - Split the image into four quadrants (Top-Left, Top-Right, Bottom-Left, Bottom-Right).
        - Generate X and Y 1D projections for each quadrant (8 profiles total).
    - **Reference Management:** Store the 8 1D profiles of the first frame (Reference Frame).
    - **Shift Measurement (PocketFFT):** 
        - Use `pocketfft` to correlate each of the 8 quadrant profiles.
        - Output 8 local shifts: $(\Delta x_{TL}, \Delta y_{TL}, \Delta x_{TR}, \Delta y_{TR}, \Delta x_{BL}, \Delta y_{BL}, \Delta x_{BR}, \Delta y_{BR})$.
    - **3-DoF Rigid Body Solver:**
        - Solve for Translation $(\Delta x_{trans}, \Delta y_{trans})$ and Field Rotation $(\Delta \theta)$ using a Least Squares fit of the 8 local shift measurements.
        - This approach natively decouples translation from rotation by observing the differential motion between quadrants.

### 3. Tiered Actuator Control Strategy
To support robust Alt-Az guiding, the corrections are split into two primary hardware control loops:
1.  **L1: Fast (Tip-Tilt / AO) - 50-100Hz:**
    - Uses direct $(\Delta x, \Delta y)$ translation components derived from the 4-axis correlation.
    - Provides zero-latency correction for high-frequency atmospheric jitter and mount tracking errors.
2.  **L2: Slow (Physical Rotator) - 0.1-1Hz:**
    - Uses the measured field rotation $(\Delta \theta)$.
    - Nudges the physical mechanical rotator to keep the sensor aligned with the sky, preventing field rotation during long exposures.

### 4. Alt-Az Specific Challenges
- **Rotation Center Offset:** The guider must handle the offset between the guide sensor center and the optical axis (rotator center).
- **Control Synchronization:** Physical rotator "nudges" must be signaled to the AO loop to prevent the high-speed loop from fighting the intentional mechanical adjustment.
- **Zenith Acceleration:** Prediction logic for field rotation rate near the Zenith to supplement the reactive DONUTS measurement.
- **Geometry Normalization:** Handle the non-uniform pixel count in centered diagonal projections to avoid edge artifacts in the FFT.

### 5. Integration into `cgmath`
- **File:** `kstars/ekos/guide/internalguide/gmath.cpp`
- **Method:** `performProcessing`
- **Changes:**
    - Detect if `Options::guideAlgorithm() == GUIDE_DONUTS`.
    - If a reference profile doesn't exist (first frame), initialize it.
    - Call `DonutsGuider::calculateTransform(imageData)`.
    - Synthesize a "virtual star position" based on the translation components $(\Delta x, \Delta y)$.
    - Provide the rotation correction $\Delta \theta$ to the `InternalGuider` to manage the physical rotator loop.

### 6. Calibration Compatibility
- The `CalibrationProcess` moves the mount and expects the "star" to move.
- Since DONUTS measures image shift, it will naturally report the movement $(\Delta x, \Delta y)$ as the image contents shift.
- Ensure the reference profile is captured at the *beginning* of the calibration process (`CAL_START`).

### 7. Optimization & Robustness
- **Pre-processing (per McCormac et al.):**
    - **Median Background Subtraction:** Subtract the global median to ensure a zero-base background for profiles.
    - **Masking:** Ignore bad/hot pixels and clip saturated stars to prevent bias in the correlation peak.
- **FFT-Specific Optimizations:**
    - **Windowing:** Apply a Hann or Hamming window to the 1D profiles before FFT to eliminate edge artifacts and spectral leakage.
    - **Geometry Normalization (All Axes):** 
        - For the X and Y axes, ensure the sums are normalized to account for any non-rectangular ROIs.
        - For the 45° (U) and 135° (V) diagonal projections centered on the frame, divide the raw sums by a "pixel count" vector. This compensates for the varying lengths of diagonal summation lines (which follow a trapezoidal or triangular distribution depending on the sensor aspect ratio).
- **Performance:** Use SIMD for projections and `pocketfft` for efficient O(N log N) correlation.

## Testing Strategy

To ensure the robustness of the 3-DoF stabilization, a three-stage testing strategy will be employed.

### Stage 1: Mathematical Verification (Unit Tests)
**File:** `Tests/ekos/guide/testdonutsguider.cpp`
- **Method:** Load a high-quality "master" FITS image (e.g., `m47_sim_stars.fits`).
- **Synthesis:** Programmatically apply sub-pixel translations $(\Delta x, \Delta y)$ and rotations $(\Delta \theta)$ to the image buffer.
- **Verification:** Feed the transformed buffers to `DonutsGuider` and assert that the recovered transform matches the input within a tolerance of $<0.01$ pixels/degrees.
- **Robustness:** Add Gaussian noise and vary global brightness (cloud flicker) to test correlation stability and SNR reporting.

### Stage 2: Functional Testing (Internal Loop)
- **Method:** Integrate `DonutsGuider` into the `cgmath` loop but bypass hardware.
- **Verification:** Ensure that the "virtual star position" generated from the DONUTS drift correctly drives the existing pulse-calculation algorithms (PID, GPG).
- **Calibration:** Verify that the `CalibrationProcess` can build a valid transformation matrix using DONUTS-measured shifts.

### Stage 3: System-Level Testing (INDI Simulators)
- **Environment:** Use the Ekos UI Test framework (`test_ekos_guide.cpp`).
- **Setup:** A profile with "CCD Simulator", "Telescope Simulator", and "Rotator Simulator".
- **Verification:** 
    - Induce drift by commanding the Telescope Simulator to slew at a non-sidereal rate.
    - Induce rotation by commanding the Rotator Simulator to a known offset.
    - Verify that the L1 (AO) and L2 (Rotator) loops correctly nullify the simulated errors and maintain a stable lock.

## Milestones

To ensure a smooth implementation and validate the math early, development will follow these phased milestones:

### Milestone 1: Mathematical Proof of Concept (Standalone Test)
- **Goal:** Prove the 4-Quadrant FFT approach can accurately recover 3-DoF transformations.
- **Action:** Create `Tests/ekos/guide/standalone_donuts_test.cpp` (No KDE/Qt dependencies).
- **Test:** Use 8 local 1D correlations to recover $(\Delta x, \Delta y, \Delta \theta)$ from a synthetic image.
- **Success:** Accuracy $<0.05$ pixels/degrees.

### Milestone 2: Basic Guiding Integration (Translation Only)
- **Goal:** Integrate DONUTS into the standard Ekos guiding loop for Equatorial mounts.
- **Action:** Add the `GuideAlgorithm` UI option. Connect the X/Y FFT shifts to `cgmath` to generate "virtual star positions".
- **Test:** Use the INDI CCD and Telescope simulators to verify that Ekos can calibrate and guide using the DONUTS translation output.
- **Success:** Stable guiding using full-image correlation instead of star centroids.

### Milestone 3: Full 3-DoF Alt-Az Integration
- **Goal:** Enable the secondary low-frequency rotation control loop.
- **Action:** Calculate $\Delta \theta$ using the U/V projections. Implement the `newRotationCorrection` signal and connect it to the INDI Rotator driver.
- **Test:** Use the simulators to induce field rotation (e.g., tracking near Zenith in Alt-Az mode) and verify the rotator counteracts the drift.
- **Success:** A complete, robust 3-DoF guiding system suitable for Alt-Az observatories.

## Success Criteria
- [ ] Guider can successfully calibrate using DONUTS.
- [ ] Guider can maintain a lock position with a precision comparable to or better than SEP multi-star.
- [ ] Algorithm handles slight transparency changes or cloud flickers via profile normalization.
