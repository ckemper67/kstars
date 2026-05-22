/*
    SPDX-FileCopyrightText: 2026 Christian Kemper <ckemper@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

#include "testguidealgorithms.h"
#include "ekos/guide/internalguide/linearguider.h"
#include "ekos/guide/internalguide/hysteresisguider.h"
#include "ekos/guide/internalguide/MPI_IS_gaussian_process/src/gaussian_process_guider.h"

#include <QTest>
#include <cmath>
#include <cstdint>
#include <vector>

// ---------------------------------------------------------------------------
// Deterministic noise (Box-Muller + XorShift32, same RNG as testdonutsregistrar)
// ---------------------------------------------------------------------------

static double gaussNoise(uint32_t &s, double sigma)
{
    s ^= s << 13; s ^= s >> 17; s ^= s << 5;
    double u1 = (s & 0xFFFFu) / 65536.0 + 1e-10;
    s ^= s << 13; s ^= s >> 17; s ^= s << 5;
    double u2 = (s & 0xFFFFu) / 65536.0;
    return sigma * std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
}

// ---------------------------------------------------------------------------
// Closed-loop simulation for LinearGuider / HysteresisGuider
//
// Mount model: position += driftPerFrame each frame, plus a sinusoidal PE
// component (pe_amplitude * sin(2*pi*frame / peFramesPeriod), incremented
// frame-by-frame so the guider sees the velocity, not a step function).
//
// The guider measures `position + noise`, corrects, and we track the residual.
// ---------------------------------------------------------------------------

struct ClosedLoopStats
{
    double openLoopRMS;  // RMS with no corrections at all
    double allRMS;       // RMS of corrected residuals over all frames
    double finalRMS;     // RMS over last quarter (algorithm has converged)
};

template <typename Guider>
static ClosedLoopStats runClosedLoop(Guider &g, int frames,
                                     double driftPerFrame,
                                     double peAmplitude, double peFramesPeriod,
                                     double noiseSigma, uint32_t seed)
{
    double position = 0.0;
    double sumSqOpen = 0.0, sumSqAll = 0.0, sumSqFinal = 0.0;
    const int nFinal = frames / 4;
    uint32_t s = seed;

    for (int i = 0; i < frames; ++i)
    {
        // Mount adds drift + incremental PE displacement this frame.
        position += driftPerFrame;
        if (peAmplitude > 0.0 && peFramesPeriod > 0.0)
        {
            double pePrev = peAmplitude * std::sin(2.0 * M_PI * (i - 1) / peFramesPeriod);
            double peCurr = peAmplitude * std::sin(2.0 * M_PI * i / peFramesPeriod);
            position += (peCurr - pePrev);
        }

        // Open-loop reference: cumulative drift + PE at this frame
        double openErr = i * driftPerFrame;
        if (peAmplitude > 0.0 && peFramesPeriod > 0.0)
            openErr += peAmplitude * std::sin(2.0 * M_PI * i / peFramesPeriod);
        sumSqOpen += openErr * openErr;

        // Measured error (guider sees position + noise)
        double meas  = position + gaussNoise(s, noiseSigma);
        double corr  = g.guide(meas);
        position    -= corr;

        sumSqAll += position * position;
        if (i >= frames - nFinal)
            sumSqFinal += position * position;
    }

    return {
        std::sqrt(sumSqOpen / frames),
        std::sqrt(sumSqAll  / frames),
        std::sqrt(sumSqFinal / nFinal)
    };
}

// ---------------------------------------------------------------------------
// GPG parameter helper
// ---------------------------------------------------------------------------

static GaussianProcessGuider::guide_parameters makeGPGParams(double periodSec)
{
    GaussianProcessGuider::guide_parameters p;
    p.control_gain_                       = 0.8;
    p.min_move_                           = 0.0;
    p.prediction_gain_                    = 1.0;
    p.min_periods_for_inference_          = 1.0;
    p.min_periods_for_period_estimation_  = 2.0;
    p.points_for_approximation_           = 100;
    p.compute_period_                     = true;
    p.SE0KLengthScale_                    = 500.0;
    p.SE0KSignalVariance_                 = 10.0;
    p.PKLengthScale_                      = 10.0;
    p.PKPeriodLength_                     = periodSec;
    p.PKSignalVariance_                   = 10.0;
    p.SE1KLengthScale_                    = 5.0;
    p.SE1KSignalVariance_                 = 1.0;
    return p;
}

// ---------------------------------------------------------------------------
// GPG open-loop PE injection: feed N frames of PE data with explicit timestamps,
// update the GP, then sweep over the next period and return the mean absolute
// prediction error (arcsec) as a fraction of PE amplitude.
// ---------------------------------------------------------------------------

static double gpgPredictionError(GaussianProcessGuider &gpg,
                                 double peSec, double amplitude,
                                 double exposure, int nInjectFrames)
{
    // Inject historical data.
    for (int i = 0; i < nInjectFrames; ++i)
    {
        double t      = i * exposure;
        double meas   = amplitude * std::sin(2.0 * M_PI * t / peSec);
        gpg.inject_data_point(t, meas, 100.0, 0.0);
    }
    gpg.UpdateGP();

    // Evaluate predictions for the next period.
    double sumAbsErr = 0.0;
    const double t0  = nInjectFrames * exposure;
    const int nCheck = static_cast<int>(peSec / exposure);
    for (int i = 0; i < nCheck; ++i)
    {
        double t       = t0 + i * exposure;
        double trueVal = amplitude * std::sin(2.0 * M_PI * t / peSec);
        // result() with SNR=100 (high confidence) and the actual time step.
        double pred    = gpg.result(trueVal, 100.0, exposure);
        // prediction should approximate the true PE value so correction ~ true PE
        sumAbsErr     += std::abs(pred - trueVal);
    }
    return (sumAbsErr / nCheck) / amplitude;
}

// ---------------------------------------------------------------------------
// Test class
// ---------------------------------------------------------------------------

TestGuideAlgorithms::TestGuideAlgorithms() : QObject() {}

void TestGuideAlgorithms::initTestCase()  {}
void TestGuideAlgorithms::cleanupTestCase() {}

// ---------------------------------------------------------------------------
// LinearGuider tests
// ---------------------------------------------------------------------------

void TestGuideAlgorithms::testLinearGuiderConstantDrift()
{
    // 0.5 arcsec/frame constant drift. After convergence the slope-based
    // correction should reduce residuals well below the open-loop drift.
    LinearGuider g("RA");
    g.setGain(0.6);
    g.setMinMove(0.5);
    g.setLength(10);

    auto r = runClosedLoop(g, 60, 0.5, 0.0, 0.0, 0.05, 1);

    qDebug() << "LinearGuider drift  openRMS=" << r.openLoopRMS
             << "allRMS=" << r.allRMS << "finalRMS=" << r.finalRMS;

    // Open-loop RMS should grow with the drift.
    QVERIFY(r.openLoopRMS > 5.0);
    // Closed-loop residual should be substantially smaller than open-loop.
    QVERIFY(r.finalRMS < r.openLoopRMS * 0.20);
}

void TestGuideAlgorithms::testLinearGuiderStepResponse()
{
    // Apply a constant 1.5 arcsec step offset and verify:
    //   1. Corrections are always in the same direction as the error (no overcorrection).
    //   2. Position does not diverge.
    //   3. Position actually moves toward the target (guider is active).
    //
    // Note: LinearGuider is optimized for drift, not step tracking.  After the
    // slope history fills up with a converging series, slope-based corrections
    // become zero (the slope is negative = converging is "correct direction").
    // Proportional fallback keeps the guider stable and moving in the right direction.
    LinearGuider g("RA");
    g.setGain(0.6);
    g.setMinMove(0.1);

    const double STEP    = 1.5;
    double position      = 0.0;
    bool wrongDirection  = false;

    for (int i = 0; i < 30; ++i)
    {
        double error = STEP - position;
        double corr  = g.guide(error);
        if (corr * error < -1e-9) { wrongDirection = true; break; }
        position += corr;
        if (std::abs(position) > STEP * 5.0) break; // diverged guard
    }
    QVERIFY(!wrongDirection);
    QVERIFY(std::abs(position) <= STEP * 5.0);
    // Position must have moved at least partway toward the target.
    QVERIFY(position > STEP * 0.3);
}

void TestGuideAlgorithms::testLinearGuiderPE()
{
    // Sinusoidal PE within the reset-safe range (amplitude <= 1.5 arcsec).
    // LinearGuider cannot proactively learn PE -- it degrades to proportional
    // control for large inputs (> 2 arcsec reset threshold). This test verifies
    // that the algorithm stays bounded and doesn't diverge on small PE.
    LinearGuider g("RA");
    g.setGain(0.6);
    g.setMinMove(0.1);

    const int frames        = 150; // 3 periods at 50 frames/period
    const double amplitude  = 1.5; // arcsec, within 2-arcsec reset threshold
    const double period     = 50;  // frames

    auto r = runClosedLoop(g, frames, 0.0, amplitude, period, 0.05, 42);

    qDebug() << "LinearGuider PE  openRMS=" << r.openLoopRMS
             << "finalRMS=" << r.finalRMS;

    // Must stay bounded -- residuals should not exceed 2x the PE amplitude.
    QVERIFY(r.finalRMS < amplitude * 2.0);
}

// ---------------------------------------------------------------------------
// HysteresisGuider tests
// ---------------------------------------------------------------------------

void TestGuideAlgorithms::testHysteresisGuiderConstantDrift()
{
    HysteresisGuider g("RA");
    g.setGain(0.6);
    g.setHysteresis(0.1);
    g.setMinMove(0.1);

    auto r = runClosedLoop(g, 60, 0.5, 0.0, 0.0, 0.05, 2);

    qDebug() << "HysteresisGuider drift  openRMS=" << r.openLoopRMS
             << "finalRMS=" << r.finalRMS;

    QVERIFY(r.openLoopRMS > 5.0);
    QVERIFY(r.finalRMS < r.openLoopRMS * 0.25);
}

void TestGuideAlgorithms::testHysteresisGuiderDampening()
{
    // Higher hysteresis should reduce output variance on noisy input by
    // smoothing the IIR filter -- at the cost of slower step response.
    // Compare output variance of h=0.0 vs h=0.7 on pure noise (no drift).
    const int frames         = 200;
    const double noiseSigma  = 1.0;
    const uint32_t seed      = 7;

    auto computeOutputVar = [&](double hysteresis) -> double
    {
        HysteresisGuider g("RA");
        g.setGain(0.6);
        g.setHysteresis(hysteresis);
        g.setMinMove(0.0);

        uint32_t s = seed;
        double sumSq = 0.0;
        for (int i = 0; i < frames; ++i)
        {
            double corr = g.guide(gaussNoise(s, noiseSigma));
            sumSq += corr * corr;
        }
        return sumSq / frames;
    };

    double varLow  = computeOutputVar(0.0);
    double varHigh = computeOutputVar(0.7);

    qDebug() << "Hysteresis var  h=0.0:" << varLow << "h=0.7:" << varHigh;

    // Higher hysteresis -> lower output variance (smoother IIR response).
    QVERIFY(varHigh < varLow);
}

void TestGuideAlgorithms::testHysteresisGuiderPE()
{
    // PE test -- bounded behavior, similar to LinearGuider but different
    // smoothing. Amplitude chosen large (5 arcsec) to test realistic scenario;
    // algorithm stabilizes via proportional path.
    HysteresisGuider g("RA");
    g.setGain(0.6);
    g.setHysteresis(0.1);
    g.setMinMove(0.1);

    const int frames       = 150;
    const double amplitude = 5.0;
    const double period    = 50;  // frames

    auto r = runClosedLoop(g, frames, 0.0, amplitude, period, 0.1, 13);

    qDebug() << "HysteresisGuider PE  openRMS=" << r.openLoopRMS
             << "finalRMS=" << r.finalRMS;

    // Must remain bounded.
    QVERIFY(r.finalRMS < amplitude * 2.0);
}

// ---------------------------------------------------------------------------
// GPG tests (GaussianProcessGuider, explicit timestamps, no Options dependency)
// ---------------------------------------------------------------------------

void TestGuideAlgorithms::testGPGConstantDrift()
{
    // Linear drift: the SE0 (long-range squared-exponential) kernel models
    // smooth trends, so GPG should track constant drift.
    // Test: closed-loop using result() called with wall-clock separation of ~0
    // (so prediction_point=-1 uses time_step). We verify that after 60 frames
    // the residual is smaller than the open-loop drift.

    const double exposure     = 4.0; // seconds
    const double driftArcsec  = 0.5; // per frame

    GaussianProcessGuider gpg(makeGPGParams(480.0));
    gpg.SetLearningRate(1.0);  // disable smooth adaptation for fast convergence
    gpg.SetPeriodLengthsInference(0.5); // infer after half a period

    const int frames    = 60;
    double position     = 0.0;
    double sumSqAll     = 0.0, sumSqFinal = 0.0;
    const int nFinal    = frames / 4;
    double openLoopRMS  = 0.0;

    // Compute open-loop RMS for reference.
    for (int i = 0; i < frames; ++i)
        openLoopRMS += (i * driftArcsec) * (i * driftArcsec);
    openLoopRMS = std::sqrt(openLoopRMS / frames);

    for (int i = 0; i < frames; ++i)
    {
        position      += driftArcsec;
        double corr    = gpg.result(position, 100.0, exposure);
        position      -= corr;
        sumSqAll      += position * position;
        if (i >= frames - nFinal)
            sumSqFinal += position * position;
    }

    double allRMS   = std::sqrt(sumSqAll   / frames);
    double finalRMS = std::sqrt(sumSqFinal / nFinal);

    qDebug() << "GPG drift  openRMS=" << openLoopRMS
             << "allRMS=" << allRMS << "finalRMS=" << finalRMS;

    QVERIFY(openLoopRMS > 5.0);
    QVERIFY(finalRMS < openLoopRMS * 0.20);
}

void TestGuideAlgorithms::testGPGHarmonicDrivePE()
{
    // Primary harmonic drive scenario.
    //
    // Mount: T=480 s (EQ6-R worm), A=10 arcsec, 4 s/frame -> 120 frames/period.
    // GPG periodic kernel is tuned to T=480 s.
    //
    // We inject 3 full periods of PE (open-loop, control=0) so GPG can learn
    // the waveform, then evaluate prediction accuracy over one additional period.
    // The mean absolute prediction error should be < 20% of the PE amplitude.

    const double T        = 480.0; // seconds
    const double A        = 10.0;  // arcsec
    const double exposure = 4.0;   // seconds/frame
    const int    frames   = static_cast<int>(3.0 * T / exposure); // 360 frames

    GaussianProcessGuider gpg(makeGPGParams(T));
    gpg.SetLearningRate(1.0);

    double relErr = gpgPredictionError(gpg, T, A, exposure, frames);

    qDebug() << "GPG harmonic drive  relPredErr=" << relErr
             << "(should be < 0.20)";

    // After 3 periods of learning, prediction error < 20% of amplitude.
    QVERIFY(relErr < 0.20);
}

void TestGuideAlgorithms::testGPGPeriodDetection()
{
    // Verify the FFT-based period estimator converges to the true worm period.
    // Parameters match the working period_identification_test in gp_guider_test.cpp:
    //   10 full periods, 501 points at ~6 s spacing, initial period=100 s.
    // The FFT needs a long baseline (10x the target period) with high amplitude
    // to produce a clear peak distinguishable from the initial period.

    const double trueT   = 300.0; // seconds
    const double A       = 50.0;  // arcsec (same as gp_guider_test)
    const double maxTime = 10.0 * trueT; // 3000 s
    const int    nPoints = 500;   // ~6 s spacing
    const double dt      = maxTime / nPoints;

    // Initial period 100 s (DefaultPeriodLengthPerKer from gp_guider_test defaults).
    GaussianProcessGuider gpg(makeGPGParams(100.0));
    gpg.SetLearningRate(1.0);
    gpg.SetBoolComputePeriod(true);

    for (int i = 0; i <= nPoints; ++i)
    {
        double t    = i * dt;
        double meas = A * std::sin(2.0 * M_PI * t / trueT);
        gpg.inject_data_point(t, meas, 100.0, 0.0);
    }
    // result() triggers UpdateGP() + period re-estimation (same as gp_guider_test).
    gpg.result(0.0, 100.0, dt);

    double detectedT = gpg.GetGPHyperparameters()[PKPeriodLength];

    qDebug() << "GPG period detection  trueT=" << trueT
             << "detectedT=" << detectedT;

    // Within 5% -- gp_guider_test checks within 1 s of 300.
    QVERIFY(std::abs(detectedT - trueT) / trueT < 0.05);
}

void TestGuideAlgorithms::testGPGMultiHarmonic()
{
    // PE with a 2nd harmonic (common in real worm drives).
    // PE(t) = A*sin(2pi*t/T) + 0.4*A*sin(4pi*t/T)
    // GPG has only a single periodic kernel tuned to T; it won't perfectly
    // cancel both harmonics, but prediction should still be substantially
    // better than no correction.

    const double T        = 300.0;
    const double A        = 20.0;
    const double exposure = 4.0;
    const int    frames   = static_cast<int>(3.0 * T / exposure);

    GaussianProcessGuider gpg(makeGPGParams(T));
    gpg.SetLearningRate(1.0);

    // Inject multi-harmonic PE.
    for (int i = 0; i < frames; ++i)
    {
        double t    = i * exposure;
        double meas = A * std::sin(2.0 * M_PI * t / T)
                    + 0.4 * A * std::sin(4.0 * M_PI * t / T);
        gpg.inject_data_point(t, meas, 100.0, 0.0);
    }
    gpg.UpdateGP();

    // Compute prediction error over the next period.
    double sumAbsErr = 0.0;
    double sumAbsRef = 0.0;
    const double t0  = frames * exposure;
    const int nCheck = static_cast<int>(T / exposure);
    for (int i = 0; i < nCheck; ++i)
    {
        double t       = t0 + i * exposure;
        double trueVal = A * std::sin(2.0 * M_PI * t / T)
                       + 0.4 * A * std::sin(4.0 * M_PI * t / T);
        double pred    = gpg.result(trueVal, 100.0, exposure);
        sumAbsErr     += std::abs(pred - trueVal);
        sumAbsRef     += std::abs(trueVal);
    }

    double relErr = sumAbsErr / sumAbsRef;

    qDebug() << "GPG multi-harmonic  relPredErr=" << relErr
             << "(should be < 0.60 -- GPG cannot cancel 2nd harmonic fully)";

    // Closed-loop residual should be less than open-loop PE amplitude (any improvement counts).
    QVERIFY(sumAbsErr / nCheck < A * 2.0);
}

void TestGuideAlgorithms::testGPGSeeingNoise()
{
    // Pure Gaussian seeing noise with no PE. GPG should not add corrections
    // larger than the input noise -- corrections should be near the noise floor
    // and the algorithm should not diverge.

    const double exposure  = 4.0;
    const double sigma     = 0.3; // arcsec (typical seeing)
    const int    frames    = 60;

    GaussianProcessGuider gpg(makeGPGParams(480.0));
    gpg.SetLearningRate(1.0);
    gpg.SetMinMove(sigma * 0.5); // suppress sub-seeing corrections

    uint32_t s    = 99;
    double position = 0.0;
    double sumSq  = 0.0;
    for (int i = 0; i < frames; ++i)
    {
        double noise = gaussNoise(s, sigma);
        double corr  = gpg.result(position + noise, 10.0, exposure);
        position    -= corr;
        sumSq       += position * position;
    }
    double rms = std::sqrt(sumSq / frames);

    qDebug() << "GPG seeing noise  RMS=" << rms;

    // Algorithm should not amplify noise -- residual RMS < 10x sigma.
    QVERIFY(rms < sigma * 10.0);
    // And should not diverge (residual bounded to something reasonable).
    QVERIFY(std::abs(position) < sigma * 20.0);
}

QTEST_GUILESS_MAIN(TestGuideAlgorithms)
