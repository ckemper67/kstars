/*
    SPDX-FileCopyrightText: 2026 Christian Kemper <ckemper@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

// Printout-style benchmark comparing LinearGuider, HysteresisGuider, and GPG.
//
// Build:  ninja -C ~/work/kstars/build-mac-tests benchmarkguidealgorithms
// Run:    build-mac-tests/Tests/ekos/guide/benchmarkguidealgorithms.app/Contents/MacOS/benchmarkguidealgorithms

#include "ekos/guide/internalguide/linearguider.h"
#include "ekos/guide/internalguide/hysteresisguider.h"
#include "ekos/guide/internalguide/MPI_IS_gaussian_process/src/gaussian_process_guider.h"

#include <QCoreApplication>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Noise: Box-Muller + XorShift32
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
// PE model
//
// Supports:
//   - Multi-harmonic:       A1*sin(2pi*t/T) + A2*sin(4pi*t/T+phi2) + A3*sin(6pi*t/T+phi3)
//   - Amplitude modulation: base *= (1 + modDepth*sin(2pi*t/modT))
//   - Directional lag:      += hystLag * sign(d(fundamental)/dt)  -- flex-spline hysteresis
//   - OU stochastic:        Ornstein-Uhlenbeck process added to above
// ---------------------------------------------------------------------------

struct PEParams
{
    double T       = 0.0;  // fundamental period (s), 0 = no PE
    double A1      = 0.0;  // fundamental amplitude (arcsec)
    double A2      = 0.0;  // 2nd harmonic amplitude
    double A3      = 0.0;  // 3rd harmonic amplitude
    double phi2    = 0.0;  // 2nd harmonic phase offset (rad)
    double phi3    = 0.0;  // 3rd harmonic phase offset (rad)
    double modT    = 0.0;  // amplitude modulation period (s), 0 = none
    double modDepth= 0.0;  // amplitude modulation depth [0, 1]
    double hystLag = 0.0;  // directional lag amplitude (arcsec), 0 = none
    double ouTau   = 0.0;  // OU correlation time (s), 0 = none
    double ouSigSS = 0.0;  // OU steady-state sigma (arcsec)
};

// Harmonic + envelope + hysteresis component only (no OU, no seed consumed).
static double absHarmonicPE(const PEParams &pe, double t)
{
    if (pe.T <= 0.0) return 0.0;
    double arg = 2.0 * M_PI * t / pe.T;
    double v = pe.A1 * std::sin(arg)
             + pe.A2 * std::sin(2.0 * arg + pe.phi2)
             + pe.A3 * std::sin(3.0 * arg + pe.phi3);
    if (pe.modT > 0.0)
        v *= 1.0 + pe.modDepth * std::sin(2.0 * M_PI * t / pe.modT);
    if (pe.hystLag > 0.0)
    {
        double vel = pe.A1 * (2.0 * M_PI / pe.T) * std::cos(arg);
        v += pe.hystLag * (vel >= 0.0 ? 1.0 : -1.0);
    }
    return v;
}

// Advance OU state by one frame and return the new value.
static double advanceOU(const PEParams &pe, double &ouX, double dt, uint32_t &s)
{
    if (pe.ouTau <= 0.0) return 0.0;
    double alpha    = std::exp(-dt / pe.ouTau);
    double sigStep  = pe.ouSigSS * std::sqrt(1.0 - alpha * alpha);
    ouX = alpha * ouX + sigStep * gaussNoise(s, 1.0);
    return ouX;
}

// ---------------------------------------------------------------------------
// Result type
// ---------------------------------------------------------------------------

struct Stats
{
    std::string          algo;
    double               openRMS;
    double               allRMS;
    double               finalRMS;
    double               reduction;
    std::vector<double>  finalPositions; // last quarter, for detrend
    std::vector<double>  allPositions;   // full history, filled when collectAll=true
};

static double detrendedRMS(const std::vector<double> &v)
{
    int n = static_cast<int>(v.size());
    if (n < 2) return 0.0;
    double mx = (n - 1) / 2.0, sy = 0.0, sxy = 0.0, sxx = 0.0;
    for (int i = 0; i < n; ++i) { sy += v[i]; sxy += (i-mx)*v[i]; sxx += (i-mx)*(i-mx); }
    double slope = sxy / sxx;
    double ic    = sy / n - slope * mx;
    double sq = 0.0;
    for (int i = 0; i < n; ++i) { double r = v[i] - (slope*i+ic); sq += r*r; }
    return std::sqrt(sq / n);
}

// ---------------------------------------------------------------------------
// Closed-loop runners
// ---------------------------------------------------------------------------

using Steps = std::vector<std::pair<int, double>>; // {frame_index, delta_arcsec}

template <typename Guider>
static Stats runCL(const std::string &name, Guider &g,
                   int frames, double exposure,
                   double driftPerFrame, const PEParams &pe,
                   double noiseSigma, uint32_t seed,
                   const Steps &steps = {}, bool collectAll = false)
{
    uint32_t sn = seed, so = seed + 1000000u;
    double pos = 0.0, ouX = 0.0;
    double prevPE = absHarmonicPE(pe, 0.0);
    double sqOpen = 0.0, sqAll = 0.0, sqFinal = 0.0;
    int nFinal = std::max(1, frames / 4);
    std::vector<double> finalPos, allPos;
    if (collectAll) allPos.reserve(frames);

    for (int i = 0; i < frames; ++i)
    {
        for (const auto &[fr, delta] : steps)
            if (fr == i) pos += delta;

        double t      = i * exposure;
        double ou     = advanceOU(pe, ouX, exposure, so);
        double currPE = absHarmonicPE(pe, t) + ou;
        pos          += driftPerFrame + (currPE - prevPE);
        prevPE        = currPE;

        double openPos = i * driftPerFrame + currPE;
        sqOpen += openPos * openPos;

        double meas = pos + gaussNoise(sn, noiseSigma);
        pos -= g.guide(meas);

        sqAll += pos * pos;
        if (i >= frames - nFinal) { sqFinal += pos*pos; finalPos.push_back(pos); }
        if (collectAll) allPos.push_back(pos);
    }

    double oRMS = std::sqrt(sqOpen / frames);
    double aRMS = std::sqrt(sqAll  / frames);
    double fRMS = std::sqrt(sqFinal / nFinal);
    double red  = oRMS > 0.0 ? (1.0 - fRMS / oRMS) * 100.0 : 0.0;
    return {name, oRMS, aRMS, fRMS, red, finalPos, allPos};
}

static Stats runCLGPG(const std::string &name, GaussianProcessGuider &gpg,
                      int frames, double exposure,
                      double driftPerFrame, const PEParams &pe,
                      double noiseSigma, uint32_t seed,
                      const Steps &steps = {}, bool collectAll = false)
{
    uint32_t sn = seed, so = seed + 1000000u;
    double pos = 0.0, ouX = 0.0;
    double prevPE = absHarmonicPE(pe, 0.0);
    double sqOpen = 0.0, sqAll = 0.0, sqFinal = 0.0;
    int nFinal = std::max(1, frames / 4);
    std::vector<double> finalPos, allPos;
    if (collectAll) allPos.reserve(frames);

    for (int i = 0; i < frames; ++i)
    {
        for (const auto &[fr, delta] : steps)
            if (fr == i) pos += delta;

        double t      = i * exposure;
        double ou     = advanceOU(pe, ouX, exposure, so);
        double currPE = absHarmonicPE(pe, t) + ou;
        pos          += driftPerFrame + (currPE - prevPE);
        prevPE        = currPE;

        double openPos = i * driftPerFrame + currPE;
        sqOpen += openPos * openPos;

        double meas = pos + gaussNoise(sn, noiseSigma);
        pos -= gpg.result(meas, 100.0, exposure);

        sqAll += pos * pos;
        if (i >= frames - nFinal) { sqFinal += pos*pos; finalPos.push_back(pos); }
        if (collectAll) allPos.push_back(pos);
    }

    double oRMS = std::sqrt(sqOpen / frames);
    double aRMS = std::sqrt(sqAll  / frames);
    double fRMS = std::sqrt(sqFinal / nFinal);
    double red  = oRMS > 0.0 ? (1.0 - fRMS / oRMS) * 100.0 : 0.0;
    return {name, oRMS, aRMS, fRMS, red, finalPos, allPos};
}

// Segment RMS helper: RMS of positions[a..b) -- assumes a < b <= positions.size()
static double segRMS(const std::vector<double> &v, int a, int b)
{
    double sq = 0.0;
    for (int i = a; i < b; ++i) sq += v[i] * v[i];
    return std::sqrt(sq / (b - a));
}

// ---------------------------------------------------------------------------
// GPG parameter factories
// ---------------------------------------------------------------------------

static GaussianProcessGuider::guide_parameters makeGPGParams(double period,
                                                              bool learnPeriod = false)
{
    GaussianProcessGuider::guide_parameters p;
    p.control_gain_                      = 0.8;
    p.min_move_                          = 0.0;
    p.prediction_gain_                   = 1.0;
    p.min_periods_for_inference_         = 1.0;
    p.min_periods_for_period_estimation_ = 2.0;
    p.points_for_approximation_          = 100;
    p.compute_period_                    = learnPeriod;
    p.SE0KLengthScale_                   = 500.0;
    p.SE0KSignalVariance_                = 10.0;
    p.PKLengthScale_                     = 10.0;
    p.PKPeriodLength_                    = period;
    p.PKSignalVariance_                  = 10.0;
    p.SE1KLengthScale_                   = 5.0;
    p.SE1KSignalVariance_                = 1.0;
    return p;
}

// ---------------------------------------------------------------------------
// Table printing
// ---------------------------------------------------------------------------

static void printHeader(const char *title, bool showDetrend = false)
{
    printf("\n%s\n", title);
    if (showDetrend)
        printf("  %-22s  %10s  %10s  %10s  %10s  %10s\n",
               "Algorithm", "Open RMS", "All RMS", "Final RMS", "Reduction", "DeTrend");
    else
        printf("  %-22s  %10s  %10s  %10s  %10s\n",
               "Algorithm", "Open RMS", "All RMS", "Final RMS", "Reduction");
    printf("  %-22s  %10s  %10s  %10s  %10s%s\n",
           "----------------------",
           "----------", "----------", "----------", "----------",
           showDetrend ? "  ----------" : "");
}

static void printRow(const Stats &r, bool showDetrend = false)
{
    if (showDetrend)
    {
        double dRMS = detrendedRMS(r.finalPositions);
        printf("  %-22s  %9.3f\"  %9.3f\"  %9.3f\"  %9.1f%%  %9.3f\"\n",
               r.algo.c_str(), r.openRMS, r.allRMS, r.finalRMS, r.reduction, dRMS);
    }
    else
    {
        printf("  %-22s  %9.3f\"  %9.3f\"  %9.3f\"  %9.1f%%\n",
               r.algo.c_str(), r.openRMS, r.allRMS, r.finalRMS, r.reduction);
    }
}

// ---------------------------------------------------------------------------
// Scenario runners
// ---------------------------------------------------------------------------

static void runScenario(const char *title,
                        int frames, double exposure,
                        double driftPerFrame,
                        double peA, double peT_sec,
                        double noiseSigma, double gpgPeriod,
                        bool gpgLearn = false)
{
    PEParams pe;
    pe.T  = peT_sec;
    pe.A1 = peA;

    const uint32_t SEED = 42;

    {
        LinearGuider g("RA"); g.setGain(0.7); g.setMinMove(0.1); g.setLength(25);
        auto r = runCL("LinearGuider", g, frames, exposure, driftPerFrame, pe, noiseSigma, SEED);
        printHeader(title);
        printRow(r);
    }
    {
        HysteresisGuider g("RA"); g.setGain(0.6); g.setHysteresis(0.1); g.setMinMove(0.1);
        auto r = runCL("HysteresisGuider", g, frames, exposure, driftPerFrame, pe, noiseSigma, SEED);
        printRow(r);
    }
    {
        GaussianProcessGuider gpg(makeGPGParams(gpgPeriod, gpgLearn));
        gpg.SetLearningRate(1.0);
        auto r = runCLGPG("GPG", gpg, frames, exposure, driftPerFrame, pe, noiseSigma, SEED);
        printRow(r);
    }
}

// H5: step disturbance scenario.
// Single step at stepFrame. Reports 3 windows:
//   Pre-step  : frames 0..stepFrame-1
//   Immediate : frames stepFrame..stepFrame+nWindow-1   (first response)
//   Late      : frames (total-nWindow)..total-1         (settled recovery)
static void runH5Scenario(const char *title,
                          int frames, double exposure,
                          const PEParams &pe, double noiseSigma,
                          double gpgInitPeriod, int stepFrame, double stepDelta,
                          int nWindow = 25)
{
    const uint32_t SEED = 42;
    const Steps steps = {{stepFrame, stepDelta}};

    auto printH5Header = [&]() {
        printf("\n%s\n", title);
        printf("  step: %+.1f\" at frame %d (t=%.0fs)  "
               "Immediate=fr%d-%d  Late=fr%d-%d\n",
               stepDelta, stepFrame, stepFrame * exposure,
               stepFrame, stepFrame + nWindow - 1,
               frames - nWindow, frames - 1);
        printf("  %-22s  %10s  %10s  %10s\n",
               "Algorithm", "Pre-step", "Immediate", "Late(recov)");
        printf("  %-22s  %10s  %10s  %10s\n",
               "----------------------", "----------", "----------", "----------");
    };

    auto printH5Row = [&](const Stats &r) {
        double preRMS  = segRMS(r.allPositions, 0,              stepFrame);
        double immRMS  = segRMS(r.allPositions, stepFrame,      stepFrame + nWindow);
        double lateRMS = segRMS(r.allPositions, frames - nWindow, frames);
        printf("  %-22s  %9.3f\"  %9.3f\"  %9.3f\"\n",
               r.algo.c_str(), preRMS, immRMS, lateRMS);
    };

    {
        LinearGuider g("RA"); g.setGain(0.7); g.setMinMove(0.1); g.setLength(25);
        auto r = runCL("LinearGuider", g, frames, exposure, 0.0, pe,
                       noiseSigma, SEED, steps, true);
        printH5Header();
        printH5Row(r);
    }
    {
        HysteresisGuider g("RA"); g.setGain(0.6); g.setHysteresis(0.1); g.setMinMove(0.1);
        auto r = runCL("HysteresisGuider", g, frames, exposure, 0.0, pe,
                       noiseSigma, SEED, steps, true);
        printH5Row(r);
    }
    {
        GaussianProcessGuider gpg(makeGPGParams(gpgInitPeriod, true));
        gpg.SetLearningRate(1.0);
        auto r = runCLGPG("GPG (learn)", gpg, frames, exposure, 0.0, pe,
                          noiseSigma, SEED, steps, true);
        printH5Row(r);
    }
}

static void runHarmonicScenario(const char *title,
                                int frames, double exposure,
                                double driftPerFrame, const PEParams &pe,
                                double noiseSigma, double gpgInitPeriod,
                                bool showDetrend = false,
                                const Steps &steps = {})
{
    const uint32_t SEED = 42;

    {
        LinearGuider g("RA"); g.setGain(0.7); g.setMinMove(0.1); g.setLength(25);
        auto r = runCL("LinearGuider", g, frames, exposure, driftPerFrame, pe,
                       noiseSigma, SEED, steps);
        printHeader(title, showDetrend);
        printRow(r, showDetrend);
    }
    {
        HysteresisGuider g("RA"); g.setGain(0.6); g.setHysteresis(0.1); g.setMinMove(0.1);
        auto r = runCL("HysteresisGuider", g, frames, exposure, driftPerFrame, pe,
                       noiseSigma, SEED, steps);
        printRow(r, showDetrend);
    }
    {
        // GPG: initial period deliberately mismatched -- learns the true period
        GaussianProcessGuider gpg(makeGPGParams(gpgInitPeriod, true));
        gpg.SetLearningRate(1.0);
        auto r = runCLGPG("GPG (learn)", gpg, frames, exposure, driftPerFrame, pe,
                          noiseSigma, SEED, steps);
        printRow(r, showDetrend);
    }
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

int main(int argc, char *argv[])
{
    QCoreApplication app(argc, argv);

    printf("=== Guide Algorithm RMSE Benchmark ===\n");
    printf("Open RMS:  uncorrected mount error over all frames\n");
    printf("All RMS:   corrected residual, all frames (includes learning transient)\n");
    printf("Final RMS: corrected residual, last quarter (algorithm should have converged)\n");
    printf("Reduction: (1 - Final/Open) x 100%%\n");
    printf("DeTrend:   final-quarter RMS after removing linear drift (H6 only)\n");

    // -----------------------------------------------------------------------
    // Original worm gear scenarios (GPG pre-tuned to true period)
    // -----------------------------------------------------------------------

    runScenario(
        "Scenario 1: Constant drift  drift=0.50\"/frame  frames=60  noise=0.00\"",
        60, 4.0, 0.5, 0.0, 0.0, 0.0, 480.0);

    runScenario(
        "Scenario 2: Fast worm PE  T=100s  A=5.0\"  frames=150 (6 periods)  noise=0.00\"",
        150, 4.0, 0.0, 5.0, 100.0, 0.0, 100.0);

    runScenario(
        "Scenario 3: Realistic worm PE  T=480s  A=10.0\"  frames=360 (3 periods)  noise=0.00\"",
        360, 4.0, 0.0, 10.0, 480.0, 0.0, 480.0);

    runScenario(
        "Scenario 4: Realistic worm PE + seeing  T=480s  A=10.0\"  frames=360  noise=0.30\"",
        360, 4.0, 0.0, 10.0, 480.0, 0.3, 480.0);

    runScenario(
        "Scenario 5: Drift + worm PE + noise  drift=0.20\"/frame  T=480s  A=5.0\"  noise=0.30\"  frames=360",
        360, 4.0, 0.2, 5.0, 480.0, 0.3, 480.0);

    runScenario(
        "Scenario 6: Pure noise  sigma=0.30\"  frames=60  (stability check -- all should worsen)",
        60, 4.0, 0.0, 0.0, 0.0, 0.3, 480.0);

    // -----------------------------------------------------------------------
    // Harmonic drive scenarios (GPG learns period from data)
    // Note: GPG initial period intentionally mismatched; Final RMS reflects
    // performance after learning. All RMS includes the learning transient.
    // -----------------------------------------------------------------------

    printf("\n\n--- Harmonic Drive Scenarios ---\n");
    printf("GPG: compute_period=ON, initial period differs from true period.\n");
    printf("Use Final RMS (last quarter) as the settled performance metric.\n");

    // H1: Classic multi-harmonic, clean
    {
        PEParams pe;
        pe.T = 30.0; pe.A1 = 4.0; pe.A2 = 1.6; pe.A3 = 0.8;
        runHarmonicScenario(
            "H1: Harmonic PE  T=30s  A1=4.0\" A2=1.6\" A3=0.8\"  frames=400 (26 periods)  noise=0.10\"\n"
            "    [clean multi-harmonic; GPG periodic kernel only covers fundamental]",
            400, 2.0, 0.0, pe, 0.10, 100.0);
    }

    // H2: Non-stationary amplitude (amplitude modulated by 10-min envelope)
    {
        PEParams pe;
        pe.T = 30.0; pe.A1 = 4.0; pe.A2 = 1.5; pe.A3 = 0.8;
        pe.modT = 600.0; pe.modDepth = 0.5;
        runHarmonicScenario(
            "H2: Non-stationary PE  T=30s  amplitude x(1+0.5*sin(t/600s))  frames=600  noise=0.15\"\n"
            "    [amplitude modulation breaks GP stationarity; linear regression adapts locally]",
            600, 2.0, 0.0, pe, 0.15, 100.0);
    }

    // H3: Stochastic non-repeatable component (OU process overlaid on harmonics)
    {
        PEParams pe;
        pe.T = 30.0; pe.A1 = 4.0; pe.A2 = 1.2; pe.A3 = 0.6;
        pe.ouTau = 15.0; pe.ouSigSS = 1.0;
        runHarmonicScenario(
            "H3: Stochastic PE  T=30s  + OU(tau=15s, sigma=1.0\")  frames=400  noise=0.20\"\n"
            "    [non-repeatable component; all algorithms converge toward same ceiling]",
            400, 2.0, 0.0, pe, 0.20, 100.0);
    }

    // H4: Directional flex-spline hysteresis
    {
        PEParams pe;
        pe.T = 30.0; pe.A1 = 4.0; pe.A2 = 1.6; pe.A3 = 0.8;
        pe.hystLag = 0.5;
        runHarmonicScenario(
            "H4: Hysteresis lag  T=30s  A1=4.0\" A2=1.6\" A3=0.8\"  hystLag=0.5\"  frames=300  noise=0.15\"\n"
            "    [sign discontinuity at each direction reversal; smooth GP kernels over-smooth]",
            300, 2.0, 0.0, pe, 0.15, 100.0);
    }

    // H5: Step disturbance -- single large step at start of last quarter.
    // 500 frames (2000s): GPG learns T=480s in first 375 frames, then step hits.
    // Pre-step = baseline PE correction quality.
    // Immediate (fr375-399) = first 100s of recovery.
    // Late (fr475-499) = last 100s -- should all converge, but rate differs.
    {
        PEParams pe;
        pe.T = 480.0; pe.A1 = 5.0;
        runH5Scenario(
            "H5: Step disturbance  T=480s  A=5.0\"  step=+5.0\" at fr375 (t=1500s)  frames=500  noise=0.20\"\n"
            "    [GPG has 375 frames to learn period before step; Immediate/Late show recovery speed]",
            500, 4.0, pe, 0.20, 200.0, 375, +5.0);
    }

    // H6: Combined drift + harmonics -- with detrended final RMS
    {
        PEParams pe;
        pe.T = 30.0; pe.A1 = 3.5; pe.A2 = 1.4; pe.A3 = 0.7;
        runHarmonicScenario(
            "H6: Drift + harmonics  drift=0.05\"/frame  T=30s  A1=3.5\" A2=1.4\" A3=0.7\"  frames=600  noise=0.25\"\n"
            "    [DeTrend removes linear drift component so harmonic correction is visible]",
            600, 2.0, 0.05, pe, 0.25, 100.0, true);
    }

    printf("\n");
    return 0;
}
