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
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstring>
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

// PE waveform shapes. Real worm-gear PE is often non-sinusoidal: sawtooth
// (slow buildup, fast snap), triangle (constant velocity reversals), or
// half-rectified (single-sided drive). These expose how well a guider
// handles waveforms whose Fourier spectrum is rich in harmonics.
constexpr int PE_HARMONIC = 0;  // A1 sin(wt) + A2 sin(2wt+phi2) + A3 sin(3wt+phi3)
constexpr int PE_SAWTOOTH = 1;  // linear ramp -A1 to +A1, then snap
constexpr int PE_TRIANGLE = 2;  // -A1 to +A1 over half period, then back
constexpr int PE_HALFRECT = 3;  // A1 * max(0, sin(wt))

struct PEParams
{
    double T       = 0.0;  // fundamental period (s), 0 = no PE
    double Tdot    = 0.0;  // period drift rate (s/s); 0 = constant period.
                           // T(t) = T + Tdot*t. Phase is the time-integral of
                           // 2*pi/T(t), so frequency = 2*pi/T(t) exactly.
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
    int    waveform= PE_HARMONIC;  // PE_HARMONIC | PE_SAWTOOTH | PE_TRIANGLE | PE_HALFRECT
};

// Time-integrated phase: integral_0^t (2*pi/T(tau)) dtau, with T(tau) = T0 + Tdot*tau.
// Linear T drift gives a logarithmic phase. With Tdot=0 this is the standard 2*pi*t/T.
static inline double driftingPhase(double T0, double Tdot, double t)
{
    if (std::abs(Tdot) < 1e-12) return 2.0 * M_PI * t / T0;
    return (2.0 * M_PI / Tdot) * std::log(1.0 + Tdot * t / T0);
}

// Harmonic + envelope + hysteresis component only (no OU, no seed consumed).
static double absHarmonicPE(const PEParams &pe, double t)
{
    if (pe.T <= 0.0) return 0.0;

    double arg = driftingPhase(pe.T, pe.Tdot, t);

    double v = 0.0;
    if (pe.waveform == PE_HARMONIC)
    {
        v = pe.A1 * std::sin(arg)
          + pe.A2 * std::sin(2.0 * arg + pe.phi2)
          + pe.A3 * std::sin(3.0 * arg + pe.phi3);
    }
    else if (pe.waveform == PE_SAWTOOTH)
    {
        // Use the cumulative phase modulo 2*pi to define the sawtooth cycle,
        // so drift propagates correctly through the wrap.
        double phase = arg / (2.0 * M_PI);
        phase -= std::floor(phase);  // [0, 1)
        v = pe.A1 * (2.0 * phase - 1.0);
    }
    else if (pe.waveform == PE_TRIANGLE)
    {
        double phase = arg / (2.0 * M_PI);
        phase -= std::floor(phase);
        v = (phase < 0.5) ? pe.A1 * (4.0 * phase - 1.0)
                          : pe.A1 * (3.0 - 4.0 * phase);
    }
    else if (pe.waveform == PE_HALFRECT)
    {
        v = pe.A1 * std::max(0.0, std::sin(arg));
    }

    if (pe.modT > 0.0)
        v *= 1.0 + pe.modDepth * std::sin(2.0 * M_PI * t / pe.modT);
    if (pe.hystLag > 0.0)
    {
        // Velocity sign of the fundamental sinusoid (for direction-lag); the
        // instantaneous frequency is 2*pi/T(t), so cos(arg)*sign matches the
        // zero-crossings even with drift.
        double vel = pe.A1 * (2.0 * M_PI / (pe.T + pe.Tdot * t)) * std::cos(arg);
        v += pe.hystLag * (vel >= 0.0 ? 1.0 : -1.0);
    }
    return v;
}

// Analytic d/dt of the harmonic component of absHarmonicPE. Uses the
// instantaneous angular frequency 2*pi/T(t) -- with T(t) = T0 + Tdot*t --
// so the rate stays consistent with the drifting phase used by
// absHarmonicPE(). Returns 0 for non-harmonic waveforms (sawtooth/triangle/
// half-rectified) because their slopes are piecewise/non-smooth; the
// 2-mass spring's velocity-damping term is dominated by motor/axis
// inertia in those scenarios, so a zero PE-rate is the safe approximation.
static inline double harmonicPERate(const PEParams &pe, double t)
{
    if (pe.T <= 0.0 || pe.waveform != PE_HARMONIC) return 0.0;
    const double Tt = pe.T + pe.Tdot * t;           // instantaneous period
    if (Tt <= 0.0) return 0.0;
    const double omega = 2.0 * M_PI / Tt;            // instantaneous freq
    const double arg = driftingPhase(pe.T, pe.Tdot, t);
    return pe.A1 * omega       * std::cos(arg)
         + pe.A2 * 2.0 * omega * std::cos(2.0 * arg + pe.phi2)
         + pe.A3 * 3.0 * omega * std::cos(3.0 * arg + pe.phi3);
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
// Measurement model: colored seeing + centroid outliers
// ---------------------------------------------------------------------------
//
// The classic noiseSigma parameter models a per-frame white Gaussian
// observation noise -- a workable stand-in for read+photon noise but a
// poor proxy for atmospheric seeing. Real seeing has correlation time
// in the 0.1-1 s range; over a multi-second guide exposure the slow
// component survives averaging and shows up as a position-domain
// disturbance whose RMS depends on the guider's smoothing length.
// Modelling it as white noise unfairly rewards algorithms with long
// smoothing windows (LinearGuider with high length, etc.).
//
// Centroid outliers handle the other measurement pathology: an
// occasional bad frame (cosmic-ray hit, saturated star, satellite
// trail) that produces a position spike several sigma above the noise
// floor. Different guide algorithms tolerate these differently.
//
// Both are opt-in (defaults leave the existing white-noise-only
// measurement model unchanged).
struct MeasParams
{
    // Colored seeing: Ornstein-Uhlenbeck process at the measurement layer.
    // seeingTau = 0 disables it.
    double seeingTau   = 0.0;   // correlation time (s)
    double seeingSigma = 0.0;   // steady-state sigma (arcsec)

    // Centroid outliers: per-frame Bernoulli + amplified Gaussian.
    // outlierProb = 0 disables it.
    double outlierProb = 0.0;        // probability per frame, [0,1]
    double outlierMult = 4.0;        // outlier amplitude as multiple of noiseSigma
};

// Per-runner measurement-noise state and helpers.
struct MeasState
{
    double seeing = 0.0;        // OU value, advanced once per frame
    uint32_t whiteSeed   = 0;   // white-noise stream
    uint32_t seeingSeed  = 0;   // seeing OU stream
    uint32_t outlierSeed = 0;   // outlier Bernoulli + amplified-Gaussian stream
};

static inline MeasState makeMeasState(uint32_t baseSeed)
{
    return MeasState{ 0.0, baseSeed,
                      baseSeed ^ 0x5EEDu,
                      baseSeed ^ 0xBADCEE7Du };
}

// Advance seeing OU by one full frame and return the new value. Reuses the
// same closed-form transition as advanceOU(); kept inline so each runner
// pays only one transcendental per frame.
static inline double advanceSeeing(MeasState &ms, const MeasParams &mp, double dt)
{
    if (mp.seeingTau <= 0.0 || mp.seeingSigma <= 0.0) { ms.seeing = 0.0; return 0.0; }
    const double a = std::exp(-dt / mp.seeingTau);
    const double sd = mp.seeingSigma * std::sqrt(1.0 - a * a);
    ms.seeing = a * ms.seeing + sd * gaussNoise(ms.seeingSeed, 1.0);
    return ms.seeing;
}

// Generate one frame of measurement noise: white + colored seeing + optional
// outlier. Outlier amplitude is in units of noiseSigma * outlierMult.
static inline double measurementNoise(MeasState &ms, const MeasParams &mp,
                                      double noiseSigma, double dt)
{
    double n = gaussNoise(ms.whiteSeed, noiseSigma);
    n += advanceSeeing(ms, mp, dt);
    if (mp.outlierProb > 0.0)
    {
        // 32-bit uniform in [0,1) from the outlier-seed stream.
        uint32_t r = ms.outlierSeed;
        r ^= r << 13; r ^= r >> 17; r ^= r << 5;
        ms.outlierSeed = r;
        const double u = (r & 0xFFFFu) / 65536.0;
        if (u < mp.outlierProb)
            n += gaussNoise(ms.outlierSeed, noiseSigma * mp.outlierMult);
    }
    return n;
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

// Forward declaration: the OU rate-disturbance step used by the DD plant
// and by the 2-mass / SW plants for the optional wind disturbance. Defined
// near the DD plant block below.
static inline double ddWindStep(double &x, double tau, double sigma_ss,
                                double dt, uint32_t &s);

// ---------------------------------------------------------------------------
// Closed-loop runners
// ---------------------------------------------------------------------------

using Steps = std::vector<std::pair<int, double>>; // {frame_index, delta_arcsec}

template <typename Guider>
static Stats runCL(const std::string &name, Guider &g,
                   int frames, double exposure,
                   double driftPerFrame, const PEParams &pe,
                   double noiseSigma, uint32_t seed,
                   const Steps &steps = {}, bool collectAll = false,
                   const MeasParams &mp = {})
{
    MeasState ms = makeMeasState(seed);
    uint32_t so = seed + 1000000u;
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

        double meas = pos + measurementNoise(ms, mp, noiseSigma, exposure);
        pos -= g.guide(meas);

        // Evaluate at the next measurement time (what the next exposure will see),
        // not the post-correction position at this frame. Predictive guiders pre-position
        // the motor for the next PE value, so post-correction is intentionally offset
        // by -delta_PE; next-meas is the actual photon-arrival residual.
        // OU/noise contributions to the next frame are stochastic and omitted here.
        double nextPE = absHarmonicPE(pe, t + exposure);
        double evalPos = pos + driftPerFrame + (nextPE - absHarmonicPE(pe, t));

        sqAll += evalPos * evalPos;
        if (i >= frames - nFinal) { sqFinal += evalPos * evalPos; finalPos.push_back(evalPos); }
        if (collectAll) allPos.push_back(evalPos);
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
                      const Steps &steps = {}, bool collectAll = false,
                      const MeasParams &mp = {})
{
    MeasState ms = makeMeasState(seed);
    uint32_t so = seed + 1000000u;
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

        double meas = pos + measurementNoise(ms, mp, noiseSigma, exposure);
        pos -= gpg.result(meas, 100.0, exposure);

        // Evaluate at next measurement time (see runCL note above).
        double nextPE = absHarmonicPE(pe, t + exposure);
        double evalPos = pos + driftPerFrame + (nextPE - absHarmonicPE(pe, t));

        sqAll += evalPos * evalPos;
        if (i >= frames - nFinal) { sqFinal += evalPos * evalPos; finalPos.push_back(evalPos); }
        if (collectAll) allPos.push_back(evalPos);
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
// Tuning mode: per-scenario parameter sweep per algorithm
// ---------------------------------------------------------------------------
// When g_tune is true, the benchmark sweeps a small parameter grid for each
// algorithm in each scenario and reports both the default-parameter result
// and the best-found result side by side. Default false (existing fast run).

static bool g_tune = false;

// Sweep grids (kept small to bound runtime; ~49 runs per scenario).
static const std::vector<double> kLinearGain    = { 0.3, 0.5, 0.7, 0.9 };
static const std::vector<int>    kLinearLength  = { 10, 25, 50 };
static const std::vector<double> kHystGain      = { 0.3, 0.5, 0.7, 0.9 };
static const std::vector<double> kHystValue     = { 0.0, 0.2, 0.5 };
static const std::vector<double> kGpgCtrlGain   = { 0.6, 0.8, 1.0 };
static const std::vector<double> kGpgPKLengthSc = { 5.0, 10.0, 20.0 };

static void printTunedHeader(const char *title)
{
    printf("\n%s\n", title);
    printf("  %-22s  %10s  %10s  %s\n",
           "Algorithm", "Default", "Best", "Best params");
    printf("  %-22s  %10s  %10s  %s\n",
           "----------------------", "----------", "----------", "------------");
}

static void printTunedRow(const char *algo, const Stats &def,
                          const Stats &best, const std::string &bestParams)
{
    printf("  %-22s  %9.3f\"  %9.3f\"  %s\n",
           algo, def.finalRMS, best.finalRMS, bestParams.c_str());
}

static inline std::string fmtParams(const char *fmt, ...)
{
    char buf[128];
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(buf, sizeof(buf), fmt, ap);
    va_end(ap);
    return std::string(buf);
}

// ---------------------------------------------------------------------------
// 2-Mass Flexible Shaft Simulation
// ---------------------------------------------------------------------------

struct Sim2MassParams {
    double Jm = 0.02;   // Motor reflected inertia
    double Ja = 0.5;    // Axis inertia
    double Bf = 0.1;    // Motor friction damping
    double Ba = 0.1;    // Axis friction damping
    double Ks = 200.0;  // Gear Stiffness
    double Bs = 2.0;    // Gear Damping
    double Kt = 1.0;
    double stiction = 1.5;
    double coulomb = 0.5;
    double backlash = 0.0; // arcsec; gear play (deadband in motor-axis coupling).
                           // 0 = rigid teeth contact; >0 = no spring force while
                           // |theta_m - theta_a + pe| < backlash.

    // Wind / balance OU rate disturbance applied to the axis. 0 = off.
    double wind_tau   = 20.0;  // s, correlation time
    double wind_sigma = 0.0;   // arcsec/s, steady-state sigma
};

// Compute the gear spring torque accounting for an optional backlash deadband.
// When |relative position| <= backlash, the gear teeth are not in contact and
// no force transmits through the gear (T_spring = 0). Beyond the deadband the
// spring engages with deflection reduced by the deadband width.
static inline double computeSpringTorque(double rel_pos, double rel_vel,
                                         double Ks, double Bs, double backlash)
{
    if (backlash <= 0.0) {
        return Ks * rel_pos + Bs * rel_vel;
    }
    if (rel_pos >  backlash) return Ks * (rel_pos - backlash) + Bs * rel_vel;
    if (rel_pos < -backlash) return Ks * (rel_pos + backlash) + Bs * rel_vel;
    return 0.0;
}

template <typename Guider>
static Stats runCL2Mass(const std::string &name, Guider &g,
                        int frames, double exposure,
                        double driftPerFrame, const PEParams &peParams,
                        double noiseSigma, uint32_t seed,
                        const Sim2MassParams &sim,
                        const Steps &steps = {}, bool collectAll = false,
                        const MeasParams &mp = {})
{
    MeasState ms = makeMeasState(seed);
    uint32_t wind_seed = seed + 4000000u;
    double wind_2m = 0.0;
    double theta_m = 0.0;
    double omega_m = 0.0;
    double theta_a = 0.0;
    double omega_a = 0.0;
    
    double sqOpen = 0.0, sqAll = 0.0, sqFinal = 0.0;
    int nFinal = std::max(1, frames / 4);
    std::vector<double> finalPos, allPos;
    if (collectAll) allPos.reserve(frames);

    double h = 0.002;
    int stepsPerFrame = static_cast<int>(exposure / h);

    for (int i = 0; i < frames; ++i)
    {
        for (const auto &[fr, delta] : steps)
            if (fr == i) theta_a += delta;

        double t_start = i * exposure;
        theta_a += driftPerFrame;
        theta_m += driftPerFrame;

        double openPos = i * driftPerFrame + absHarmonicPE(peParams, t_start);
        sqOpen += openPos * openPos;

        double meas = theta_a + measurementNoise(ms, mp, noiseSigma, exposure);
        double corr = g.guide(meas);
        theta_m -= corr;

        // Multi-rate physical integration loop for compliance/vibrations/friction
        for (int s = 0; s < stepsPerFrame; ++s) {
            double current_t = t_start + s * h;

            double pe = absHarmonicPE(peParams, current_t);
            double pe_rate = harmonicPERate(peParams, current_t);

            double T_spring = computeSpringTorque(theta_m - theta_a + pe,
                                                  omega_m - omega_a + pe_rate,
                                                  sim.Ks, sim.Bs, sim.backlash);

            double omega_m_dot = (-sim.Bf * omega_m - T_spring) / sim.Jm;
            theta_m += h * omega_m;
            omega_m += h * omega_m_dot;

            double T_friction = 0.0;
            if (std::abs(omega_a) < 1e-4) {
                if (std::abs(T_spring) < sim.stiction) {
                    T_friction = T_spring; // stuck
                    omega_a = 0.0;
                } else {
                    T_friction = sim.coulomb * (T_spring > 0.0 ? 1.0 : -1.0);
                }
            } else {
                T_friction = sim.coulomb * (omega_a > 0.0 ? 1.0 : -1.0);
            }

            double omega_a_dot = (T_spring - sim.Ba * omega_a - T_friction) / sim.Ja;
            theta_a += h * omega_a;
            omega_a += h * omega_a_dot;

            // Wind: OU rate disturbance directly on the axis position.
            if (sim.wind_sigma > 0.0)
            {
                ddWindStep(wind_2m, sim.wind_tau, sim.wind_sigma, h, wind_seed);
                theta_a += wind_2m * h;
            }
        }

        sqAll += theta_a * theta_a;
        if (i >= frames - nFinal) {
            sqFinal += theta_a * theta_a;
            finalPos.push_back(theta_a);
        }
        if (collectAll) allPos.push_back(theta_a);
    }

    double oRMS = std::sqrt(sqOpen / frames);
    double aRMS = std::sqrt(sqAll  / frames);
    double fRMS = std::sqrt(sqFinal / nFinal);
    double red  = oRMS > 0.0 ? (1.0 - fRMS / oRMS) * 100.0 : 0.0;
    return {name, oRMS, aRMS, fRMS, red, finalPos, allPos};
}

static Stats runCLGPG2Mass(const std::string &name, GaussianProcessGuider &gpg,
                           int frames, double exposure,
                           double driftPerFrame, const PEParams &peParams,
                           double noiseSigma, uint32_t seed,
                           const Sim2MassParams &sim,
                           const Steps &steps = {}, bool collectAll = false,
                           const MeasParams &mp = {})
{
    MeasState ms = makeMeasState(seed);
    uint32_t wind_seed = seed + 4000000u;
    double wind_2m = 0.0;
    double theta_m = 0.0;
    double omega_m = 0.0;
    double theta_a = 0.0;
    double omega_a = 0.0;
    
    double sqOpen = 0.0, sqAll = 0.0, sqFinal = 0.0;
    int nFinal = std::max(1, frames / 4);
    std::vector<double> finalPos, allPos;
    if (collectAll) allPos.reserve(frames);

    double h = 0.002;
    int stepsPerFrame = static_cast<int>(exposure / h);

    for (int i = 0; i < frames; ++i)
    {
        for (const auto &[fr, delta] : steps)
            if (fr == i) theta_a += delta;

        double t_start = i * exposure;
        theta_a += driftPerFrame;
        theta_m += driftPerFrame;

        double openPos = i * driftPerFrame + absHarmonicPE(peParams, t_start);
        sqOpen += openPos * openPos;

        double meas = theta_a + measurementNoise(ms, mp, noiseSigma, exposure);
        double corr = gpg.result(meas, 100.0, exposure);
        theta_m -= corr;

        // Multi-rate physical integration loop for compliance/vibrations/friction
        for (int s = 0; s < stepsPerFrame; ++s) {
            double current_t = t_start + s * h;
            
            double pe = absHarmonicPE(peParams, current_t);
            double pe_rate = harmonicPERate(peParams, current_t);
            
            double T_spring = computeSpringTorque(theta_m - theta_a + pe,
                                                  omega_m - omega_a + pe_rate,
                                                  sim.Ks, sim.Bs, sim.backlash);
            
            double omega_m_dot = (-sim.Bf * omega_m - T_spring) / sim.Jm;
            theta_m += h * omega_m;
            omega_m += h * omega_m_dot;
            
            double T_friction = 0.0;
            if (std::abs(omega_a) < 1e-4) {
                if (std::abs(T_spring) < sim.stiction) {
                    T_friction = T_spring; // stuck
                    omega_a = 0.0;
                } else {
                    T_friction = sim.coulomb * (T_spring > 0.0 ? 1.0 : -1.0);
                }
            } else {
                T_friction = sim.coulomb * (omega_a > 0.0 ? 1.0 : -1.0);
            }
            
            double omega_a_dot = (T_spring - sim.Ba * omega_a - T_friction) / sim.Ja;
            theta_a += h * omega_a;
            omega_a += h * omega_a_dot;

            // Wind: OU rate disturbance directly on the axis position.
            if (sim.wind_sigma > 0.0)
            {
                ddWindStep(wind_2m, sim.wind_tau, sim.wind_sigma, h, wind_seed);
                theta_a += wind_2m * h;
            }
        }

        sqAll += theta_a * theta_a;
        if (i >= frames - nFinal) {
            sqFinal += theta_a * theta_a;
            finalPos.push_back(theta_a);
        }
        if (collectAll) allPos.push_back(theta_a);
    }

    double oRMS = std::sqrt(sqOpen / frames);
    double aRMS = std::sqrt(sqAll  / frames);
    double fRMS = std::sqrt(sqFinal / nFinal);
    double red  = oRMS > 0.0 ? (1.0 - fRMS / oRMS) * 100.0 : 0.0;
    return {name, oRMS, aRMS, fRMS, red, finalPos, allPos};
}

// ---------------------------------------------------------------------------
// Direct-Drive Mount Simulation
// ---------------------------------------------------------------------------
//
// Models a direct-drive (DD) mount: motor rigidly coupled to the axis with
// no gearbox, so no torsional compliance, no backlash, no worm PE. The
// dominant disturbances are cogging/torque ripple (a torque, integrated to
// position by the inner servo), wind gusts (modeled as an Ornstein-Uhlenbeck
// rate disturbance), and the inner velocity servo's finite tracking
// bandwidth.
//
// State: theta_a (arcsec) and axis_rate (arcsec/s). The guider issues a
// position pulse interpreted by the inner servo as a commanded rate
// (rate_cmd = -correction / exposure over the upcoming frame). The actual
// axis rate slews toward rate_cmd with first-order time constant tau_servo.
// Cogging enters as a rate disturbance (faithful: cogging is a torque, the
// servo integrates it to position). Wind enters as a rate disturbance via
// an OU process.

struct SimDirectDriveParams {
    // Inner closed-loop velocity servo seen from the guider:
    double tau_servo = 0.3;     // s, first-order rate-tracking time constant.

    // Cogging / torque ripple, injected as a rate disturbance (arcsec/s).
    // Period is in seconds (mount-internal, not tied to motor angle here --
    // the strain-wave follow-up will index by motor angle).
    double ripple_T  = 30.0;    // s, fundamental period
    double ripple_A1 = 0.0;     // arcsec/s, fundamental rate amplitude
    double ripple_A2 = 0.0;     // arcsec/s, 2nd harmonic
    double ripple_A3 = 0.0;     // arcsec/s, 3rd harmonic

    // Wind / balance, OU rate disturbance:
    double wind_tau   = 20.0;   // s, OU correlation time
    double wind_sigma = 0.0;    // arcsec/s, OU steady-state sigma

    // Optional under-damped structural mode (mount/OTA cantilever), aliased
    // into the band at guide cadence. Off by default.
    double mode_hz   = 0.0;
    double mode_zeta = 0.05;
    double mode_A    = 0.0;     // arcsec, displacement amplitude

    // Inner servo noise floor (arcsec/s rms, added once per micro-step).
    double servo_noise = 0.0;

    // Pulse-train rate model. If > 0, the guider's position correction is
    // executed as a fixed-rate command of magnitude pulse_rate (arcsec/s)
    // for duration |corr|/pulse_rate, then the rate command returns to 0
    // for the remainder of the exposure. This matches real direct-drive
    // mounts (typical guide-pulse rate is ~0.5x sidereal, so even a 2"
    // correction completes in ~0.27 s of a 4 s exposure -- nothing like
    // the "smear -corr/exposure across the whole frame" behavior the
    // legacy model assumes). If 0 (the default), the legacy "smear"
    // interpretation is used.
    double pulse_rate = 0.0;        // arcsec/s; 0 = legacy smear behavior
};

static inline double ddRipple(const SimDirectDriveParams &dd, double t)
{
    if (dd.ripple_T <= 0.0) return 0.0;
    const double w = 2.0 * M_PI / dd.ripple_T;
    return dd.ripple_A1 * std::sin(w * t)
         + dd.ripple_A2 * std::sin(2.0 * w * t)
         + dd.ripple_A3 * std::sin(3.0 * w * t);
}

static inline double ddMode(const SimDirectDriveParams &dd, double t)
{
    if (dd.mode_hz <= 0.0 || dd.mode_A <= 0.0) return 0.0;
    const double w = 2.0 * M_PI * dd.mode_hz;
    return dd.mode_A * std::exp(-dd.mode_zeta * w * t) * std::sin(w * t);
}

// OU step: dx/dt = -x/tau + sigma_ss * sqrt(2/tau) * dW
static inline double ddWindStep(double &x, double tau, double sigma_ss,
                                double dt, uint32_t &s)
{
    if (tau <= 0.0 || sigma_ss <= 0.0) { x = 0.0; return 0.0; }
    const double a = std::exp(-dt / tau);
    const double sd = sigma_ss * std::sqrt(1.0 - a * a);
    x = a * x + sd * (gaussNoise(s, 1.0));
    return x;
}

template <typename Guider>
static Stats runCLDirectDrive(const std::string &name, Guider &g,
                              int frames, double exposure,
                              double driftPerFrame,
                              double noiseSigma, uint32_t seed,
                              const SimDirectDriveParams &dd,
                              const Steps &steps = {}, bool collectAll = false,
                              const MeasParams &mp = {})
{
    MeasState ms = makeMeasState(seed);
    uint32_t sw = seed + 2000000u, sv = seed + 3000000u;
    double theta_a = 0.0;
    double axis_rate = 0.0;
    double wind = 0.0;
    double open_pos = 0.0;  // disturbance accumulation, no servo correction

    double sqOpen = 0.0, sqAll = 0.0, sqFinal = 0.0;
    int nFinal = std::max(1, frames / 4);
    std::vector<double> finalPos, allPos;
    if (collectAll) allPos.reserve(frames);

    const double h = 0.002;
    const int stepsPerFrame = static_cast<int>(exposure / h);

    for (int i = 0; i < frames; ++i)
    {
        for (const auto &[fr, delta] : steps)
            if (fr == i) { theta_a += delta; open_pos += delta; }

        const double t_start = i * exposure;
        theta_a += driftPerFrame;
        open_pos += driftPerFrame;

        sqOpen += open_pos * open_pos;

        const double meas = theta_a + measurementNoise(ms, mp, noiseSigma, exposure);
        const double corr = g.guide(meas);

        // Interpret the guider's position pulse. The legacy model
        // (pulse_rate == 0) smears the entire correction across the whole
        // exposure: rate_cmd = -corr/exposure. The pulse-train model
        // (pulse_rate > 0) executes the correction at a fixed rate
        // -sign(corr)*pulse_rate for duration |corr|/pulse_rate, then the
        // command returns to 0 for the rest of the exposure. Real DD mounts
        // match the latter (typical guide rate ~0.5x sidereal -> even a 2"
        // correction completes in ~0.3 s of a 4 s frame).
        const double cmd_legacy   = -corr / exposure;
        const double cmd_pulse    = (corr >= 0.0) ? -dd.pulse_rate : dd.pulse_rate;
        const double pulse_dur    = (dd.pulse_rate > 0.0)
                                    ? std::min(std::abs(corr) / dd.pulse_rate, exposure)
                                    : exposure;

        for (int s = 0; s < stepsPerFrame; ++s) {
            const double current_t = t_start + s * h;
            const double active = (dd.pulse_rate > 0.0)
                                  ? ((current_t - t_start) < pulse_dur ? cmd_pulse : 0.0)
                                  : cmd_legacy;
            const double rate_cmd = active;

            // Inner servo: first-order rate tracking.
            const double tau = std::max(dd.tau_servo, 1e-6);
            axis_rate += (rate_cmd - axis_rate) * (h / tau);

            // Wind OU disturbance (slow rate term).
            ddWindStep(wind, dd.wind_tau, dd.wind_sigma, h, sw);

            // Cogging ripple, structural mode, servo noise.
            const double ripple = ddRipple(dd, current_t);
            const double mode   = ddMode(dd, current_t);
            const double vnoise = (dd.servo_noise > 0.0)
                                  ? gaussNoise(sv, dd.servo_noise) : 0.0;

            theta_a += (axis_rate + ripple + wind + mode + vnoise) * h;
            // Open-loop reference: same disturbances, no servo correction.
            // Reuses the same wind/noise stream so open and closed see the
            // same disturbance realization (apples-to-apples reduction).
            open_pos += (ripple + wind + mode) * h;
        }

        sqAll += theta_a * theta_a;
        if (i >= frames - nFinal) {
            sqFinal += theta_a * theta_a;
            finalPos.push_back(theta_a);
        }
        if (collectAll) allPos.push_back(theta_a);
    }

    const double oRMS = std::sqrt(sqOpen / frames);
    const double aRMS = std::sqrt(sqAll  / frames);
    const double fRMS = std::sqrt(sqFinal / nFinal);
    const double red  = oRMS > 0.0 ? (1.0 - fRMS / oRMS) * 100.0 : 0.0;
    return {name, oRMS, aRMS, fRMS, red, finalPos, allPos};
}

static Stats runCLGPGDirectDrive(const std::string &name,
                                 GaussianProcessGuider &gpg,
                                 int frames, double exposure,
                                 double driftPerFrame,
                                 double noiseSigma, uint32_t seed,
                                 const SimDirectDriveParams &dd,
                                 const Steps &steps = {}, bool collectAll = false,
                                 const MeasParams &mp = {})
{
    MeasState ms = makeMeasState(seed);
    uint32_t sw = seed + 2000000u, sv = seed + 3000000u;
    double theta_a = 0.0;
    double axis_rate = 0.0;
    double wind = 0.0;
    double open_pos = 0.0;

    double sqOpen = 0.0, sqAll = 0.0, sqFinal = 0.0;
    int nFinal = std::max(1, frames / 4);
    std::vector<double> finalPos, allPos;
    if (collectAll) allPos.reserve(frames);

    const double h = 0.002;
    const int stepsPerFrame = static_cast<int>(exposure / h);

    for (int i = 0; i < frames; ++i)
    {
        for (const auto &[fr, delta] : steps)
            if (fr == i) { theta_a += delta; open_pos += delta; }

        const double t_start = i * exposure;
        theta_a += driftPerFrame;
        open_pos += driftPerFrame;

        sqOpen += open_pos * open_pos;

        const double meas = theta_a + measurementNoise(ms, mp, noiseSigma, exposure);
        const double corr = gpg.result(meas, 100.0, exposure);

        const double rate_cmd = -corr / exposure;

        for (int s = 0; s < stepsPerFrame; ++s) {
            const double current_t = t_start + s * h;

            const double tau = std::max(dd.tau_servo, 1e-6);
            axis_rate += (rate_cmd - axis_rate) * (h / tau);

            ddWindStep(wind, dd.wind_tau, dd.wind_sigma, h, sw);

            const double ripple = ddRipple(dd, current_t);
            const double mode   = ddMode(dd, current_t);
            const double vnoise = (dd.servo_noise > 0.0)
                                  ? gaussNoise(sv, dd.servo_noise) : 0.0;

            theta_a += (axis_rate + ripple + wind + mode + vnoise) * h;
            open_pos += (ripple + wind + mode) * h;
        }

        sqAll += theta_a * theta_a;
        if (i >= frames - nFinal) {
            sqFinal += theta_a * theta_a;
            finalPos.push_back(theta_a);
        }
        if (collectAll) allPos.push_back(theta_a);
    }

    const double oRMS = std::sqrt(sqOpen / frames);
    const double aRMS = std::sqrt(sqAll  / frames);
    const double fRMS = std::sqrt(sqFinal / nFinal);
    const double red  = oRMS > 0.0 ? (1.0 - fRMS / oRMS) * 100.0 : 0.0;
    return {name, oRMS, aRMS, fRMS, red, finalPos, allPos};
}

// ---------------------------------------------------------------------------
// Strain-Wave (Harmonic Drive) Mount Simulation
// ---------------------------------------------------------------------------
//
// Models a strain-wave drive mount: ZWO HEM/AM3/AM5, iOptron HEM, etc.
// Three SW-specific features that the 2-mass model gets wrong:
//
//   1. Kinematic error (KE) is motor-angle-indexed, not time-indexed. The
//      wave generator is a 2-lobe elliptical cam so the dominant KE is at
//      2x per motor revolution. T_KE_seconds below is a *calibration at
//      sidereal rate*, NOT a wall-clock period. During dithers/large
//      corrections, KE phase advances proportionally to motor travel; if
//      you replace this with a t-based PE the dither behavior breaks
//      silently.
//   2. Torsional stiffness is torque-dependent (soft-zone). Low stiffness
//      near zero load, high stiffness under load -- visible in real guide
//      logs as "stiffens after a few corrections."
//   3. Hysteresis is rate-independent at direction reversals (Dahl-style
//      one-state model). Captures the torque-vs-angle loop that a
//      first-order time lag would miss.
//
// No deadband backlash -- real SW drives are near-zero-backlash by design.
// Tooth-skip is modeled via the existing Steps mechanism if needed.

struct SimStrainWaveParams {
    // Motor + axis (rigid motor side; axis carries the compliance / hysteresis):
    double Jm = 0.02;
    double Ja = 0.5;
    double Bf = 0.1;
    double Ba = 0.1;
    double Kt = 1.0;
    double stiction = 0.5;     // preload friction is meaningful in SW drives
    double coulomb = 0.2;

    // KE: motor-angle-indexed. T_KE is calibration at sidereal rate (see
    // file-header comment). KE_arcsec_per_motor_rev is what determines the
    // angle-to-phase scaling; we set it so a full KE cycle at sidereal rate
    // (15.04 arcsec/s) takes T_KE seconds.
    double T_KE = 60.0;        // s, KE period as seen at sidereal rate
    double A0 = 0.0;           // arcsec, bearing-eccentricity sub-harmonic
    int    N_bearing = 30;     // bearing-eccentricity period (motor revs)
    double A1 = 3.0;           // arcsec, fundamental at 2x motor rev
    double A2 = 0.8;           // arcsec, 4x motor rev
    double A3 = 0.2;           // arcsec, 6x motor rev
    double phi2 = 0.0;
    double phi3 = 0.0;

    // Soft-zone torque-dependent stiffness:
    double Ks_low  = 50.0;
    double Ks_high = 800.0;
    double T_knee  = 1.0;
    double Bs      = 4.0;

    // Dahl hysteresis (one state):
    double h_max = 0.2;        // arcsec, saturation displacement
    double Kh    = 10.0;       // torque-domain gain on h

    // Wind / balance OU rate disturbance applied to the axis. 0 = off.
    double wind_tau   = 20.0;  // s, correlation time
    double wind_sigma = 0.0;   // arcsec/s, steady-state sigma
};

// Motor-angle-indexed kinematic error. phi is the wave-generator angle in
// radians (theta_m mapped to motor-rev phase). Sidereal-rate calibration:
// at constant sidereal motion (15.04 arcsec/s) one motor rev takes
// T_KE seconds of guiding time, so KE_arcsec_per_motor_rev = 15.04 * T_KE.
static inline double swKE(const SimStrainWaveParams &sw, double theta_m_arcsec)
{
    constexpr double SIDEREAL = 15.04108; // arcsec/s
    const double KE_arcsec_per_motor_rev = SIDEREAL * sw.T_KE;
    if (KE_arcsec_per_motor_rev <= 0.0) return 0.0;
    const double phi = 2.0 * M_PI * theta_m_arcsec / KE_arcsec_per_motor_rev;
    double v = sw.A1 * std::sin(2.0 * phi)
             + sw.A2 * std::sin(4.0 * phi + sw.phi2)
             + sw.A3 * std::sin(6.0 * phi + sw.phi3);
    if (sw.A0 > 0.0 && sw.N_bearing > 0)
        v += sw.A0 * std::sin(phi / static_cast<double>(sw.N_bearing));
    return v;
}

// Torque-dependent stiffness (lagged by one micro-step to avoid an implicit
// solve; one-step phase error is benign because tanh saturates slowly).
static inline double swKsEff(const SimStrainWaveParams &sw, double T_prev)
{
    const double t = std::abs(T_prev) / std::max(sw.T_knee, 1e-9);
    return sw.Ks_low + (sw.Ks_high - sw.Ks_low) * std::tanh(t);
}

// Dahl hysteresis update. Drives h toward sign(d_theta) * h_max with rate
// proportional to |d_theta|; rate-independent, saturates at +/-h_max.
static inline void swDahlUpdate(double &h, double d_theta, double h_max)
{
    if (h_max <= 0.0 || d_theta == 0.0) return;
    const double s = (d_theta > 0.0) ? 1.0 : -1.0;
    h += (h_max - s * h) * s * std::abs(d_theta) / h_max;
    if (h >  h_max) h =  h_max;
    if (h < -h_max) h = -h_max;
}

template <typename Guider>
static Stats runCLStrainWave(const std::string &name, Guider &g,
                             int frames, double exposure,
                             double driftPerFrame,
                             double noiseSigma, uint32_t seed,
                             const SimStrainWaveParams &sw,
                             const Steps &steps = {}, bool collectAll = false,
                             const MeasParams &mp = {})
{
    MeasState ms = makeMeasState(seed);
    uint32_t wind_seed = seed + 4000000u;
    double wind_2m = 0.0;
    double theta_m = 0.0;
    double omega_m = 0.0;
    double theta_a = 0.0;
    double omega_a = 0.0;
    double h_dahl  = 0.0;
    double T_prev  = 0.0;
    double open_pos = 0.0;

    double sqOpen = 0.0, sqAll = 0.0, sqFinal = 0.0;
    int nFinal = std::max(1, frames / 4);
    std::vector<double> finalPos, allPos;
    if (collectAll) allPos.reserve(frames);

    const double h_dt = 0.002;
    const int stepsPerFrame = static_cast<int>(exposure / h_dt);

    // Motor advances at sidereal even when residuals are near zero. KE is
    // indexed against absolute motor angle: motor_abs = i*sidereal*exposure
    // + theta_m (residual). The residual contribution is small vs sidereal
    // (~1" vs ~60"/frame at 4s exposure) so the within-frame variation of
    // motor_abs is dominated by sidereal; we use the frame-start value as
    // the KE phase for the whole inner loop (KE cycles cleanly).
    constexpr double SIDEREAL = 15.04108;  // arcsec/s

    for (int i = 0; i < frames; ++i)
    {
        for (const auto &[fr, delta] : steps)
            if (fr == i) { theta_a += delta; open_pos += delta; }

        theta_a += driftPerFrame;
        theta_m += driftPerFrame;
        open_pos += driftPerFrame;

        // Absolute motor angle for KE indexing: sidereal accumulation +
        // closed-loop residual. open_motor is the uncorrected analogue.
        const double motor_abs  = i * SIDEREAL * exposure + theta_m;
        const double open_motor = i * SIDEREAL * exposure;
        const double frame_KE = swKE(sw, motor_abs);
        // Open-loop reference: an uncorrected mount would lag the motor by
        // KE(open_motor); residual axis = -KE(open_motor). Drift adds on top.
        const double open_residual = open_pos - swKE(sw, open_motor);
        sqOpen += open_residual * open_residual;

        const double meas = theta_a + measurementNoise(ms, mp, noiseSigma, exposure);
        const double corr = g.guide(meas);
        theta_m -= corr;

        for (int s = 0; s < stepsPerFrame; ++s) {
            const double rel_pos = theta_m - theta_a + frame_KE;
            const double rel_vel = omega_m - omega_a;
            const double Ks_eff = swKsEff(sw, T_prev);
            const double T_spring = Ks_eff * rel_pos + sw.Bs * rel_vel;
            const double T_total  = T_spring + sw.Kh * h_dahl;

            // Motor side: rigid, low friction.
            const double omega_m_dot = (-sw.Bf * omega_m - T_total) / sw.Jm;
            theta_m += h_dt * omega_m;
            omega_m += h_dt * omega_m_dot;

            // Axis side with stiction/Coulomb (preload friction matters in SW):
            double T_friction = 0.0;
            if (std::abs(omega_a) < 1e-4) {
                if (std::abs(T_total) < sw.stiction) {
                    T_friction = T_total;
                    omega_a = 0.0;
                } else {
                    T_friction = sw.coulomb * (T_total > 0.0 ? 1.0 : -1.0);
                }
            } else {
                T_friction = sw.coulomb * (omega_a > 0.0 ? 1.0 : -1.0);
            }

            const double omega_a_dot = (T_total - sw.Ba * omega_a - T_friction) / sw.Ja;
            const double d_theta_a = h_dt * omega_a;
            theta_a += d_theta_a;
            omega_a += h_dt * omega_a_dot;

            // Dahl hysteresis driven by axis-side displacement (the side
            // that flexes against load).
            swDahlUpdate(h_dahl, d_theta_a, sw.h_max);

            T_prev = T_spring;

            // Wind: OU rate disturbance directly on the axis position.
            if (sw.wind_sigma > 0.0)
            {
                ddWindStep(wind_2m, sw.wind_tau, sw.wind_sigma, h_dt, wind_seed);
                theta_a += wind_2m * h_dt;
            }
        }

        sqAll += theta_a * theta_a;
        if (i >= frames - nFinal) {
            sqFinal += theta_a * theta_a;
            finalPos.push_back(theta_a);
        }
        if (collectAll) allPos.push_back(theta_a);
    }

    const double oRMS = std::sqrt(sqOpen / frames);
    const double aRMS = std::sqrt(sqAll  / frames);
    const double fRMS = std::sqrt(sqFinal / nFinal);
    const double red  = oRMS > 0.0 ? (1.0 - fRMS / oRMS) * 100.0 : 0.0;
    return {name, oRMS, aRMS, fRMS, red, finalPos, allPos};
}

static Stats runCLGPGStrainWave(const std::string &name,
                                GaussianProcessGuider &gpg,
                                int frames, double exposure,
                                double driftPerFrame,
                                double noiseSigma, uint32_t seed,
                                const SimStrainWaveParams &sw,
                                const Steps &steps = {}, bool collectAll = false,
                                const MeasParams &mp = {})
{
    MeasState ms = makeMeasState(seed);
    uint32_t wind_seed = seed + 4000000u;
    double wind_2m = 0.0;
    double theta_m = 0.0;
    double omega_m = 0.0;
    double theta_a = 0.0;
    double omega_a = 0.0;
    double h_dahl  = 0.0;
    double T_prev  = 0.0;
    double open_pos = 0.0;

    double sqOpen = 0.0, sqAll = 0.0, sqFinal = 0.0;
    int nFinal = std::max(1, frames / 4);
    std::vector<double> finalPos, allPos;
    if (collectAll) allPos.reserve(frames);

    const double h_dt = 0.002;
    const int stepsPerFrame = static_cast<int>(exposure / h_dt);

    // See runCLStrainWave for sidereal-indexing rationale.
    constexpr double SIDEREAL = 15.04108;

    for (int i = 0; i < frames; ++i)
    {
        for (const auto &[fr, delta] : steps)
            if (fr == i) { theta_a += delta; open_pos += delta; }

        theta_a += driftPerFrame;
        theta_m += driftPerFrame;
        open_pos += driftPerFrame;

        const double motor_abs  = i * SIDEREAL * exposure + theta_m;
        const double open_motor = i * SIDEREAL * exposure;
        const double frame_KE = swKE(sw, motor_abs);
        const double open_residual = open_pos - swKE(sw, open_motor);
        sqOpen += open_residual * open_residual;

        const double meas = theta_a + measurementNoise(ms, mp, noiseSigma, exposure);
        const double corr = gpg.result(meas, 100.0, exposure);
        theta_m -= corr;

        for (int s = 0; s < stepsPerFrame; ++s) {
            const double rel_pos = theta_m - theta_a + frame_KE;
            const double rel_vel = omega_m - omega_a;
            const double Ks_eff = swKsEff(sw, T_prev);
            const double T_spring = Ks_eff * rel_pos + sw.Bs * rel_vel;
            const double T_total  = T_spring + sw.Kh * h_dahl;

            const double omega_m_dot = (-sw.Bf * omega_m - T_total) / sw.Jm;
            theta_m += h_dt * omega_m;
            omega_m += h_dt * omega_m_dot;

            double T_friction = 0.0;
            if (std::abs(omega_a) < 1e-4) {
                if (std::abs(T_total) < sw.stiction) {
                    T_friction = T_total;
                    omega_a = 0.0;
                } else {
                    T_friction = sw.coulomb * (T_total > 0.0 ? 1.0 : -1.0);
                }
            } else {
                T_friction = sw.coulomb * (omega_a > 0.0 ? 1.0 : -1.0);
            }

            const double omega_a_dot = (T_total - sw.Ba * omega_a - T_friction) / sw.Ja;
            const double d_theta_a = h_dt * omega_a;
            theta_a += d_theta_a;
            omega_a += h_dt * omega_a_dot;

            swDahlUpdate(h_dahl, d_theta_a, sw.h_max);

            T_prev = T_spring;

            // Wind: OU rate disturbance directly on the axis position.
            if (sw.wind_sigma > 0.0)
            {
                ddWindStep(wind_2m, sw.wind_tau, sw.wind_sigma, h_dt, wind_seed);
                theta_a += wind_2m * h_dt;
            }
        }

        sqAll += theta_a * theta_a;
        if (i >= frames - nFinal) {
            sqFinal += theta_a * theta_a;
            finalPos.push_back(theta_a);
        }
        if (collectAll) allPos.push_back(theta_a);
    }

    const double oRMS = std::sqrt(sqOpen / frames);
    const double aRMS = std::sqrt(sqAll  / frames);
    const double fRMS = std::sqrt(sqFinal / nFinal);
    const double red  = oRMS > 0.0 ? (1.0 - fRMS / oRMS) * 100.0 : 0.0;
    return {name, oRMS, aRMS, fRMS, red, finalPos, allPos};
}

static void runComplianceScenario(const char *title,
                                  int frames, double exposure,
                                  double driftPerFrame, const PEParams &pe,
                                  double noiseSigma, double gpgInitPeriod,
                                  const Sim2MassParams &sim,
                                  bool gpgLearn = true,
                                  bool declareCompliance = false)
{
    const uint32_t SEED = 42;

    if (g_tune)
    {
        printTunedHeader(title);
        // LinearGuider
        {
            LinearGuider gd("RA"); gd.setGain(0.7); gd.setMinMove(0.1); gd.setLength(25);
            Stats sDef = runCL2Mass("LinearGuider", gd, frames, exposure, driftPerFrame, pe, noiseSigma, SEED, sim);
            Stats sBest = sDef; std::string bestP = "gain=0.7 len=25";
            for (double gain : kLinearGain) for (int length : kLinearLength) {
                LinearGuider g("RA"); g.setGain(gain); g.setMinMove(0.1); g.setLength(length);
                Stats s = runCL2Mass("LinearGuider", g, frames, exposure, driftPerFrame, pe, noiseSigma, SEED, sim);
                if (s.finalRMS < sBest.finalRMS) { sBest = s; bestP = fmtParams("gain=%.1f len=%d", gain, length); }
            }
            printTunedRow("LinearGuider", sDef, sBest, bestP);
        }
        // HysteresisGuider
        {
            HysteresisGuider gd("RA"); gd.setGain(0.6); gd.setHysteresis(0.1); gd.setMinMove(0.1);
            Stats sDef = runCL2Mass("HysteresisGuider", gd, frames, exposure, driftPerFrame, pe, noiseSigma, SEED, sim);
            Stats sBest = sDef; std::string bestP = "gain=0.6 hyst=0.1";
            for (double gain : kHystGain) for (double hyst : kHystValue) {
                HysteresisGuider g("RA"); g.setGain(gain); g.setHysteresis(hyst); g.setMinMove(0.1);
                Stats s = runCL2Mass("HysteresisGuider", g, frames, exposure, driftPerFrame, pe, noiseSigma, SEED, sim);
                if (s.finalRMS < sBest.finalRMS) { sBest = s; bestP = fmtParams("gain=%.1f hyst=%.1f", gain, hyst); }
            }
            printTunedRow("HysteresisGuider", sDef, sBest, bestP);
        }
        // GPG
        {
            GaussianProcessGuider gd(makeGPGParams(gpgInitPeriod, gpgLearn));
            gd.SetLearningRate(1.0);
            Stats sDef = runCLGPG2Mass("GPG (learn)", gd, frames, exposure, driftPerFrame, pe, noiseSigma, SEED, sim);
            Stats sBest = sDef; std::string bestP = "cg=0.8 pkls=10";
            for (double cg : kGpgCtrlGain) for (double pkls : kGpgPKLengthSc) {
                auto p = makeGPGParams(gpgInitPeriod, gpgLearn);
                p.control_gain_ = cg; p.PKLengthScale_ = pkls;
                GaussianProcessGuider gpg(p); gpg.SetLearningRate(1.0);
                Stats s = runCLGPG2Mass("GPG (learn)", gpg, frames, exposure, driftPerFrame, pe, noiseSigma, SEED, sim);
                if (s.finalRMS < sBest.finalRMS) { sBest = s; bestP = fmtParams("cg=%.1f pkls=%.0f", cg, pkls); }
            }
            printTunedRow("GPG (learn)", sDef, sBest, bestP);
        }
        return;
    }

    {
        LinearGuider g("RA"); g.setGain(0.7); g.setMinMove(0.1); g.setLength(25);
        auto r = runCL2Mass("LinearGuider", g, frames, exposure, driftPerFrame, pe, noiseSigma, SEED, sim);
        printHeader(title);
        printRow(r);
    }
    {
        HysteresisGuider g("RA"); g.setGain(0.6); g.setHysteresis(0.1); g.setMinMove(0.1);
        auto r = runCL2Mass("HysteresisGuider", g, frames, exposure, driftPerFrame, pe, noiseSigma, SEED, sim);
        printRow(r);
    }
    {
        GaussianProcessGuider gpg(makeGPGParams(gpgInitPeriod, gpgLearn));
        gpg.SetLearningRate(1.0);
        auto r = runCLGPG2Mass("GPG (learn)", gpg, frames, exposure, driftPerFrame, pe, noiseSigma, SEED, sim);
        printRow(r);
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

    if (g_tune)
    {
        printTunedHeader(title);
        // LinearGuider sweep
        {
            LinearGuider gd("RA"); gd.setGain(0.7); gd.setMinMove(0.1); gd.setLength(25);
            Stats sDef = runCL("LinearGuider", gd, frames, exposure, driftPerFrame, pe, noiseSigma, SEED);
            Stats sBest = sDef; std::string bestP = "gain=0.7 len=25";
            for (double gain : kLinearGain) for (int length : kLinearLength) {
                LinearGuider g("RA"); g.setGain(gain); g.setMinMove(0.1); g.setLength(length);
                Stats s = runCL("LinearGuider", g, frames, exposure, driftPerFrame, pe, noiseSigma, SEED);
                if (s.finalRMS < sBest.finalRMS) { sBest = s; bestP = fmtParams("gain=%.1f len=%d", gain, length); }
            }
            printTunedRow("LinearGuider", sDef, sBest, bestP);
        }
        // HysteresisGuider sweep
        {
            HysteresisGuider gd("RA"); gd.setGain(0.6); gd.setHysteresis(0.1); gd.setMinMove(0.1);
            Stats sDef = runCL("HysteresisGuider", gd, frames, exposure, driftPerFrame, pe, noiseSigma, SEED);
            Stats sBest = sDef; std::string bestP = "gain=0.6 hyst=0.1";
            for (double gain : kHystGain) for (double hyst : kHystValue) {
                HysteresisGuider g("RA"); g.setGain(gain); g.setHysteresis(hyst); g.setMinMove(0.1);
                Stats s = runCL("HysteresisGuider", g, frames, exposure, driftPerFrame, pe, noiseSigma, SEED);
                if (s.finalRMS < sBest.finalRMS) { sBest = s; bestP = fmtParams("gain=%.1f hyst=%.1f", gain, hyst); }
            }
            printTunedRow("HysteresisGuider", sDef, sBest, bestP);
        }
        // GPG sweep (over control_gain x PKLengthScale; PKSignalVariance, period kept at defaults)
        {
            GaussianProcessGuider gd(makeGPGParams(gpgPeriod, gpgLearn));
            gd.SetLearningRate(1.0);
            Stats sDef = runCLGPG("GPG", gd, frames, exposure, driftPerFrame, pe, noiseSigma, SEED);
            Stats sBest = sDef; std::string bestP = "cg=0.8 pkls=10";
            for (double cg : kGpgCtrlGain) for (double pkls : kGpgPKLengthSc) {
                auto p = makeGPGParams(gpgPeriod, gpgLearn);
                p.control_gain_ = cg; p.PKLengthScale_ = pkls;
                GaussianProcessGuider gpg(p); gpg.SetLearningRate(1.0);
                Stats s = runCLGPG("GPG", gpg, frames, exposure, driftPerFrame, pe, noiseSigma, SEED);
                if (s.finalRMS < sBest.finalRMS) { sBest = s; bestP = fmtParams("cg=%.1f pkls=%.0f", cg, pkls); }
            }
            printTunedRow("GPG", sDef, sBest, bestP);
        }
        return;
    }

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

// WG5: step disturbance scenario.
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

    // For WG5-style scenarios we sweep on the "Late recovery" segment RMS,
    // since that captures settled performance after the step disturbance.
    auto lateRMSof = [&](const Stats &r) {
        return segRMS(r.allPositions, frames - nWindow, frames);
    };

    if (g_tune)
    {
        printf("\n%s\n", title);
        printf("  step: %+.1f\" at frame %d (t=%.0fs)  swept on Late(recov)\n",
               stepDelta, stepFrame, stepFrame * exposure);
        printf("  %-22s  %10s  %10s  %s\n",
               "Algorithm", "Default", "Best", "Best params");
        printf("  %-22s  %10s  %10s  %s\n",
               "----------------------", "----------", "----------", "------------");

        auto printH5Tuned = [&](const char *algo, double defLate, double bestLate, const std::string &p) {
            printf("  %-22s  %9.3f\"  %9.3f\"  %s\n", algo, defLate, bestLate, p.c_str());
        };

        // LinearGuider
        {
            LinearGuider gd("RA"); gd.setGain(0.7); gd.setMinMove(0.1); gd.setLength(25);
            Stats sDef = runCL("LinearGuider", gd, frames, exposure, 0.0, pe, noiseSigma, SEED, steps, true);
            double defLate = lateRMSof(sDef), bestLate = defLate;
            std::string bestP = "gain=0.7 len=25";
            for (double gain : kLinearGain) for (int length : kLinearLength) {
                LinearGuider g("RA"); g.setGain(gain); g.setMinMove(0.1); g.setLength(length);
                Stats s = runCL("LinearGuider", g, frames, exposure, 0.0, pe, noiseSigma, SEED, steps, true);
                double late = lateRMSof(s);
                if (late < bestLate) { bestLate = late; bestP = fmtParams("gain=%.1f len=%d", gain, length); }
            }
            printH5Tuned("LinearGuider", defLate, bestLate, bestP);
        }
        // HysteresisGuider
        {
            HysteresisGuider gd("RA"); gd.setGain(0.6); gd.setHysteresis(0.1); gd.setMinMove(0.1);
            Stats sDef = runCL("HysteresisGuider", gd, frames, exposure, 0.0, pe, noiseSigma, SEED, steps, true);
            double defLate = lateRMSof(sDef), bestLate = defLate;
            std::string bestP = "gain=0.6 hyst=0.1";
            for (double gain : kHystGain) for (double hyst : kHystValue) {
                HysteresisGuider g("RA"); g.setGain(gain); g.setHysteresis(hyst); g.setMinMove(0.1);
                Stats s = runCL("HysteresisGuider", g, frames, exposure, 0.0, pe, noiseSigma, SEED, steps, true);
                double late = lateRMSof(s);
                if (late < bestLate) { bestLate = late; bestP = fmtParams("gain=%.1f hyst=%.1f", gain, hyst); }
            }
            printH5Tuned("HysteresisGuider", defLate, bestLate, bestP);
        }
        // GPG
        {
            GaussianProcessGuider gd(makeGPGParams(gpgInitPeriod, true));
            gd.SetLearningRate(1.0);
            Stats sDef = runCLGPG("GPG (learn)", gd, frames, exposure, 0.0, pe, noiseSigma, SEED, steps, true);
            double defLate = lateRMSof(sDef), bestLate = defLate;
            std::string bestP = "cg=0.8 pkls=10";
            for (double cg : kGpgCtrlGain) for (double pkls : kGpgPKLengthSc) {
                auto p = makeGPGParams(gpgInitPeriod, true);
                p.control_gain_ = cg; p.PKLengthScale_ = pkls;
                GaussianProcessGuider gpg(p); gpg.SetLearningRate(1.0);
                Stats s = runCLGPG("GPG (learn)", gpg, frames, exposure, 0.0, pe, noiseSigma, SEED, steps, true);
                double late = lateRMSof(s);
                if (late < bestLate) { bestLate = late; bestP = fmtParams("cg=%.1f pkls=%.0f", cg, pkls); }
            }
            printH5Tuned("GPG (learn)", defLate, bestLate, bestP);
        }
        return;
    }

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

    if (g_tune)
    {
        printTunedHeader(title);
        // LinearGuider
        {
            LinearGuider gd("RA"); gd.setGain(0.7); gd.setMinMove(0.1); gd.setLength(25);
            Stats sDef = runCL("LinearGuider", gd, frames, exposure, driftPerFrame, pe, noiseSigma, SEED, steps);
            Stats sBest = sDef; std::string bestP = "gain=0.7 len=25";
            for (double gain : kLinearGain) for (int length : kLinearLength) {
                LinearGuider g("RA"); g.setGain(gain); g.setMinMove(0.1); g.setLength(length);
                Stats s = runCL("LinearGuider", g, frames, exposure, driftPerFrame, pe, noiseSigma, SEED, steps);
                if (s.finalRMS < sBest.finalRMS) { sBest = s; bestP = fmtParams("gain=%.1f len=%d", gain, length); }
            }
            printTunedRow("LinearGuider", sDef, sBest, bestP);
        }
        // HysteresisGuider
        {
            HysteresisGuider gd("RA"); gd.setGain(0.6); gd.setHysteresis(0.1); gd.setMinMove(0.1);
            Stats sDef = runCL("HysteresisGuider", gd, frames, exposure, driftPerFrame, pe, noiseSigma, SEED, steps);
            Stats sBest = sDef; std::string bestP = "gain=0.6 hyst=0.1";
            for (double gain : kHystGain) for (double hyst : kHystValue) {
                HysteresisGuider g("RA"); g.setGain(gain); g.setHysteresis(hyst); g.setMinMove(0.1);
                Stats s = runCL("HysteresisGuider", g, frames, exposure, driftPerFrame, pe, noiseSigma, SEED, steps);
                if (s.finalRMS < sBest.finalRMS) { sBest = s; bestP = fmtParams("gain=%.1f hyst=%.1f", gain, hyst); }
            }
            printTunedRow("HysteresisGuider", sDef, sBest, bestP);
        }
        // GPG (learn=true)
        {
            GaussianProcessGuider gd(makeGPGParams(gpgInitPeriod, true));
            gd.SetLearningRate(1.0);
            Stats sDef = runCLGPG("GPG (learn)", gd, frames, exposure, driftPerFrame, pe, noiseSigma, SEED, steps);
            Stats sBest = sDef; std::string bestP = "cg=0.8 pkls=10";
            for (double cg : kGpgCtrlGain) for (double pkls : kGpgPKLengthSc) {
                auto p = makeGPGParams(gpgInitPeriod, true);
                p.control_gain_ = cg; p.PKLengthScale_ = pkls;
                GaussianProcessGuider gpg(p); gpg.SetLearningRate(1.0);
                Stats s = runCLGPG("GPG (learn)", gpg, frames, exposure, driftPerFrame, pe, noiseSigma, SEED, steps);
                if (s.finalRMS < sBest.finalRMS) { sBest = s; bestP = fmtParams("cg=%.1f pkls=%.0f", cg, pkls); }
            }
            printTunedRow("GPG (learn)", sDef, sBest, bestP);
        }
        return;
    }

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

static void runDirectDriveScenario(const char *title,
                                   int frames, double exposure,
                                   double driftPerFrame,
                                   double noiseSigma,
                                   const SimDirectDriveParams &dd,
                                   const Steps &steps = {})
{
    const uint32_t SEED = 42;

    // GPG given a nominal period equal to the ripple period when present,
    // so it can learn the repeatable cogging component.
    const double gpgInitPeriod = (dd.ripple_T > 0.0) ? dd.ripple_T : 100.0;

    if (g_tune)
    {
        printTunedHeader(title);
        {
            LinearGuider gd("RA"); gd.setGain(0.7); gd.setMinMove(0.1); gd.setLength(25);
            Stats sDef = runCLDirectDrive("LinearGuider", gd, frames, exposure, driftPerFrame, noiseSigma, SEED, dd, steps);
            Stats sBest = sDef; std::string bestP = "gain=0.7 len=25";
            for (double gain : kLinearGain) for (int length : kLinearLength) {
                LinearGuider g("RA"); g.setGain(gain); g.setMinMove(0.1); g.setLength(length);
                Stats s = runCLDirectDrive("LinearGuider", g, frames, exposure, driftPerFrame, noiseSigma, SEED, dd, steps);
                if (s.finalRMS < sBest.finalRMS) { sBest = s; bestP = fmtParams("gain=%.1f len=%d", gain, length); }
            }
            printTunedRow("LinearGuider", sDef, sBest, bestP);
        }
        {
            HysteresisGuider gd("RA"); gd.setGain(0.6); gd.setHysteresis(0.1); gd.setMinMove(0.1);
            Stats sDef = runCLDirectDrive("HysteresisGuider", gd, frames, exposure, driftPerFrame, noiseSigma, SEED, dd, steps);
            Stats sBest = sDef; std::string bestP = "gain=0.6 hyst=0.1";
            for (double gain : kHystGain) for (double hyst : kHystValue) {
                HysteresisGuider g("RA"); g.setGain(gain); g.setHysteresis(hyst); g.setMinMove(0.1);
                Stats s = runCLDirectDrive("HysteresisGuider", g, frames, exposure, driftPerFrame, noiseSigma, SEED, dd, steps);
                if (s.finalRMS < sBest.finalRMS) { sBest = s; bestP = fmtParams("gain=%.1f hyst=%.1f", gain, hyst); }
            }
            printTunedRow("HysteresisGuider", sDef, sBest, bestP);
        }
        {
            GaussianProcessGuider gd(makeGPGParams(gpgInitPeriod, true));
            gd.SetLearningRate(1.0);
            Stats sDef = runCLGPGDirectDrive("GPG (learn)", gd, frames, exposure, driftPerFrame, noiseSigma, SEED, dd, steps);
            Stats sBest = sDef; std::string bestP = "cg=0.8 pkls=10";
            for (double cg : kGpgCtrlGain) for (double pkls : kGpgPKLengthSc) {
                auto p = makeGPGParams(gpgInitPeriod, true);
                p.control_gain_ = cg; p.PKLengthScale_ = pkls;
                GaussianProcessGuider gpg(p); gpg.SetLearningRate(1.0);
                Stats s = runCLGPGDirectDrive("GPG (learn)", gpg, frames, exposure, driftPerFrame, noiseSigma, SEED, dd, steps);
                if (s.finalRMS < sBest.finalRMS) { sBest = s; bestP = fmtParams("cg=%.1f pkls=%.0f", cg, pkls); }
            }
            printTunedRow("GPG (learn)", sDef, sBest, bestP);
        }
        return;
    }

    {
        LinearGuider g("RA"); g.setGain(0.7); g.setMinMove(0.1); g.setLength(25);
        auto r = runCLDirectDrive("LinearGuider", g, frames, exposure, driftPerFrame, noiseSigma, SEED, dd, steps);
        printHeader(title);
        printRow(r);
    }
    {
        HysteresisGuider g("RA"); g.setGain(0.6); g.setHysteresis(0.1); g.setMinMove(0.1);
        auto r = runCLDirectDrive("HysteresisGuider", g, frames, exposure, driftPerFrame, noiseSigma, SEED, dd, steps);
        printRow(r);
    }
    {
        GaussianProcessGuider gpg(makeGPGParams(gpgInitPeriod, true));
        gpg.SetLearningRate(1.0);
        auto r = runCLGPGDirectDrive("GPG (learn)", gpg, frames, exposure, driftPerFrame, noiseSigma, SEED, dd, steps);
        printRow(r);
    }
}

static void runStrainWaveScenario(const char *title,
                                  int frames, double exposure,
                                  double driftPerFrame,
                                  double noiseSigma,
                                  const SimStrainWaveParams &sw,
                                  const Steps &steps = {})
{
    const uint32_t SEED = 42;
    // GPG given the KE period at sidereal rate as a learnable period.
    const double gpgInitPeriod = (sw.T_KE > 0.0) ? sw.T_KE : 100.0;

    if (g_tune)
    {
        printTunedHeader(title);
        {
            LinearGuider gd("RA"); gd.setGain(0.7); gd.setMinMove(0.1); gd.setLength(25);
            Stats sDef = runCLStrainWave("LinearGuider", gd, frames, exposure, driftPerFrame, noiseSigma, SEED, sw, steps);
            Stats sBest = sDef; std::string bestP = "gain=0.7 len=25";
            for (double gain : kLinearGain) for (int length : kLinearLength) {
                LinearGuider g("RA"); g.setGain(gain); g.setMinMove(0.1); g.setLength(length);
                Stats s = runCLStrainWave("LinearGuider", g, frames, exposure, driftPerFrame, noiseSigma, SEED, sw, steps);
                if (s.finalRMS < sBest.finalRMS) { sBest = s; bestP = fmtParams("gain=%.1f len=%d", gain, length); }
            }
            printTunedRow("LinearGuider", sDef, sBest, bestP);
        }
        {
            HysteresisGuider gd("RA"); gd.setGain(0.6); gd.setHysteresis(0.1); gd.setMinMove(0.1);
            Stats sDef = runCLStrainWave("HysteresisGuider", gd, frames, exposure, driftPerFrame, noiseSigma, SEED, sw, steps);
            Stats sBest = sDef; std::string bestP = "gain=0.6 hyst=0.1";
            for (double gain : kHystGain) for (double hyst : kHystValue) {
                HysteresisGuider g("RA"); g.setGain(gain); g.setHysteresis(hyst); g.setMinMove(0.1);
                Stats s = runCLStrainWave("HysteresisGuider", g, frames, exposure, driftPerFrame, noiseSigma, SEED, sw, steps);
                if (s.finalRMS < sBest.finalRMS) { sBest = s; bestP = fmtParams("gain=%.1f hyst=%.1f", gain, hyst); }
            }
            printTunedRow("HysteresisGuider", sDef, sBest, bestP);
        }
        {
            GaussianProcessGuider gd(makeGPGParams(gpgInitPeriod, true));
            gd.SetLearningRate(1.0);
            Stats sDef = runCLGPGStrainWave("GPG (learn)", gd, frames, exposure, driftPerFrame, noiseSigma, SEED, sw, steps);
            Stats sBest = sDef; std::string bestP = "cg=0.8 pkls=10";
            for (double cg : kGpgCtrlGain) for (double pkls : kGpgPKLengthSc) {
                auto p = makeGPGParams(gpgInitPeriod, true);
                p.control_gain_ = cg; p.PKLengthScale_ = pkls;
                GaussianProcessGuider gpg(p); gpg.SetLearningRate(1.0);
                Stats s = runCLGPGStrainWave("GPG (learn)", gpg, frames, exposure, driftPerFrame, noiseSigma, SEED, sw, steps);
                if (s.finalRMS < sBest.finalRMS) { sBest = s; bestP = fmtParams("cg=%.1f pkls=%.0f", cg, pkls); }
            }
            printTunedRow("GPG (learn)", sDef, sBest, bestP);
        }
        return;
    }

    {
        LinearGuider g("RA"); g.setGain(0.7); g.setMinMove(0.1); g.setLength(25);
        auto r = runCLStrainWave("LinearGuider", g, frames, exposure, driftPerFrame, noiseSigma, SEED, sw, steps);
        printHeader(title);
        printRow(r);
    }
    {
        HysteresisGuider g("RA"); g.setGain(0.6); g.setHysteresis(0.1); g.setMinMove(0.1);
        auto r = runCLStrainWave("HysteresisGuider", g, frames, exposure, driftPerFrame, noiseSigma, SEED, sw, steps);
        printRow(r);
    }
    {
        GaussianProcessGuider gpg(makeGPGParams(gpgInitPeriod, true));
        gpg.SetLearningRate(1.0);
        auto r = runCLGPGStrainWave("GPG (learn)", gpg, frames, exposure, driftPerFrame, noiseSigma, SEED, sw, steps);
        printRow(r);
    }
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

int main(int argc, char *argv[])
{
    QCoreApplication app(argc, argv);

    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--tune") == 0) g_tune = true;
    }

    printf("=== Guide Algorithm RMSE Benchmark ===\n");
    if (g_tune) {
        printf("Mode: TUNING (per-scenario parameter sweep, ~10 min)\n");
        printf("Each row shows the default-param result and the best from a small\n");
        printf("sweep. Sweeps: LinearGuider gain*length (4*3=12), HysteresisGuider\n");
        printf("gain*hysteresis (4*3=12), GPG control_gain*PKLengthScale (3*3=9).\n");
        printf("Pass without --tune for the fast benchmark.\n\n");
    } else {
        printf("Open RMS:  uncorrected mount error over all frames\n");
        printf("All RMS:   corrected residual, all frames (includes learning transient)\n");
        printf("Final RMS: corrected residual, last quarter (algorithm should have converged)\n");
        printf("Reduction: (1 - Final/Open) x 100%%\n");
        printf("DeTrend:   final-quarter RMS after removing linear drift (WG6 only)\n");
        printf("(Pass --tune to run per-scenario parameter sweeps for each algorithm.)\n");
    }

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

    printf("\n\n--- Worm-Gear Mount Scenarios ---\n");
    printf("GPG: compute_period=ON, initial period differs from true period.\n");
    printf("Use Final RMS (last quarter) as the settled performance metric.\n");

    // WG1: Classic multi-harmonic, clean
    {
        PEParams pe;
        pe.T = 30.0; pe.A1 = 4.0; pe.A2 = 1.6; pe.A3 = 0.8;
        runHarmonicScenario(
            "WG1: Harmonic PE  T=30s  A1=4.0\" A2=1.6\" A3=0.8\"  frames=400 (26 periods)  noise=0.10\"\n"
            "    [clean multi-harmonic; GPG periodic kernel only covers fundamental]",
            400, 2.0, 0.0, pe, 0.10, 100.0);
    }

    // WG2: Non-stationary amplitude (amplitude modulated by 10-min envelope)
    {
        PEParams pe;
        pe.T = 30.0; pe.A1 = 4.0; pe.A2 = 1.5; pe.A3 = 0.8;
        pe.modT = 600.0; pe.modDepth = 0.5;
        runHarmonicScenario(
            "WG2: Non-stationary PE  T=30s  amplitude x(1+0.5*sin(t/600s))  frames=600  noise=0.15\"\n"
            "    [amplitude modulation breaks GP stationarity; linear regression adapts locally]",
            600, 2.0, 0.0, pe, 0.15, 100.0);
    }

    // WG3: Stochastic non-repeatable component (OU process overlaid on harmonics)
    {
        PEParams pe;
        pe.T = 30.0; pe.A1 = 4.0; pe.A2 = 1.2; pe.A3 = 0.6;
        pe.ouTau = 15.0; pe.ouSigSS = 1.0;
        runHarmonicScenario(
            "WG3: Stochastic PE  T=30s  + OU(tau=15s, sigma=1.0\")  frames=400  noise=0.20\"\n"
            "    [non-repeatable component; all algorithms converge toward same ceiling]",
            400, 2.0, 0.0, pe, 0.20, 100.0);
    }

    // WG4: Directional flex-spline hysteresis
    {
        PEParams pe;
        pe.T = 30.0; pe.A1 = 4.0; pe.A2 = 1.6; pe.A3 = 0.8;
        pe.hystLag = 0.5;
        runHarmonicScenario(
            "WG4: Hysteresis lag  T=30s  A1=4.0\" A2=1.6\" A3=0.8\"  hystLag=0.5\"  frames=300  noise=0.15\"\n"
            "    [sign discontinuity at each direction reversal; smooth GP kernels over-smooth]",
            300, 2.0, 0.0, pe, 0.15, 100.0);
    }

    // WG5: Step disturbance -- single large step at start of last quarter.
    // 500 frames (2000s): GPG learns T=480s in first 375 frames, then step hits.
    // Pre-step = baseline PE correction quality.
    // Immediate (fr375-399) = first 100s of recovery.
    // Late (fr475-499) = last 100s -- should all converge, but rate differs.
    {
        PEParams pe;
        pe.T = 480.0; pe.A1 = 5.0;
        runH5Scenario(
            "WG5: Step disturbance  T=480s  A=5.0\"  step=+5.0\" at fr375 (t=1500s)  frames=500  noise=0.20\"\n"
            "    [GPG has 375 frames to learn period before step; Immediate/Late show recovery speed]",
            500, 4.0, pe, 0.20, 200.0, 375, +5.0);
    }

    // WG6: Combined drift + harmonics -- with detrended final RMS
    {
        PEParams pe;
        pe.T = 30.0; pe.A1 = 3.5; pe.A2 = 1.4; pe.A3 = 0.7;
        runHarmonicScenario(
            "WG6: Drift + harmonics  drift=0.05\"/frame  T=30s  A1=3.5\" A2=1.4\" A3=0.7\"  frames=600  noise=0.25\"\n"
            "    [DeTrend removes linear drift component so harmonic correction is visible]",
            600, 2.0, 0.05, pe, 0.25, 100.0, true);
    }

    // WG7: Torsional Compliance & Resonance
    {
        PEParams pe;
        pe.T = 30.0; pe.A1 = 4.0; pe.A2 = 1.6; pe.A3 = 0.8;
        Sim2MassParams sim;
        sim.Ks = 250.0; // Moderate stiffness
        sim.Bs = 2.0;   // Moderate damping
        sim.stiction = 0.0; // Clean compliance
        sim.coulomb = 0.0;
        runComplianceScenario(
            "WG7: Torsional Compliance & Resonance  T=30s  A1=4.0\" A2=1.6\" A3=0.8\"  Ks=250.0  Bs=2.0\n"
            "    [Simulates torsional wind-up & resonant mount vibrations; tests predictive cancellation]",
            400, 2.0, 0.0, pe, 0.10, 100.0, sim);
    }

    // WG8: Static Friction + Compliance Wind-Up
    {
        PEParams pe;
        pe.T = 30.0; pe.A1 = 4.0; pe.A2 = 1.6; pe.A3 = 0.8;
        Sim2MassParams sim;
        sim.Ks = 200.0; // More flexible shaft
        sim.Bs = 2.0;
        sim.stiction = 1.5; // High stiction
        sim.coulomb = 0.5;   // High Coulomb friction
        runComplianceScenario(
            "WG8: Static Friction + Compliance Wind-Up  T=30s  A1=4.0\" A2=1.6\" A3=0.8\"  Ks=200  stiction=1.5Nm\n"
            "    [Simulates stiction stick-slip / deadband play; tests deadband-aware recovery]",
            400, 2.0, 0.0, pe, 0.10, 100.0, sim);
    }

    // WG9: Fully Adaptive PE Learning
    {
        PEParams pe;
        pe.T = 40.0; pe.A1 = 4.0; pe.A2 = 1.6; pe.A3 = 0.8; // Unknown frequency T=40s
        runHarmonicScenario(
            "WG9: Fully Adaptive PE Learning  T=40s  A1=4.0\" A2=1.6\" A3=0.8\"  (Unknown period)\n"
            "    [Verifies active frequency estimator converges to unknown period T=40s and schedules IMP cancellation]",
            500, 2.0, 0.0, pe, 0.10, 100.0);
    }

    // WG10: Sawtooth PE -- worm with slow drift, fast snap-back
    // Real worms with single-tooth-pitch error look like sawtooth: linear
    // build-up then sharp return at the tooth-engagement transition.
    // FFT picks up the fundamental at T plus a Fourier comb 2T, 3T, ...;
    // Demonstrates how algorithms handle harmonic content beyond their model.
    {
        PEParams pe;
        pe.T = 30.0; pe.A1 = 4.0; pe.waveform = PE_SAWTOOTH;
        runHarmonicScenario(
            "WG10: Sawtooth PE  T=30s  A=4.0\"  frames=400  noise=0.10\"\n"
            "    [non-sinusoidal worm error; harmonic-rich spectrum tests how many modes IMP captures]",
            400, 2.0, 0.0, pe, 0.10, 100.0);
    }

    // WG11: Triangle PE -- constant-velocity reversals
    // Common in mounts with friction-limited motor velocity; the worm tracks
    // linearly then reverses cleanly. Discontinuity in velocity (not position).
    {
        PEParams pe;
        pe.T = 30.0; pe.A1 = 4.0; pe.waveform = PE_TRIANGLE;
        runHarmonicScenario(
            "WG11: Triangle PE  T=30s  A=4.0\"  frames=400  noise=0.10\"\n"
            "    [constant-velocity ramps with mid-period sign reversal; tests transient response at direction change]",
            400, 2.0, 0.0, pe, 0.10, 100.0);
    }

    // WG12: Half-rectified PE -- single-sided drive train error
    // Half-period zero, half-period sine. Strong DC + harmonic content. The
    // Tests how algorithms respond to a DC bias overlaid on sinusoidal PE.
    // component without help from the motor state.
    {
        PEParams pe;
        pe.T = 30.0; pe.A1 = 4.0; pe.waveform = PE_HALFRECT;
        runHarmonicScenario(
            "WG12: Half-rectified PE  T=30s  A=4.0\"  frames=400  noise=0.10\"\n"
            "    [single-sided worm-tooth error; large DC bias plus harmonic spectrum]",
            400, 2.0, 0.0, pe, 0.10, 100.0);
    }

    // WG13: True backlash (direction-dependent deadband)
    // Worm-and-wheel gear play: when motor reverses, axis stays put until
    // the +/- backlash gap is traversed. Different from WG4's hysteresis lag
    // Backlash creates a signed asymmetry in the gear-mesh deflection
    // path that triggers on commanded direction changes -- this scenario tests
    // it. Compliance is mild (Ks=400) so the backlash effect is isolated.
    {
        PEParams pe;
        pe.T = 30.0; pe.A1 = 4.0; pe.A2 = 1.6; pe.A3 = 0.8;
        Sim2MassParams sim;
        sim.Ks = 400.0;     // Stiffer than WG7/WG8 to isolate backlash effect
        sim.Bs = 4.0;       // More damped
        sim.stiction = 0.0; // No stiction
        sim.coulomb = 0.0;
        sim.backlash = 1.0; // 1 arcsec of gear play (typical for hobby mounts)
        runComplianceScenario(
            "WG13: True backlash  T=30s  A1=4.0\" A2=1.6\" A3=0.8\"  Ks=400.0  backlash=1.0\"\n"
            "    [direction-dependent deadband at each PE zero-crossing; exercises punch-through compensation]",
            400, 2.0, 0.0, pe, 0.10, 100.0, sim);
    }

    // WG14: Period drift -- worm period changes slowly over the run
    // Real harmonic-drive/worm-gear mounts drift in period with temperature.
    // ~1-5% per hour is typical. Test uses 10% drift over the run to make
    // adaptation visible across the run; GPG re-fits its kernel,
    // gradient-based period inference. Both should adapt; the question is
    // how cleanly.
    // T0=30s, Tdot=0.0025 s/s -> T grows from 30s to ~33s over 1200s = 600 frames.
    {
        PEParams pe;
        pe.T = 30.0; pe.Tdot = 0.0025; pe.A1 = 4.0;
        runHarmonicScenario(
            "WG14: Period drift  T0=30s -> ~33s over 600 frames  A=4.0\"  noise=0.10\"\n"
            "    [worm period drifts 10% over the run; tests period re-tracking across algorithms]",
            600, 2.0, 0.0, pe, 0.10, 30.0);
    }

    // WG15: Mid-run plant change -- step disturbance at frame 200 on top of
    // ongoing PE. Simulates focus shift, refraction change, or worm contact
    // zone change. Tests adaptive recovery: pre-step settled performance vs
    // post-step transient.
    {
        PEParams pe;
        pe.T = 30.0; pe.A1 = 5.0;
        runH5Scenario(
            "WG15: Mid-run amplitude step  T=30s  A1=5.0\"  +3\" step at fr200\n"
            "    [adaptive recovery test: sudden plant change requiring re-learning]",
            400, 2.0, pe, 0.10, 30.0, 200, +3.0);
    }

    // WG16: Slow worm PE + backlash. Long worm period (T=480s) means current_u
    // reverses only ~once per ~60 frames. The smarter punch-through trigger
    // requires kPunchSustained=30 frames of sustained direction before firing,
    // so it stays dormant on WG13's fast PE (T=30s, reversal every ~7 frames)
    // but engages here where actual backlash-gap traversal happens.
    {
        PEParams pe;
        pe.T = 480.0; pe.A1 = 5.0;
        Sim2MassParams sim;
        sim.Ks = 400.0;
        sim.Bs = 4.0;
        sim.stiction = 0.0;
        sim.coulomb = 0.0;
        sim.backlash = 1.0;
        runComplianceScenario(
            "WG16: Slow worm + backlash  T=480s  A=5.0\"  Ks=400.0  backlash=1.0\"\n"
            "    [long-period PE so sustained direction holds; punch-through engages on each true reversal;\n"
            "     simulates 2-mass flexible plant with declared compliance]",
            360, 4.0, 0.0, pe, 0.10, 480.0, sim,
            /*gpgLearn=*/true,
            /*declareCompliance=*/true);
    }

    // WG17: Slow-response motor + larger backlash. Models a direct-drive or
    // low-reduction mount where the motor itself has significant inertia
    // (Jm=0.5 vs the default 0.02 -> ~25x heavier). Motor time constant
    // tau = Jm/Bf = 0.5/0.05 = 10s -> several frames at dt=4s to settle.
    // Backlash=2" is double WG16 to make gap traversal a multi-frame event.
    // This is the scenario where backlash punch-through should actually help:
    // the motor takes time to cross the gap, so pre-positioning it via
    // punch-through is a real benefit.
    {
        PEParams pe;
        pe.T = 480.0; pe.A1 = 5.0;
        Sim2MassParams sim;
        sim.Jm = 0.5;       // Slow motor (vs default 0.02)
        sim.Bf = 0.05;      // Less damping
        sim.Ks = 300.0;
        sim.Bs = 3.0;
        sim.stiction = 0.0;
        sim.coulomb = 0.0;
        sim.backlash = 2.0; // Larger gap
        runComplianceScenario(
            "WG17: Slow-response motor + backlash  T=480s  A=5.0\"  Jm=0.5  backlash=2.0\"\n"
            "    [direct-drive style mount: motor takes several frames to settle;\n"
            "     gap traversal is a multi-frame event where punch-through should help;\n"
            "     simulates 2-mass flexible plant with declared compliance]",
            360, 4.0, 0.0, pe, 0.10, 480.0, sim,
            /*gpgLearn=*/true,
            /*declareCompliance=*/true);
    }

    // ----- Direct-Drive Scenarios -----
    //
    // A real direct-drive mount has no gearbox, no backlash, no torsional
    // compliance -- the motor is the axis. The 2-mass scenarios above
    // (including WG17) cannot model this faithfully. These scenarios use a
    // dedicated single-state plant with a closed inner velocity servo;
    // cogging is injected as a rate disturbance (physically a torque, the
    // servo integrates it to position), wind as an OU rate disturbance.
    printf("\n\n--- Direct Drive Scenarios ---\n");

    // DD1: Clean fast servo + tiny ripple. Easiest scenario in the suite;
    // all algorithms should track to near the noise floor.
    {
        SimDirectDriveParams dd;
        dd.tau_servo = 0.3;
        dd.ripple_T  = 30.0;
        dd.ripple_A1 = 0.02;     // arcsec/s -- ~0.1" peak position-equivalent
        runDirectDriveScenario(
            "DD1: Clean DD  tau_servo=0.3s  ripple_T=30s  A1=0.02\"/s  noise=0.10\"\n"
            "    [fast inner servo + tiny cogging; baseline DD performance]",
            300, 4.0, 0.0, 0.10, dd);
    }

    // DD2: Cogging harmonics. Fundamental + 2nd/3rd at a short period
    // (well above the per-frame Nyquist for repeatability). Tests how each
    // algorithm handles a fast repeatable disturbance very different from
    // a 480s worm period.
    {
        SimDirectDriveParams dd;
        dd.tau_servo = 0.3;
        dd.ripple_T  = 30.0;
        dd.ripple_A1 = 0.08;
        dd.ripple_A2 = 0.03;
        dd.ripple_A3 = 0.015;
        runDirectDriveScenario(
            "DD2: DD + cogging harmonics  ripple_T=30s  A1=0.08 A2=0.03 A3=0.015\"/s  noise=0.10\"\n"
            "    [repeatable torque ripple with 2nd/3rd harmonics; GPG should learn it]",
            300, 4.0, 0.0, 0.10, dd);
    }

    // DD3: Wind OU dominant, no repeatable structure. GPG cannot learn
    // anything; pure disturbance-rejection test.
    {
        SimDirectDriveParams dd;
        dd.tau_servo = 0.3;
        dd.ripple_T  = 0.0;
        dd.wind_tau   = 20.0;
        dd.wind_sigma = 0.15;   // arcsec/s rate disturbance
        runDirectDriveScenario(
            "DD3: DD + wind OU  tau=20s  sigma=0.15\"/s  noise=0.10\"\n"
            "    [no repeatable structure; pure disturbance-rejection test, GPG cannot learn]",
            300, 4.0, 0.0, 0.10, dd);
    }

    // DD4: Sluggish inner servo (tau_servo = 1.5s vs default 0.3s). Models
    // a cheap or under-tuned DD servo. Corrections lag the command;
    // aggressive guiders can excite oscillation.
    {
        SimDirectDriveParams dd;
        dd.tau_servo = 1.5;
        dd.ripple_T  = 30.0;
        dd.ripple_A1 = 0.05;
        runDirectDriveScenario(
            "DD4: Sluggish DD servo  tau_servo=1.5s  ripple_T=30s  A1=0.05\"/s  noise=0.10\"\n"
            "    [slow inner loop: correction takes several seconds to settle;\n"
            "     aggressive gains may oscillate]",
            300, 4.0, 0.0, 0.10, dd);
    }

    // DD5: Structural mode aliased into the band. Small under-damped mount
    // resonance; sanity check that algorithms tolerate aliased noise.
    {
        SimDirectDriveParams dd;
        dd.tau_servo = 0.3;
        dd.ripple_T  = 30.0;
        dd.ripple_A1 = 0.02;
        dd.mode_hz   = 6.0;    // well above guide Nyquist -> aliased
        dd.mode_zeta = 0.05;
        dd.mode_A    = 0.05;
        runDirectDriveScenario(
            "DD5: DD + structural mode  mode=6Hz zeta=0.05 A=0.05\"  noise=0.10\"\n"
            "    [under-damped mount resonance aliased into the guide band; nothing should blow up]",
            300, 4.0, 0.0, 0.10, dd);
    }

    // ----- Strain-Wave (Harmonic Drive) Scenarios -----
    //
    // Models strain-wave mounts (ZWO HEM/AM3/AM5, iOptron HEM, etc.). Three
    // SW-specific features that the 2-mass and DD plants get wrong:
    //   1. Motor-angle-indexed KE with dominant 2x harmonic
    //   2. Soft-zone torque-dependent stiffness (Ks_eff = tanh-saturated)
    //   3. Rate-independent (Dahl) hysteresis at direction reversals
    // No deadband backlash.
    //
    // SW1-SW3 each isolate one feature; SW4 combines all three to expose
    // the soft-zone x hysteresis interaction that none of SW1-SW3 catches.
    printf("\n\n--- Strain-Wave Drive Scenarios ---\n");

    // SW1: Clean SW. Moderate KE, high enough load that the system sits in
    // the stiff regime. The 2x-per-motor-rev signature should dominate the
    // residual; soft-zone and hysteresis stay quiet.
    {
        SimStrainWaveParams sw;
        sw.T_KE = 60.0;
        sw.A1 = 3.0; sw.A2 = 0.8; sw.A3 = 0.2;
        sw.A0 = 0.3;  // small bearing-eccentricity sub-harmonic
        sw.h_max = 0.02;  // tiny: keep hysteresis dormant in SW1
        sw.Kh = 10.0;
        runStrainWaveScenario(
            "SW1: Clean SW  T_KE=60s  A1=3.0\" A2=0.8\" A3=0.2\" A0=0.3\"  noise=0.10\"\n"
            "    [stiff regime; 2x-per-motor-rev KE signature dominates residual]",
            300, 4.0, 0.0, 0.10, sw);
    }

    // SW2: Soft-zone exposed. Small KE so spring torque stays low and the
    // stiffness sits near Ks_low. Any guider expecting Ks_high will
    // overshoot. Short exposures / low noise expose this most.
    {
        SimStrainWaveParams sw;
        sw.T_KE = 60.0;
        sw.A1 = 0.8; sw.A2 = 0.2; sw.A3 = 0.05;
        sw.Ks_low = 30.0;     // softer low regime
        sw.Ks_high = 600.0;
        sw.T_knee = 1.5;      // higher knee so we stay in soft zone
        sw.h_max = 0.02;
        sw.Kh = 10.0;
        runStrainWaveScenario(
            "SW2: SW soft-zone exposed  T_KE=60s  A1=0.8\"  Ks_low=30 Ks_high=600 T_knee=1.5  noise=0.05\"\n"
            "    [low load -> system stays in low-stiffness regime; tests overshoot under reduced Ks]",
            300, 4.0, 0.0, 0.05, sw);
    }

    // SW3: Hysteresis exercised. Inject dither steps to force direction
    // reversals; each one traverses the Dahl loop. Residual should show
    // step-like offsets at reversal frames, not clean recovery.
    {
        SimStrainWaveParams sw;
        sw.T_KE = 60.0;
        sw.A1 = 1.0; sw.A2 = 0.3; sw.A3 = 0.1;
        sw.h_max = 0.3;       // saturation displacement
        sw.Kh = 12.0;         // hysteresis ~10-30% of spring torque
        // Dither pattern: 5 reversals at 50-frame intervals
        Steps dither = {{50, 2.0}, {100, -2.0}, {150, 2.0}, {200, -2.0}, {250, 2.0}};
        runStrainWaveScenario(
            "SW3: SW + hysteresis at reversals  T_KE=60s  h_max=0.3\" Kh=12  +/-2\" dithers every 50fr  noise=0.10\"\n"
            "    [each dither reversal traverses the Dahl loop; expect step-like residual offsets]",
            300, 4.0, 0.0, 0.10, sw, dither);
    }

    // SW4: Realistic combined. All three features active simultaneously
    // at moderate load -- closest analogue to a real session. Specifically
    // catches the soft-zone x hysteresis interaction: at low load, the
    // hysteresis offset is a much larger fraction of total restoring
    // torque, so reversals near zero load behave qualitatively differently
    // from reversals at high load. SW1-SW3 isolate; SW4 integrates.
    {
        SimStrainWaveParams sw;
        sw.T_KE = 60.0;
        sw.A1 = 2.0; sw.A2 = 0.5; sw.A3 = 0.15;
        sw.A0 = 0.2;
        sw.Ks_low = 40.0;
        sw.Ks_high = 800.0;
        sw.T_knee = 1.0;
        sw.h_max = 0.2;
        sw.Kh = 10.0;
        Steps dither = {{75, 1.5}, {150, -1.5}, {225, 1.5}};
        runStrainWaveScenario(
            "SW4: Realistic combined SW  KE+soft-zone+hysteresis  dithers at fr75/150/225  noise=0.10\"\n"
            "    [all three SW features active together; catches soft-zone x hysteresis interaction at reversals]",
            300, 4.0, 0.0, 0.10, sw, dither);
    }

    // ----- Belt-Drive Scenarios -----
    //
    // Belt-driven mounts (Skywatcher belt-mod HEQ5/EQ6, AZ-EQ6, Celestron
    // CGX-L) couple the motor to the worm shaft (or axis directly) through
    // a toothed belt. Belt physics = single-stage compliance:
    //   - Ks ~50-100 (softer than worm-gear; belt is the compliant element)
    //   - Bs ~5-10 (belt damps more than steel gear teeth)
    //   - backlash ~0 (pre-tensioned belts)
    // No SW-style soft-zone, no DD-style servo. Worm PE downstream of the
    // belt is unchanged, so the existing time-indexed PE applies.
    //
    // Belt physics is a parameter regime on the 2-mass plant, not a new
    // physical model -- so these scenarios reuse runComplianceScenario with
    // belt-tuned Sim2MassParams rather than duplicating ~400 lines for a
    // dedicated SimBeltParams. Tandem belt+worm (EQ6 belt-mod) is a 3-mass
    // problem deferred to a follow-up; see plans/belt-mount-simulator.md.
    printf("\n\n--- Belt Drive Scenarios ---\n");

    // B1: Clean belt-drive. Long worm PE through a soft belt coupling.
    // Tests how each algorithm handles a softer plant than worm-gear at
    // the same PE signature. Compare with WG6 (long worm PE, stiff plant).
    {
        PEParams pe;
        pe.T = 480.0; pe.A1 = 5.0;
        Sim2MassParams sim;
        sim.Ks = 80.0;       // belt-soft
        sim.Bs = 8.0;        // belt-damped
        sim.backlash = 0.0;  // pre-tensioned
        sim.stiction = 0.5;
        sim.coulomb = 0.2;
        runComplianceScenario(
            "B1: Clean belt drive  T=480s  A=5.0\"  Ks=80 Bs=8 (soft belt)  noise=0.10\"\n"
            "    [worm PE through a soft belt coupling; tests algorithms against softer plant than WG6]",
            360, 4.0, 0.0, pe, 0.10, 480.0, sim,
            /*gpgLearn=*/true,
            /*declareCompliance=*/false);
    }

    // B2: Belt + drift. Same belt regime with constant drift added. Tests
    // drift rejection on a soft plant.
    {
        PEParams pe;
        pe.T = 480.0; pe.A1 = 5.0;
        Sim2MassParams sim;
        sim.Ks = 80.0;
        sim.Bs = 8.0;
        sim.backlash = 0.0;
        sim.stiction = 0.5;
        sim.coulomb = 0.2;
        runComplianceScenario(
            "B2: Belt + drift  T=480s  A=5.0\"  drift=0.05\"/frame  Ks=80 Bs=8  noise=0.10\"\n"
            "    [soft belt plant with constant drift; tests drift rejection under belt compliance]",
            360, 4.0, 0.05, pe, 0.10, 480.0, sim,
            /*gpgLearn=*/true,
            /*declareCompliance=*/false);
    }

    // B3: Belt + slip event. A single mid-run step disturbance models a
    // brief belt slip or tensioner give. Tests step recovery on a soft
    // plant; compare with WG5 (worm-gear mid-run step).
    //
    // runComplianceScenario does not accept a Steps argument, so B3 is
    // wired inline using runCL2Mass / runCLGPG2Mass directly to inject the
    // step.
    {
        PEParams pe;
        pe.T = 480.0; pe.A1 = 5.0;
        Sim2MassParams sim;
        sim.Ks = 80.0;
        sim.Bs = 8.0;
        sim.backlash = 0.0;
        sim.stiction = 0.5;
        sim.coulomb = 0.2;
        Steps slip = {{180, 2.5}};   // single 2.5" slip at mid-run

        const char *title =
            "B3: Belt + slip event  T=480s  A=5.0\"  +2.5\" slip at fr180  Ks=80 Bs=8  noise=0.10\"\n"
            "    [soft belt plant with mid-run slip step; tests step recovery vs WG5 (stiff-plant step)]";
        const uint32_t SEED = 42;
        {
            LinearGuider g("RA"); g.setGain(0.7); g.setMinMove(0.1); g.setLength(25);
            auto r = runCL2Mass("LinearGuider", g, 360, 4.0, 0.0, pe, 0.10, SEED, sim, slip);
            printHeader(title);
            printRow(r);
        }
        {
            HysteresisGuider g("RA"); g.setGain(0.6); g.setHysteresis(0.1); g.setMinMove(0.1);
            auto r = runCL2Mass("HysteresisGuider", g, 360, 4.0, 0.0, pe, 0.10, SEED, sim, slip);
            printRow(r);
        }
        {
            GaussianProcessGuider gpg(makeGPGParams(480.0, true));
            gpg.SetLearningRate(1.0);
            auto r = runCLGPG2Mass("GPG (learn)", gpg, 360, 4.0, 0.0, pe, 0.10, SEED, sim, slip);
            printRow(r);
        }
    }

    // -----------------------------------------------------------------------
    // Realistic seeing scenarios -- colored noise + occasional centroid outliers
    // -----------------------------------------------------------------------
    //
    // The other scenarios use white Gaussian measurement noise, which an
    // algorithm with a long smoothing window (LinearGuider with large length,
    // GPG with broad SE kernel) can average away cleanly. Real atmospheric
    // seeing has a colored spectrum: a slow component on top of the white
    // floor. The slow component does not average out, so the apparent
    // "noise-rejection" advantage of long smoothing is partially fictional.
    //
    // These scenarios opt into a MeasParams that models colored seeing via an
    // OU process (tau=8s, sigma=0.4") on top of the existing white floor, and
    // adds a 4% chance of a 4x-sigma centroid outlier per frame (bright-star
    // saturation, satellite trail). Compare the algorithms' Final RMS in
    // these scenarios against their counterparts above to see how robust each
    // is to realistic measurement pathology.
    {
        PEParams pe;
        pe.T = 480.0; pe.A1 = 10.0;
        const uint32_t SEED = 42;
        MeasParams seeing;
        seeing.seeingTau   = 8.0;   // slow-component correlation time (s)
        seeing.seeingSigma = 0.4;   // slow-component sigma (arcsec)
        seeing.outlierProb = 0.04;  // 4% chance per frame
        seeing.outlierMult = 4.0;   // outlier amplitude = 4x white sigma
        const char *title =
            "S1: Realistic worm PE + seeing  T=480s  A=10.0\"  white=0.30\"  "
            "seeing OU(tau=8s,sigma=0.4\")  outliers=4%\n"
            "    [colored seeing exposes algorithms whose smoothing assumes white noise]";
        {
            LinearGuider g("RA"); g.setGain(0.7); g.setMinMove(0.1); g.setLength(25);
            auto r = runCL("LinearGuider", g, 360, 4.0, 0.0, pe, 0.30, SEED, {}, false, seeing);
            printHeader(title);
            printRow(r);
        }
        {
            HysteresisGuider g("RA"); g.setGain(0.6); g.setHysteresis(0.1); g.setMinMove(0.1);
            auto r = runCL("HysteresisGuider", g, 360, 4.0, 0.0, pe, 0.30, SEED, {}, false, seeing);
            printRow(r);
        }
        {
            GaussianProcessGuider gpg(makeGPGParams(480.0, true));
            gpg.SetLearningRate(1.0);
            auto r = runCLGPG("GPG (learn)", gpg, 360, 4.0, 0.0, pe, 0.30, SEED, {}, false, seeing);
            printRow(r);
        }
    }

    // S2: Same plant compliance as WG7 (resonant 2-mass) but with realistic
    // seeing on top, to see how compliance + colored noise interact.
    {
        PEParams pe;
        pe.T = 30.0; pe.A1 = 4.0; pe.A2 = 1.6; pe.A3 = 0.8;
        Sim2MassParams sim;
        sim.Ks = 250.0; sim.Bs = 2.0;
        const uint32_t SEED = 42;
        MeasParams seeing;
        seeing.seeingTau   = 8.0;
        seeing.seeingSigma = 0.4;
        seeing.outlierProb = 0.04;
        seeing.outlierMult = 4.0;
        const char *title =
            "S2: 2-mass resonant + seeing  T=30s  A1=4.0\" A2=1.6\" A3=0.8\"  "
            "Ks=250  white=0.10\"  seeing OU(tau=8s,sigma=0.4\")  outliers=4%\n"
            "    [compliance + colored seeing: long-window smoothers pay more than under white noise]";
        {
            LinearGuider g("RA"); g.setGain(0.7); g.setMinMove(0.1); g.setLength(25);
            auto r = runCL2Mass("LinearGuider", g, 400, 2.0, 0.0, pe, 0.10, SEED, sim, {}, false, seeing);
            printHeader(title);
            printRow(r);
        }
        {
            HysteresisGuider g("RA"); g.setGain(0.6); g.setHysteresis(0.1); g.setMinMove(0.1);
            auto r = runCL2Mass("HysteresisGuider", g, 400, 2.0, 0.0, pe, 0.10, SEED, sim, {}, false, seeing);
            printRow(r);
        }
        {
            GaussianProcessGuider gpg(makeGPGParams(30.0, true));
            gpg.SetLearningRate(1.0);
            auto r = runCLGPG2Mass("GPG (learn)", gpg, 400, 2.0, 0.0, pe, 0.10, SEED, sim, {}, false, seeing);
            printRow(r);
        }
    }

    printf("\n");
    return 0;
}
