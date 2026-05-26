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
#include "ekos/guide/internalguide/mpcguider.h"
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
static const std::vector<double> kMpcQ          = { 5.0, 10.0, 20.0, 50.0 };
static const std::vector<double> kMpcR          = { 0.05, 0.1, 0.5, 1.0 };

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
                        const Steps &steps = {}, bool collectAll = false)
{
    uint32_t sn = seed;
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

        double meas = theta_a + gaussNoise(sn, noiseSigma);
        double corr = g.guide(meas);
        theta_m -= corr;

        // Multi-rate physical integration loop for compliance/vibrations/friction
        for (int s = 0; s < stepsPerFrame; ++s) {
            double current_t = t_start + s * h;
            
            double pe = absHarmonicPE(peParams, current_t);
            double arg = 2.0 * M_PI * current_t / peParams.T;
            double pe_rate = (peParams.T > 0.0) ? (peParams.A1 * (2.0 * M_PI / peParams.T) * std::cos(arg)
                             + peParams.A2 * (4.0 * M_PI / peParams.T) * std::cos(2.0 * arg + peParams.phi2)
                             + peParams.A3 * (6.0 * M_PI / peParams.T) * std::cos(3.0 * arg + peParams.phi3)) : 0.0;
            
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
                           const Steps &steps = {}, bool collectAll = false)
{
    uint32_t sn = seed;
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

        double meas = theta_a + gaussNoise(sn, noiseSigma);
        double corr = gpg.result(meas, 100.0, exposure);
        theta_m -= corr;

        // Multi-rate physical integration loop for compliance/vibrations/friction
        for (int s = 0; s < stepsPerFrame; ++s) {
            double current_t = t_start + s * h;
            
            double pe = absHarmonicPE(peParams, current_t);
            double arg = 2.0 * M_PI * current_t / peParams.T;
            double pe_rate = (peParams.T > 0.0) ? (peParams.A1 * (2.0 * M_PI / peParams.T) * std::cos(arg)
                             + peParams.A2 * (4.0 * M_PI / peParams.T) * std::cos(2.0 * arg + peParams.phi2)
                             + peParams.A3 * (6.0 * M_PI / peParams.T) * std::cos(3.0 * arg + peParams.phi3)) : 0.0;
            
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

static void runComplianceScenario(const char *title,
                                  int frames, double exposure,
                                  double driftPerFrame, const PEParams &pe,
                                  double noiseSigma, double gpgInitPeriod,
                                  const Sim2MassParams &sim,
                                  bool gpgLearn = true)
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
        // MPCGuider
        {
            MPCGuider gd("RA"); gd.setParameters(10.0, 0.1); gd.setMinMove(0.1);
            Stats sDef = runCL2Mass("MPCGuider", gd, frames, exposure, driftPerFrame, pe, noiseSigma, SEED, sim);
            Stats sBest = sDef; std::string bestP = "Q=10 R=0.10";
            for (double Q : kMpcQ) for (double R : kMpcR) {
                MPCGuider g("RA"); g.setParameters(Q, R); g.setMinMove(0.1);
                Stats s = runCL2Mass("MPCGuider", g, frames, exposure, driftPerFrame, pe, noiseSigma, SEED, sim);
                if (s.finalRMS < sBest.finalRMS) { sBest = s; bestP = fmtParams("Q=%.0f R=%.2f", Q, R); }
            }
            printTunedRow("MPCGuider", sDef, sBest, bestP);
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
    {
        MPCGuider g("RA"); g.setParameters(10.0, 0.1); g.setMinMove(0.1);
        // Note: do NOT pass sim.backlash through setBacklash() here. The
        // current MPCSolver punch-through triggers on EVERY commanded
        // current_u zero-crossing, which fires constantly under sinusoidal
        // PE and over-compensates (H13 1.19" -> 1.40" when enabled).
        // The punch-through model needs gating (e.g., only after a long
        // direction hold, or scaled by predicted motion) before it can be
        // wired up safely. The API on MPCGuider is in place for future use.
        auto r = runCL2Mass("MPCGuider", g, frames, exposure, driftPerFrame, pe, noiseSigma, SEED, sim);
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
        // MPCGuider sweep
        {
            MPCGuider gd("RA"); gd.setParameters(10.0, 0.1); gd.setMinMove(0.1);
            Stats sDef = runCL("MPCGuider", gd, frames, exposure, driftPerFrame, pe, noiseSigma, SEED);
            Stats sBest = sDef; std::string bestP = "Q=10 R=0.10";
            for (double Q : kMpcQ) for (double R : kMpcR) {
                MPCGuider g("RA"); g.setParameters(Q, R); g.setMinMove(0.1);
                Stats s = runCL("MPCGuider", g, frames, exposure, driftPerFrame, pe, noiseSigma, SEED);
                if (s.finalRMS < sBest.finalRMS) { sBest = s; bestP = fmtParams("Q=%.0f R=%.2f", Q, R); }
            }
            printTunedRow("MPCGuider", sDef, sBest, bestP);
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
    {
        MPCGuider g("RA"); g.setParameters(10.0, 0.1); g.setMinMove(0.1);
        auto r = runCL("MPCGuider", g, frames, exposure, driftPerFrame, pe, noiseSigma, SEED);
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

    // For H5-style scenarios we sweep on the "Late recovery" segment RMS,
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
        // MPCGuider
        {
            MPCGuider gd("RA"); gd.setParameters(10.0, 0.1); gd.setMinMove(0.1);
            Stats sDef = runCL("MPCGuider", gd, frames, exposure, 0.0, pe, noiseSigma, SEED, steps, true);
            double defLate = lateRMSof(sDef), bestLate = defLate;
            std::string bestP = "Q=10 R=0.10";
            for (double Q : kMpcQ) for (double R : kMpcR) {
                MPCGuider g("RA"); g.setParameters(Q, R); g.setMinMove(0.1);
                Stats s = runCL("MPCGuider", g, frames, exposure, 0.0, pe, noiseSigma, SEED, steps, true);
                double late = lateRMSof(s);
                if (late < bestLate) { bestLate = late; bestP = fmtParams("Q=%.0f R=%.2f", Q, R); }
            }
            printH5Tuned("MPCGuider", defLate, bestLate, bestP);
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
    {
        MPCGuider g("RA"); g.setParameters(10.0, 0.1); g.setMinMove(0.1);
        auto r = runCL("MPCGuider", g, frames, exposure, 0.0, pe,
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
        // MPCGuider
        {
            MPCGuider gd("RA"); gd.setParameters(10.0, 0.1); gd.setMinMove(0.1);
            Stats sDef = runCL("MPCGuider", gd, frames, exposure, driftPerFrame, pe, noiseSigma, SEED, steps);
            Stats sBest = sDef; std::string bestP = "Q=10 R=0.10";
            for (double Q : kMpcQ) for (double R : kMpcR) {
                MPCGuider g("RA"); g.setParameters(Q, R); g.setMinMove(0.1);
                Stats s = runCL("MPCGuider", g, frames, exposure, driftPerFrame, pe, noiseSigma, SEED, steps);
                if (s.finalRMS < sBest.finalRMS) { sBest = s; bestP = fmtParams("Q=%.0f R=%.2f", Q, R); }
            }
            printTunedRow("MPCGuider", sDef, sBest, bestP);
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
    {
        MPCGuider g("RA"); g.setParameters(10.0, 0.1); g.setMinMove(0.1);
        auto r = runCL("MPCGuider", g, frames, exposure, driftPerFrame, pe,
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

    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--tune") == 0) g_tune = true;
    }

    printf("=== Guide Algorithm RMSE Benchmark ===\n");
    if (g_tune) {
        printf("Mode: TUNING (per-scenario parameter sweep, ~10 min)\n");
        printf("Each row shows the default-param result and the best from a small\n");
        printf("sweep. Sweeps: LinearGuider gain*length (4*3=12), HysteresisGuider\n");
        printf("gain*hysteresis (4*3=12), GPG control_gain*PKLengthScale (3*3=9),\n");
        printf("MPCGuider Q*R (4*4=16). Pass without --tune for the fast benchmark.\n\n");
    } else {
        printf("Open RMS:  uncorrected mount error over all frames\n");
        printf("All RMS:   corrected residual, all frames (includes learning transient)\n");
        printf("Final RMS: corrected residual, last quarter (algorithm should have converged)\n");
        printf("Reduction: (1 - Final/Open) x 100%%\n");
        printf("DeTrend:   final-quarter RMS after removing linear drift (H6 only)\n");
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

    // H7: Torsional Compliance & Resonance
    {
        PEParams pe;
        pe.T = 30.0; pe.A1 = 4.0; pe.A2 = 1.6; pe.A3 = 0.8;
        Sim2MassParams sim;
        sim.Ks = 250.0; // Moderate stiffness
        sim.Bs = 2.0;   // Moderate damping
        sim.stiction = 0.0; // Clean compliance
        sim.coulomb = 0.0;
        runComplianceScenario(
            "H7: Torsional Compliance & Resonance  T=30s  A1=4.0\" A2=1.6\" A3=0.8\"  Ks=250.0  Bs=2.0\n"
            "    [Simulates torsional wind-up & resonant mount vibrations; MPC cancels pre-emptively]",
            400, 2.0, 0.0, pe, 0.10, 100.0, sim);
    }

    // H8: Static Friction + Compliance Wind-Up
    {
        PEParams pe;
        pe.T = 30.0; pe.A1 = 4.0; pe.A2 = 1.6; pe.A3 = 0.8;
        Sim2MassParams sim;
        sim.Ks = 200.0; // More flexible shaft
        sim.Bs = 2.0;
        sim.stiction = 1.5; // High stiction
        sim.coulomb = 0.5;   // High Coulomb friction
        runComplianceScenario(
            "H8: Static Friction + Compliance Wind-Up  T=30s  A1=4.0\" A2=1.6\" A3=0.8\"  Ks=200  stiction=1.5Nm\n"
            "    [Simulates stiction stick-slip / deadband play; MPC deadband punch-through suppresses wind-up]",
            400, 2.0, 0.0, pe, 0.10, 100.0, sim);
    }

    // H9: Fully Adaptive PE Learning
    {
        PEParams pe;
        pe.T = 40.0; pe.A1 = 4.0; pe.A2 = 1.6; pe.A3 = 0.8; // Unknown frequency T=40s
        runHarmonicScenario(
            "H9: Fully Adaptive PE Learning  T=40s  A1=4.0\" A2=1.6\" A3=0.8\"  (Unknown period)\n"
            "    [Verifies active frequency estimator converges to unknown period T=40s and schedules IMP cancellation]",
            500, 2.0, 0.0, pe, 0.10, 100.0);
    }

    // H10: Sawtooth PE -- worm with slow drift, fast snap-back
    // Real worms with single-tooth-pitch error look like sawtooth: linear
    // build-up then sharp return at the tooth-engagement transition.
    // FFT picks up the fundamental at T plus a Fourier comb 2T, 3T, ...;
    // MPC's 2-oscillator IMP only models two of those harmonics.
    {
        PEParams pe;
        pe.T = 30.0; pe.A1 = 4.0; pe.waveform = PE_SAWTOOTH;
        runHarmonicScenario(
            "H10: Sawtooth PE  T=30s  A=4.0\"  frames=400  noise=0.10\"\n"
            "    [non-sinusoidal worm error; harmonic-rich spectrum tests how many modes IMP captures]",
            400, 2.0, 0.0, pe, 0.10, 100.0);
    }

    // H11: Triangle PE -- constant-velocity reversals
    // Common in mounts with friction-limited motor velocity; the worm tracks
    // linearly then reverses cleanly. Discontinuity in velocity (not position).
    {
        PEParams pe;
        pe.T = 30.0; pe.A1 = 4.0; pe.waveform = PE_TRIANGLE;
        runHarmonicScenario(
            "H11: Triangle PE  T=30s  A=4.0\"  frames=400  noise=0.10\"\n"
            "    [constant-velocity ramps with mid-period sign reversal; tests transient response at direction change]",
            400, 2.0, 0.0, pe, 0.10, 100.0);
    }

    // H12: Half-rectified PE -- single-sided drive train error
    // Half-period zero, half-period sine. Strong DC + harmonic content. The
    // MPC's IMP, which models pure sinusoidal d1/d2, cannot match the DC
    // component without help from the motor state.
    {
        PEParams pe;
        pe.T = 30.0; pe.A1 = 4.0; pe.waveform = PE_HALFRECT;
        runHarmonicScenario(
            "H12: Half-rectified PE  T=30s  A=4.0\"  frames=400  noise=0.10\"\n"
            "    [single-sided worm-tooth error; large DC bias plus harmonic spectrum]",
            400, 2.0, 0.0, pe, 0.10, 100.0);
    }

    // H13: True backlash (direction-dependent deadband)
    // Worm-and-wheel gear play: when motor reverses, axis stays put until
    // the +/- backlash gap is traversed. Different from H4's hysteresis lag
    // (which is a signed force offset). MPCSolver has a backlash punch-through
    // path that triggers on commanded direction changes -- this scenario tests
    // it. Compliance is mild (Ks=400) so the backlash effect is isolated.
    {
        PEParams pe;
        pe.T = 30.0; pe.A1 = 4.0; pe.A2 = 1.6; pe.A3 = 0.8;
        Sim2MassParams sim;
        sim.Ks = 400.0;     // Stiffer than H7/H8 to isolate backlash effect
        sim.Bs = 4.0;       // More damped
        sim.stiction = 0.0; // No stiction
        sim.coulomb = 0.0;
        sim.backlash = 1.0; // 1 arcsec of gear play (typical for hobby mounts)
        runComplianceScenario(
            "H13: True backlash  T=30s  A1=4.0\" A2=1.6\" A3=0.8\"  Ks=400.0  backlash=1.0\"\n"
            "    [direction-dependent deadband at each PE zero-crossing; exercises punch-through compensation]",
            400, 2.0, 0.0, pe, 0.10, 100.0, sim);
    }

    // H14: Period drift -- worm period changes slowly over the run
    // Real harmonic-drive/worm-gear mounts drift in period with temperature.
    // ~1-5% per hour is typical. Test uses 10% drift over the run to make
    // adaptation visible. MPC's FFT updates every 10 frames; GPG runs
    // gradient-based period inference. Both should adapt; the question is
    // how cleanly.
    // T0=30s, Tdot=0.0025 s/s -> T grows from 30s to ~33s over 1200s = 600 frames.
    {
        PEParams pe;
        pe.T = 30.0; pe.Tdot = 0.0025; pe.A1 = 4.0;
        runHarmonicScenario(
            "H14: Period drift  T0=30s -> ~33s over 600 frames  A=4.0\"  noise=0.10\"\n"
            "    [worm period drifts 10% over the run; tests FFT (MPC) vs gradient (GPG) re-tracking]",
            600, 2.0, 0.0, pe, 0.10, 30.0);
    }

    // H15: Mid-run plant change -- step disturbance at frame 200 on top of
    // ongoing PE. Simulates focus shift, refraction change, or worm contact
    // zone change. Tests adaptive recovery: pre-step settled performance vs
    // post-step transient.
    {
        PEParams pe;
        pe.T = 30.0; pe.A1 = 5.0;
        runH5Scenario(
            "H15: Mid-run amplitude step  T=30s  A1=5.0\"  +3\" step at fr200\n"
            "    [adaptive recovery test: sudden plant change requiring re-learning]",
            400, 2.0, pe, 0.10, 30.0, 200, +3.0);
    }

    printf("\n");
    return 0;
}
