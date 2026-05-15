/*
    SPDX-FileCopyrightText: 2026 Christian Kemper <ckemper@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

#include "donuts.h"

// pocketfft is a private implementation detail; it never appears in donuts.h.
#include "../internalguide/pocketfft_hdronly.h"

#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <utility>
#include <cstddef>

namespace Donuts
{

// ---------------------------------------------------------------------------
// Internal types
// ---------------------------------------------------------------------------

struct Profiles
{
    std::vector<double> xProf;  // sum of bright pixels along each column -> x shift
    std::vector<double> yProf;  // sum of bright pixels along each row    -> y shift
};

struct FrameStats
{
    double median { 0 };
    double stddev { 0 };
    double clip   { 0 };
};

// ---------------------------------------------------------------------------
// Statistics
// ---------------------------------------------------------------------------

static FrameStats computeStats(const double *buf, int n)
{
    if (n <= 0) return {};

    // Welford one-pass mean + variance for stddev; track max for clip.
    double mean = 0, m2 = 0, maxVal = 0;
    for (int i = 0; i < n; ++i)
    {
        double delta = buf[i] - mean;
        mean += delta / (i + 1);
        m2   += delta * (buf[i] - mean);
        if (buf[i] > maxVal) maxVal = buf[i];
    }
    double stddev = (n > 1) ? std::sqrt(m2 / (n - 1)) : 0.0;

    // Median via nth_element on a subsample (capped at 64 K for speed).
    const int SAMPLE = std::min(n, 65536);
    std::vector<double> sample(SAMPLE);
    double step = static_cast<double>(n) / SAMPLE;
    for (int i = 0; i < SAMPLE; ++i)
        sample[i] = buf[static_cast<int>(i * step)];
    auto mid = sample.begin() + SAMPLE / 2;
    std::nth_element(sample.begin(), mid, sample.end());
    double median = *mid;

    return { median, stddev, maxVal * 0.95 };
}

// ---------------------------------------------------------------------------
// Algorithm helpers
// ---------------------------------------------------------------------------

static void filterSpikes(std::vector<double> &prof, double spikeRatio)
{
    if (prof.size() < 3) return;
    std::vector<double> clean = prof;
    for (std::size_t i = 1; i + 1 < prof.size(); ++i)
    {
        double avg = (prof[i - 1] + prof[i + 1]) / 2.0;
        if (prof[i] > avg * spikeRatio)
            clean[i] = avg;
    }
    prof = clean;
}

// FFT-based 1-D phase correlation with parabolic sub-pixel refinement.
// Returns {shift_pixels, correlation_SNR}.
static std::pair<double, double> correlate(
    const std::vector<double> &ref,
    const std::vector<double> &curr,
    double tukeyAlpha, double lpCutoff)
{
    const std::size_t n = ref.size();
    if (n < 4) return {0.0, 0.0};

    // Tukey window to reduce spectral leakage.
    std::vector<double> r(n), c(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        double w = 1.0;
        double t = static_cast<double>(i) / (n - 1);
        if (t < tukeyAlpha / 2.0)
            w = 0.5 * (1.0 + std::cos(M_PI * (2.0 * t / tukeyAlpha - 1.0)));
        else if (t > 1.0 - tukeyAlpha / 2.0)
            w = 0.5 * (1.0 + std::cos(M_PI * (2.0 * t / tukeyAlpha - 2.0 / tukeyAlpha + 1.0)));
        r[i] = ref[i]  * w;
        c[i] = curr[i] * w;
    }

    using namespace pocketfft;
    shape_t  shape{n};
    stride_t si{sizeof(double)}, sc{sizeof(std::complex<double>)};
    std::vector<std::complex<double>> rf(n / 2 + 1), cf(n / 2 + 1), cp(n / 2 + 1);

    r2c(shape, si, sc, 0, FORWARD, r.data(), rf.data(), 1.0);
    r2c(shape, si, sc, 0, FORWARD, c.data(), cf.data(), 1.0);

    for (std::size_t i = 0; i < rf.size(); ++i)
    {
        cp[i] = cf[i] * std::conj(rf[i]);
        double m = std::abs(cp[i]);
        const double k_full = rf.size() * lpCutoff * 0.9;
        const double k_cut  = rf.size() * lpCutoff;
        double ki = static_cast<double>(i);
        double filter;
        if (ki <= k_full)
            filter = 1.0;
        else if (ki >= k_cut)
            filter = 0.0;
        else
            filter = 0.5 * (1.0 + std::cos(M_PI * (ki - k_full) / (k_cut - k_full)));
        if (m > 1e-9) cp[i] = (cp[i] / m) * filter;
    }

    std::vector<double> out(n);
    c2r(shape, sc, si, 0, BACKWARD, cp.data(), out.data(), 1.0 / n);

    auto   max_it = std::max_element(out.begin(), out.end());
    int    peak   = static_cast<int>(std::distance(out.begin(), max_it));
    double pval   = *max_it;

    double avg = 0.0, var = 0.0;
    for (double v : out) avg += v;
    avg /= n;
    for (double v : out) var += (v - avg) * (v - avg);
    double snr = (pval - avg) / (std::sqrt(var / n) + 1e-9);

    double shift = peak;
    {
        // Use wrap-around neighbors so boundary peaks (e.g. index 0 for a small
        // negative shift) are refined correctly.
        int    pm1 = (peak == 0) ? static_cast<int>(n) - 1 : peak - 1;
        int    pp1 = (peak == static_cast<int>(n) - 1) ? 0 : peak + 1;
        double y1  = out[pm1], y2 = out[peak], y3 = out[pp1];
        double den = y1 - 2.0 * y2 + y3;
        if (std::abs(den) > 1e-6)
            shift += 0.5 * (y1 - y3) / den;
    }
    if (shift > static_cast<double>(n) / 2.0) shift -= n;

    return {shift, snr};
}

// ---------------------------------------------------------------------------
// Profile building
// ---------------------------------------------------------------------------

static Profiles buildProfiles(
    const double *buf, int w, int h,
    const FrameStats &stats, const Config &cfg)
{
    const double threshold = stats.median + cfg.sigmaThreshold * stats.stddev;

    Profiles p;
    p.xProf.assign(w, 0.0);
    p.yProf.assign(h, 0.0);

    for (int y = 0; y < h; ++y)
    {
        for (int x = 0; x < w; ++x)
        {
            double val = buf[y * w + x];
            if (val < threshold) continue;
            if (val > stats.clip) val = stats.clip;
            val -= stats.median;
            p.xProf[x] += val;
            p.yProf[y] += val;
        }
    }

    filterSpikes(p.xProf, cfg.spikeRatio);
    filterSpikes(p.yProf, cfg.spikeRatio);
    return p;
}

// ---------------------------------------------------------------------------
// Guider implementation
// ---------------------------------------------------------------------------

struct Guider::Impl
{
    Config     cfg;
    Profiles   refProfiles;
    FrameStats refStats;
    bool       hasRef { false };
};

Guider::Guider(Config cfg)
    : m_impl(std::make_unique<Impl>())
{
    m_impl->cfg = cfg;
}

Guider::~Guider() = default;

void Guider::setReference(const double *pixels, int width, int height)
{
    m_impl->refStats    = computeStats(pixels, width * height);
    m_impl->refProfiles = buildProfiles(pixels, width, height, m_impl->refStats, m_impl->cfg);
    m_impl->hasRef = true;
}

Transform Guider::measure(const double *pixels, int width, int height)
{
    if (!m_impl->hasRef) return {};
    auto curr = buildProfiles(pixels, width, height, m_impl->refStats, m_impl->cfg);

    auto rx = correlate(m_impl->refProfiles.xProf, curr.xProf,
                        m_impl->cfg.tukeyAlpha, m_impl->cfg.lpCutoff);
    auto ry = correlate(m_impl->refProfiles.yProf, curr.yProf,
                        m_impl->cfg.tukeyAlpha, m_impl->cfg.lpCutoff);

    Transform t;
    t.dx  = rx.first;
    t.dy  = ry.first;
    t.snr = std::min(rx.second, ry.second);
    return t;
}

void Guider::reset()
{
    m_impl->refProfiles = Profiles{};
    m_impl->hasRef = false;
}

bool Guider::hasReference() const
{
    return m_impl->hasRef;
}

} // namespace Donuts
