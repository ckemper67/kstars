/*
    SPDX-FileCopyrightText: 2024 Christian Kemper <ckemper@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

#pragma once

// Pure algorithmic core of the 4-quadrant DONUTS guider.
// No Qt, no FITSData -- only STL and PocketFFT.
// Included by both donutsguider.cpp (Qt/production) and standalone tests.

#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <limits>
#include <utility>

#include "pocketfft_hdronly.h"

namespace DonutsCore
{

// Algorithm tuning constants -- change here to affect all call sites.
static constexpr double TUKEY_ALPHA      = 0.15;  // fraction of window that is tapered
static constexpr double LP_CUTOFF        = 0.75;  // fraction of Nyquist kept (low-pass)
static constexpr double PHASE_MIN_NORM   = 1e-9;  // avoid divide-by-zero in phase normalisation
static constexpr double SUBPIX_MIN_DEN   = 1e-6;  // parabolic peak must have |denom| > this
static constexpr double SNR_EPSILON      = 1e-9;  // floor added to SNR stddev denominator
static constexpr double SIGMA_THRESHOLD  = 3.0;   // pixel must be median + N*sigma to count
static constexpr double QUADRANT_OVERLAP = 0.15;  // fractional overlap between adjacent quadrants
static constexpr double SPIKE_RATIO      = 4.0;   // neighbour average multiplier for spike test
static constexpr double TIKHONOV_THETA   = 1e-3;  // Tikhonov regularisation on rotation DoF
static constexpr double DET_MIN          = 1e-9;  // degenerate system guard

struct Quadrant
{
    std::vector<double> xProf, yProf;
    double xc { 0 };
    double yc { 0 };
};

struct ProfileSet
{
    Quadrant q[4];
    int width  { 0 };
    int height { 0 };
};

struct Transform
{
    double dx     { 0 };
    double dy     { 0 };
    double dtheta { 0 };
    double snr    { 0 };
};

// Remove obvious hot-pixel spikes from a 1D projection profile.
inline void filterSpikes(std::vector<double> &prof)
{
    if (prof.size() < 3) return;
    std::vector<double> clean = prof;
    for (size_t i = 1; i + 1 < prof.size(); ++i)
    {
        double avg = (prof[i - 1] + prof[i + 1]) / 2.0;
        if (prof[i] > avg * SPIKE_RATIO)
            clean[i] = avg;
    }
    prof = clean;
}

// FFT-based 1D phase correlation with parabolic sub-pixel interpolation.
// Returns {shift_pixels, correlation_SNR}.
inline std::pair<double, double> correlate(
    const std::vector<double> &ref, const std::vector<double> &curr)
{
    size_t n = ref.size();
    if (n < 4) return {0.0, 0.0};

    // Tukey window to reduce spectral leakage.
    std::vector<double> r(n), c(n);
    for (size_t i = 0; i < n; ++i)
    {
        double w = 1.0;
        if (i < TUKEY_ALPHA * (n - 1) / 2.0)
            w = 0.5 * (1.0 + std::cos(M_PI * (2.0 * i / (TUKEY_ALPHA * (n - 1)) - 1.0)));
        else if (i > (n - 1) * (1.0 - TUKEY_ALPHA / 2.0))
            w = 0.5 * (1.0 + std::cos(M_PI * (2.0 * i / (TUKEY_ALPHA * (n - 1)) - 2.0 / TUKEY_ALPHA + 1.0)));
        r[i] = ref[i]  * w;
        c[i] = curr[i] * w;
    }

    using namespace pocketfft;
    shape_t  shape{n};
    stride_t si{sizeof(double)}, sc{sizeof(std::complex<double>)};
    std::vector<std::complex<double>> rf(n / 2 + 1), cf(n / 2 + 1), cp(n / 2 + 1);

    r2c(shape, si, sc, 0, FORWARD, r.data(), rf.data(), 1.0);
    r2c(shape, si, sc, 0, FORWARD, c.data(), cf.data(), 1.0);

    for (size_t i = 0; i < rf.size(); ++i)
    {
        cp[i] = cf[i] * std::conj(rf[i]);
        double m = std::abs(cp[i]);
        // Low-pass rectangular cutoff to suppress high-frequency interpolation noise.
        double filter = (i < rf.size() * LP_CUTOFF) ? 1.0 : 0.0;
        if (m > PHASE_MIN_NORM) cp[i] = (cp[i] / m) * filter;
    }

    std::vector<double> out(n);
    c2r(shape, sc, si, 0, BACKWARD, cp.data(), out.data(), 1.0 / n);

    auto   max_it = std::max_element(out.begin(), out.end());
    int    peak   = std::distance(out.begin(), max_it);

    double avg = 0.0;
    for (double v : out) avg += v;
    avg /= n;
    double var = 0.0;
    for (double v : out) var += (v - avg) * (v - avg);
    double snr = (*max_it - avg) / (std::sqrt(var / n) + SNR_EPSILON);

    // Parabolic sub-pixel refinement -- works even when side-lobes are negative.
    double shift = peak;
    if (peak > 0 && peak < (int)n - 1)
    {
        double y1 = out[peak - 1], y2 = out[peak], y3 = out[peak + 1];
        double den = y1 - 2.0 * y2 + y3;
        if (std::abs(den) > SUBPIX_MIN_DEN)
            shift += 0.5 * (y1 - y3) / den;
    }
    if (shift > (double)n / 2.0) shift -= n;
    return {shift, snr};
}

// Build quadrant projection profiles from a row-major double buffer.
// median/stddev/clip are pre-computed by the caller (allows bit-depth-aware clipping).
inline ProfileSet buildProfiles(
    const double *buf, int w, int h, double median, double stddev, double clip)
{
    ProfileSet p;
    p.width  = w;
    p.height = h;

    const double threshold = median + SIGMA_THRESHOLD * stddev;
    const int    qw        = static_cast<int>(w * (0.5 + QUADRANT_OVERLAP / 2.0));
    const int    qh        = static_cast<int>(h * (0.5 + QUADRANT_OVERLAP / 2.0));

    struct Range { int x0, y0, x1, y1; } qr[4] = {
        {0,      0,      qw, qh},  // TL
        {w - qw, 0,      w,  qh},  // TR
        {0,      h - qh, qw, h},   // BL
        {w - qw, h - qh, w,  h}    // BR
    };

    for (int i = 0; i < 4; ++i)
    {
        p.q[i].xProf.assign(qr[i].x1 - qr[i].x0, 0.0);
        p.q[i].yProf.assign(qr[i].y1 - qr[i].y0, 0.0);
    }

    for (int y = 0; y < h; ++y)
    {
        for (int x = 0; x < w; ++x)
        {
            double val = buf[y * w + x];
            if (val < threshold) continue;
            if (val > clip) val = clip;
            val -= median;

            for (int i = 0; i < 4; ++i)
            {
                if (x >= qr[i].x0 && x < qr[i].x1 && y >= qr[i].y0 && y < qr[i].y1)
                {
                    p.q[i].xProf[x - qr[i].x0] += val;
                    p.q[i].yProf[y - qr[i].y0] += val;
                }
            }
        }
    }

    // Filter spikes first, then compute centroids from cleaned profiles so that
    // hot pixels do not bias the rotation lever arms used by solveTransform.
    const double mx = w / 2.0, my = h / 2.0;
    for (int i = 0; i < 4; ++i)
    {
        filterSpikes(p.q[i].xProf);
        filterSpikes(p.q[i].yProf);

        double xFlux = 0.0, xMom = 0.0;
        for (int k = 0; k < (int)p.q[i].xProf.size(); ++k)
        {
            xFlux += p.q[i].xProf[k];
            xMom  += (qr[i].x0 + k - mx) * p.q[i].xProf[k];
        }
        if (xFlux > 0.0) p.q[i].xc = xMom / xFlux;

        double yFlux = 0.0, yMom = 0.0;
        for (int k = 0; k < (int)p.q[i].yProf.size(); ++k)
        {
            yFlux += p.q[i].yProf[k];
            yMom  += (qr[i].y0 + k - my) * p.q[i].yProf[k];
        }
        if (yFlux > 0.0) p.q[i].yc = yMom / yFlux;
    }
    return p;
}

// Weighted Least Squares 3-DoF rigid-body solver (dx, dy, dtheta).
// Transform::snr holds the worst per-quadrant correlation SNR.
// Returns a zero Transform (snr==0) if the system is degenerate.
inline Transform solveTransform(const ProfileSet &refP, const ProfileSet &currP)
{
    Transform result;

    double m11 = 0, m13 = 0, m22 = 0, m23 = 0, m33 = TIKHONOV_THETA;
    double b1  = 0, b2  = 0, b3  = 0;
    double minSNR = std::numeric_limits<double>::max();

    for (int i = 0; i < 4; ++i)
    {
        auto rx = correlate(refP.q[i].xProf, currP.q[i].xProf);
        auto ry = correlate(refP.q[i].yProf, currP.q[i].yProf);

        minSNR = std::min({minSNR, rx.second, ry.second});

        const double xi = refP.q[i].xc;
        const double yi = refP.q[i].yc;
        const double wx = rx.second * rx.second;
        const double wy = ry.second * ry.second;

        m11 += wx;
        m13 += -yi * wx;
        m22 += wy;
        m23 += xi * wy;
        m33 += (yi * yi * wx) + (xi * xi * wy);
        b1  += rx.first * wx;
        b2  += ry.first * wy;
        b3  += (-yi * rx.first * wx) + (xi * ry.first * wy);
    }

    result.snr = minSNR;

    const double det = m11 * (m22 * m33 - m23 * m23) - m13 * (m22 * m13);
    if (std::abs(det) < DET_MIN) return result;

    result.dx     = (b1 * (m22 * m33 - m23 * m23) + m13 * (b2 * m23 - b3 * m22)) / det;
    result.dy     = (b2 * (m11 * m33 - m13 * m13) + m23 * (b1 * m13 - b3 * m11)) / det;
    result.dtheta = (m11 * (m22 * b3 - b2 * m23) + m13 * (0.0 - b1 * m22))       / det;
    return result;
}

} // namespace DonutsCore
