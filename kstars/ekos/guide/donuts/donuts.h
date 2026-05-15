/*
    SPDX-FileCopyrightText: 2026 Christian Kemper <ckemper@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

#pragma once
#include <memory>

// 4-quadrant guider inspired by DONUTS (McCormac et al. 2013, PASP 125, 548).
// Departs from the original in three ways:
//   - 4-quadrant profile split instead of single full-frame projections
//   - Phase-only correlation (Kuglin & Hines 1975) instead of standard cross-correlation
//   - Weighted least-squares solver for rotation (dtheta) in addition to translation
// No Qt, no FITSData, no pocketfft in this header.

namespace Donuts
{

struct Config
{
    double tukeyAlpha      = 0.15;  // spectral leakage taper (fraction of profile)
    double lpCutoff        = 0.75;  // low-pass cutoff as fraction of Nyquist
    double sigmaThreshold  = 3.0;   // pixel selection: median + N * stddev
    double quadrantOverlap = 0.15;  // fractional overlap between adjacent quadrants
    double tikhonov        = 1e-3;  // Tikhonov regularisation on rotation DoF
    double spikeRatio      = 4.0;   // hot-pixel suppression: value > ratio * neighbour avg
};

struct Transform
{
    double dx     = 0;  // x translation in pixels (positive = right)
    double dy     = 0;  // y translation in pixels (positive = down)
    double dtheta = 0;  // rotation in radians; positive = CCW in math coords (y-up),
                        // which is visually CW on screen (y-down image coords).
                        // Wire to a rotator with the appropriate sign for your mount.
    double snr    = 0;  // minimum per-quadrant correlation SNR (< 3 = unreliable)

    bool valid() const { return snr >= 3.0; }
};

// 3-DoF rigid-body guider using 4-quadrant 1-D phase-only correlation.
// Not thread-safe; use external locking when sharing across threads.
class Guider
{
public:
    explicit Guider(Config cfg = {});
    ~Guider();

    // Store frame as guiding reference.
    // pixels: row-major double buffer, width * height elements.
    // Any intensity scale is accepted; the library computes its own statistics.
    void setReference(const double *pixels, int width, int height);

    // Measure the 3-DoF transform of pixels relative to the stored reference.
    // Returns a Transform with snr < 3.0 when the result is unreliable.
    Transform measure(const double *pixels, int width, int height);

    void reset();
    bool hasReference() const;

private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};

} // namespace Donuts
