/*
    SPDX-FileCopyrightText: 2026 Christian Kemper <ckemper@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

#pragma once
#include <memory>

// Translation registrar inspired by DONUTS (McCormac et al. 2013, PASP 125, 548).
// Departs from the original in one way: uses phase-only correlation
// (Kuglin & Hines 1975) instead of standard cross-correlation.
// Otherwise follows the full-frame 1-D projection approach of the paper.
// No Qt, no FITSData, no pocketfft in this header.

namespace Donuts
{

struct Config
{
    double tukeyAlpha           = 0.15;  // spectral leakage taper (fraction of profile)
    double lpCutoff             = 0.75;  // low-pass cutoff as fraction of Nyquist
    double sigmaThreshold       = 3.0;   // pixel selection: median + N * stddev
    double spikeRatio           = 4.0;   // hot-pixel suppression: value > ratio * neighbour avg
    double streakRidgeThreshold = 8.0;   // peak/surround ratio above which a profile bin is
                                         // zeroed as a satellite streak; 0 disables detection
};

struct Transform
{
    double dx  = 0;   // x translation in pixels (positive = right)
    double dy  = 0;   // y translation in pixels (positive = down)
    double snr = 0;   // minimum of x and y correlation SNR (< 3 = unreliable)

    bool valid() const { return snr >= 3.0; }
};

// 2-DoF image registrar using full-frame 1-D phase-only correlation.
// Not thread-safe; use external locking when sharing across threads.
class Registrar
{
public:
    explicit Registrar(Config cfg = {});
    ~Registrar();

    // Store frame as guiding reference.
    // pixels: row-major double buffer, width * height elements.
    // Any intensity scale is accepted; the library computes its own statistics.
    void setReference(const double *pixels, int width, int height);

    // Measure the 2-DoF translation of pixels relative to the stored reference.
    // Returns a Transform with snr < 3.0 when the result is unreliable.
    Transform measure(const double *pixels, int width, int height);

    void reset();
    bool hasReference() const;

private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};

} // namespace Donuts
