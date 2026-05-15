/*
    SPDX-FileCopyrightText: 2024 Christian Kemper <ckemper@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

#pragma once
#include <memory>

// Public API for the DONUTS translation guiding algorithm.
// No Qt, no FITSData, no pocketfft in this header.

namespace Donuts
{

struct Config
{
    double tukeyAlpha     = 0.15;  // spectral leakage taper (fraction of profile)
    double lpCutoff       = 0.75;  // low-pass cutoff as fraction of Nyquist
    double sigmaThreshold = 3.0;   // pixel selection: median + N * stddev
    double spikeRatio     = 4.0;   // hot-pixel suppression: value > ratio * neighbour avg
    double clipSigmas     = 20.0;  // saturation clip: median + N * stddev
};

struct Transform
{
    double dx     = 0;  // x translation in pixels (positive = right)
    double dy     = 0;  // y translation in pixels (positive = down)
    double dtheta = 0;  // unused in translation-only mode; always 0
    double snr    = 0;  // minimum of x and y correlation SNR (< 3 = unreliable)

    bool valid() const { return snr >= 3.0; }
};

// 2-DoF translation guider using global 1-D phase correlation.
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
