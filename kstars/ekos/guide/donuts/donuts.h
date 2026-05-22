/*
    SPDX-FileCopyrightText: 2026 Christian Kemper <ckemper@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

#pragma once
#include <memory>

// 4-quadrant image registrar inspired by DONUTS (McCormac et al. 2013, PASP 125, 548).
// Departs from the original in three ways:
//   - 4-quadrant profile split instead of single full-frame projections
//   - Phase-only correlation (Kuglin & Hines 1975) instead of standard cross-correlation
//   - WLS solver for rotation (dtheta) and optionally isotropic scale in addition to translation
// No Qt, no FITSData, no pocketfft in this header.

namespace Donuts
{

struct Config
{
    double tukeyAlpha      = 0.15;   // spectral leakage taper (fraction of profile)
    double lpCutoff        = 0.75;   // low-pass cutoff as fraction of Nyquist
    double sigmaThreshold  = 3.0;    // pixel selection: median + N * stddev
    double quadrantOverlap = 0.15;   // fractional overlap between adjacent quadrants
    double tikhonov        = 1e-3;   // Tikhonov regularisation on rotation/scale DoF
    double spikeRatio            = 4.0;   // hot-pixel suppression: value > ratio * neighbour avg
    double streakRidgeThreshold  = 8.0;   // streak suppression: peak > ratio * max-outside-guard-band
    bool   detectScale     = false;  // add isotropic scale as a 4th DoF (for focus tracking)
};

// 2x3 affine matrix in row-major order.
// Maps destination pixel (xd, yd) to source pixel (xs, ys):
//   xs = a*xd + b*yd + tx
//   ys = c*xd + d*yd + ty
struct AffineMatrix
{
    double a = 1, b = 0, tx = 0;
    double c = 0, d = 1, ty = 0;

    static AffineMatrix identity() { return {}; }
};

struct Transform
{
    double dx     = 0;    // x translation in pixels (positive = right)
    double dy     = 0;    // y translation in pixels (positive = down)
    double dtheta = 0;    // rotation in radians; positive = CCW in math coords (y-up),
                          // which is visually CW on screen (y-down image coords).
    double scale  = 1.0;  // isotropic scale factor relative to reference (1.0 = unchanged).
                          // Only populated when Config::detectScale is true.
                          // Focus drift: delta_focus_um = (scale - 1.0) * focal_length_mm * 1e3
    double snr    = 0;    // minimum per-quadrant correlation SNR (< 3 = unreliable)

    bool valid() const { return snr >= 3.0; }
};

// 4-DoF image registrar using 4-quadrant 1-D phase-only correlation.
// Solves for translation (dx, dy), rotation (dtheta), and optionally scale.
// Not thread-safe; use external locking when sharing across threads.
class Registrar
{
public:
    explicit Registrar(Config cfg = {});
    ~Registrar();

    // Store frame as registration reference.
    // pixels: row-major double buffer, width * height elements.
    // Any intensity scale is accepted; the library computes its own statistics.
    // Resets the accumulated field-rotation to zero.
    void setReference(const double *pixels, int width, int height);

    // Measure the transform of pixels relative to the stored reference.
    // Solves for dx, dy, dtheta always; also solves for scale when Config::detectScale is true.
    // Returns a Transform with snr < 3.0 when the result is unreliable.
    //
    // Field rotation is tracked internally: each call de-rotates the frame by the accumulated
    // rotation so the algorithm always sees only a small residual.  Transform::dtheta on return
    // is the total rotation from the original reference frame.
    Transform measure(const double *pixels, int width, int height);

    // Reset registration state (clears reference and accumulated rotation).
    void reset();
    bool hasReference() const;

private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};

} // namespace Donuts
