/*
    SPDX-FileCopyrightText: 2026 Christian Kemper <ckemper@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

// Isolated translation unit for SIMD image warping.
// Kept separate from donuts.cpp so arm_neon.h / immintrin.h cannot affect
// the floating-point contraction mode of the double-precision phase correlation.

#include "simd_warp.h"

#include <algorithm>
#include <cmath>
#include <utility>

#ifdef __ARM_NEON
#include <arm_neon.h>
static inline float hsum_neon(float32x4_t v)
{
#ifdef __aarch64__
    return vaddvq_f32(v);
#else
    float32x2_t t = vadd_f32(vget_high_f32(v), vget_low_f32(v));
    return vget_lane_f32(vpadd_f32(t, t), 0);
#endif
}
#endif

#ifdef __SSE2__
#include <immintrin.h>
static inline float hsum_ps(__m128 v)
{
#ifdef __SSE3__
    v = _mm_hadd_ps(v, v);
    v = _mm_hadd_ps(v, v);
#else
    __m128 shuf = _mm_shuffle_ps(v, v, _MM_SHUFFLE(2, 3, 0, 1));
    v           = _mm_add_ps(v, shuf);
    shuf        = _mm_movehl_ps(shuf, v);
    v           = _mm_add_ss(v, shuf);
#endif
    return _mm_cvtss_f32(v);
}
#endif

namespace Donuts
{

// Catmull-Rom weights for fractional position t in [0,1).
// w[0]=weight for ix-1, w[1]=ix, w[2]=ix+1, w[3]=ix+2.
static inline void cubicWeights4(float t, float w[4])
{
    float t2 = t * t, t3 = t2 * t;
    w[0] = -0.5f*t3 + t2 - 0.5f*t;
    w[1] =  1.5f*t3 - 2.5f*t2 + 1.0f;
    w[2] = -1.5f*t3 + 2.0f*t2 + 0.5f*t;
    w[3] =  0.5f*t3 - 0.5f*t2;
}

// Returns the x range [lo, hi] (inclusive) where all 16 source pixels in the
// 4x4 Catmull-Rom kernel are in-bounds. Returns {0, -1} when no such range exists.
// Source mapping is linear in x: sx(x) = m.a*x + sx0, sy(x) = m.c*x + sy0.
static std::pair<int,int> interiorXRange(double sx0, double sy0,
                                         double dsx_dx, double dsy_dx,
                                         int w, int h)
{
    double lo = 0.0, hi = static_cast<double>(w - 1);

    if (std::abs(dsx_dx) > 1e-12)
    {
        double a = (1.0             - sx0) / dsx_dx;
        double b = (w - 2.0 - 1e-9 - sx0) / dsx_dx;
        if (dsx_dx > 0) { lo = std::max(lo, a); hi = std::min(hi, b); }
        else            { lo = std::max(lo, b); hi = std::min(hi, a); }
    }
    else if (sx0 < 1.0 || sx0 >= w - 2.0) return {0, -1};

    if (std::abs(dsy_dx) > 1e-12)
    {
        double a = (1.0             - sy0) / dsy_dx;
        double b = (h - 2.0 - 1e-9 - sy0) / dsy_dx;
        if (dsy_dx > 0) { lo = std::max(lo, a); hi = std::min(hi, b); }
        else            { lo = std::max(lo, b); hi = std::min(hi, a); }
    }
    else if (sy0 < 1.0 || sy0 >= h - 2.0) return {0, -1};

    int ilo = std::max(0,     static_cast<int>(std::ceil(lo)));
    int ihi = std::min(w - 1, static_cast<int>(hi));
    return (ilo <= ihi) ? std::make_pair(ilo, ihi) : std::make_pair(0, -1);
}

void warpRows1(const std::vector<float> &src, int w, int h,
               const AffineMatrix &m,
               std::vector<float> &dst,
               int ylo, int yhi)
{
    for (int y = ylo; y < yhi; ++y)
    {
        const double sx0 = m.b * y + m.tx;
        const double sy0 = m.d * y + m.ty;

        auto [xInLo, xInHi] = interiorXRange(sx0, sy0, m.a, m.c, w, h);

        // Slow (clamped) path for border pixels where the 4x4 kernel clips.
        auto slowPixel = [&](int x)
        {
            const double sx = m.a * x + sx0;
            const double sy = m.c * x + sy0;
            const std::size_t k = static_cast<std::size_t>(y) * w + x;
            if (sx < 0.0 || sx >= w || sy < 0.0 || sy >= h) { dst[k] = 0.0f; return; }
            const int   ix = static_cast<int>(sx), iy = static_cast<int>(sy);
            const float fx = static_cast<float>(sx - ix), fy = static_cast<float>(sy - iy);
            float wx[4], wy[4];
            cubicWeights4(fx, wx);
            cubicWeights4(fy, wy);
            float v = 0.0f;
            for (int dr = -1; dr <= 2; ++dr)
            {
                const int   py  = std::max(0, std::min(h - 1, iy + dr));
                const float wyd = wy[dr + 1];
                for (int dc = -1; dc <= 2; ++dc)
                {
                    const int px = std::max(0, std::min(w - 1, ix + dc));
                    v += wx[dc + 1] * wyd * src[static_cast<std::size_t>(py) * w + px];
                }
            }
            dst[k] = v;
        };

        for (int x = 0; x < xInLo; ++x) slowPixel(x);

        // Fast interior path: no bounds clamping, SIMD 4-wide accumulation.
        for (int x = xInLo; x <= xInHi; ++x)
        {
            const double      sx  = m.a * x + sx0;
            const double      sy  = m.c * x + sy0;
            const std::size_t k   = static_cast<std::size_t>(y) * w + x;
            const int         ix  = static_cast<int>(sx), iy = static_cast<int>(sy);
            const float       fx  = static_cast<float>(sx - ix), fy = static_cast<float>(sy - iy);
            float wx[4], wy[4];
            cubicWeights4(fx, wx);
            cubicWeights4(fy, wy);
#ifdef __ARM_NEON
            const float32x4_t wxv = vld1q_f32(wx);
            float32x4_t accV = vdupq_n_f32(0.0f);
            for (int dr = 0; dr < 4; ++dr)
            {
                const std::size_t base = static_cast<std::size_t>(iy + dr - 1) * w + (ix - 1);
                const float32x4_t wyd  = vdupq_n_f32(wy[dr]);
                accV = vmlaq_f32(accV, vmulq_f32(vld1q_f32(src.data() + base), wxv), wyd);
            }
            dst[k] = hsum_neon(accV);
#elif defined(__SSE2__)
            const __m128 wxv = _mm_loadu_ps(wx);
            __m128 accV = _mm_setzero_ps();
            for (int dr = 0; dr < 4; ++dr)
            {
                const std::size_t base = static_cast<std::size_t>(iy + dr - 1) * w + (ix - 1);
                const __m128 wxw  = _mm_mul_ps(wxv, _mm_set1_ps(wy[dr]));
                accV = _mm_add_ps(accV, _mm_mul_ps(_mm_loadu_ps(src.data() + base), wxw));
            }
            dst[k] = hsum_ps(accV);
#else
            float v = 0.0f;
            for (int dr = 0; dr < 4; ++dr)
            {
                const std::size_t base = static_cast<std::size_t>(iy + dr - 1) * w + (ix - 1);
                const float       wyd  = wy[dr];
                for (int dc = 0; dc < 4; ++dc)
                    v += wx[dc] * wyd * src[base + dc];
            }
            dst[k] = v;
#endif
        }

        for (int x = xInHi + 1; x < w; ++x) slowPixel(x);
    }
}

} // namespace Donuts
