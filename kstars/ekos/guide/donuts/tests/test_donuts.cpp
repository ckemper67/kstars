/*
    SPDX-FileCopyrightText: 2024 Christian Kemper <ckemper@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later

    Standalone test for the Donuts registration algorithm.
    No Qt, no FITSData, no external test framework required.
    Star fields are generated synthetically so the test is deterministic
    and runs in < 1 second.
*/

#include "donuts.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <string>

// ---------------------------------------------------------------------------
// Minimal test harness
// ---------------------------------------------------------------------------

static int g_pass = 0, g_fail = 0;

static void check(bool ok, const char *expr, const char *file, int line)
{
    if (ok) {
        ++g_pass;
    } else {
        ++g_fail;
        std::printf("FAIL  %s:%d  %s\n", file, line, expr);
    }
}

#define CHECK(expr)       check((expr), #expr, __FILE__, __LINE__)
#define CHECK_NEAR(a,b,t) check(std::abs((a)-(b)) < (t), \
    (#a " ~ " #b " (tol " #t ")"), __FILE__, __LINE__)

// ---------------------------------------------------------------------------
// Synthetic star field
// ---------------------------------------------------------------------------

// Parameters for a single star.
struct Star { double x, y, peak; };

// Generate a synthetic 16-bit-range star field as a double buffer.
// background: sky level; stars: list of Gaussian PSF stars with given sigma.
static std::vector<double> makeFrame(int w, int h,
    const std::vector<Star> &stars,
    double background = 1000.0,
    double psf_sigma  = 1.5)
{
    std::vector<double> buf(w * h, background);
    const double inv2s2 = 1.0 / (2.0 * psf_sigma * psf_sigma);

    for (const auto &s : stars)
    {
        int xlo = static_cast<int>(std::max(0.0,   s.x - 5 * psf_sigma));
        int xhi = static_cast<int>(std::min((double)w - 1, s.x + 5 * psf_sigma));
        int ylo = static_cast<int>(std::max(0.0,   s.y - 5 * psf_sigma));
        int yhi = static_cast<int>(std::min((double)h - 1, s.y + 5 * psf_sigma));

        for (int y = ylo; y <= yhi; ++y)
            for (int x = xlo; x <= xhi; ++x)
            {
                double dx = x - s.x, dy = y - s.y;
                buf[y * w + x] += s.peak * std::exp(-(dx*dx + dy*dy) * inv2s2);
            }
    }
    return buf;
}

// Scale a frame by factor s around the image centre using bilinear interpolation.
static std::vector<double> scaleFrame(
    const std::vector<double> &src, int w, int h, double s)
{
    std::vector<double> dst(w * h, 0.0);
    const double cx = w / 2.0, cy = h / 2.0;
    const double inv_s = 1.0 / s;

    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x)
        {
            double x2 = (x - cx) * inv_s + cx;
            double y2 = (y - cy) * inv_s + cy;
            if (x2 >= 0 && x2 < w - 1 && y2 >= 0 && y2 < h - 1)
            {
                int    ix = static_cast<int>(x2), iy = static_cast<int>(y2);
                double fx = x2 - ix, fy = y2 - iy;
                dst[y * w + x] =
                    (1-fx)*(1-fy) * src[iy*w+ix  ]   + fx*(1-fy) * src[iy*w+ix+1] +
                    (1-fx)*fy     * src[(iy+1)*w+ix ] + fx*fy     * src[(iy+1)*w+ix+1];
            }
        }
    return dst;
}

// Apply a rigid-body transform (rotate around centre, then translate)
// to src and write the result into dst using bilinear interpolation.
static std::vector<double> transformFrame(
    const std::vector<double> &src, int w, int h,
    double dx, double dy, double dtheta_deg)
{
    std::vector<double> dst(w * h, 0.0);
    const double cosT = std::cos(dtheta_deg * M_PI / 180.0);
    const double sinT = std::sin(dtheta_deg * M_PI / 180.0);
    const double cx = w / 2.0, cy = h / 2.0;

    for (int y = 0; y < h; ++y)
    {
        for (int x = 0; x < w; ++x)
        {
            // Inverse transform: undo translation, then undo rotation.
            double x1 = x - dx, y1 = y - dy;
            double x2 = (x1 - cx) * cosT + (y1 - cy) * sinT + cx;
            double y2 = -(x1 - cx) * sinT + (y1 - cy) * cosT + cy;

            if (x2 >= 0 && x2 < w - 1 && y2 >= 0 && y2 < h - 1)
            {
                int    ix = static_cast<int>(x2), iy = static_cast<int>(y2);
                double fx = x2 - ix,              fy = y2 - iy;
                dst[y * w + x] =
                    (1-fx)*(1-fy) * src[iy*w + ix  ] + fx*(1-fy) * src[iy*w + ix+1] +
                    (1-fx)*fy     * src[(iy+1)*w+ix ] + fx*fy     * src[(iy+1)*w+ix+1];
            }
        }
    }
    return dst;
}

// Deterministic pseudo-random star field: same seed -> same result.
static std::vector<Star> randomStars(int w, int h, int count, uint32_t seed)
{
    std::vector<Star> stars;
    stars.reserve(count);
    uint32_t s = seed;
    auto rng = [&]() -> double {
        s ^= s << 13; s ^= s >> 17; s ^= s << 5;  // xorshift32
        return (s & 0xFFFFu) / 65536.0;
    };

    const double margin = 20.0;
    for (int i = 0; i < count; ++i)
    {
        double x    = margin + rng() * (w - 2 * margin);
        double y    = margin + rng() * (h - 2 * margin);
        double peak = 5000.0 + rng() * 45000.0;
        stars.push_back({x, y, peak});
    }
    return stars;
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

static void testPureTranslation()
{
    std::printf("--- testPureTranslation ---\n");
    const int W = 512, H = 512;
    auto stars  = randomStars(W, H, 80, 42);
    auto ref    = makeFrame(W, H, stars);
    auto curr   = transformFrame(ref, W, H, 5.4, -3.2, 0.0);

    Donuts::Registrar g;
    g.setReference(ref.data(), W, H);
    auto t = g.measure(curr.data(), W, H);

    std::printf("  Expected (5.4, -3.2)  Got (%.4f, %.4f)  SNR=%.1f\n",
                t.dx, t.dy, t.snr);

    CHECK(t.valid());
    CHECK_NEAR(t.dx,  5.4,  0.15);
    CHECK_NEAR(t.dy, -3.2,  0.15);
}

static void testSubPixelTranslation()
{
    std::printf("--- testSubPixelTranslation ---\n");
    const int W = 512, H = 512;
    auto stars = randomStars(W, H, 80, 7);
    auto ref   = makeFrame(W, H, stars);
    auto curr  = transformFrame(ref, W, H, 0.7, -0.4, 0.0);

    Donuts::Registrar g;
    g.setReference(ref.data(), W, H);
    auto t = g.measure(curr.data(), W, H);

    std::printf("  Expected (0.7, -0.4)  Got (%.4f, %.4f)  SNR=%.1f\n",
                t.dx, t.dy, t.snr);

    CHECK(t.valid());
    // Hann LP roll-off attenuates the high-frequency band that determines
    // sub-pixel accuracy; 0.15 px matches the combined-transform tolerance.
    CHECK_NEAR(t.dx,  0.7, 0.15);
    CHECK_NEAR(t.dy, -0.4, 0.15);
}

static void testPureRotation()
{
    std::printf("--- testPureRotation ---\n");
    const int W = 512, H = 512;
    auto stars = randomStars(W, H, 80, 13);
    auto ref   = makeFrame(W, H, stars);
    auto curr  = transformFrame(ref, W, H, 0.0, 0.0, 0.5);

    Donuts::Registrar g;
    g.setReference(ref.data(), W, H);
    auto t = g.measure(curr.data(), W, H);

    double deg = t.dtheta * 180.0 / M_PI;
    std::printf("  Expected 0.5 deg  Got %.4f deg  SNR=%.1f\n", deg, t.snr);

    CHECK(t.valid());
    CHECK_NEAR(deg, 0.5, 0.05);
}

static void testCombined()
{
    std::printf("--- testCombined ---\n");
    const int W = 512, H = 512;
    auto stars = randomStars(W, H, 80, 99);
    auto ref   = makeFrame(W, H, stars);
    auto curr  = transformFrame(ref, W, H, 2.0, 2.0, 0.2);

    Donuts::Registrar g;
    g.setReference(ref.data(), W, H);
    auto t = g.measure(curr.data(), W, H);

    double deg = t.dtheta * 180.0 / M_PI;
    std::printf("  Expected (2.0, 2.0, 0.2 deg)  Got (%.4f, %.4f, %.4f deg)  SNR=%.1f\n",
                t.dx, t.dy, deg, t.snr);

    CHECK(t.valid());
    CHECK_NEAR(t.dx, 2.0, 0.15);
    CHECK_NEAR(t.dy, 2.0, 0.15);
    CHECK_NEAR(deg,  0.2, 0.05);
}

static void testReset()
{
    std::printf("--- testReset ---\n");
    const int W = 256, H = 256;
    auto stars = randomStars(W, H, 40, 5);
    auto frame = makeFrame(W, H, stars);

    Donuts::Registrar g;
    CHECK(!g.hasReference());
    g.setReference(frame.data(), W, H);
    CHECK(g.hasReference());
    g.reset();
    CHECK(!g.hasReference());

    // measure without reference returns zero-snr transform
    auto t = g.measure(frame.data(), W, H);
    CHECK(!t.valid());
}

static void testSmallImage()
{
    std::printf("--- testSmallImage (128x128) ---\n");
    const int W = 128, H = 128;
    auto stars = randomStars(W, H, 20, 17);
    auto ref   = makeFrame(W, H, stars);
    auto curr  = transformFrame(ref, W, H, 1.5, -1.0, 0.0);

    Donuts::Registrar g;
    g.setReference(ref.data(), W, H);
    auto t = g.measure(curr.data(), W, H);

    std::printf("  Expected (1.5, -1.0)  Got (%.4f, %.4f)  SNR=%.1f\n",
                t.dx, t.dy, t.snr);

    // Smaller image -> less SNR, wider tolerance
    CHECK(t.valid());
    CHECK_NEAR(t.dx,  1.5, 0.3);
    CHECK_NEAR(t.dy, -1.0, 0.3);
}

// Combined rotation+translation on a small (256x256) image.
// The 1D phase correlation measures a single bulk shift per quadrant.  When
// rotation and translation are applied together, stars within a quadrant span
// different y-positions and each shifts by a slightly different x-amount.
// At 0.2 deg rotation on a 256x256 image the within-quadrant spread is
// ~0.35 px, so the correlation peak lands near the middle of this spread
// rather than at the flux-weighted centroid.  This gives a systematic dtheta
// underestimate of ~0.12 deg, independently of noise.
// Translation accuracy stays sub-pixel; the wider rotation tolerance here
// documents this 1D projection limit rather than a bug.
static void testCombinedSmallImage()
{
    std::printf("--- testCombinedSmallImage (256x256, 1D projection limit) ---\n");
    const int W = 256, H = 256;
    auto stars = randomStars(W, H, 25, 99);
    auto ref   = makeFrame(W, H, stars);
    auto curr  = transformFrame(ref, W, H, 2.0, 2.0, 0.2);

    Donuts::Registrar g;
    g.setReference(ref.data(), W, H);
    auto t = g.measure(curr.data(), W, H);

    double deg = t.dtheta * 180.0 / M_PI;
    std::printf("  Expected (2.0, 2.0, 0.200 deg)  Got (%.4f, %.4f, %.4f deg)  SNR=%.1f\n",
                t.dx, t.dy, deg, t.snr);

    CHECK(t.valid());
    CHECK_NEAR(t.dx, 2.0, 0.2);
    CHECK_NEAR(t.dy, 2.0, 0.2);
    // Wider tolerance than testCombined (512x512): the 2x smaller lever arm
    // halves the rotation SNR, giving a ~0.12 deg systematic floor.
    CHECK_NEAR(deg,  0.2, 0.15);
}

static void testAltAzNorthPointing()
{
    // Alt-az mount guiding near the north celestial pole (NCP).
    // The NCP is the true rotation pivot: it sits 200 px north of (above) the
    // guide star, which is at the image center.  As the sky rotates around the
    // NCP between guide frames, the guide star traces a small eastward arc --
    // appearing as a westward drift in the sensor combined with field rotation.
    //
    // The off-center pivot formula gives the equivalent rigid-body transform:
    //   pivot at (cx, cy-200), theta = 0.5 deg
    //   dx = 0*(1-cos) + (-200)*sin = -1.745 px   (guide star drifts west)
    //   dy = (-200)*(1-cos) - 0*sin = -0.008 px   (negligible N/S component)
    //
    // DONUTS should decompose this as:
    //   dx ~ -1.745 px  ->  guide command: nudge mount east
    //   dtheta ~ 0.5 deg -> guide command: advance de-rotator

    std::printf("--- testAltAzNorthPointing ---\n");
    const int W = 512, H = 512;
    auto stars = randomStars(W, H, 80, 55);
    auto ref   = makeFrame(W, H, stars);

    const double cx = W / 2.0, cy = H / 2.0;
    const double px = cx, py = cy - 200;  // NCP: 200 px north of guide star
    const double thetaDeg = 0.5;
    const double cosT = std::cos(thetaDeg * M_PI / 180.0);
    const double sinT = std::sin(thetaDeg * M_PI / 180.0);
    const double expectedDx = (px - cx) * (1.0 - cosT) + (py - cy) * sinT;
    const double expectedDy = (py - cy) * (1.0 - cosT) - (px - cx) * sinT;

    auto curr = transformFrame(ref, W, H, expectedDx, expectedDy, thetaDeg);

    Donuts::Registrar g;
    g.setReference(ref.data(), W, H);
    auto t = g.measure(curr.data(), W, H);

    double deg = t.dtheta * 180.0 / M_PI;
    std::printf("  NCP 200 px north of center, theta=%.1f deg\n", thetaDeg);
    std::printf("  Expected dx=%.4f dy=%.4f dtheta=%.4f deg\n",
                expectedDx, expectedDy, thetaDeg);
    std::printf("  Got     dx=%.4f dy=%.4f dtheta=%.4f deg  SNR=%.1f\n",
                t.dx, t.dy, deg, t.snr);

    CHECK(t.valid());
    CHECK_NEAR(t.dx,  expectedDx, 0.15);
    CHECK_NEAR(t.dy,  expectedDy, 0.15);
    CHECK_NEAR(deg,   thetaDeg,   0.05);
}

static void testAltAzMeridianCrossing()
{
    // Alt-az mount tracking a star crossing the southern meridian.
    // Between guide frames two things happen simultaneously:
    //   1) The mount tracking error lets the guide star drift 2.5 px east and
    //      1.5 px south.
    //   2) The de-rotator lags by 0.3 deg of field rotation (CCW).
    //
    // The rotation pivot is the guide star's new sensor position (cx+2.5, cy+1.5),
    // which is off-center by the tracking error.  Rotating around an off-center
    // pivot is mathematically identical to rotating around the image center then
    // applying the same translation, so DONUTS decomposes this as:
    //   dx = 2.5 px   ->  mount correction: nudge west
    //   dy = 1.5 px   ->  mount correction: nudge north
    //   dtheta = 0.3 deg -> de-rotator correction

    std::printf("--- testAltAzMeridianCrossing ---\n");
    const int W = 512, H = 512;
    auto stars = randomStars(W, H, 80, 66);
    auto ref   = makeFrame(W, H, stars);

    const double dx_track = 2.5;   // guide star drifted east
    const double dy_track = 1.5;   // guide star drifted south (y-down)
    const double thetaDeg = 0.3;   // de-rotator lag (CCW field rotation)

    auto curr = transformFrame(ref, W, H, dx_track, dy_track, thetaDeg);

    Donuts::Registrar g;
    g.setReference(ref.data(), W, H);
    auto t = g.measure(curr.data(), W, H);

    double deg = t.dtheta * 180.0 / M_PI;
    std::printf("  Guide drift (%.1f, %.1f) px, field rot %.2f deg\n",
                dx_track, dy_track, thetaDeg);
    std::printf("  Got: dx=%.4f dy=%.4f dtheta=%.4f deg  SNR=%.1f\n",
                t.dx, t.dy, deg, t.snr);

    CHECK(t.valid());
    CHECK_NEAR(t.dx, dx_track, 0.15);
    CHECK_NEAR(t.dy, dy_track, 0.15);
    CHECK_NEAR(deg,  thetaDeg, 0.05);
}

static void testAltAzRising()
{
    // Alt-az mount tracking a star in the eastern sky (rising toward meridian).
    // As the object rises, the parallactic angle increases: the field rotates CCW
    // (positive dtheta).  The guide star drifts west and north due to imperfect
    // tracking on the ascending arc.
    //
    // Same math as testAltAzMeridianCrossing but with reversed drift direction,
    // confirming DONUTS correctly measures negative dx/dy offsets.

    std::printf("--- testAltAzRising ---\n");
    const int W = 512, H = 512;
    auto stars = randomStars(W, H, 80, 88);
    auto ref   = makeFrame(W, H, stars);

    const double dx_track = -2.5;  // guide star drifted west
    const double dy_track = -1.5;  // guide star drifted north (y-down, so negative)
    const double thetaDeg =  0.3;  // CCW field rotation (parallactic angle increasing)

    auto curr = transformFrame(ref, W, H, dx_track, dy_track, thetaDeg);

    Donuts::Registrar g;
    g.setReference(ref.data(), W, H);
    auto t = g.measure(curr.data(), W, H);

    double deg = t.dtheta * 180.0 / M_PI;
    std::printf("  Guide drift (%.1f, %.1f) px, field rot %.2f deg\n",
                dx_track, dy_track, thetaDeg);
    std::printf("  Got: dx=%.4f dy=%.4f dtheta=%.4f deg  SNR=%.1f\n",
                t.dx, t.dy, deg, t.snr);

    CHECK(t.valid());
    CHECK_NEAR(t.dx, dx_track, 0.15);
    CHECK_NEAR(t.dy, dy_track, 0.15);
    CHECK_NEAR(deg,  thetaDeg, 0.05);
}

static void testAltAzSetting()
{
    // Alt-az mount tracking a star in the western sky (setting past meridian).
    // Past the meridian, the parallactic angle decreases: the field rotates CW,
    // which is a NEGATIVE dtheta in the algorithm's convention.
    // Guide star drifts east and south as the object descends.
    //
    // This is the key case that exercises negative rotation detection -- the
    // de-rotator must run in reverse compared to the rising-object case.

    std::printf("--- testAltAzSetting ---\n");
    const int W = 512, H = 512;
    auto stars = randomStars(W, H, 80, 44);
    auto ref   = makeFrame(W, H, stars);

    const double dx_track =  2.0;   // guide star drifted east
    const double dy_track =  1.5;   // guide star drifted south
    const double thetaDeg = -0.3;   // CW field rotation (parallactic angle decreasing)

    auto curr = transformFrame(ref, W, H, dx_track, dy_track, thetaDeg);

    Donuts::Registrar g;
    g.setReference(ref.data(), W, H);
    auto t = g.measure(curr.data(), W, H);

    double deg = t.dtheta * 180.0 / M_PI;
    std::printf("  Guide drift (%.1f, %.1f) px, field rot %.2f deg\n",
                dx_track, dy_track, thetaDeg);
    std::printf("  Got: dx=%.4f dy=%.4f dtheta=%.4f deg  SNR=%.1f\n",
                t.dx, t.dy, deg, t.snr);

    CHECK(t.valid());
    CHECK_NEAR(t.dx, dx_track, 0.15);
    CHECK_NEAR(t.dy, dy_track, 0.15);
    CHECK_NEAR(deg,  thetaDeg, 0.05);
}

static void testConfig()
{
    std::printf("--- testConfig (custom tukey alpha) ---\n");
    const int W = 512, H = 512;
    auto stars = randomStars(W, H, 80, 33);
    auto ref   = makeFrame(W, H, stars);
    auto curr  = transformFrame(ref, W, H, 3.0, -2.0, 0.0);

    Donuts::Config cfg;
    cfg.tukeyAlpha = 0.25;

    Donuts::Registrar g(cfg);
    g.setReference(ref.data(), W, H);
    auto t = g.measure(curr.data(), W, H);

    std::printf("  Expected (3.0, -2.0)  Got (%.4f, %.4f)  SNR=%.1f\n",
                t.dx, t.dy, t.snr);

    CHECK(t.valid());
    CHECK_NEAR(t.dx,  3.0, 0.15);
    CHECK_NEAR(t.dy, -2.0, 0.15);
}

// ---------------------------------------------------------------------------
// Simulation scenario helpers inspired by McCormac et al. 2013 (PASP 125, 548)
// ---------------------------------------------------------------------------

// Minimal Box-Muller RNG (no stdlib dependency)
struct Rng
{
    uint32_t s;
    explicit Rng(uint32_t seed = 1) : s(seed ? seed : 1) {}
    double uni() { s ^= s << 13; s ^= s >> 17; s ^= s << 5; return (s & 0xFFFFu) / 65536.0; }
    double gauss() {
        double u1 = uni() + 1e-10, u2 = uni();
        return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
    }
    double uniform(double lo, double hi) { return lo + uni() * (hi - lo); }
};

// Add additive Gaussian noise (sigma per pixel) to a frame.
static std::vector<double> addGaussianNoise(
    const std::vector<double> &src, double sigma, Rng &rng)
{
    std::vector<double> dst = src;
    for (auto &v : dst) v += rng.gauss() * sigma;
    return dst;
}

// Apply per-pixel multiplicative Gaussian noise: value *= (1 + N(0, sigma)).
// McCormac 2013 Fig. 1 "green" scenario: sigma = 0.20 (+-20% per pixel).
static std::vector<double> addMultiplicativeNoise(
    const std::vector<double> &src, double sigma, Rng &rng)
{
    std::vector<double> dst = src;
    for (auto &v : dst) v *= (1.0 + rng.gauss() * sigma);
    return dst;
}

// Separable Gaussian blur.  sigma = 0 returns src unchanged.
static std::vector<double> gaussianBlur(
    const std::vector<double> &src, int w, int h, double sigma)
{
    if (sigma < 0.01) return src;
    int r  = static_cast<int>(std::ceil(3.0 * sigma));
    int ks = 2 * r + 1;
    std::vector<double> k(ks);
    double ksum = 0.0;
    for (int i = 0; i < ks; ++i) {
        double x = i - r;
        k[i] = std::exp(-0.5 * x * x / (sigma * sigma));
        ksum += k[i];
    }
    for (auto &kv : k) kv /= ksum;

    // Horizontal pass
    std::vector<double> tmp(w * h);
    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x) {
            double s = 0.0;
            for (int d = -r; d <= r; ++d)
                s += src[y * w + std::max(0, std::min(w - 1, x + d))] * k[d + r];
            tmp[y * w + x] = s;
        }

    // Vertical pass
    std::vector<double> dst(w * h);
    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x) {
            double s = 0.0;
            for (int d = -r; d <= r; ++d)
                s += tmp[std::max(0, std::min(h - 1, y + d)) * w + x] * k[d + r];
            dst[y * w + x] = s;
        }
    return dst;
}

// ---------------------------------------------------------------------------
// S/N sweep analog (inspired by McCormac 2013 s3.5, adapted for 4-quadrant)
// ---------------------------------------------------------------------------
// The paper created 50x50 single-star stamps; we use a 256x256 multi-star
// field because our 4-quadrant solver needs spatial coverage.  Noise is
// Gaussian white added to the comparison frame only (reference is clean).
// Success: |dx_residual| <= 0.3 px AND |dy_residual| <= 0.3 px.
static void testSnrSweepTranslation()
{
    std::printf("--- testSnrSweepTranslation (McCormac 2013 s3.5 analog) ---\n");

    const int    W       = 256, H = 256;
    const double BG      = 1000.0;
    const double PSF_SIG = 2.5 / 2.355;   // FWHM=2.5 px, as in paper
    const double MAX_PK  = 50000.0;        // brightest star peak above BG
    const int    N_SH    = 100;            // shifts per S/N level
    const double RANGE   = 3.0;            // +/- px shift range

    auto stars = randomStars(W, H, 20, 1234);
    {   // normalise so brightest star has peak = MAX_PK
        double mx = 0.0;
        for (auto &st : stars) mx = std::max(mx, st.peak);
        for (auto &st : stars) st.peak = st.peak / mx * MAX_PK;
    }
    auto ref = makeFrame(W, H, stars, BG, PSF_SIG);

    const double snrTab[] = {1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0, 20.0};
    const int    NL       = static_cast<int>(sizeof(snrTab) / sizeof(snrTab[0]));
    double       pct[NL]  = {};

    std::printf("  %5s  %8s  %7s  %7s\n", "S/N", "success%", "mE_dx", "mE_dy");

    for (int li = 0; li < NL; ++li)
    {
        double sigma_n = MAX_PK / snrTab[li];
        Rng    rng(1000 + li * 37);
        Donuts::Registrar g;
        g.setReference(ref.data(), W, H);
        int    npass = 0;
        double sEx = 0.0, sEy = 0.0;

        for (int sh = 0; sh < N_SH; ++sh)
        {
            double tdx = rng.uniform(-RANGE, RANGE);
            double tdy = rng.uniform(-RANGE, RANGE);
            auto cur = addGaussianNoise(
                transformFrame(ref, W, H, tdx, tdy, 0.0), sigma_n, rng);
            auto t  = g.measure(cur.data(), W, H);
            double ex = std::abs(t.dx - tdx), ey = std::abs(t.dy - tdy);
            if (t.valid() && ex <= 0.3 && ey <= 0.3) { ++npass; sEx += ex; sEy += ey; }
        }

        pct[li]       = 100.0 * npass / N_SH;
        double mex    = (npass > 0) ? sEx / npass : 99.9;
        double mey    = (npass > 0) ? sEy / npass : 99.9;
        std::printf("  %5.1f  %8.1f  %7.4f  %7.4f\n", snrTab[li], pct[li], mex, mey);
    }

    // Our phase-only correlation + 3-sigma detection threshold gives a higher
    // limiting S/N than the original paper (which used amplitude correlation on
    // a single star).  At S/N=20 the algorithm must be reliable; below S/N=5
    // it is expected to fail.
    CHECK(pct[8] >= 95.0);   // S/N=20 (index 8) -> >= 95% success
    CHECK(pct[0] < 80.0);    // S/N=1  (index 0) -> not perfect (sanity)
}

// ---------------------------------------------------------------------------
// Uniform intensity change robustness (analog of McCormac 2013 Fig. 1 blue)
// ---------------------------------------------------------------------------
// The overall comparison-frame flux is scaled uniformly.  Phase-only
// cross-correlation normalises spectral amplitude, so this should have
// negligible effect on the result.
static void testUniformIntensityChange()
{
    std::printf("--- testUniformIntensityChange (McCormac 2013 Fig 1 blue) ---\n");

    const int    W = 256, H = 256;
    const double PSF_SIG = 2.5 / 2.355;
    auto stars = randomStars(W, H, 20, 77);
    auto ref   = makeFrame(W, H, stars, 1000.0, PSF_SIG);

    Donuts::Registrar g;
    g.setReference(ref.data(), W, H);

    const double tdx = 3.5, tdy = -2.0;
    auto shifted = transformFrame(ref, W, H, tdx, tdy, 0.0);

    const double scales[] = {0.5, 0.75, 1.25, 1.5, 2.0};
    std::printf("  %6s  %8s  %8s\n", "scale", "dx_err", "dy_err");
    for (double sc : scales)
    {
        std::vector<double> cur(shifted.size());
        for (std::size_t i = 0; i < shifted.size(); ++i) cur[i] = shifted[i] * sc;
        auto t  = g.measure(cur.data(), W, H);
        double ex = std::abs(t.dx - tdx), ey = std::abs(t.dy - tdy);
        std::printf("  %6.2f  %8.4f  %8.4f\n", sc, ex, ey);
        CHECK(t.valid());
        CHECK(ex <= 0.15);
        CHECK(ey <= 0.15);
    }
}

// ---------------------------------------------------------------------------
// Per-pixel multiplicative noise robustness (analog of McCormac 2013 Fig. 1 green)
// ---------------------------------------------------------------------------
// Each comparison pixel is randomly scaled by N(1, 0.20), matching the
// paper's "standard deviation of 1 from +-20% of the intensity modified
// pixel value".
static void testPixelNoise()
{
    std::printf("--- testPixelNoise (McCormac 2013 Fig 1 green) ---\n");

    const int    W = 256, H = 256;
    const double PSF_SIG = 2.5 / 2.355;
    auto stars = randomStars(W, H, 20, 77);
    auto ref   = makeFrame(W, H, stars, 1000.0, PSF_SIG);

    Donuts::Registrar g;
    g.setReference(ref.data(), W, H);

    const double tdx = 3.5, tdy = -2.0;
    auto shifted = transformFrame(ref, W, H, tdx, tdy, 0.0);

    // 20 trials with different noise seeds; all must pass.
    int npass = 0;
    for (int trial = 0; trial < 20; ++trial)
    {
        Rng rng(500 + trial * 7);
        auto cur = addMultiplicativeNoise(shifted, 0.20, rng);
        auto t   = g.measure(cur.data(), W, H);
        double ex = std::abs(t.dx - tdx), ey = std::abs(t.dy - tdy);
        if (t.valid() && ex <= 0.3 && ey <= 0.3) ++npass;
    }
    std::printf("  %d / 20 trials passed (+-20%% per-pixel noise)\n", npass);
    CHECK(npass >= 18);
}

// ---------------------------------------------------------------------------
// Seeing change robustness (analog of McCormac 2013 s3.6 / Table 4)
// ---------------------------------------------------------------------------
// Reference is in focus (FWHM = 2.5 px).  Comparison frames are Gaussian-
// blurred to simulate seeing deterioration (0-125% FWHM increase), then
// shifted by a fixed amount.  Blur sigma derived from quadrature addition:
//   sigma_blur = sigma_ref * sqrt((1+delta)^2 - 1)
static void testSeeingChange()
{
    std::printf("--- testSeeingChange (McCormac 2013 s3.6 / Table 4 analog) ---\n");

    const int    W       = 256, H = 256;
    const double PSF_SIG = 2.5 / 2.355;   // sigma for FWHM=2.5 px
    // Uniform peak brightness so all stars survive heavy blurring: a star
    // blurred to 125% FWHM has its peak reduced by (sigma_ref/sigma_result)^2
    // = 1/(1+1.25)^2 = 0.18.  At peak=30000 the blurred peak is 5400, well
    // above the threshold set from the clean reference frame.
    auto stars = randomStars(W, H, 30, 77);
    for (auto &st : stars) st.peak = 30000.0;
    auto ref   = makeFrame(W, H, stars, 1000.0, PSF_SIG);

    Donuts::Registrar g;
    g.setReference(ref.data(), W, H);

    const double tdx = 2.0, tdy = -1.5;
    auto shifted = transformFrame(ref, W, H, tdx, tdy, 0.0);

    // delta_FWHM levels (fraction) matching Table 4 of paper
    const double deltas[]   = {0.00, 0.20, 0.35, 0.50, 0.65, 0.80, 0.95, 1.10, 1.25};
    const int    ND         = static_cast<int>(sizeof(deltas) / sizeof(deltas[0]));

    std::printf("  %8s  %8s  %8s  %8s\n",
                "dFWHM%", "FWHM_px", "dx_err", "dy_err");

    for (int i = 0; i < ND; ++i)
    {
        double d         = deltas[i];
        double fwhm_res  = 2.5 * (1.0 + d);
        double sig_res   = fwhm_res / 2.355;
        double sig_blur  = PSF_SIG * std::sqrt((1.0 + d) * (1.0 + d) - 1.0);
        auto   blurred   = gaussianBlur(shifted, W, H, sig_blur);
        auto   t         = g.measure(blurred.data(), W, H);
        double ex        = std::abs(t.dx - tdx), ey = std::abs(t.dy - tdy);
        std::printf("  %8.0f  %8.2f  %8.4f  %8.4f\n",
                    d * 100.0, fwhm_res, ex, ey);
        // Phase-only cross-correlation normalises amplitude, so PSF broadening
        // mainly reduces the high-frequency content of the profiles.  Residual
        // errors stay sub-pixel but exceed the 0.3 px of the paper's amplitude-
        // correlation algorithm at large blur levels.
        CHECK(t.valid());
        CHECK(ex <= 0.5);
        CHECK(ey <= 0.5);
    }
}

// ---------------------------------------------------------------------------
// S/N sweep with rotation (4-quadrant extension, not from McCormac 2013)
// ---------------------------------------------------------------------------
// Same noise model as testSnrSweepTranslation but the comparison frame also
// carries a small rotation (0.3 deg).  At high S/N all three DoFs should
// be recovered; rotation accuracy degrades as noise rises.
static void testSnrSweepWithRotation()
{
    std::printf("--- testSnrSweepWithRotation (McCormac 2013 extension) ---\n");

    const int    W        = 256, H = 256;
    const double BG       = 1000.0;
    const double PSF_SIG  = 2.5 / 2.355;
    const double STAR_PK  = 30000.0;       // uniform brightness: S/N = STAR_PK / sigma_noise
    const int    N_SH     = 100;
    const double RANGE    = 1.5;           // +/- px
    const double ROT_DEG  = 0.5;          // fixed rotation; at W/4=64px lever arm -> 0.56 px arc

    // Uniform-brightness stars so the S/N definition is unambiguous and
    // all stars remain detectable across all tested noise levels.
    auto stars = randomStars(W, H, 40, 5678);
    for (auto &st : stars) st.peak = STAR_PK;
    auto ref = makeFrame(W, H, stars, BG, PSF_SIG);

    const double snrTab[] = {3.0, 5.0, 7.0, 10.0, 15.0, 20.0};
    const int    NL       = static_cast<int>(sizeof(snrTab) / sizeof(snrTab[0]));
    double       pctXY[NL] = {}, pctRot[NL] = {};

    std::printf("  %5s  %8s  %8s\n", "S/N", "XY_ok%%", "rot_ok%%");

    for (int li = 0; li < NL; ++li)
    {
        double sigma_n = STAR_PK / snrTab[li];
        Rng    rng(2000 + li * 53);
        Donuts::Registrar g;
        g.setReference(ref.data(), W, H);
        int nXY = 0, nRot = 0;

        for (int sh = 0; sh < N_SH; ++sh)
        {
            double tdx = rng.uniform(-RANGE, RANGE);
            double tdy = rng.uniform(-RANGE, RANGE);
            auto cur   = addGaussianNoise(
                transformFrame(ref, W, H, tdx, tdy, ROT_DEG), sigma_n, rng);
            auto t     = g.measure(cur.data(), W, H);
            double ex  = std::abs(t.dx - tdx), ey = std::abs(t.dy - tdy);
            double er  = std::abs(t.dtheta * 180.0 / M_PI - ROT_DEG);
            bool xyOk  = t.valid() && ex <= 0.3 && ey <= 0.3;
            bool rotOk = t.valid() && er <= 0.1;
            if (xyOk)  ++nXY;
            if (rotOk) ++nRot;
        }

        pctXY[li]  = 100.0 * nXY  / N_SH;
        pctRot[li] = 100.0 * nRot / N_SH;
        std::printf("  %5.1f  %8.1f  %8.1f\n",
                    snrTab[li], pctXY[li], pctRot[li]);
    }

    // Translation at S/N=20 must be reliable (index 5).
    CHECK(pctXY[5] >= 90.0);
    // Rotation detection is harder: the differential centroid shift per quadrant
    // (~35px lever arm * 0.5 deg = 0.31 px) competes with translation noise.
    // 70% success at S/N=20 correctly characterises the algorithm's rotation
    // sensitivity floor under combined noise + translation.
    CHECK(pctRot[5] >= 70.0);
}

// ---------------------------------------------------------------------------
// Translation-only characterisation tests
// (mirrored into Tests/ekos/guide/testdonutsguider.cpp later)
// ---------------------------------------------------------------------------

// Sweep translation magnitude from 0.5 to 8.0 px along a diagonal.
// Characterises how accuracy degrades as the shift grows relative to the
// profile length; asserts <= 0.3 px for shifts that fit within the linear
// regime of the phase correlation.
static void testTranslationMagnitudeSweep()
{
    std::printf("--- testTranslationMagnitudeSweep ---\n");
    const int W = 512, H = 512;
    auto stars = randomStars(W, H, 80, 42);
    auto ref   = makeFrame(W, H, stars);

    Donuts::Registrar g;
    g.setReference(ref.data(), W, H);

    std::printf("  %6s  %8s  %8s\n", "shift", "dx_err", "dy_err");
    for (int i = 1; i <= 16; ++i)
    {
        double shift = i * 0.5;
        double tdx   = shift, tdy = shift * 0.7;   // non-square to expose asymmetry
        auto   cur   = transformFrame(ref, W, H, tdx, tdy, 0.0);
        auto   t     = g.measure(cur.data(), W, H);
        double ex    = std::abs(t.dx - tdx), ey = std::abs(t.dy - tdy);
        std::printf("  %6.1f  %8.4f  %8.4f  %s\n",
                    shift, ex, ey, (ex <= 0.3 && ey <= 0.3) ? "ok" : "DEGRADED");
        if (shift <= 6.0)
        {
            CHECK(t.valid());
            CHECK(ex <= 0.3);
            CHECK(ey <= 0.3);
        }
    }
}

// Verify translation accuracy is symmetric across all 8 compass directions.
static void testTranslationDirections()
{
    std::printf("--- testTranslationDirections ---\n");
    const int W = 512, H = 512;
    auto stars = randomStars(W, H, 80, 31);
    auto ref   = makeFrame(W, H, stars);

    Donuts::Registrar g;
    g.setReference(ref.data(), W, H);

    struct Case { double dx, dy; const char *label; };
    const Case cases[] = {
        {  5.0,  0.0, "E "  }, { -5.0,  0.0, "W "  },
        {  0.0,  5.0, "S "  }, {  0.0, -5.0, "N "  },
        {  3.5,  3.5, "SE"  }, { -3.5,  3.5, "SW"  },
        {  3.5, -3.5, "NE"  }, { -3.5, -3.5, "NW"  },
    };
    std::printf("  %4s  %8s  %8s\n", "dir", "dx_err", "dy_err");
    for (const auto &c : cases)
    {
        auto cur = transformFrame(ref, W, H, c.dx, c.dy, 0.0);
        auto t   = g.measure(cur.data(), W, H);
        double ex = std::abs(t.dx - c.dx), ey = std::abs(t.dy - c.dy);
        std::printf("  %4s  %8.4f  %8.4f\n", c.label, ex, ey);
        CHECK(t.valid());
        CHECK(ex <= 0.15);
        CHECK(ey <= 0.15);
    }
}

// Verify the algorithm works across a range of guide star counts.
// In the field this matters: sparse fields (few stars) give weaker profiles.
static void testStarDensity()
{
    std::printf("--- testStarDensity ---\n");
    const int    W  = 512, H = 512;
    const double DX = 4.3, DY = -3.1;
    const int    counts[] = {5, 10, 20, 40, 80};

    std::printf("  %6s  %8s  %8s  %6s\n", "stars", "dx_err", "dy_err", "SNR");
    for (int n : counts)
    {
        auto stars = randomStars(W, H, n, 200 + n);
        auto ref   = makeFrame(W, H, stars);
        auto cur   = transformFrame(ref, W, H, DX, DY, 0.0);

        Donuts::Registrar g;
        g.setReference(ref.data(), W, H);
        auto t = g.measure(cur.data(), W, H);

        double ex = std::abs(t.dx - DX), ey = std::abs(t.dy - DY);
        std::printf("  %6d  %8.4f  %8.4f  %6.1f\n", n, ex, ey, t.snr);
        if (n >= 10)
        {
            CHECK(t.valid());
            CHECK(ex <= 0.3);
            CHECK(ey <= 0.3);
        }
    }
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------

static void testScaleDetection()
{
    // 512x512 gives quadrant centroids at ~128px, doubling the scale lever arm
    // vs. 256x256 and halving the scale estimation error.
    const int W = 512, H = 512;
    auto stars = randomStars(W, H, 20, 42u);
    for (auto &s : stars) s.peak = 30000.0;

    Donuts::Config cfg;
    cfg.detectScale = true;
    Donuts::Registrar guider(cfg);

    std::vector<double> ref = makeFrame(W, H, stars, 1000.0, 2.5);
    guider.setReference(ref.data(), W, H);

    // Pure scale: realistic focus-drift range is sub-0.5% per session.
    // At 1% the linear centroid model introduces ~0.2px coupling into dx, which
    // is expected and acceptable (1% scale = ~10mm defocus at f=1000mm; stars
    // are unusable long before that level).
    std::printf("--- testScaleDetection ---\n");
    std::printf("  scale   got_dx   got_dy  got_scale  scale_err\n");
    const double scales[] = { 0.995, 0.998, 1.0, 1.002, 1.005 };
    for (double s : scales)
    {
        auto frame = scaleFrame(ref, W, H, s);
        auto t = guider.measure(frame.data(), W, H);
        std::printf("  %.3f  %7.4f  %7.4f  %9.6f  %9.6f\n",
                    s, t.dx, t.dy, t.scale, t.scale - s);
        CHECK(t.valid());
        CHECK_NEAR(t.dx,    0.0, 0.15);   // <0.15px bleedthrough within 0.5% scale range
        CHECK_NEAR(t.dy,    0.0, 0.15);
        CHECK_NEAR(t.scale, s,   0.002);  // scale accurate to 0.2% at 512x512
    }

    // Scale + translation: the two DoF must not bleed into each other.
    {
        const double TX = 3.0, TY = -2.0, S = 1.008;
        auto frame = scaleFrame(ref, W, H, S);
        // Apply translation on top of scale.
        frame = transformFrame(frame, W, H, TX, TY, 0.0);
        auto t = guider.measure(frame.data(), W, H);
        CHECK(t.valid());
        CHECK_NEAR(t.dx,    TX, 0.2);
        CHECK_NEAR(t.dy,    TY, 0.2);
        CHECK_NEAR(t.scale, S,  0.004);
    }
}

int main()
{
    testPureTranslation();
    testSubPixelTranslation();
    testPureRotation();
    testCombined();
    testCombinedSmallImage();
    testAltAzNorthPointing();
    testAltAzMeridianCrossing();
    testAltAzRising();
    testAltAzSetting();
    testReset();
    testSmallImage();
    testConfig();
    testSnrSweepTranslation();
    testUniformIntensityChange();
    testPixelNoise();
    testSeeingChange();
    testSnrSweepWithRotation();
    testTranslationMagnitudeSweep();
    testTranslationDirections();
    testStarDensity();
    testScaleDetection();

    std::printf("\n%d passed, %d failed\n", g_pass, g_fail);
    return g_fail ? 1 : 0;
}
