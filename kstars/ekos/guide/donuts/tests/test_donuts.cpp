/*
    SPDX-FileCopyrightText: 2024 Christian Kemper <ckemper@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later

    Standalone test for the Donuts guiding algorithm.
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

    Donuts::Guider g;
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

    Donuts::Guider g;
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

    Donuts::Guider g;
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

    Donuts::Guider g;
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

    Donuts::Guider g;
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

    Donuts::Guider g;
    g.setReference(ref.data(), W, H);
    auto t = g.measure(curr.data(), W, H);

    std::printf("  Expected (1.5, -1.0)  Got (%.4f, %.4f)  SNR=%.1f\n",
                t.dx, t.dy, t.snr);

    // Smaller image -> less SNR, wider tolerance
    CHECK(t.valid());
    CHECK_NEAR(t.dx,  1.5, 0.3);
    CHECK_NEAR(t.dy, -1.0, 0.3);
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

    Donuts::Guider g;
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

    Donuts::Guider g;
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

    Donuts::Guider g;
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

    Donuts::Guider g;
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

    Donuts::Guider g(cfg);
    g.setReference(ref.data(), W, H);
    auto t = g.measure(curr.data(), W, H);

    std::printf("  Expected (3.0, -2.0)  Got (%.4f, %.4f)  SNR=%.1f\n",
                t.dx, t.dy, t.snr);

    CHECK(t.valid());
    CHECK_NEAR(t.dx,  3.0, 0.15);
    CHECK_NEAR(t.dy, -2.0, 0.15);
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------

int main()
{
    testPureTranslation();
    testSubPixelTranslation();
    testPureRotation();
    testCombined();
    testAltAzNorthPointing();
    testAltAzMeridianCrossing();
    testAltAzRising();
    testAltAzSetting();
    testReset();
    testSmallImage();
    testConfig();

    std::printf("\n%d passed, %d failed\n", g_pass, g_fail);
    return g_fail ? 1 : 0;
}
