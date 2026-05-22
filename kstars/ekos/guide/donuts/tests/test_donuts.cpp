/*
    SPDX-FileCopyrightText: 2026 Christian Kemper <ckemper@gmail.com>

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

struct Star { double x, y, peak; };

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

// Apply a rigid-body transform (translate then rotate around centre) using
// bilinear interpolation.  dtheta_deg = 0 for pure translation.
static std::vector<double> transformFrame(
    const std::vector<double> &src, int w, int h,
    double dx, double dy, double dtheta_deg = 0.0)
{
    std::vector<double> dst(w * h, 0.0);
    const double cosT = std::cos(dtheta_deg * M_PI / 180.0);
    const double sinT = std::sin(dtheta_deg * M_PI / 180.0);
    const double cx = w / 2.0, cy = h / 2.0;

    for (int y = 0; y < h; ++y)
    {
        for (int x = 0; x < w; ++x)
        {
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

static std::vector<Star> randomStars(int w, int h, int count, uint32_t seed)
{
    std::vector<Star> stars;
    stars.reserve(count);
    uint32_t s = seed;
    auto rng = [&]() -> double {
        s ^= s << 13; s ^= s >> 17; s ^= s << 5;
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
// Noise / perturbation helpers (for McCormac 2013 scenario tests)
// ---------------------------------------------------------------------------

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

static std::vector<double> addGaussianNoise(
    const std::vector<double> &src, double sigma, Rng &rng)
{
    std::vector<double> dst = src;
    for (auto &v : dst) v += rng.gauss() * sigma;
    return dst;
}

static std::vector<double> addMultiplicativeNoise(
    const std::vector<double> &src, double sigma, Rng &rng)
{
    std::vector<double> dst = src;
    for (auto &v : dst) v *= (1.0 + rng.gauss() * sigma);
    return dst;
}

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

    std::vector<double> tmp(w * h);
    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x) {
            double s = 0.0;
            for (int d = -r; d <= r; ++d)
                s += src[y * w + std::max(0, std::min(w - 1, x + d))] * k[d + r];
            tmp[y * w + x] = s;
        }

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
// Tests
// ---------------------------------------------------------------------------

static void testPureTranslation()
{
    std::printf("--- testPureTranslation ---\n");
    const int W = 512, H = 512;
    auto stars = randomStars(W, H, 80, 42);
    auto ref   = makeFrame(W, H, stars);
    auto curr  = transformFrame(ref, W, H, 5.4, -3.2);

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
    auto curr  = transformFrame(ref, W, H, 0.7, -0.4);

    Donuts::Registrar g;
    g.setReference(ref.data(), W, H);
    auto t = g.measure(curr.data(), W, H);

    std::printf("  Expected (0.7, -0.4)  Got (%.4f, %.4f)  SNR=%.1f\n",
                t.dx, t.dy, t.snr);

    CHECK(t.valid());
    CHECK_NEAR(t.dx,  0.7, 0.15);
    CHECK_NEAR(t.dy, -0.4, 0.15);
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

    auto t = g.measure(frame.data(), W, H);
    CHECK(!t.valid());
}

static void testSmallImage()
{
    std::printf("--- testSmallImage (128x128) ---\n");
    const int W = 128, H = 128;
    auto stars = randomStars(W, H, 20, 17);
    auto ref   = makeFrame(W, H, stars);
    auto curr  = transformFrame(ref, W, H, 1.5, -1.0);

    Donuts::Registrar g;
    g.setReference(ref.data(), W, H);
    auto t = g.measure(curr.data(), W, H);

    std::printf("  Expected (1.5, -1.0)  Got (%.4f, %.4f)  SNR=%.1f\n",
                t.dx, t.dy, t.snr);

    CHECK(t.valid());
    CHECK_NEAR(t.dx,  1.5, 0.3);
    CHECK_NEAR(t.dy, -1.0, 0.3);
}

static void testConfig()
{
    std::printf("--- testConfig (custom tukeyAlpha) ---\n");
    const int W = 512, H = 512;
    auto stars = randomStars(W, H, 80, 33);
    auto ref   = makeFrame(W, H, stars);
    auto curr  = transformFrame(ref, W, H, 3.0, -2.0);

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

static void testSnrSweepTranslation()
{
    std::printf("--- testSnrSweepTranslation (McCormac 2013 s3.5 analog) ---\n");

    const int    W       = 256, H = 256;
    const double BG      = 1000.0;
    const double PSF_SIG = 2.5 / 2.355;
    const double MAX_PK  = 50000.0;
    const int    N_SH    = 100;
    const double RANGE   = 3.0;

    auto stars = randomStars(W, H, 20, 1234);
    {
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
                transformFrame(ref, W, H, tdx, tdy), sigma_n, rng);
            auto t  = g.measure(cur.data(), W, H);
            double ex = std::abs(t.dx - tdx), ey = std::abs(t.dy - tdy);
            if (t.valid() && ex <= 0.3 && ey <= 0.3) { ++npass; sEx += ex; sEy += ey; }
        }

        pct[li]    = 100.0 * npass / N_SH;
        double mex = (npass > 0) ? sEx / npass : 99.9;
        double mey = (npass > 0) ? sEy / npass : 99.9;
        std::printf("  %5.1f  %8.1f  %7.4f  %7.4f\n", snrTab[li], pct[li], mex, mey);
    }

    CHECK(pct[8] >= 95.0);   // S/N=20 -> >= 95% success
    CHECK(pct[0] < 80.0);    // S/N=1  -> not perfect (sanity)
}

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
    auto shifted = transformFrame(ref, W, H, tdx, tdy);

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
    auto shifted = transformFrame(ref, W, H, tdx, tdy);

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

static void testSeeingChange()
{
    std::printf("--- testSeeingChange (McCormac 2013 s3.6 / Table 4 analog) ---\n");

    const int    W       = 256, H = 256;
    const double PSF_SIG = 2.5 / 2.355;
    auto stars = randomStars(W, H, 30, 77);
    for (auto &st : stars) st.peak = 30000.0;
    auto ref   = makeFrame(W, H, stars, 1000.0, PSF_SIG);

    Donuts::Registrar g;
    g.setReference(ref.data(), W, H);

    const double tdx = 2.0, tdy = -1.5;
    auto shifted = transformFrame(ref, W, H, tdx, tdy);

    const double deltas[] = {0.00, 0.20, 0.35, 0.50, 0.65, 0.80, 0.95, 1.10, 1.25};
    const int    ND       = static_cast<int>(sizeof(deltas) / sizeof(deltas[0]));

    std::printf("  %8s  %8s  %8s  %8s\n", "dFWHM%", "FWHM_px", "dx_err", "dy_err");

    for (int i = 0; i < ND; ++i)
    {
        double d        = deltas[i];
        double fwhm_res = 2.5 * (1.0 + d);
        double sig_blur = PSF_SIG * std::sqrt((1.0 + d) * (1.0 + d) - 1.0);
        auto   blurred  = gaussianBlur(shifted, W, H, sig_blur);
        auto   t        = g.measure(blurred.data(), W, H);
        double ex       = std::abs(t.dx - tdx), ey = std::abs(t.dy - tdy);
        std::printf("  %8.0f  %8.2f  %8.4f  %8.4f\n",
                    d * 100.0, fwhm_res, ex, ey);
        CHECK(t.valid());
        CHECK(ex <= 0.5);
        CHECK(ey <= 0.5);
    }
}

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
        double tdx   = shift, tdy = shift * 0.7;
        auto   cur   = transformFrame(ref, W, H, tdx, tdy);
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
        auto cur = transformFrame(ref, W, H, c.dx, c.dy);
        auto t   = g.measure(cur.data(), W, H);
        double ex = std::abs(t.dx - c.dx), ey = std::abs(t.dy - c.dy);
        std::printf("  %4s  %8.4f  %8.4f\n", c.label, ex, ey);
        CHECK(t.valid());
        CHECK(ex <= 0.15);
        CHECK(ey <= 0.15);
    }
}

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
        auto cur   = transformFrame(ref, W, H, DX, DY);

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

int main()
{
    testPureTranslation();
    testSubPixelTranslation();
    testReset();
    testSmallImage();
    testConfig();
    testSnrSweepTranslation();
    testUniformIntensityChange();
    testPixelNoise();
    testSeeingChange();
    testTranslationMagnitudeSweep();
    testTranslationDirections();
    testStarDensity();

    std::printf("\n%d passed, %d failed\n", g_pass, g_fail);
    return g_fail ? 1 : 0;
}
