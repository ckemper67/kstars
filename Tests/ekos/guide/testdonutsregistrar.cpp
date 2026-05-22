/*
    SPDX-FileCopyrightText: 2026 Christian Kemper <ckemper@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

#include "testdonutsregistrar.h"
#include "ekos/guide/donuts/donuts.h"
#include <QTest>
#include <cmath>
#include <cstdint>

TestDonutsRegistrar::TestDonutsRegistrar() : QObject()
{
}

void TestDonutsRegistrar::initTestCase() {}
void TestDonutsRegistrar::cleanupTestCase() {}

void TestDonutsRegistrar::testMathRecovery()
{
    const int W = 256, H = 256;
    Buffer master = makeStarField(W, H);

    Donuts::Registrar guider;
    guider.setReference(master.data(), W, H);

    // Pure translation
    {
        double tx = 5.4, ty = -3.2;
        Buffer shifted = transformBuffer(master, W, H, tx, ty);
        Donuts::Transform res = guider.measure(shifted.data(), W, H);
        qDebug() << "Translation: Expected (" << tx << "," << ty
                 << ") Got (" << res.dx << "," << res.dy << ")";
        QVERIFY(std::abs(res.dx - tx) < 0.15);
        QVERIFY(std::abs(res.dy - ty) < 0.15);
    }

    // Sub-pixel translation
    {
        double tx = 0.7, ty = -0.4;
        Buffer shifted = transformBuffer(master, W, H, tx, ty);
        Donuts::Transform res = guider.measure(shifted.data(), W, H);
        qDebug() << "Sub-pixel: Expected (" << tx << "," << ty
                 << ") Got (" << res.dx << "," << res.dy << ")";
        QVERIFY(std::abs(res.dx - tx) < 0.2);
        QVERIFY(std::abs(res.dy - ty) < 0.2);
    }
}

void TestDonutsRegistrar::testTranslationMagnitudeSweep()
{
    const int W = 256, H = 256;
    Buffer ref = makeStarField(W, H);

    Donuts::Registrar g;
    g.setReference(ref.data(), W, H);

    for (int i = 1; i <= 12; ++i)
    {
        double shift = i * 0.5;
        double tdx   = shift, tdy = shift * 0.7;
        Buffer cur   = transformBuffer(ref, W, H, tdx, tdy);
        Donuts::Transform t = g.measure(cur.data(), W, H);
        QVERIFY(t.valid());
        QVERIFY(std::abs(t.dx - tdx) < 0.3);
        QVERIFY(std::abs(t.dy - tdy) < 0.3);
    }
}

void TestDonutsRegistrar::testTranslationDirections()
{
    const int W = 256, H = 256;
    Buffer ref = makeStarField(W, H, 31);

    Donuts::Registrar g;
    g.setReference(ref.data(), W, H);

    struct Case { double dx, dy; };
    const Case cases[] = {
        { 5.0, 0.0}, {-5.0, 0.0},
        { 0.0, 5.0}, { 0.0,-5.0},
        { 3.5, 3.5}, {-3.5, 3.5},
        { 3.5,-3.5}, {-3.5,-3.5},
    };
    for (const auto &c : cases)
    {
        Buffer cur = transformBuffer(ref, W, H, c.dx, c.dy);
        Donuts::Transform t = g.measure(cur.data(), W, H);
        QVERIFY(t.valid());
        QVERIFY(std::abs(t.dx - c.dx) < 0.15);
        QVERIFY(std::abs(t.dy - c.dy) < 0.15);
    }
}

void TestDonutsRegistrar::testStarDensity()
{
    const int    W = 256, H = 256;
    const double DX = 4.3, DY = -3.1;
    const int    counts[] = {10, 20, 40, 80};

    for (int n : counts)
    {
        // Simple deterministic star field using makeStarField as a base
        // (makeRandomField inlined for brevity)
        std::vector<double> buf(W * H, 100.0);
        uint32_t s = 200u + static_cast<uint32_t>(n);
        auto rng = [&]() -> double {
            s ^= s << 13; s ^= s >> 17; s ^= s << 5;
            return (s & 0xFFFFu) / 65536.0;
        };
        const double margin = 20.0, psf_sig = 2.5;
        const double inv2s2 = 1.0 / (2.0 * psf_sig * psf_sig);
        const int    radius = static_cast<int>(4 * psf_sig) + 1;
        for (int i = 0; i < n; ++i)
        {
            double sx   = margin + rng() * (W - 2 * margin);
            double sy   = margin + rng() * (H - 2 * margin);
            double peak = 5000.0 + rng() * 45000.0;
            int x0 = static_cast<int>(sx), y0 = static_cast<int>(sy);
            for (int dy = -radius; dy <= radius; ++dy)
                for (int dx = -radius; dx <= radius; ++dx)
                {
                    int xi = x0 + dx, yi = y0 + dy;
                    if (xi < 0 || xi >= W || yi < 0 || yi >= H) continue;
                    double d2 = (sx - xi) * (sx - xi) + (sy - yi) * (sy - yi);
                    buf[yi * W + xi] += peak * std::exp(-d2 * inv2s2);
                }
        }
        Buffer cur = transformBuffer(buf, W, H, DX, DY);
        Donuts::Registrar g;
        g.setReference(buf.data(), W, H);
        Donuts::Transform t = g.measure(cur.data(), W, H);
        QVERIFY(t.valid());
        QVERIFY(std::abs(t.dx - DX) < 0.3);
        QVERIFY(std::abs(t.dy - DY) < 0.3);
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

TestDonutsRegistrar::Buffer TestDonutsRegistrar::makeStarField(int w, int h, uint32_t seed)
{
    Buffer img(w * h, 1000.0);

    uint32_t s = seed;
    auto rng = [&]() -> double {
        s ^= s << 13; s ^= s >> 17; s ^= s << 5;
        return (s & 0xFFFFu) / 65536.0;
    };

    constexpr double margin  = 20.0;
    constexpr double sigma   = 1.5;
    constexpr double inv2s2  = 1.0 / (2.0 * sigma * sigma);
    constexpr int    radius  = static_cast<int>(4.0 * sigma) + 1;

    for (int i = 0; i < 80; ++i)
    {
        double sx   = margin + rng() * (w - 2 * margin);
        double sy   = margin + rng() * (h - 2 * margin);
        double peak = 5000.0 + rng() * 45000.0;
        int x0 = static_cast<int>(sx), y0 = static_cast<int>(sy);
        for (int dy = -radius; dy <= radius; ++dy)
            for (int dx = -radius; dx <= radius; ++dx)
            {
                int xi = x0 + dx, yi = y0 + dy;
                if (xi < 0 || xi >= w || yi < 0 || yi >= h) continue;
                double d2 = (sx - xi) * (sx - xi) + (sy - yi) * (sy - yi);
                img[yi * w + xi] += peak * std::exp(-d2 * inv2s2);
            }
    }
    return img;
}

TestDonutsRegistrar::Buffer TestDonutsRegistrar::transformBuffer(
    const Buffer &src, int w, int h, double dx, double dy, double dthetaDeg)
{
    Buffer dst(w * h, 0.0);
    const double cosT = std::cos(dthetaDeg * M_PI / 180.0);
    const double sinT = std::sin(dthetaDeg * M_PI / 180.0);
    const double cx = w / 2.0, cy = h / 2.0;

    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x)
        {
            double x1 = x - dx, y1 = y - dy;
            double x2 = (x1 - cx) * cosT + (y1 - cy) * sinT + cx;
            double y2 = -(x1 - cx) * sinT + (y1 - cy) * cosT + cy;
            if (x2 >= 0 && x2 < w - 1 && y2 >= 0 && y2 < h - 1)
            {
                int    ix = static_cast<int>(x2), iy = static_cast<int>(y2);
                double fx = x2 - ix, fy = y2 - iy;
                dst[y * w + x] =
                    (1-fx)*(1-fy)*src[iy*w+ix]     + fx*(1-fy)*src[iy*w+ix+1] +
                    (1-fx)*fy   *src[(iy+1)*w+ix]  + fx*fy   *src[(iy+1)*w+ix+1];
            }
        }
    return dst;
}

QTEST_GUILESS_MAIN(TestDonutsRegistrar)
