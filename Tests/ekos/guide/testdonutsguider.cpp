/*
    SPDX-FileCopyrightText: 2024 Christian Kemper <ckemper@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

#include "testdonutsguider.h"
#include "ekos/guide/donuts/donuts.h"
#include <QTest>
#include <cmath>
#include <cstdint>

TestDonutsGuider::TestDonutsGuider() : QObject()
{
}

void TestDonutsGuider::initTestCase() {}
void TestDonutsGuider::cleanupTestCase() {}

void TestDonutsGuider::testMathRecovery()
{
    const int W = 256, H = 256;
    Buffer master = makeStarField(W, H);

    Donuts::Guider guider;
    guider.setReference(master.data(), W, H);

    // Test 1: Pure translation
    {
        double tx = 5.4, ty = -3.2;
        Buffer shifted = transformBuffer(master, W, H, tx, ty, 0);
        Donuts::Transform res = guider.measure(shifted.data(), W, H);
        qDebug() << "Test 1 (Translation): Expected (" << tx << "," << ty
                 << ") Got (" << res.dx << "," << res.dy << ")";
        QVERIFY(std::abs(res.dx - tx) < 0.1);
        QVERIFY(std::abs(res.dy - ty) < 0.1);
    }

    // Test 2: Pure rotation
    {
        double rotDeg = 0.5;
        Buffer rotated = transformBuffer(master, W, H, 0, 0, rotDeg);
        Donuts::Transform res = guider.measure(rotated.data(), W, H);
        double gotDeg = res.dtheta * 180.0 / M_PI;
        qDebug() << "Test 2 (Rotation): Expected" << rotDeg << "deg Got" << gotDeg << "deg";
        QVERIFY(std::abs(gotDeg - rotDeg) < 0.1);
    }

    // Test 3: Combined translation + rotation
    {
        Buffer combined = transformBuffer(master, W, H, 2.0, 2.0, 0.2);
        Donuts::Transform res = guider.measure(combined.data(), W, H);
        double gotDeg = res.dtheta * 180.0 / M_PI;
        qDebug() << "Test 3 (Combined): Expected (2.0, 2.0, 0.2) Got ("
                 << res.dx << "," << res.dy << "," << gotDeg << ")";
        // 4-quadrant algorithm trades ~0.15px accuracy for rotation detection under combined
        // rotation+translation (stars near quadrant boundaries shift across them).
        QVERIFY(std::abs(res.dx - 2.0) < 0.2);
        QVERIFY(std::abs(res.dy - 2.0) < 0.2);
        QVERIFY(std::abs(gotDeg - 0.2) < 0.05);
    }
}

void TestDonutsGuider::testOffCenterRotation()
{
    // Rotating around an off-center pivot is equivalent to rotating around the image
    // center plus a residual translation.  Given pivot (px, py) and image center (cx, cy):
    //   dx = (px-cx)*(1-cos(theta)) + (py-cy)*sin(theta)
    //   dy = (py-cy)*(1-cos(theta)) - (px-cx)*sin(theta)
    // DONUTS should recover this (dx, dy, theta) regardless of the pivot.

    const int W = 256, H = 256;
    Buffer master = makeStarField(W, H);

    Donuts::Guider guider;
    guider.setReference(master.data(), W, H);

    const double cx = W / 2.0;
    const double cy = H / 2.0;
    const double px = cx + 100;  // pivot 100px right of center
    const double py = cy + 50;   // pivot 50px below center
    const double thetaDeg = 1.0;

    const double cosT = std::cos(thetaDeg * M_PI / 180.0);
    const double sinT = std::sin(thetaDeg * M_PI / 180.0);
    const double expectedDx = (px - cx) * (1.0 - cosT) + (py - cy) * sinT;
    const double expectedDy = (py - cy) * (1.0 - cosT) - (px - cx) * sinT;

    Buffer rotated = transformBuffer(master, W, H, expectedDx, expectedDy, thetaDeg);
    Donuts::Transform res = guider.measure(rotated.data(), W, H);
    double gotDeg = res.dtheta * 180.0 / M_PI;

    qDebug() << "Off-center rotation: pivot=(" << px << "," << py << ") theta=" << thetaDeg
             << "deg; expected dx=" << expectedDx << "dy=" << expectedDy
             << "; got dx=" << res.dx << "dy=" << res.dy << "dtheta=" << gotDeg << "deg";

    QVERIFY(std::abs(res.dx - expectedDx) < 0.2);
    QVERIFY(std::abs(res.dy - expectedDy) < 0.2);
    QVERIFY(std::abs(gotDeg - thetaDeg) < 0.05);
}

// ---------------------------------------------------------------------------
// Synthetic double-buffer helpers (no FITSData dependency)
// ---------------------------------------------------------------------------

TestDonutsGuider::Buffer TestDonutsGuider::makeStarField(int w, int h)
{
    Buffer img(w * h, 100.0);

    struct Star { double x, y, peak; };
    const Star stars[] = {
        {40, 30, 8000}, {80, 50, 7000}, {120, 80, 9000}, {160, 40, 6000}, {200, 70, 8500},
        {30, 100, 7500}, {70, 130, 8000}, {110, 110, 9500}, {150, 140, 7000}, {190, 120, 8200},
        {50, 170, 8800}, {90, 190, 7200}, {130, 160, 9000}, {170, 180, 6500}, {210, 200, 8000},
        {20, 210, 7800}, {60, 220, 8300}, {100, 240, 7600}, {140, 230, 9100}, {180, 210, 7400},
        {220, 50, 8600}, {230, 150, 7900}, {240, 230, 8100}, {10, 50, 7700}, {10, 150, 8400},
    };

    constexpr double sigma  = 2.5;
    constexpr double inv2s2 = 1.0 / (2.0 * sigma * sigma);
    constexpr int    radius = static_cast<int>(4 * sigma) + 1;

    for (const auto &s : stars)
    {
        int x0 = static_cast<int>(s.x), y0 = static_cast<int>(s.y);
        for (int dy = -radius; dy <= radius; ++dy)
            for (int dx = -radius; dx <= radius; ++dx)
            {
                int xi = x0 + dx, yi = y0 + dy;
                if (xi < 0 || xi >= w || yi < 0 || yi >= h) continue;
                double d2 = (s.x - xi) * (s.x - xi) + (s.y - yi) * (s.y - yi);
                img[yi * w + xi] += s.peak * std::exp(-d2 * inv2s2);
            }
    }
    return img;
}

TestDonutsGuider::Buffer TestDonutsGuider::transformBuffer(
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

// ---------------------------------------------------------------------------
// Translation-only characterisation tests
// Ported from kstars/ekos/guide/donuts/tests/test_donuts.cpp
// ---------------------------------------------------------------------------

// Generate a frame with N Gaussian stars at pseudo-random positions.
// Simpler than makeStarField: no fixed layout, seed controls reproducibility.
static std::vector<double> makeRandomField(int w, int h, int count, uint32_t seed)
{
    std::vector<double> buf(w * h, 100.0);
    uint32_t s = seed ? seed : 1u;
    auto rng = [&]() -> double {
        s ^= s << 13; s ^= s >> 17; s ^= s << 5;
        return (s & 0xFFFFu) / 65536.0;
    };

    const double margin   = 20.0;
    const double psf_sig  = 2.5;
    const double inv2s2   = 1.0 / (2.0 * psf_sig * psf_sig);
    const int    radius   = static_cast<int>(4 * psf_sig) + 1;

    for (int i = 0; i < count; ++i)
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
                buf[yi * w + xi] += peak * std::exp(-d2 * inv2s2);
            }
    }
    return buf;
}

void TestDonutsGuider::testTranslationMagnitudeSweep()
{
    const int W = 256, H = 256;
    Buffer ref = makeStarField(W, H);

    Donuts::Guider g;
    g.setReference(ref.data(), W, H);

    for (int i = 1; i <= 12; ++i)
    {
        double shift = i * 0.5;
        double tdx   = shift, tdy = shift * 0.7;
        Buffer cur   = transformBuffer(ref, W, H, tdx, tdy, 0.0);
        Donuts::Transform t = g.measure(cur.data(), W, H);
        QVERIFY(t.valid());
        QVERIFY(std::abs(t.dx - tdx) < 0.3);
        QVERIFY(std::abs(t.dy - tdy) < 0.3);
    }
}

void TestDonutsGuider::testTranslationDirections()
{
    const int W = 256, H = 256;
    Buffer ref = makeStarField(W, H);

    Donuts::Guider g;
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
        Buffer cur = transformBuffer(ref, W, H, c.dx, c.dy, 0.0);
        Donuts::Transform t = g.measure(cur.data(), W, H);
        QVERIFY(t.valid());
        QVERIFY(std::abs(t.dx - c.dx) < 0.15);
        QVERIFY(std::abs(t.dy - c.dy) < 0.15);
    }
}

void TestDonutsGuider::testStarDensity()
{
    const int    W = 256, H = 256;
    const double DX = 4.3, DY = -3.1;
    const int    counts[] = {10, 20, 40, 80};

    for (int n : counts)
    {
        auto buf = makeRandomField(W, H, n, 200u + static_cast<uint32_t>(n));
        auto cur = TestDonutsGuider::transformBuffer(
            buf, W, H, DX, DY, 0.0);
        Donuts::Guider g;
        g.setReference(buf.data(), W, H);
        Donuts::Transform t = g.measure(cur.data(), W, H);
        QVERIFY(t.valid());
        QVERIFY(std::abs(t.dx - DX) < 0.3);
        QVERIFY(std::abs(t.dy - DY) < 0.3);
    }
}

QTEST_GUILESS_MAIN(TestDonutsGuider)
