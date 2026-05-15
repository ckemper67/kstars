/*
 * SPDX-FileCopyrightText: 2026 Christian Kemper <ckemper@gmail.com>
 *
 * SPDX-License-Identifier: GPL-2.0-or-later
 *
 * Standalone test suite for the DONUTS guider algorithm.
 * Uses the production Donuts::Guider API -- same code path as the Qt build.
 * No Qt, no FITSData required.
 *
 * Build (from repo root):
 *   c++ -std=c++17 -O2 \
 *       -I kstars/ekos/guide/donuts \
 *       -I kstars/ekos/guide/internalguide \
 *       -I /opt/homebrew/include \
 *       Tests/ekos/guide/standalone_donuts_test.cpp \
 *       kstars/ekos/guide/donuts/donuts.cpp \
 *       /opt/homebrew/lib/libcfitsio.dylib \
 *       -o Tests/ekos/guide/standalone_donuts_test
 *
 * Run:
 *   ./Tests/ekos/guide/standalone_donuts_test Tests/fitsviewer/m47_sim_stars.fits
 */

#include "donuts.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <fitsio.h>

using namespace std;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static vector<double> loadFITS(const string &filename, int &w, int &h)
{
    fitsfile *fptr = nullptr;
    int s = 0;
    if (fits_open_file(&fptr, filename.c_str(), READONLY, &s)) return {};
    long naxes[2] = {};
    fits_get_img_size(fptr, 2, naxes, &s);
    w = naxes[0]; h = naxes[1];
    vector<double> buf(w * h);
    fits_read_img(fptr, TDOUBLE, 1, (long)w * h, nullptr, buf.data(), nullptr, &s);
    fits_close_file(fptr, &s);
    return buf;
}

// Bilinear-interpolation rigid-body transform.
static vector<double> transformImage(
    const vector<double> &src, int w, int h, double dx, double dy, double dtheta)
{
    vector<double> dst(w * h, 0.0);
    const double rad = dtheta * M_PI / 180.0;
    const double co  = cos(rad), si = sin(rad);
    const double mx  = w / 2.0,  my = h / 2.0;

    for (int y = 0; y < h; ++y)
    {
        for (int x = 0; x < w; ++x)
        {
            double x1 = x - dx - mx, y1 = y - dy - my;
            double sx  = x1 * co + y1 * si + mx;
            double sy  = -x1 * si + y1 * co + my;
            if (sx >= 0 && sx < w - 1 && sy >= 0 && sy < h - 1)
            {
                int    ix = (int)sx, iy = (int)sy;
                double fx = sx - ix,  fy = sy - iy;
                dst[y * w + x] =
                    (1 - fx) * (1 - fy) * src[iy * w + ix]          +
                    fx       * (1 - fy) * src[iy * w + ix + 1]       +
                    (1 - fx) * fy       * src[(iy + 1) * w + ix]     +
                    fx       * fy       * src[(iy + 1) * w + ix + 1];
            }
        }
    }
    return dst;
}

// ---------------------------------------------------------------------------
// Test cases
// ---------------------------------------------------------------------------

static bool testHotPixelResilience(const vector<double> &clean, int w, int h)
{
    cout << "\n-- Test: Hot Pixel Resilience (dx=1.5, dy=-1.2, rot=0.3 deg) --\n";
    const double tx = 1.5, ty = -1.2, tr = 0.3;

    vector<double> refBuf = clean;
    for (int i = 0; i < 10; ++i)
        refBuf[(100 + i * 50) * w + (100 + i * 40)] = 65000.0;

    vector<double> curBuf = transformImage(refBuf, w, h, tx, ty, tr);
    // Hot pixels are static -- same detector defects in both frames.
    for (int i = 0; i < 10; ++i)
        curBuf[(100 + i * 50) * w + (100 + i * 40)] = 65000.0;

    Donuts::Guider g;
    g.setReference(refBuf.data(), w, h);
    auto result = g.measure(curBuf.data(), w, h);

    cout << "  Recovered: dx=" << result.dx << " dy=" << result.dy
         << "  SNR=" << result.snr << "\n";

    bool ok = fabs(result.dx - tx) < 0.2 && fabs(result.dy - ty) < 0.2;
    cout << "  Result: " << (ok ? "PASS" : "FAIL") << "\n";
    return ok;
}

static bool testLargeShift(const vector<double> &clean, int w, int h)
{
    cout << "\n-- Test: Large Shift (dx=20, dy=-15, rot=0) --\n";
    const double tx = 20.0, ty = -15.0;

    Donuts::Guider g;
    g.setReference(clean.data(), w, h);

    auto cur = transformImage(clean, w, h, tx, ty, 0.0);
    auto res = g.measure(cur.data(), w, h);

    cout << "  Recovered: dx=" << res.dx << " dy=" << res.dy
         << "  SNR=" << res.snr << "\n";

    bool ok = fabs(res.dx - tx) < 0.5 && fabs(res.dy - ty) < 0.5 && res.snr > 3.0;
    cout << "  Result: " << (ok ? "PASS" : "FAIL") << "\n";
    return ok;
}

static bool testNoiseInjection(const vector<double> &clean, int w, int h)
{
    cout << "\n-- Test: Noise Injection (dx=3, dy=-2, noise=stddev*0.5) --\n";
    const double tx = 3.0, ty = -2.0;

    // Estimate stddev for noise scaling.
    vector<double> sorted = clean;
    sort(sorted.begin(), sorted.end());
    double median = sorted[sorted.size() / 2];
    double var = 0.0;
    for (double v : clean) var += (v - median) * (v - median);
    double stddev = sqrt(var / clean.size());

    vector<double> noisy = transformImage(clean, w, h, tx, ty, 0.0);
    uint32_t seed = 42;
    for (double &v : noisy)
    {
        seed = seed * 1664525u + 1013904223u;
        double noise = (static_cast<double>(seed) / 0xFFFFFFFFu - 0.5) * stddev * 0.5;
        v = max(0.0, v + noise);
    }

    Donuts::Guider g;
    g.setReference(clean.data(), w, h);
    auto res = g.measure(noisy.data(), w, h);

    cout << "  Recovered: dx=" << res.dx << " dy=" << res.dy
         << "  SNR=" << res.snr << "\n";

    bool ok = fabs(res.dx - tx) < 0.5 && fabs(res.dy - ty) < 0.5 && res.snr > 3.0;
    cout << "  Result: " << (ok ? "PASS" : "FAIL") << "\n";
    return ok;
}

static bool testLostStar(int w, int h)
{
    cout << "\n-- Test: Lost Star / Low SNR (blank frame -> snr < 3) --\n";

    // Uniform grey reference -- no stars, no structure.
    vector<double> grey(w * h, 1000.0);

    Donuts::Guider g;
    g.setReference(grey.data(), w, h);
    auto res = g.measure(grey.data(), w, h);

    cout << "  SNR=" << res.snr << " (must be < 3)\n";

    bool ok = res.snr < 3.0;
    cout << "  Result: " << (ok ? "PASS" : "FAIL") << "\n";
    return ok;
}

// ---------------------------------------------------------------------------

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        cerr << "Usage: " << argv[0] << " <fits_file>\n";
        return 1;
    }

    int w = 0, h = 0;
    vector<double> clean = loadFITS(argv[1], w, h);
    if (clean.empty()) { cerr << "Failed to load FITS: " << argv[1] << "\n"; return 1; }

    int passed = 0, total = 0;
    auto run = [&](bool result) { ++total; if (result) ++passed; };

    run(testHotPixelResilience(clean, w, h));
    run(testLargeShift        (clean, w, h));
    run(testNoiseInjection    (clean, w, h));
    run(testLostStar          (      w, h));

    cout << "\n========================================\n";
    cout << "Results: " << passed << "/" << total << " tests passed\n";
    if (passed == total)
        cout << "ALL TESTS PASSED\n";
    else
        cout << "SOME TESTS FAILED\n";

    return (passed == total) ? 0 : 1;
}
