/*
 * Standalone test suite for the DONUTS guider core algorithm.
 * Uses DonutsCore (donutscore.h) directly -- no Qt, no FITSData required.
 *
 * Build (from repo root):
 *   c++ -std=c++17 -O2 \
 *       -I kstars/ekos/guide/internalguide \
 *       -I /opt/homebrew/include \
 *       Tests/ekos/guide/standalone_donuts_test.cpp \
 *       /opt/homebrew/lib/libcfitsio.dylib \
 *       -o Tests/ekos/guide/standalone_donuts_test
 *
 * Run:
 *   ./Tests/ekos/guide/standalone_donuts_test Tests/fitsviewer/m47_sim_stars.fits
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <fitsio.h>

#include "donutscore.h"

using namespace std;
using namespace DonutsCore;

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

// Bilinear-interpolation rigid-body transform (same as the Qt test helper).
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

static void computeStats(const vector<double> &buf, double &median, double &stddev)
{
    vector<double> sorted = buf;
    sort(sorted.begin(), sorted.end());
    median = sorted[sorted.size() / 2];
    double var = 0.0;
    for (double v : buf) var += (v - median) * (v - median);
    stddev = sqrt(var / buf.size());
}

// ---------------------------------------------------------------------------
// Test cases
// ---------------------------------------------------------------------------

static bool testHotPixelResilience(const vector<double> &clean, int w, int h,
                                   double median, double stddev)
{
    cout << "\n-- Test: Hot Pixel Resilience (dx=1.5, dy=-1.2, rot=0.3 deg) --\n";
    const double clip = 65535.0 * 0.95;
    const double tx = 1.5, ty = -1.2, tr = 0.3;

    vector<double> refBuf = clean;
    // Inject 10 static hot pixels identical in both frames.
    for (int i = 0; i < 10; ++i)
        refBuf[(100 + i * 50) * w + (100 + i * 40)] = 65000.0;

    ProfileSet ref = buildProfiles(refBuf.data(), w, h, median, stddev, clip);

    vector<double> curBuf = transformImage(refBuf, w, h, tx, ty, tr);
    for (int i = 0; i < 10; ++i)
        curBuf[(100 + i * 50) * w + (100 + i * 40)] = 65000.0;

    ProfileSet curr   = buildProfiles(curBuf.data(), w, h, median, stddev, clip);
    Transform  result = solveTransform(ref, curr);

    cout << "  Recovered: dx=" << result.dx << " dy=" << result.dy
         << " rot=" << result.dtheta * 180.0 / M_PI << " deg  SNR=" << result.snr << "\n";

    bool ok = fabs(result.dx - tx) < 0.2 && fabs(result.dy - ty) < 0.2;
    cout << "  Result: " << (ok ? "PASS" : "FAIL") << "\n";
    return ok;
}

static bool testLargeShift(const vector<double> &clean, int w, int h,
                           double median, double stddev)
{
    cout << "\n-- Test: Large Shift (dx=20, dy=-15, rot=0) --\n";
    const double clip  = 65535.0 * 0.95;
    const double tx = 20.0, ty = -15.0;

    ProfileSet ref  = buildProfiles(clean.data(), w, h, median, stddev, clip);
    vector<double> cur = transformImage(clean, w, h, tx, ty, 0.0);
    ProfileSet curr = buildProfiles(cur.data(), w, h, median, stddev, clip);
    Transform  res  = solveTransform(ref, curr);

    cout << "  Recovered: dx=" << res.dx << " dy=" << res.dy
         << "  SNR=" << res.snr << "\n";

    bool ok = fabs(res.dx - tx) < 0.5 && fabs(res.dy - ty) < 0.5 && res.snr > 3.0;
    cout << "  Result: " << (ok ? "PASS" : "FAIL") << "\n";
    return ok;
}

static bool testNoiseInjection(const vector<double> &clean, int w, int h,
                               double median, double stddev)
{
    cout << "\n-- Test: Noise Injection (dx=3, dy=-2, noise=stddev*0.5) --\n";
    const double clip = 65535.0 * 0.95;
    const double tx = 3.0, ty = -2.0;

    // Deterministic LCG noise so the test is reproducible.
    vector<double> noisy = transformImage(clean, w, h, tx, ty, 0.0);
    uint32_t seed = 42;
    for (double &v : noisy)
    {
        seed = seed * 1664525u + 1013904223u;
        double noise = (static_cast<double>(seed) / 0xFFFFFFFFu - 0.5) * stddev * 0.5;
        v = max(0.0, v + noise);
    }

    ProfileSet ref  = buildProfiles(clean.data(), w, h, median, stddev, clip);
    ProfileSet curr = buildProfiles(noisy.data(), w, h, median, stddev, clip);
    Transform  res  = solveTransform(ref, curr);

    cout << "  Recovered: dx=" << res.dx << " dy=" << res.dy
         << "  SNR=" << res.snr << "\n";

    bool ok = fabs(res.dx - tx) < 0.5 && fabs(res.dy - ty) < 0.5 && res.snr > 3.0;
    cout << "  Result: " << (ok ? "PASS" : "FAIL") << "\n";
    return ok;
}

static bool testLostStar(int w, int h, double median, double stddev)
{
    cout << "\n-- Test: Lost Star / Low SNR (blank frame -> snr < 3) --\n";
    const double clip = 65535.0 * 0.95;

    // Uniform grey frame -- no stars, no gradient.
    vector<double> refBuf(w * h, median);
    vector<double> blankBuf(w * h, median);

    ProfileSet ref  = buildProfiles(refBuf.data(),   w, h, median, stddev, clip);
    ProfileSet curr = buildProfiles(blankBuf.data(), w, h, median, stddev, clip);
    Transform  res  = solveTransform(ref, curr);

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

    double median = 0.0, stddev = 0.0;
    computeStats(clean, median, stddev);

    int passed = 0, total = 0;

    auto run = [&](bool result) { ++total; if (result) ++passed; };

    run(testHotPixelResilience(clean, w, h, median, stddev));
    run(testLargeShift        (clean, w, h, median, stddev));
    run(testNoiseInjection    (clean, w, h, median, stddev));
    run(testLostStar          (      w, h, median, stddev));

    cout << "\n========================================\n";
    cout << "Results: " << passed << "/" << total << " tests passed\n";
    if (passed == total)
        cout << "ALL TESTS PASSED\n";
    else
        cout << "SOME TESTS FAILED\n";

    return (passed == total) ? 0 : 1;
}
