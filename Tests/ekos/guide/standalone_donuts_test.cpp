/*
 * SPDX-FileCopyrightText: 2026 Christian Kemper <ckemper@gmail.com>
 *
 * SPDX-License-Identifier: GPL-2.0-or-later
 *
 * Standalone test suite for the DONUTS registration algorithm.
 * Uses the production Donuts::Registrar API -- same code path as the Qt build.
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
// Image manipulation helpers
// ---------------------------------------------------------------------------

// Paint a satellite streak from (x0,y0) to (x1,y1) with a Gaussian cross-section.
// peak   -- ADU added at the centreline per pixel
// sigma  -- half-width of the streak in pixels (typically 0.5-1.5)
static void addStreak(vector<double> &buf, int w, int h,
                      double x0, double y0, double x1, double y1,
                      double peak, double sigma = 1.0)
{
    double lx = x1 - x0, ly = y1 - y0;
    double len2 = lx * lx + ly * ly;
    if (len2 == 0.0) return;

    int margin = static_cast<int>(4.0 * sigma) + 2;
    int xlo = max(0,   (int)min(x0, x1) - margin);
    int xhi = min(w-1, (int)max(x0, x1) + margin);
    int ylo = max(0,   (int)min(y0, y1) - margin);
    int yhi = min(h-1, (int)max(y0, y1) + margin);

    for (int y = ylo; y <= yhi; ++y)
        for (int x = xlo; x <= xhi; ++x)
        {
            // perpendicular distance from pixel centre to the line segment
            double t  = ((x - x0) * lx + (y - y0) * ly) / len2;
            t = max(0.0, min(1.0, t));
            double px = x0 + t * lx, py = y0 + t * ly;
            double d2 = (x - px) * (x - px) + (y - py) * (y - py);
            buf[y * w + x] += peak * exp(-d2 / (2.0 * sigma * sigma));
        }
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

    Donuts::Registrar g;
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

    Donuts::Registrar g;
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

    Donuts::Registrar g;
    g.setReference(clean.data(), w, h);
    auto res = g.measure(noisy.data(), w, h);

    cout << "  Recovered: dx=" << res.dx << " dy=" << res.dy
         << "  SNR=" << res.snr << "\n";

    bool ok = fabs(res.dx - tx) < 0.5 && fabs(res.dy - ty) < 0.5 && res.snr > 3.0;
    cout << "  Result: " << (ok ? "PASS" : "FAIL") << "\n";
    return ok;
}

static bool testLostStar(const vector<double> &clean, int w, int h)
{
    cout << "\n-- Test: Lost Star / Low SNR (blank current frame -> snr < 3) --\n";

    // Real star-field reference; blank current (all zeros -- every pixel below
    // the stellar threshold so all profiles are zero).  Zero profiles produce a
    // flat correlation output and SNR must be < 3.
    Donuts::Registrar g;
    g.setReference(clean.data(), w, h);

    vector<double> blank(w * h, 0.0);
    auto res = g.measure(blank.data(), w, h);

    cout << "  SNR=" << res.snr << " (must be < 3)\n";

    bool ok = res.snr < 3.0;
    cout << "  Result: " << (ok ? "PASS" : "FAIL") << "\n";
    return ok;
}

// Sweep satellite streak brightness levels across one diagonal angle.
// The test is intentionally diagnostic: it prints per-level results so you can
// see exactly where the algorithm starts to struggle, rather than just PASS/FAIL.
static bool testSatelliteStreak(const vector<double> &clean, int w, int h)
{
    cout << "\n-- Test: Satellite Streak (dx=3.0, dy=-2.5) --\n";
    const double tx = 3.0, ty = -2.5;
    const double tol = 0.5;

    // Diagonal streak across the frame (lower-left to upper-right).
    const double sx0 = 0,   sy0 = h * 0.75;
    const double sx1 = w-1, sy1 = h * 0.25;
    const double streakSigma = 1.0;  // ~1px half-width, sharp satellite

    // Reference frame has no streak.
    Donuts::Registrar g;
    g.setReference(clean.data(), w, h);

    // ADU/pixel levels: dim debris -> typical LEO -> bright Starlink -> saturated flare.
    // expectDegrade: diagonal flare at 60K ADU deposits flux into every row, so the
    // 1D projection gets a broad hump rather than a spike -- peak/surround stays ~1,
    // below the detection threshold.  This is a known limit of the projection approach.
    struct Level { double peak; const char *label; bool expectDegrade; };
    const Level levels[] = {
        {    500.0, "dim debris    ~mag6",  false},
        {   3000.0, "typical LEO   ~mag4",  false},
        {  15000.0, "bright Starlink~mag2", false},
        {  60000.0, "saturated flare     ", true },
    };

    bool allPass = true;
    for (const auto &lv : levels)
    {
        vector<double> cur = transformImage(clean, w, h, tx, ty, 0.0);
        addStreak(cur, w, h, sx0, sy0, sx1, sy1, lv.peak, streakSigma);

        auto res = g.measure(cur.data(), w, h);
        double ex = fabs(res.dx - tx), ey = fabs(res.dy - ty);
        bool ok = ex < tol && ey < tol && res.valid();

        cout << "  " << lv.label
             << "  peak=" << (int)lv.peak
             << "  dx_err=" << ex << " dy_err=" << ey
             << "  SNR=" << res.snr
             << "  -> " << (ok ? "OK" : (lv.expectDegrade ? "DEGRADED (expected)" : "DEGRADED")) << "\n";

        if (!ok && !lv.expectDegrade) allPass = false;
    }

    // Horizontal streak: worst case for the Y-projection profile.
    // Every pixel in the streak row sums into one yProf bin -> enormous spike.
    {
        vector<double> cur = transformImage(clean, w, h, tx, ty, 0.0);
        addStreak(cur, w, h, 0, h/2, w-1, h/2, 15000.0, streakSigma);

        auto res = g.measure(cur.data(), w, h);
        double ex = fabs(res.dx - tx), ey = fabs(res.dy - ty);
        bool ok = ex < tol && ey < tol && res.valid();
        cout << "  horizontal streak peak=15000"
             << "  dx_err=" << ex << " dy_err=" << ey
             << "  SNR=" << res.snr
             << "  -> " << (ok ? "OK" : "DEGRADED") << "\n";
        if (!ok) allPass = false;
    }

    // Vertical streak: symmetric case -- spike lands in xProf instead.
    {
        vector<double> cur = transformImage(clean, w, h, tx, ty, 0.0);
        addStreak(cur, w, h, w/2, 0, w/2, h-1, 15000.0, streakSigma);

        auto res = g.measure(cur.data(), w, h);
        double ex = fabs(res.dx - tx), ey = fabs(res.dy - ty);
        bool ok = ex < tol && ey < tol && res.valid();
        cout << "  vertical   streak peak=15000"
             << "  dx_err=" << ex << " dy_err=" << ey
             << "  SNR=" << res.snr
             << "  -> " << (ok ? "OK" : "DEGRADED") << "\n";
        if (!ok) allPass = false;
    }

    // Note: saturated diagonal flare (60 000 ADU) degrades to ~2 px -- known
    // limitation.  Diagonal streaks distribute flux across every row so the
    // per-row spike stays comparable to stellar peaks and is not detectable.

    cout << "  Result: " << (allPass ? "PASS" : "SOME LEVELS DEGRADED") << "\n";
    return allPass;
}

// Model blinking airplane navigation lights.
//
// Physics (typical commercial flight, 45-deg elevation, 5 arcsec/px guide scale):
//   Angular velocity  ~500 px/s  ->  25 px per 50 ms strobe flash.
//   The airplane crosses the guide FOV in ~0.3 s, leaving 0-1 flashes.
//   Here we place 3 flashes (slow/far airplane, or multiple exposures stacked)
//   along a diagonal path to stress-test the algorithm.
//
// Key difference from a satellite streak: each flash is a SHORT dash (~25 px)
// rather than a line spanning the full image width.  The per-row profile spike
// is therefore only  flash_length * peak  (not  image_width * peak), which is
// much smaller relative to stellar peaks and may not trigger the ridge suppressor.
static bool testAirplaneLights(const vector<double> &clean, int w, int h)
{
    cout << "\n-- Test: Airplane Strobe Lights (dx=3.0, dy=-2.5, 3 flashes diagonal) --\n";
    const double tx = 3.0, ty = -2.5;
    const double tol = 0.5;
    const double flashLen = 25.0;   // pixels; 50 ms flash at ~500 px/s
    const double flashSigma = 1.5;  // slightly wider PSF than a point star

    // Three flashes spread along a diagonal path across the frame.
    struct Flash { double x0, y0, x1, y1; };
    const Flash flashes[] = {
        { w*0.20, h*0.30, w*0.20 + flashLen*0.71, h*0.30 + flashLen*0.71 },
        { w*0.50, h*0.50, w*0.50 + flashLen*0.71, h*0.50 + flashLen*0.71 },
        { w*0.75, h*0.68, w*0.75 + flashLen*0.71, h*0.68 + flashLen*0.71 },
    };

    Donuts::Registrar g;
    g.setReference(clean.data(), w, h);

    struct Level { double peak; const char *label; };
    const Level levels[] = {
        {   5000.0, "distant/dim      "},
        {  15000.0, "Starlink-equiv   "},
        {  60000.0, "bright strobe    "},
        { 200000.0, "very bright flare"},
    };

    bool allPass = true;
    for (const auto &lv : levels)
    {
        vector<double> cur = transformImage(clean, w, h, tx, ty, 0.0);
        for (const auto &f : flashes)
            addStreak(cur, w, h, f.x0, f.y0, f.x1, f.y1, lv.peak, flashSigma);

        auto res = g.measure(cur.data(), w, h);
        double ex = fabs(res.dx - tx), ey = fabs(res.dy - ty);
        bool ok = ex < tol && ey < tol && res.valid();

        cout << "  " << lv.label
             << "  peak=" << (int)lv.peak
             << "  dx_err=" << ex << " dy_err=" << ey
             << "  SNR=" << res.snr
             << "  -> " << (ok ? "OK" : "DEGRADED") << "\n";
        if (!ok) allPass = false;
    }

    cout << "  Result: " << (allPass ? "PASS" : "SOME LEVELS DEGRADED") << "\n";
    return allPass;
}

// Verify that incremental field-rotation de-rotation handles large accumulated
// alt-az field rotation.  Without pre-rotation the algorithm fails silently
// above ~2 deg (dx_err > 30 px at 3 deg, 80 px at 10 deg).  With 0.5-deg
// steps the per-frame residual stays in the sub-degree regime where the
// phase-correlation and WLS solver remain accurate.
static bool testRotationTracking(const vector<double> &clean, int w, int h)
{
    cout << "\n-- Test: Rotation Tracking (0.5-deg steps to 15 deg, dx=2.0 dy=-1.5) --\n";
    const double tx = 2.0, ty = -1.5;
    // 1.0px tolerance: at 15 deg accumulated rotation, sub-pixel accuracy is a
    // bonus.  Without de-rotation the error would exceed 80 px at 10 deg.
    const double tol = 1.0;

    Donuts::Registrar g;
    g.setReference(clean.data(), w, h);

    // The Registrar tracks accumulated rotation internally; no external state needed.
    bool allPass = true;

    // 0.5-deg steps simulate alt-az guiding at 15 deg/hr with a 2-min cadence.
    // Only a subset of steps is printed to keep output concise.
    const int printAt[] = { 1, 2, 4, 6, 10, 15, 20, 25, 30, -1 };
    int printIdx = 0;

    for (int i = 1; i <= 30; ++i)
    {
        double rot_deg = i * 0.5;
        auto cur = transformImage(clean, w, h, tx, ty, rot_deg);
        auto res = g.measure(cur.data(), w, h);

        double dx_err     = fabs(res.dx - tx);
        double dy_err     = fabs(res.dy - ty);
        double dtheta_err = fabs(res.dtheta * 180.0 / M_PI - rot_deg);
        bool ok = dx_err < tol && dy_err < tol;
        if (!ok) allPass = false;

        if (printAt[printIdx] == i)
        {
            printf("  rot=%5.1f deg: dx_err=%.3f  dy_err=%.3f"
                   "  dtheta_err=%.3f deg  SNR=%.1f  -> %s\n",
                   rot_deg, dx_err, dy_err, dtheta_err, res.snr,
                   ok ? "OK" : "FAIL");
            ++printIdx;
        }
    }

    cout << "  Result: " << (allPass ? "PASS" : "FAIL") << "\n";
    return allPass;
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
    run(testLostStar          (clean, w, h));
    run(testSatelliteStreak   (clean, w, h));
    run(testAirplaneLights    (clean, w, h));
    run(testRotationTracking  (clean, w, h));

    cout << "\n========================================\n";
    cout << "Results: " << passed << "/" << total << " tests passed\n";
    if (passed == total)
        cout << "ALL TESTS PASSED\n";
    else
        cout << "SOME TESTS FAILED\n";

    return (passed == total) ? 0 : 1;
}
