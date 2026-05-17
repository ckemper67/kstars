/*
 * SPDX-FileCopyrightText: 2026 Christian Kemper <ckemper@gmail.com>
 *
 * SPDX-License-Identifier: GPL-2.0-or-later
 *
 * Standalone FITS stacker using the DONUTS rotation/translation solver.
 *
 * Monochrome: aligns and stacks directly.
 * Bayer (BAYERPAT=RGGB): extracts green at half resolution for DONUTS
 * registration, then 2x2 debayers and applies the rigid-body transform to
 * each colour channel before stacking. Output is a 3-plane FITS (R/G/B) at
 * half resolution.
 *
 * Alignment: Catmull-Rom bicubic with clamped boundaries.
 *
 * Frame rejection (two-pass):
 *   1. Hard bounds: |dtheta| <= 5 deg, |dx|/|dy| <= 150 half-res px.
 *   2. MAD sigma-clip (3-sigma) on dx, dy, dtheta of hard-bound survivors.
 *
 * Per-pixel coverage count so zero-padded border regions don't dilute edges.
 *
 * Performance optimizations vs v1:
 *   - float pixel buffers (halves memory, enables 4-wide ARM NEON SIMD)
 *   - 3-channel fused bicubic warp: coordinate math amortized over R/G/B
 *   - ARM NEON fast path for interior pixels (no per-lane bounds clamping)
 *   - Per-row interior x-range precomputed analytically to split fast/slow paths
 *   - Thread-parallel warp, accumulate, and reduce via std::thread (-j N)
 *   - Per-thread scratch buffer eliminates per-pixel malloc in median/sigmaclip
 *   - Per-pixel median+MAD sigma-clip (robust to satellites/airplane trails)
 *   - uint8_t mask (no bit-pack overhead of vector<bool>)
 *   - Welford streaming mode (-m wstream): O(pixels) memory regardless of
 *     frame count; per-pixel online rejection gate (kappa*sigma) using
 *     Welford's numerically stable recurrence. Ideal for long EAA sessions.
 *
 * Build (from repo root):
 *   c++ -std=c++17 -O2 -march=native \
 *       -I kstars/ekos/guide/donuts \
 *       -I /opt/homebrew/include \
 *       Tests/ekos/guide/fits_stack.cpp \
 *       kstars/ekos/guide/donuts/donuts.cpp \
 *       /opt/homebrew/lib/libcfitsio.dylib \
 *       /opt/homebrew/lib/libpng.dylib \
 *       -o Tests/ekos/guide/fits_stack
 *
 * Usage:
 *   ./fits_stack [-o output.fits] [-p output.png]
 *               [-m mean|median|sigmaclip|wstream] [-k kappa] [-j threads]
 *               frame1.fits frame2.fits ...
 *
 *   Inputs are sorted alphabetically; the first is the reference.
 *   Default output: stacked.fits.  Default mode: mean.
 *   -j 0 (or omitted) uses std::thread::hardware_concurrency().
 */

#include "donuts.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fitsio.h>
#include <functional>
#include <iostream>
#include <numeric>
#include <png.h>
#include <string>
#include <thread>
#include <vector>

#ifdef __ARM_NEON
#include <arm_neon.h>
#endif

using namespace std;

using Pix = float;

// ---------------------------------------------------------------------------
// Rejection thresholds
// ---------------------------------------------------------------------------

constexpr double MAX_ROTATION_RAD   = 5.0 * M_PI / 180.0;
constexpr double MAX_TRANSLATION_PX = 150.0;
constexpr double MAD_NSIGMA         = 3.0;

// ---------------------------------------------------------------------------
// Stacking mode
// ---------------------------------------------------------------------------

enum class StackMode { Mean, Median, SigmaClip, WelfordStream };

// ---------------------------------------------------------------------------
// Thread pool: partition [0, n) into nThreads contiguous blocks.
// ---------------------------------------------------------------------------

static void parallelFor(int n, int nThreads, function<void(int, int)> fn)
{
    if (nThreads <= 1 || n <= nThreads)
    {
        fn(0, n);
        return;
    }
    vector<thread> threads;
    threads.reserve(nThreads);
    int chunk = (n + nThreads - 1) / nThreads;
    for (int t = 0; t < nThreads; ++t)
    {
        int lo = t * chunk, hi = min(lo + chunk, n);
        if (lo >= n) break;
        threads.emplace_back(fn, lo, hi);
    }
    for (auto &t : threads) t.join();
}

// ---------------------------------------------------------------------------
// FITS I/O
// ---------------------------------------------------------------------------

static vector<Pix> loadFITS(const string &path, int &w, int &h)
{
    fitsfile *fptr = nullptr;
    int s = 0;
    if (fits_open_file(&fptr, path.c_str(), READONLY, &s))
    {
        cerr << "Cannot open " << path << " (cfitsio status " << s << ")\n";
        return {};
    }
    long naxes[2] = {};
    fits_get_img_size(fptr, 2, naxes, &s);
    w = (int)naxes[0];
    h = (int)naxes[1];
    // Read as double for correct BZERO/BSCALE handling, then narrow to float.
    vector<double> tmp((size_t)w * h);
    fits_read_img(fptr, TDOUBLE, 1, (long)w * h, nullptr, tmp.data(), nullptr, &s);
    fits_close_file(fptr, &s);
    if (s) { cerr << "Read error on " << path << " (status " << s << ")\n"; return {}; }
    vector<Pix> buf(tmp.size());
    for (size_t i = 0; i < tmp.size(); ++i) buf[i] = (Pix)tmp[i];
    return buf;
}

static bool isBayerRGGB(const string &path)
{
    fitsfile *fptr = nullptr;
    int s = 0;
    if (fits_open_file(&fptr, path.c_str(), READONLY, &s)) return false;
    char value[FLEN_VALUE] = {};
    fits_read_key(fptr, TSTRING, "BAYERPAT", value, nullptr, &s);
    fits_close_file(fptr, &s);
    if (s != 0) return false;
    string pat(value);
    pat.erase(0, pat.find_first_not_of(" '\""));
    pat.erase(pat.find_last_not_of(" '\"") + 1);
    return pat == "RGGB";
}

static bool writeFITS(const string &path, const vector<Pix> &pixels, int w, int h)
{
    fitsfile *fptr = nullptr;
    int s = 0;
    fits_create_file(&fptr, ("!" + path).c_str(), &s);
    long naxes[2] = { (long)w, (long)h };
    fits_create_img(fptr, FLOAT_IMG, 2, naxes, &s);
    fits_write_img(fptr, TFLOAT, 1, (long)w * h,
                   const_cast<Pix *>(pixels.data()), &s);
    fits_close_file(fptr, &s);
    if (s) { cerr << "Write error on " << path << " (status " << s << ")\n"; return false; }
    return true;
}

static bool writeFITS3(const string &path,
                       const vector<Pix> &r, const vector<Pix> &g, const vector<Pix> &b,
                       int w, int h)
{
    fitsfile *fptr = nullptr;
    int s = 0;
    fits_create_file(&fptr, ("!" + path).c_str(), &s);
    long naxes[3] = { (long)w, (long)h, 3 };
    fits_create_img(fptr, FLOAT_IMG, 3, naxes, &s);
    long npix = (long)w * h;
    fits_write_img(fptr, TFLOAT,           1, npix, const_cast<Pix *>(r.data()), &s);
    fits_write_img(fptr, TFLOAT,   npix + 1, npix, const_cast<Pix *>(g.data()), &s);
    fits_write_img(fptr, TFLOAT, 2*npix + 1, npix, const_cast<Pix *>(b.data()), &s);
    fits_close_file(fptr, &s);
    if (s) { cerr << "Write error on " << path << " (status " << s << ")\n"; return false; }
    return true;
}

// ---------------------------------------------------------------------------
// PNG output: median background subtraction + linked luminance asinh stretch
// ---------------------------------------------------------------------------

static vector<Pix> subtractMedian(vector<Pix> ch)
{
    vector<Pix> s(ch);
    std::sort(s.begin(), s.end());
    Pix med = s[s.size() / 2];
    for (Pix &v : ch) v -= med;
    return ch;
}

static vector<uint8_t> asinhStretch(const vector<Pix> &ch, double lo, double hi)
{
    double softening = 0.05 * (hi - lo);
    if (softening < 1.0) softening = 1.0;
    double norm = std::asinh((hi - lo) / softening);
    vector<uint8_t> out(ch.size());
    for (size_t i = 0; i < ch.size(); ++i)
    {
        double v = std::asinh(((double)ch[i] - lo) / softening) / norm * 255.0;
        out[i] = (uint8_t)std::max(0.0, std::min(255.0, v));
    }
    return out;
}

static bool writePNG(const string &path,
                     const vector<Pix> &r, const vector<Pix> &g, const vector<Pix> &b,
                     int w, int h)
{
    auto rn = subtractMedian(r);
    auto gn = subtractMedian(g);
    auto bn = subtractMedian(b);

    vector<Pix> lum(rn.size());
    for (size_t i = 0; i < lum.size(); ++i)
        lum[i] = 0.299f * rn[i] + 0.587f * gn[i] + 0.114f * bn[i];

    vector<Pix> slum(lum);
    std::sort(slum.begin(), slum.end());
    double lo = slum[(size_t)(0.005 * (slum.size() - 1))];
    double hi = slum[(size_t)(0.995 * (slum.size() - 1))];

    auto r8 = asinhStretch(rn, lo, hi);
    auto g8 = asinhStretch(gn, lo, hi);
    auto b8 = asinhStretch(bn, lo, hi);

    FILE *fp = fopen(path.c_str(), "wb");
    if (!fp) { cerr << "Cannot write PNG: " << path << "\n"; return false; }

    png_structp png  = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    png_infop   info = png_create_info_struct(png);
    if (setjmp(png_jmpbuf(png))) { fclose(fp); return false; }

    png_init_io(png, fp);
    png_set_IHDR(png, info, (png_uint_32)w, (png_uint_32)h,
                 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    png_write_info(png, info);

    vector<uint8_t> row((size_t)w * 3);
    for (int y = 0; y < h; ++y)
    {
        for (int x = 0; x < w; ++x)
        {
            size_t i = (size_t)y * w + x;
            row[x * 3]     = r8[i];
            row[x * 3 + 1] = g8[i];
            row[x * 3 + 2] = b8[i];
        }
        png_write_row(png, row.data());
    }
    png_write_end(png, nullptr);
    png_destroy_write_struct(&png, &info);
    fclose(fp);
    return true;
}

// Donuts::Guider takes const double*; convert float buffer for registration.
static vector<double> toDouble(const vector<Pix> &v)
{
    vector<double> d(v.size());
    for (size_t i = 0; i < v.size(); ++i) d[i] = v[i];
    return d;
}

// ---------------------------------------------------------------------------
// Bayer helpers (RGGB, XBAYROFF=0 YBAYROFF=0)
// ---------------------------------------------------------------------------

static vector<Pix> extractGreen(const vector<Pix> &bayer, int bw, int bh, int &gw, int &gh)
{
    gw = bw / 2;
    gh = bh / 2;
    vector<Pix> green((size_t)gw * gh);
    for (int r = 0; r < gh; ++r)
        for (int c = 0; c < gw; ++c)
        {
            Pix g1 = bayer[(size_t)(2*r)     * bw + (2*c + 1)];
            Pix g2 = bayer[(size_t)(2*r + 1)  * bw + (2*c)];
            green[(size_t)r * gw + c] = (g1 + g2) * 0.5f;
        }
    return green;
}

static void debayerRGGB(const vector<Pix> &bayer, int bw, int bh,
                        vector<Pix> &r, vector<Pix> &g, vector<Pix> &b,
                        int &w, int &h)
{
    w = bw / 2;
    h = bh / 2;
    size_t npix = (size_t)w * h;
    r.resize(npix); g.resize(npix); b.resize(npix);
    for (int row = 0; row < h; ++row)
        for (int col = 0; col < w; ++col)
        {
            size_t i = (size_t)row * w + col;
            r[i] = bayer[(size_t)(2*row)     * bw + (2*col)];
            g[i] = (bayer[(size_t)(2*row)     * bw + (2*col + 1)] +
                    bayer[(size_t)(2*row + 1)  * bw + (2*col)]) * 0.5f;
            b[i] = bayer[(size_t)(2*row + 1)  * bw + (2*col + 1)];
        }
}

// ---------------------------------------------------------------------------
// Bicubic helpers
// ---------------------------------------------------------------------------

// Catmull-Rom weights for fractional position t in [0,1).
// w[0]=weight for ix-1, w[1]=ix+0, w[2]=ix+1, w[3]=ix+2.
static inline void cubicWeights4(float t, float w[4])
{
    float t2 = t * t, t3 = t2 * t;
    w[0] = -0.5f*t3 + t2 - 0.5f*t;
    w[1] =  1.5f*t3 - 2.5f*t2 + 1.0f;
    w[2] = -1.5f*t3 + 2.0f*t2 + 0.5f*t;
    w[3] =  0.5f*t3 - 0.5f*t2;
}

// Per-row x range where all 16 source pixels in the 4x4 kernel are in-bounds.
// Source mapping: sx(x) = dsx_dx*x + sx0,  sy(x) = dsy_dx*x + sy0  (linear in x).
// These are just the first column of the affine matrix (m.a, m.c) and the
// per-row offsets (m.b*y + m.tx, m.d*y + m.ty).
// Safe condition: ix in [1, w-3] and iy in [1, h-3].
// Returns [xlo, xhi] inclusive; xhi < xlo means the whole row needs clamping.
static pair<int,int> interiorXRange(double sx0, double sy0,
                                    double dsx_dx, double dsy_dx,
                                    int w, int h)
{
    double lo = 0.0, hi = (double)(w - 1);

    // sx constraint: dsx_dx*x + sx0 in [1, w-2)
    if (std::abs(dsx_dx) > 1e-12)
    {
        double a = (1.0   - sx0) / dsx_dx;
        double b = (w-2.0 - sx0) / dsx_dx - 1e-9;
        if (dsx_dx > 0) { lo = std::max(lo, a); hi = std::min(hi, b); }
        else            { lo = std::max(lo, b); hi = std::min(hi, a); }
    }
    else if (sx0 < 1.0 || sx0 >= w - 2.0) return {0, -1};

    // sy constraint: dsy_dx*x + sy0 in [1, h-2)
    if (std::abs(dsy_dx) > 1e-12)
    {
        double a = (1.0   - sy0) / dsy_dx;
        double b = (h-2.0 - sy0) / dsy_dx - 1e-9;
        if (dsy_dx > 0) { lo = std::max(lo, a); hi = std::min(hi, b); }
        else            { lo = std::max(lo, b); hi = std::min(hi, a); }
    }
    else if (sy0 < 1.0 || sy0 >= h - 2.0) return {0, -1};

    int ilo = std::max(0,   (int)std::ceil(lo));
    int ihi = std::min(w-1, (int)hi);
    return (ilo <= ihi) ? make_pair(ilo, ihi) : make_pair(0, -1);
}

// ---------------------------------------------------------------------------
// Fused 3-channel bicubic warp using a precomputed AffineMatrix.
// Processes dst rows [ylo, yhi) for parallel dispatch.
// Pass t.alignmentMatrix(w, h) to undo a measured DONUTS transform.
// ---------------------------------------------------------------------------

static void alignFrame3rows(
    const vector<Pix> &srcR, const vector<Pix> &srcG, const vector<Pix> &srcB,
    int w, int h,
    const Donuts::AffineMatrix &m,
    vector<Pix> &dstR, vector<Pix> &dstG, vector<Pix> &dstB,
    vector<uint8_t> &mask,
    int ylo, int yhi)
{
    for (int y = ylo; y < yhi; ++y)
    {
        // Per-row offsets: sx(x) = m.a*x + sx0,  sy(x) = m.c*x + sy0
        const double sx0 = m.b * y + m.tx;
        const double sy0 = m.d * y + m.ty;

        auto [xInLo, xInHi] = interiorXRange(sx0, sy0, m.a, m.c, w, h);

        for (int x = 0; x < w; ++x)
        {
            const double sx = m.a * x + sx0;
            const double sy = m.c * x + sy0;
            const size_t k  = (size_t)y * w + x;

            if (sx < 0.0 || sx >= w || sy < 0.0 || sy >= h)
            {
                mask[k] = 0;
                continue;
            }

            const int   ix = (int)sx, iy = (int)sy;
            const float fx = (float)(sx - ix), fy = (float)(sy - iy);
            float wx[4], wy[4];
            cubicWeights4(fx, wx);
            cubicWeights4(fy, wy);

            mask[k] = 1;

            if (x >= xInLo && x <= xInHi)
            {
                // Fast path: all 16 neighbors in-bounds, no clamping needed.
#ifdef __ARM_NEON
                const float32x4_t wxv = vld1q_f32(wx);
                float accR = 0.0f, accG = 0.0f, accB = 0.0f;
                for (int dr = 0; dr < 4; ++dr)
                {
                    const size_t base = (size_t)(iy + dr - 1) * w + (ix - 1);
                    const float  wyd  = wy[dr];
                    accR += wyd * vaddvq_f32(vmulq_f32(vld1q_f32(srcR.data() + base), wxv));
                    accG += wyd * vaddvq_f32(vmulq_f32(vld1q_f32(srcG.data() + base), wxv));
                    accB += wyd * vaddvq_f32(vmulq_f32(vld1q_f32(srcB.data() + base), wxv));
                }
                dstR[k] = accR;
                dstG[k] = accG;
                dstB[k] = accB;
#else
                float vR = 0.0f, vG = 0.0f, vB = 0.0f;
                for (int dr = 0; dr < 4; ++dr)
                {
                    const size_t base = (size_t)(iy + dr - 1) * w + (ix - 1);
                    const float  wyd  = wy[dr];
                    for (int dc = 0; dc < 4; ++dc)
                    {
                        const float  wdc = wx[dc] * wyd;
                        const size_t s   = base + dc;
                        vR += wdc * srcR[s];
                        vG += wdc * srcG[s];
                        vB += wdc * srcB[s];
                    }
                }
                dstR[k] = vR;
                dstG[k] = vG;
                dstB[k] = vB;
#endif
            }
            else
            {
                // Slow path: clamp each of the 16 source coordinates.
                float vR = 0.0f, vG = 0.0f, vB = 0.0f;
                for (int dr = -1; dr <= 2; ++dr)
                {
                    const int   py  = std::max(0, std::min(h - 1, iy + dr));
                    const float wyd = wy[dr + 1];
                    for (int dc = -1; dc <= 2; ++dc)
                    {
                        const int   px  = std::max(0, std::min(w - 1, ix + dc));
                        const float wdc = wx[dc + 1] * wyd;
                        const size_t s  = (size_t)py * w + px;
                        vR += wdc * srcR[s];
                        vG += wdc * srcG[s];
                        vB += wdc * srcB[s];
                    }
                }
                dstR[k] = vR;
                dstG[k] = vG;
                dstB[k] = vB;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Per-pixel combine helpers -- caller provides scratch[depth]; no malloc here.
// ---------------------------------------------------------------------------

// O(n) median via nth_element.
static Pix pixelMedian(const Pix *vals, Pix *scratch, int n)
{
    std::copy(vals, vals + n, scratch);
    std::nth_element(scratch, scratch + n / 2, scratch + n);
    return scratch[n / 2];
}

// Iterative median+MAD sigma-clip with in-place compaction.
//
// Uses median as the centre and 1.4826*MAD as the scale estimate.
// Mean+stddev fails on satellite/airplane trails: a single pixel at 950
// among background ~12 inflates mean to ~169 and sigma to ~383, so the
// outlier is inside the 3-sigma band (hi = 1317) and survives every pass.
// Median = 12.5 and MAD-sigma = 1.48 clips 950 immediately (hi = 16.9).
//
// devbuf must be caller-allocated with capacity >= n (same size as scratch).
static Pix pixelSigmaClip(const Pix *vals, Pix *scratch, Pix *devbuf,
                           int n, double kappa)
{
    std::copy(vals, vals + n, scratch);
    int cnt = n;
    for (int iter = 0; iter < 5 && cnt > 2; ++iter)
    {
        // Robust centre: median via nth_element (O(n) average).
        std::nth_element(scratch, scratch + cnt/2, scratch + cnt);
        double med = scratch[cnt/2];

        // Robust scale: MAD = median(|x - median|), sigma = 1.4826 * MAD.
        for (int i = 0; i < cnt; ++i) devbuf[i] = std::abs(scratch[i] - (Pix)med);
        std::nth_element(devbuf, devbuf + cnt/2, devbuf + cnt);
        double sigma = 1.4826 * devbuf[cnt/2];
        if (sigma < 1e-10) break;

        double lo = med - kappa * sigma, hi = med + kappa * sigma;
        int keep = 0;
        for (int i = 0; i < cnt; ++i)
            if (scratch[i] >= lo && scratch[i] <= hi)
                scratch[keep++] = scratch[i];
        if (keep == cnt) break;
        cnt = keep;
    }
    double sum = 0.0;
    for (int i = 0; i < cnt; ++i) sum += scratch[i];
    return cnt > 0 ? (Pix)(sum / cnt) : (Pix)vals[0];
}

// Reduce flat per-pixel stacks to final images.
// stk layout: stk[pixelIndex * depth + frameIndex].
// Each thread allocates one scratch buffer of size depth -- no per-pixel malloc.
static void reduceStack3(
    const vector<Pix> &stkR, const vector<Pix> &stkG, const vector<Pix> &stkB,
    const vector<int> &scount, int depth,
    StackMode mode, double kappa, int nThreads,
    vector<Pix> &outR, vector<Pix> &outG, vector<Pix> &outB)
{
    size_t npix = scount.size();
    outR.resize(npix); outG.resize(npix); outB.resize(npix);

    parallelFor((int)npix, nThreads, [&](int klo, int khi)
    {
        vector<Pix> scratch(depth), devbuf(depth);
        for (int k = klo; k < khi; ++k)
        {
            int n = scount[k];
            if (n == 0) { outR[k] = outG[k] = outB[k] = 0.0f; continue; }
            const Pix *vr = stkR.data() + (size_t)k * depth;
            const Pix *vg = stkG.data() + (size_t)k * depth;
            const Pix *vb = stkB.data() + (size_t)k * depth;
            if (mode == StackMode::Median)
            {
                outR[k] = pixelMedian   (vr, scratch.data(), n);
                outG[k] = pixelMedian   (vg, scratch.data(), n);
                outB[k] = pixelMedian   (vb, scratch.data(), n);
            }
            else
            {
                outR[k] = pixelSigmaClip(vr, scratch.data(), devbuf.data(), n, kappa);
                outG[k] = pixelSigmaClip(vg, scratch.data(), devbuf.data(), n, kappa);
                outB[k] = pixelSigmaClip(vb, scratch.data(), devbuf.data(), n, kappa);
            }
        }
    });
}

// ---------------------------------------------------------------------------
// Statistics helpers for per-frame MAD sigma-clip
// ---------------------------------------------------------------------------

static double medianOf(vector<double> v)
{
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

static double madSigmaOf(const vector<double> &v, double med)
{
    vector<double> dev(v.size());
    for (size_t i = 0; i < v.size(); ++i) dev[i] = std::abs(v[i] - med);
    return 1.4826 * medianOf(dev);
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        cerr << "Usage: " << argv[0]
             << " [-o output.fits] [-p output.png]"
                " [-m mean|median|sigmaclip] [-k kappa] [-j threads]"
                " frame1.fits ...\n";
        return 1;
    }

    string    outputPath = "stacked.fits";
    string    pngPath;
    StackMode mode     = StackMode::Mean;
    double    kappa    = 3.0;
    int       nThreads = (int)thread::hardware_concurrency();
    if (nThreads < 1) nThreads = 1;

    vector<string> inputs;

    for (int i = 1; i < argc; ++i)
    {
        string a = argv[i];
        if      (a == "-o" && i+1 < argc) outputPath = argv[++i];
        else if (a == "-p" && i+1 < argc) pngPath    = argv[++i];
        else if (a == "-k" && i+1 < argc) kappa      = stod(argv[++i]);
        else if (a == "-j" && i+1 < argc)
        {
            nThreads = stoi(argv[++i]);
            if (nThreads <= 0) nThreads = (int)thread::hardware_concurrency();
            if (nThreads <  1) nThreads = 1;
        }
        else if (a == "-m" && i+1 < argc)
        {
            string m = argv[++i];
            if      (m == "median")    mode = StackMode::Median;
            else if (m == "sigmaclip") mode = StackMode::SigmaClip;
            else if (m == "wstream")   mode = StackMode::WelfordStream;
            else if (m != "mean") { cerr << "Unknown mode: " << m << "\n"; return 1; }
        }
        else inputs.push_back(a);
    }
    if (inputs.empty()) { cerr << "No input files.\n"; return 1; }
    std::sort(inputs.begin(), inputs.end());

    const bool bayer = isBayerRGGB(inputs[0]);

    // ------------------------------------------------------------------
    // Load reference; set DONUTS reference on green channel or full frame.
    // ------------------------------------------------------------------

    int bw = 0, bh = 0;
    vector<Pix> refRaw = loadFITS(inputs[0], bw, bh);
    if (refRaw.empty()) return 1;

    const char *modeNames[] = { "mean", "median", "sigmaclip", "wstream" };
    cout << "Reference: " << inputs[0] << "  (" << bw << "x" << bh << ")";
    if (bayer) cout << "  [RGGB Bayer]";
    cout << "  mode=" << modeNames[(int)mode];
    if (mode == StackMode::SigmaClip || mode == StackMode::WelfordStream)
        cout << " kappa=" << kappa;
    cout << "  threads=" << nThreads << "\n";

    Donuts::Guider guider;
    if (bayer)
    {
        int gw = 0, gh = 0;
        auto refGreen = extractGreen(refRaw, bw, bh, gw, gh);
        auto refGreenD = toDouble(refGreen);
        guider.setReference(refGreenD.data(), gw, gh);
    }
    else { auto d = toDouble(refRaw); guider.setReference(d.data(), bw, bh); }

    // ------------------------------------------------------------------
    // Pass 1: measure transforms for all non-reference frames.
    // ------------------------------------------------------------------

    struct FrameInfo { string path; Donuts::Transform t; bool accepted=true; string reason; };
    vector<FrameInfo> frames;
    frames.reserve(inputs.size() - 1);

    for (size_t i = 1; i < inputs.size(); ++i)
    {
        FrameInfo fi;
        fi.path = inputs[i];
        int fw = 0, fh = 0;
        vector<Pix> raw = loadFITS(fi.path, fw, fh);
        if (raw.empty())          { fi.accepted=false; fi.reason="load failed";   frames.push_back(fi); continue; }
        if (fw != bw || fh != bh) { fi.accepted=false; fi.reason="size mismatch"; frames.push_back(fi); continue; }

        if (bayer)
        {
            int gw=0, gh=0;
            auto g  = extractGreen(raw, fw, fh, gw, gh);
            auto gd = toDouble(g);
            fi.t    = guider.measure(gd.data(), gw, gh);
        }
        else { auto d = toDouble(raw); fi.t = guider.measure(d.data(), fw, fh); }

        if      (!fi.t.valid())
            { fi.accepted=false; fi.reason="SNR < 3"; }
        else if (std::abs(fi.t.dtheta) > MAX_ROTATION_RAD)
            { fi.accepted=false; fi.reason="rotation bound"; }
        else if (std::abs(fi.t.dx) > MAX_TRANSLATION_PX || std::abs(fi.t.dy) > MAX_TRANSLATION_PX)
            { fi.accepted=false; fi.reason="translation bound"; }

        frames.push_back(fi);
    }

    // Pass 1b: MAD sigma-clip on accepted transforms.
    {
        vector<double> dxV, dyV, dtV;
        for (const auto &fi : frames)
            if (fi.accepted) { dxV.push_back(fi.t.dx); dyV.push_back(fi.t.dy); dtV.push_back(fi.t.dtheta); }
        if (dxV.size() >= 3)
        {
            double mDx=medianOf(dxV), sDx=std::max(madSigmaOf(dxV,mDx), 0.5);
            double mDy=medianOf(dyV), sDy=std::max(madSigmaOf(dyV,mDy), 0.5);
            double mDt=medianOf(dtV), sDt=std::max(madSigmaOf(dtV,mDt), 0.5*M_PI/180.0);
            for (auto &fi : frames)
                if (fi.accepted &&
                    (std::abs(fi.t.dx     - mDx) > MAD_NSIGMA * sDx ||
                     std::abs(fi.t.dy     - mDy) > MAD_NSIGMA * sDy ||
                     std::abs(fi.t.dtheta - mDt) > MAD_NSIGMA * sDt))
                    { fi.accepted=false; fi.reason="sigma-clip"; }
        }
    }

    for (size_t i = 0; i < frames.size(); ++i)
    {
        const auto &fi = frames[i];
        cout << "  [" << (i+1) << "] " << fi.path
             << "  dx=" << fi.t.dx << "  dy=" << fi.t.dy
             << "  dtheta=" << fi.t.dtheta * 180.0 / M_PI << " deg"
             << "  snr=" << fi.t.snr;
        cout << (fi.accepted ? "  -> OK\n" : ("  -> SKIP (" + fi.reason + ")\n"));
    }

    // ------------------------------------------------------------------
    // Allocate accumulators.
    // ------------------------------------------------------------------

    int dw = bayer ? bw/2 : bw;
    int dh = bayer ? bh/2 : bh;
    size_t npix = (size_t)dw * dh;

    int nAccepted = 1;
    for (const auto &fi : frames) if (fi.accepted) ++nAccepted;
    int depth = nAccepted;

    vector<Pix> accumR(npix, 0.0f), accumG(npix, 0.0f), accumB(npix, 0.0f);
    vector<int> coverage(npix, 0);

    vector<Pix> stkR, stkG, stkB;
    vector<int> scount;
    if (mode == StackMode::Median || mode == StackMode::SigmaClip)
    {
        stkR.assign(npix * depth, 0.0f);
        stkG.assign(npix * depth, 0.0f);
        stkB.assign(npix * depth, 0.0f);
        scount.assign(npix, 0);
        size_t mb = 3 * npix * depth * sizeof(Pix) / (1024*1024);
        cout << "Allocated " << mb << " MB for pixel stacks.\n";
    }

    // Welford streaming state: {mean, M2, count} per pixel -- O(pixels) memory
    // regardless of frame count. mean[] is the live stacked image at all times.
    vector<Pix>      wMeanR, wMeanG, wMeanB;
    vector<Pix>      wM2R,   wM2G,   wM2B;
    vector<uint32_t> wCount;
    if (mode == StackMode::WelfordStream)
    {
        wMeanR.assign(npix, 0.0f); wMeanG.assign(npix, 0.0f); wMeanB.assign(npix, 0.0f);
        wM2R.assign(npix, 0.0f);   wM2G.assign(npix, 0.0f);   wM2B.assign(npix, 0.0f);
        wCount.assign(npix, 0u);
        size_t mb = (6 * sizeof(Pix) + sizeof(uint32_t)) * npix / (1024*1024);
        cout << "Welford streaming: " << mb << " MB (constant, frame-count-independent).\n";
    }

    // Threaded accumulate: pixel ranges are disjoint across threads -- no locking.
    auto accumulateFrame = [&](const vector<Pix> &r, const vector<Pix> &g,
                               const vector<Pix> &b, const vector<uint8_t> *maskPtr)
    {
        parallelFor((int)npix, nThreads, [&](int klo, int khi)
        {
            for (int k = klo; k < khi; ++k)
            {
                if (maskPtr && !(*maskPtr)[k]) continue;
                if (mode == StackMode::Mean)
                {
                    accumR[k] += r[k]; accumG[k] += g[k]; accumB[k] += b[k];
                    coverage[k]++;
                }
                else
                {
                    int f = scount[k]++;
                    stkR[(size_t)k*depth + f] = r[k];
                    stkG[(size_t)k*depth + f] = g[k];
                    stkB[(size_t)k*depth + f] = b[k];
                }
            }
        });
    };

    // Welford online update for streaming sigma-clip.
    // bootstrap=true for the reference frame: always accept, no gate.
    // bootstrap=false for subsequent frames: gate with kappa*sigma before update.
    // Gate rejects a pixel if ANY channel exceeds kappa*sigma -- avoids color
    // contamination from satellite trails that hit only one channel (Bayer).
    // Requires n >= 2 to have a meaningful variance estimate; earlier samples
    // are always accepted to seed the statistics.
    auto welfordFrame = [&](const vector<Pix> &r, const vector<Pix> &g,
                            const vector<Pix> &b, const vector<uint8_t> *maskPtr,
                            bool bootstrap)
    {
        const float kf = (float)kappa;
        parallelFor((int)npix, nThreads, [&](int klo, int khi)
        {
            for (int k = klo; k < khi; ++k)
            {
                if (maskPtr && !(*maskPtr)[k]) continue;
                uint32_t n = wCount[k];

                // Rejection gate: skip pixel if any channel is an outlier.
                if (!bootstrap && n >= 2)
                {
                    const float nf1 = (float)(n - 1);
                    auto outlier = [&](float val, float mean, float m2) -> bool {
                        if (m2 <= 0.0f) return false;
                        float sig = std::sqrt(m2 / nf1);
                        return std::abs(val - mean) > kf * sig;
                    };
                    if (outlier(r[k], wMeanR[k], wM2R[k]) ||
                        outlier(g[k], wMeanG[k], wM2G[k]) ||
                        outlier(b[k], wMeanB[k], wM2B[k]))
                        continue;
                }

                // Welford recurrence (numerically stable online mean + M2).
                ++n;
                const float nf = (float)n;
                float d;
                d = r[k] - wMeanR[k]; wMeanR[k] += d / nf; wM2R[k] += d * (r[k] - wMeanR[k]);
                d = g[k] - wMeanG[k]; wMeanG[k] += d / nf; wM2G[k] += d * (g[k] - wMeanG[k]);
                d = b[k] - wMeanB[k]; wMeanB[k] += d / nf; wM2B[k] += d * (b[k] - wMeanB[k]);
                wCount[k] = n;
            }
        });
    };

    // ------------------------------------------------------------------
    // Add reference frame (no transform; all pixels valid).
    // ------------------------------------------------------------------

    {
        vector<Pix> r, g, b;
        if (bayer) { int tw=0,th=0; debayerRGGB(refRaw, bw, bh, r, g, b, tw, th); }
        else       { r = refRaw; g = refRaw; b = refRaw; }
        if (mode == StackMode::WelfordStream) welfordFrame(r, g, b, nullptr, true);
        else                                  accumulateFrame(r, g, b, nullptr);
    }

    // ------------------------------------------------------------------
    // Pass 2: load, align, accumulate each accepted frame.
    // ------------------------------------------------------------------

    int stackCount = 1;
    for (const auto &fi : frames)
    {
        if (!fi.accepted) continue;
        int fw=0, fh=0;
        vector<Pix> raw = loadFITS(fi.path, fw, fh);
        if (raw.empty()) continue;

        vector<Pix> r, g, b;
        if (bayer) { int tw=0,th=0; debayerRGGB(raw, fw, fh, r, g, b, tw, th); }
        else       { r = raw; g = raw; b = raw; }

        vector<Pix>     dstR(npix), dstG(npix), dstB(npix);
        vector<uint8_t> mask(npix, 0);

        const auto M = fi.t.alignmentMatrix(dw, dh);
        parallelFor(dh, nThreads, [&](int ylo, int yhi)
        {
            alignFrame3rows(r, g, b, dw, dh, M,
                            dstR, dstG, dstB, mask, ylo, yhi);
        });

        if (mode == StackMode::WelfordStream) welfordFrame(dstR, dstG, dstB, &mask, false);
        else                                  accumulateFrame(dstR, dstG, dstB, &mask);
        ++stackCount;
    }

    cout << "\nStacked " << stackCount << "/" << (int)inputs.size()
         << " frames -> " << outputPath << "\n";

    // ------------------------------------------------------------------
    // Reduce to final image.
    // ------------------------------------------------------------------

    vector<Pix> finalR(npix), finalG(npix), finalB(npix);

    if (mode == StackMode::WelfordStream)
    {
        // mean[] is already the final stacked image -- no separate reduce pass needed.
        finalR = wMeanR; finalG = wMeanG; finalB = wMeanB;
    }
    else if (mode == StackMode::Mean)
    {
        parallelFor((int)npix, nThreads, [&](int klo, int khi)
        {
            for (int k = klo; k < khi; ++k)
            {
                float inv = coverage[k] > 0 ? 1.0f / coverage[k] : 0.0f;
                finalR[k] = accumR[k] * inv;
                finalG[k] = accumG[k] * inv;
                finalB[k] = accumB[k] * inv;
            }
        });
    }
    else
    {
        reduceStack3(stkR, stkG, stkB, scount, depth, mode, kappa, nThreads,
                     finalR, finalG, finalB);
    }

    // ------------------------------------------------------------------
    // Write outputs.
    // ------------------------------------------------------------------

    if (bayer)
    {
        if (!writeFITS3(outputPath, finalR, finalG, finalB, dw, dh)) return 1;
        cout << "Output: " << dw << "x" << dh << " 3-plane FITS (R/G/B)\n";
        if (!pngPath.empty())
        {
            if (!writePNG(pngPath, finalR, finalG, finalB, dw, dh)) return 1;
            cout << "PNG:    " << pngPath << "\n";
        }
    }
    else
    {
        if (!writeFITS(outputPath, finalR, dw, dh)) return 1;
        if (!pngPath.empty())
        {
            if (!writePNG(pngPath, finalR, finalR, finalR, dw, dh)) return 1;
            cout << "PNG:    " << pngPath << "\n";
        }
    }

    return 0;
}
