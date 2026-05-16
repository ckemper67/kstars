/*
 * SPDX-FileCopyrightText: 2026 Christian Kemper <ckemper@gmail.com>
 *
 * SPDX-License-Identifier: GPL-2.0-or-later
 *
 * Standalone FITS stacker using the DONUTS rotation/translation solver.
 *
 * Monochrome: aligns and mean-combines directly.
 * Bayer (BAYERPAT header): extracts the green channel at half resolution for
 * DONUTS registration, then does a 2x2 RGGB debayer and applies the full
 * rigid-body transform (dx, dy, dtheta) to each color channel before stacking.
 * Output for Bayer input is a 3-plane FITS (R/G/B) at half the input resolution.
 * Only RGGB with XBAYROFF=0 YBAYROFF=0 is supported.
 *
 * Alignment uses Catmull-Rom bicubic interpolation with clamped boundaries.
 *
 * Frame rejection (two-pass):
 *   1. Hard bounds: |dtheta| <= 5 deg, |dx|/|dy| <= 150 half-res px.
 *   2. MAD sigma-clip (3-sigma) on dx, dy, dtheta of hard-bound survivors.
 *
 * Accumulation uses a per-pixel coverage count (mean mode) or per-pixel
 * value store (median and sigmaclip modes) so zero-padded border regions
 * from shifted frames don't dilute edge pixel values.
 *
 * PNG stretch: per-channel median subtraction, then linked luminance-based
 * asinh stretch (same tone curve applied to all three channels).
 *
 * Build (from repo root):
 *   c++ -std=c++17 -O2 \
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
 *               [-m mean|median|sigmaclip] [-k kappa]
 *               frame1.fits frame2.fits ...
 *
 *   Inputs are sorted alphabetically. The first file (after sorting) is the
 *   reference frame.  Default output: stacked.fits.  Default mode: mean.
 *   -k sets the sigma-clip rejection threshold (default 3.0).
 */

#include "donuts.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fitsio.h>
#include <iostream>
#include <numeric>
#include <png.h>
#include <string>
#include <vector>

using namespace std;

// ---------------------------------------------------------------------------
// Rejection thresholds
// ---------------------------------------------------------------------------

constexpr double MAX_ROTATION_RAD   = 5.0 * M_PI / 180.0;
constexpr double MAX_TRANSLATION_PX = 150.0;
constexpr double MAD_NSIGMA         = 3.0;

// ---------------------------------------------------------------------------
// Stacking mode
// ---------------------------------------------------------------------------

enum class StackMode { Mean, Median, SigmaClip };

// ---------------------------------------------------------------------------
// FITS I/O
// ---------------------------------------------------------------------------

static vector<double> loadFITS(const string &path, int &w, int &h)
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
    vector<double> buf((size_t)w * h);
    fits_read_img(fptr, TDOUBLE, 1, (long)w * h, nullptr, buf.data(), nullptr, &s);
    fits_close_file(fptr, &s);
    if (s) { cerr << "Read error on " << path << " (status " << s << ")\n"; return {}; }
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

static bool writeFITS(const string &path, const vector<double> &pixels, int w, int h)
{
    fitsfile *fptr = nullptr;
    int s = 0;
    fits_create_file(&fptr, ("!" + path).c_str(), &s);
    long naxes[2] = { (long)w, (long)h };
    fits_create_img(fptr, DOUBLE_IMG, 2, naxes, &s);
    fits_write_img(fptr, TDOUBLE, 1, (long)w * h,
                   const_cast<double *>(pixels.data()), &s);
    fits_close_file(fptr, &s);
    if (s) { cerr << "Write error on " << path << " (status " << s << ")\n"; return false; }
    return true;
}

static bool writeFITS3(const string &path,
                       const vector<double> &r, const vector<double> &g, const vector<double> &b,
                       int w, int h)
{
    fitsfile *fptr = nullptr;
    int s = 0;
    fits_create_file(&fptr, ("!" + path).c_str(), &s);
    long naxes[3] = { (long)w, (long)h, 3 };
    fits_create_img(fptr, DOUBLE_IMG, 3, naxes, &s);
    long npix = (long)w * h;
    fits_write_img(fptr, TDOUBLE,           1, npix, const_cast<double *>(r.data()), &s);
    fits_write_img(fptr, TDOUBLE,   npix + 1, npix, const_cast<double *>(g.data()), &s);
    fits_write_img(fptr, TDOUBLE, 2*npix + 1, npix, const_cast<double *>(b.data()), &s);
    fits_close_file(fptr, &s);
    if (s) { cerr << "Write error on " << path << " (status " << s << ")\n"; return false; }
    return true;
}

// ---------------------------------------------------------------------------
// PNG output: median background subtraction + linked luminance asinh stretch
// ---------------------------------------------------------------------------

static vector<double> subtractMedian(vector<double> ch)
{
    vector<double> s(ch);
    std::sort(s.begin(), s.end());
    double med = s[s.size() / 2];
    for (double &v : ch) v -= med;
    return ch;
}

static vector<uint8_t> asinhStretch(const vector<double> &ch, double lo, double hi)
{
    const double softening = 0.05 * (hi - lo);
    const double norm      = std::asinh((hi - lo) / softening);
    vector<uint8_t> out(ch.size());
    for (size_t i = 0; i < ch.size(); ++i)
    {
        double v = std::asinh((ch[i] - lo) / softening) / norm * 255.0;
        out[i] = (uint8_t)std::max(0.0, std::min(255.0, v));
    }
    return out;
}

static bool writePNG(const string &path,
                     const vector<double> &r, const vector<double> &g, const vector<double> &b,
                     int w, int h)
{
    auto rn = subtractMedian(r);
    auto gn = subtractMedian(g);
    auto bn = subtractMedian(b);

    vector<double> lum(rn.size());
    for (size_t i = 0; i < lum.size(); ++i)
        lum[i] = 0.299 * rn[i] + 0.587 * gn[i] + 0.114 * bn[i];

    vector<double> slum(lum);
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

// ---------------------------------------------------------------------------
// Bayer helpers (RGGB, XBAYROFF=0 YBAYROFF=0)
// ---------------------------------------------------------------------------

static vector<double> extractGreen(const vector<double> &bayer, int bw, int bh, int &gw, int &gh)
{
    gw = bw / 2;
    gh = bh / 2;
    vector<double> green((size_t)gw * gh);
    for (int r = 0; r < gh; ++r)
        for (int c = 0; c < gw; ++c)
        {
            double g1 = bayer[(size_t)(2*r)     * bw + (2*c + 1)];
            double g2 = bayer[(size_t)(2*r + 1)  * bw + (2*c)];
            green[(size_t)r * gw + c] = (g1 + g2) * 0.5;
        }
    return green;
}

static void debayerRGGB(const vector<double> &bayer, int bw, int bh,
                        vector<double> &r, vector<double> &g, vector<double> &b,
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
                    bayer[(size_t)(2*row + 1)  * bw + (2*col)]) * 0.5;
            b[i] = bayer[(size_t)(2*row + 1)  * bw + (2*col + 1)];
        }
}

// ---------------------------------------------------------------------------
// Bicubic (Catmull-Rom) rigid-body alignment, clamped boundaries.
// mask (if provided) is set true for pixels with a source sample in-frame.
// dtheta in radians; to undo a measured transform pass (-dx, -dy, -dtheta).
// ---------------------------------------------------------------------------

static double cubicWeight(double x)
{
    x = std::abs(x);
    if (x < 1.0) return  1.5*x*x*x - 2.5*x*x + 1.0;
    if (x < 2.0) return -0.5*x*x*x + 2.5*x*x - 4.0*x + 2.0;
    return 0.0;
}

static vector<double> alignFrame(
    const vector<double> &src, int w, int h,
    double dx, double dy, double dtheta,
    vector<bool> *mask = nullptr)
{
    if (mask) mask->assign((size_t)w * h, false);
    vector<double> dst((size_t)w * h, 0.0);
    const double co = cos(dtheta), si = sin(dtheta);
    const double mx = w / 2.0, my = h / 2.0;

    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x)
        {
            double x1 = x - dx - mx, y1 = y - dy - my;
            double sx  =  x1 * co + y1 * si + mx;
            double sy  = -x1 * si + y1 * co + my;

            if (sx < 0.0 || sx >= w || sy < 0.0 || sy >= h) continue;

            int    ix = (int)sx, iy = (int)sy;
            double fx = sx - ix,  fy = sy - iy;
            double val = 0.0;

            for (int dr = -1; dr <= 2; ++dr)
            {
                double wy = cubicWeight(fy - dr);
                int py = std::max(0, std::min(h - 1, iy + dr));
                for (int dc = -1; dc <= 2; ++dc)
                {
                    double wx = cubicWeight(fx - dc);
                    int px = std::max(0, std::min(w - 1, ix + dc));
                    val += wx * wy * src[(size_t)py * w + px];
                }
            }

            dst[(size_t)y * w + x] = val;
            if (mask) (*mask)[(size_t)y * w + x] = true;
        }
    return dst;
}

// ---------------------------------------------------------------------------
// Per-pixel combine helpers
// ---------------------------------------------------------------------------

static double pixelMedian(double *vals, int n)
{
    vector<double> tmp(vals, vals + n);
    std::sort(tmp.begin(), tmp.end());
    return tmp[n / 2];
}

static double pixelSigmaClip(double *vals, int n, double kappa)
{
    vector<bool> keep(n, true);
    for (int iter = 0; iter < 5; ++iter)
    {
        double sum = 0.0; int cnt = 0;
        for (int i = 0; i < n; ++i) if (keep[i]) { sum += vals[i]; ++cnt; }
        if (cnt == 0) break;
        double mean = sum / cnt;
        double var  = 0.0;
        for (int i = 0; i < n; ++i) if (keep[i]) var += (vals[i] - mean) * (vals[i] - mean);
        double sigma = std::sqrt(var / cnt);
        bool changed = false;
        for (int i = 0; i < n; ++i)
            if (keep[i] && std::abs(vals[i] - mean) > kappa * sigma)
                { keep[i] = false; changed = true; }
        if (!changed) break;
    }
    double sum = 0.0; int cnt = 0;
    for (int i = 0; i < n; ++i) if (keep[i]) { sum += vals[i]; ++cnt; }
    return cnt > 0 ? sum / cnt : vals[0];
}

// Reduce per-pixel stacks to final image according to mode.
// stk layout: stk[pixelIndex * depth + frameIndex], count[pixelIndex] = valid frames.
static vector<double> reduceStack(const vector<double> &stk, const vector<int> &count,
                                  int depth, StackMode mode, double kappa)
{
    size_t npix = count.size();
    vector<double> out(npix, 0.0);
    for (size_t k = 0; k < npix; ++k)
    {
        int n = count[k];
        if (n == 0) continue;
        double *vals = const_cast<double *>(&stk[k * depth]);
        if      (mode == StackMode::Median)    out[k] = pixelMedian   (vals, n);
        else if (mode == StackMode::SigmaClip) out[k] = pixelSigmaClip(vals, n, kappa);
        else { double s = 0; for (int i = 0; i < n; ++i) s += vals[i]; out[k] = s / n; }
    }
    return out;
}

// ---------------------------------------------------------------------------
// Statistics helpers for frame sigma-clipping
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
                " [-m mean|median|sigmaclip] [-k kappa]"
                " frame1.fits ...\n";
        return 1;
    }

    string    outputPath = "stacked.fits";
    string    pngPath;
    StackMode mode  = StackMode::Mean;
    double    kappa = 3.0;
    vector<string> inputs;

    for (int i = 1; i < argc; ++i)
    {
        string a = argv[i];
        if      (a == "-o" && i+1 < argc) outputPath = argv[++i];
        else if (a == "-p" && i+1 < argc) pngPath    = argv[++i];
        else if (a == "-k" && i+1 < argc) kappa      = stod(argv[++i]);
        else if (a == "-m" && i+1 < argc)
        {
            string m = argv[++i];
            if      (m == "median")    mode = StackMode::Median;
            else if (m == "sigmaclip") mode = StackMode::SigmaClip;
            else if (m != "mean") { cerr << "Unknown mode: " << m << "\n"; return 1; }
        }
        else inputs.push_back(a);
    }
    if (inputs.empty()) { cerr << "No input files.\n"; return 1; }
    std::sort(inputs.begin(), inputs.end());

    const bool bayer = isBayerRGGB(inputs[0]);

    // ------------------------------------------------------------------
    // Load reference and set DONUTS reference.
    // ------------------------------------------------------------------
    int bw = 0, bh = 0;
    vector<double> refRaw = loadFITS(inputs[0], bw, bh);
    if (refRaw.empty()) return 1;

    cout << "Reference: " << inputs[0] << "  (" << bw << "x" << bh << ")";
    if (bayer) cout << "  [RGGB Bayer]";
    const char *modeName[] = { "mean", "median", "sigmaclip" };
    cout << "  mode=" << modeName[(int)mode];
    if (mode == StackMode::SigmaClip) cout << " kappa=" << kappa;
    cout << "\n";

    Donuts::Guider guider;
    if (bayer)
    {
        int gw = 0, gh = 0;
        auto refGreen = extractGreen(refRaw, bw, bh, gw, gh);
        guider.setReference(refGreen.data(), gw, gh);
    }
    else guider.setReference(refRaw.data(), bw, bh);

    // ------------------------------------------------------------------
    // Pass 1: measure transforms for all non-reference frames.
    // ------------------------------------------------------------------
    struct FrameInfo { string path; Donuts::Transform t; bool accepted = true; string reason; };
    vector<FrameInfo> frames;
    frames.reserve(inputs.size() - 1);

    for (size_t i = 1; i < inputs.size(); ++i)
    {
        FrameInfo fi;
        fi.path = inputs[i];
        int fw = 0, fh = 0;
        vector<double> raw = loadFITS(fi.path, fw, fh);
        if (raw.empty())            { fi.accepted = false; fi.reason = "load failed";    frames.push_back(fi); continue; }
        if (fw != bw || fh != bh)   { fi.accepted = false; fi.reason = "size mismatch";  frames.push_back(fi); continue; }

        if (bayer) { int gw=0,gh=0; auto g=extractGreen(raw,fw,fh,gw,gh); fi.t=guider.measure(g.data(),gw,gh); }
        else       { fi.t = guider.measure(raw.data(), fw, fh); }

        if      (!fi.t.valid())                                                      { fi.accepted=false; fi.reason="SNR < 3"; }
        else if (std::abs(fi.t.dtheta) > MAX_ROTATION_RAD)                          { fi.accepted=false; fi.reason="rotation bound"; }
        else if (std::abs(fi.t.dx) > MAX_TRANSLATION_PX || std::abs(fi.t.dy) > MAX_TRANSLATION_PX) { fi.accepted=false; fi.reason="translation bound"; }
        frames.push_back(fi);
    }

    // Pass 1b: MAD sigma-clip on accepted transforms.
    {
        vector<double> dxV, dyV, dtV;
        for (const auto &fi : frames) if (fi.accepted) { dxV.push_back(fi.t.dx); dyV.push_back(fi.t.dy); dtV.push_back(fi.t.dtheta); }
        if (dxV.size() >= 3)
        {
            double mDx=medianOf(dxV), sDx=std::max(madSigmaOf(dxV,mDx), 0.5);
            double mDy=medianOf(dyV), sDy=std::max(madSigmaOf(dyV,mDy), 0.5);
            double mDt=medianOf(dtV), sDt=std::max(madSigmaOf(dtV,mDt), 0.5*M_PI/180.0);
            for (auto &fi : frames)
                if (fi.accepted &&
                    (std::abs(fi.t.dx-mDx)>MAD_NSIGMA*sDx ||
                     std::abs(fi.t.dy-mDy)>MAD_NSIGMA*sDy ||
                     std::abs(fi.t.dtheta-mDt)>MAD_NSIGMA*sDt))
                    { fi.accepted=false; fi.reason="sigma-clip"; }
        }
    }

    // Print summary.
    for (size_t i = 0; i < frames.size(); ++i)
    {
        const auto &fi = frames[i];
        cout << "  [" << (i+1) << "] " << fi.path
             << "  dx=" << fi.t.dx << "  dy=" << fi.t.dy
             << "  dtheta=" << fi.t.dtheta*180.0/M_PI << " deg"
             << "  snr=" << fi.t.snr;
        cout << (fi.accepted ? "  -> OK\n" : ("  -> SKIP (" + fi.reason + ")\n"));
    }

    // ------------------------------------------------------------------
    // Determine accepted count; allocate accumulators.
    // ------------------------------------------------------------------
    int dw = bayer ? bw/2 : bw;
    int dh = bayer ? bh/2 : bh;
    size_t npix = (size_t)dw * dh;

    int nAccepted = 1; // reference always counts
    for (const auto &fi : frames) if (fi.accepted) ++nAccepted;
    int depth = nAccepted; // max frames per pixel in stack

    // Mean uses running sum + coverage. Median/SigmaClip stores all values.
    vector<double> accumR(npix,0.0), accumG(npix,0.0), accumB(npix,0.0);
    vector<int>    coverage(npix, 0);

    vector<double> stackR, stackG, stackB; // only for median/sigmaclip
    vector<int>    scount;                 // per-pixel valid frame count
    if (mode != StackMode::Mean)
    {
        stackR.assign(npix * depth, 0.0);
        stackG.assign(npix * depth, 0.0);
        stackB.assign(npix * depth, 0.0);
        scount.assign(npix, 0);
        cout << "Allocated " << (3*npix*depth*8/1024/1024) << " MB for pixel stacks.\n";
    }

    // Lambda to add a debayered/mono frame into the appropriate accumulator.
    auto accumulate = [&](const vector<double> &r, const vector<double> &g,
                          const vector<double> &b, const vector<bool> *mask)
    {
        for (size_t k = 0; k < npix; ++k)
        {
            if (mask && !(*mask)[k]) continue;
            if (mode == StackMode::Mean)
            {
                accumR[k] += r[k]; accumG[k] += g[k]; accumB[k] += b[k];
                coverage[k]++;
            }
            else
            {
                int f = scount[k]++;
                stackR[k*depth + f] = r[k];
                stackG[k*depth + f] = g[k];
                stackB[k*depth + f] = b[k];
            }
        }
    };

    // ------------------------------------------------------------------
    // Add reference (no transform; all pixels valid).
    // ------------------------------------------------------------------
    {
        vector<double> r, g, b;
        if (bayer) { int tw=0,th=0; debayerRGGB(refRaw,bw,bh,r,g,b,tw,th); }
        else       { r = refRaw; g = refRaw; b = refRaw; }
        accumulate(r, g, b, nullptr); // null mask = all valid
    }

    // ------------------------------------------------------------------
    // Pass 2: load, debayer, align, accumulate accepted frames.
    // ------------------------------------------------------------------
    int stackCount = 1;
    for (const auto &fi : frames)
    {
        if (!fi.accepted) continue;
        int fw=0,fh=0;
        vector<double> raw = loadFITS(fi.path, fw, fh);
        if (raw.empty()) continue;

        vector<double> r, g, b;
        if (bayer) { int tw=0,th=0; debayerRGGB(raw,fw,fh,r,g,b,tw,th); }
        else       { r = raw; g = raw; b = raw; }

        vector<bool> mask;
        auto ar = alignFrame(r, dw, dh, -fi.t.dx, -fi.t.dy, -fi.t.dtheta, &mask);
        auto ag = alignFrame(g, dw, dh, -fi.t.dx, -fi.t.dy, -fi.t.dtheta);
        auto ab = alignFrame(b, dw, dh, -fi.t.dx, -fi.t.dy, -fi.t.dtheta);

        accumulate(ar, ag, ab, &mask);
        ++stackCount;
    }

    cout << "\nStacked " << stackCount << "/" << (int)inputs.size() << " frames -> " << outputPath << "\n";

    // ------------------------------------------------------------------
    // Reduce to final image.
    // ------------------------------------------------------------------
    vector<double> finalR(npix), finalG(npix), finalB(npix);
    if (mode == StackMode::Mean)
    {
        for (size_t k = 0; k < npix; ++k)
        {
            double inv = coverage[k] > 0 ? 1.0 / coverage[k] : 0.0;
            finalR[k] = accumR[k] * inv;
            finalG[k] = accumG[k] * inv;
            finalB[k] = accumB[k] * inv;
        }
    }
    else
    {
        finalR = reduceStack(stackR, scount, depth, mode, kappa);
        finalG = reduceStack(stackG, scount, depth, mode, kappa);
        finalB = reduceStack(stackB, scount, depth, mode, kappa);
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
        if (!writeFITS(outputPath, finalR, bw, bh)) return 1;
        if (!pngPath.empty())
        {
            if (!writePNG(pngPath, finalR, finalG, finalB, dw, dh)) return 1;
            cout << "PNG:    " << pngPath << "\n";
        }
    }

    return 0;
}
