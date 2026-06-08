/*
    SPDX-FileCopyrightText: 2025 Christian Kemper <ckemper@gmail.com>

    Adapted from INDI SkyRenderer by Jasem Mutlaq / Christian Kemper.

    SPDX-License-Identifier: GPL-2.0-or-later
*/

#include "starfieldrenderer.h"
#include "bright_stars_catalog.h"

#include <fitsio.h>

#include <QCoreApplication>
#include <QDir>
#include <QStandardPaths>

#include <algorithm>
#include <clocale>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static double range360(double d)
{
    while (d < 0.0)   d += 360.0;
    while (d >= 360.0) d -= 360.0;
    return d;
}

static double rangeDec(double d)
{
    return std::max(-90.0, std::min(90.0, d));
}

// ---------------------------------------------------------------------------
// Frame buffer (replaces INDI::CCDChip)
// ---------------------------------------------------------------------------

struct FrameBuffer
{
    int width  = 0;
    int height = 0;
    std::vector<uint16_t> pixels;

    void resize(int w, int h)
    {
        width  = w;
        height = h;
        pixels.assign(static_cast<size_t>(w) * h, 0);
    }
};

// ---------------------------------------------------------------------------
// Rendering helpers adapted from INDI SkyRenderer
// ---------------------------------------------------------------------------

static double flux(float mag, float limitingMag, float saturationMag, int maxVal)
{
    if (limitingMag == saturationMag)
        return 1.0;
    double const k = 2.5 * log10(maxVal) / (limitingMag - saturationMag);
    return pow(10.0, (limitingMag - mag) * k / 2.5);
}

static int addToPixel(FrameBuffer &fb, int x, int y, int val)
{
    if (x < 0 || x >= fb.width || y < 0 || y >= fb.height)
        return 0;
    int idx = y * fb.width + x;
    int newval = static_cast<int>(fb.pixels[idx]) + val;
    if (newval > 65535)
        newval = 65535;
    fb.pixels[idx] = static_cast<uint16_t>(newval);
    return 1;
}

static int drawStar(FrameBuffer &fb, float mag, float x, float y,
                    float exp_s, float imageScaleX, float imageScaleY,
                    const RenderParams &cfg)
{
    if (imageScaleX <= 0.0f || imageScaleY <= 0.0f)
        return 0;

    static constexpr float kMaxStarInfluencePx = 100.0f;
    if (x < -kMaxStarInfluencePx || x > fb.width + kMaxStarInfluencePx ||
            y < -kMaxStarInfluencePx || y > fb.height + kMaxStarInfluencePx)
        return 0;

    float const totalFlux = static_cast<float>(
                                flux(mag, cfg.limitingMag, cfg.saturationMag, cfg.maxVal) * exp_s);

    float const pixel_area = imageScaleX * imageScaleY;
    if (pixel_area <= 0.0f)
        return 0;

    float const beta     = 2.5f;
    float const betaTerm = 4.0f * (std::pow(2.0f, 1.0f / beta) - 1.0f);

    float fwhm2 = cfg.seeingArcsec * cfg.seeingArcsec;

    float const minFwhm  = std::min(imageScaleX, imageScaleY);
    float const minFwhm2 = minFwhm * minFwhm;
    if (fwhm2 < minFwhm2)
        fwhm2 = minFwhm2;

    float const alpha2     = fwhm2 / betaTerm;
    float const moffatNorm = (beta - 1.0f) / (static_cast<float>(M_PI) * alpha2);

    float const threshold = std::max(moffatNorm * totalFlux * pixel_area * 2.0f, 1.0f);
    float const r2_max    = alpha2 * (std::pow(threshold, 1.0f / beta) - 1.0f);
    int   const moffatBox = (r2_max > 0.0f)
                            ? static_cast<int>(std::sqrt(r2_max) / std::min(imageScaleX, imageScaleY)) + 1
                            : static_cast<int>(3.0f * cfg.seeingArcsec / std::min(imageScaleX, imageScaleY)) + 1;
    static constexpr int maxRenderBox = 40;
    int const boxsize = std::min(moffatBox, maxRenderBox);

    float const fracX = x - std::floor(x);
    float const fracY = y - std::floor(y);
    int   const ix    = static_cast<int>(std::floor(x));
    int   const iy    = static_cast<int>(std::floor(y));

    int drew = 0;
    int const box2 = boxsize * boxsize;
    for (int sy = -boxsize; sy <= boxsize; sy++)
    {
        for (int sx = -boxsize; sx <= boxsize; sx++)
        {
            if (sx * sx + sy * sy > box2)
                continue;
            float const dx  = (sx - fracX) * imageScaleX;
            float const dy  = (sy - fracY) * imageScaleY;
            float const dc2 = dx * dx + dy * dy;
            float const moffat = moffatNorm * std::pow(1.0f + dc2 / alpha2, -beta);
            int   const fp     = static_cast<int>(moffat * totalFlux * pixel_area);
            if (fp > 0)
            {
                if (addToPixel(fb, ix + sx, iy + sy, fp) != 0)
                    drew = 1;
            }
        }
    }
    return drew;
}

static void drawSkyGlow(FrameBuffer &fb, float exp_s, float imageScaleX,
                        float imageScaleY, const RenderParams &cfg)
{
    float const skyflux = static_cast<float>(
                              flux(cfg.skyGlowMag, cfg.limitingMag, cfg.saturationMag, cfg.maxVal)) * exp_s;

    float const vig = std::min(fb.width, fb.height) * imageScaleX;
    float const invVig2 = 1.0f / (vig * vig);

    for (int y = 0; y < fb.height; y++)
    {
        float const sy = fb.height / 2.0f - y;
        for (int x = 0; x < fb.width; x++)
        {
            float const sx = fb.width / 2.0f - x;
            float const dc2 = sx * sx * imageScaleX * imageScaleX
                              + sy * sy * imageScaleY * imageScaleY;
            float const fa = std::exp(-2.0f * 0.7f * dc2 * invVig2);

            int idx = y * fb.width + x;
            float fp = (fb.pixels[idx] + skyflux) * fa;
            if (fp > cfg.maxVal) fp = static_cast<float>(cfg.maxVal);
            if (fp < fb.pixels[idx]) fp = fb.pixels[idx];
            fb.pixels[idx] = static_cast<uint16_t>(fp);
        }
    }
}

static void applyReadNoise(FrameBuffer &fb, int bias, int readNoise, int maxVal)
{
    if (readNoise <= 0 && bias <= 0)
        return;
    for (size_t i = 0; i < fb.pixels.size(); i++)
    {
        int newval = static_cast<int>(fb.pixels[i]) + bias;
        if (readNoise > 0)
            newval += rand() % readNoise;
        if (newval > maxVal)
            newval = maxVal;
        if (newval < 0)
            newval = 0;
        fb.pixels[i] = static_cast<uint16_t>(newval);
    }
}

// ---------------------------------------------------------------------------
// FITS writing with WCS headers
// ---------------------------------------------------------------------------

static bool writeFITS(const FrameBuffer &fb, const RenderParams &params,
                      double pixscale, const QString &outputPath)
{
    fitsfile *fptr = nullptr;
    int status = 0;
    long naxes[2] = { params.width, params.height };

    QString bangPath = "!" + outputPath;
    fits_create_file(&fptr, bangPath.toLocal8Bit().constData(), &status);
    if (status)
        return false;

    fits_create_img(fptr, USHORT_IMG, 2, naxes, &status);

    char ctype1[] = "RA---TAN";
    char ctype2[] = "DEC--TAN";
    char radecsys[] = "FK5";
    double crval1  = params.ra;
    double crval2  = params.dec;
    // The renderer treats integer C-array indices as pixel centers, so a star
    // at the optical axis lands at C-coord (width/2, height/2) which is the
    // center of FITS pixel (width/2 + 1, height/2 + 1) in 1-indexed convention.
    double crpix1  = params.width  / 2.0 + 1.0;
    double crpix2  = params.height / 2.0 + 1.0;
    double degpix  = pixscale / 3600.0;
    double cdelt1  = -degpix;
    double cdelt2  =  degpix;
    double rotRad  = params.rotationDeg * (M_PI / 180.0);
    double cosRot  = std::cos(rotRad);
    double sinRot  = std::sin(rotRad);
    // Standard Calabretta & Greisen CDELT/CROTA2 -> CD matrix relation:
    double cd1_1   =  cdelt1 * cosRot;
    double cd1_2   = -cdelt2 * sinRot;
    double cd2_1   =  cdelt1 * sinRot;
    double cd2_2   =  cdelt2 * cosRot;
    double equinox = 2000.0;
    double focal   = params.focalLengthMM;
    double pxsz    = params.pixelSizeUM;

    fits_write_key(fptr, TSTRING, "CTYPE1",   ctype1,   nullptr, &status);
    fits_write_key(fptr, TSTRING, "CTYPE2",   ctype2,   nullptr, &status);
    fits_write_key(fptr, TSTRING, "RADECSYS", radecsys, nullptr, &status);
    fits_write_key(fptr, TDOUBLE, "EQUINOX",  &equinox, nullptr, &status);
    fits_write_key(fptr, TDOUBLE, "CRVAL1",   &crval1,  nullptr, &status);
    fits_write_key(fptr, TDOUBLE, "CRVAL2",   &crval2,  nullptr, &status);
    fits_write_key(fptr, TDOUBLE, "CRPIX1",   &crpix1,  nullptr, &status);
    fits_write_key(fptr, TDOUBLE, "CRPIX2",   &crpix2,  nullptr, &status);
    fits_write_key(fptr, TDOUBLE, "CD1_1",    &cd1_1,   nullptr, &status);
    fits_write_key(fptr, TDOUBLE, "CD1_2",    &cd1_2,   nullptr, &status);
    fits_write_key(fptr, TDOUBLE, "CD2_1",    &cd2_1,   nullptr, &status);
    fits_write_key(fptr, TDOUBLE, "CD2_2",    &cd2_2,   nullptr, &status);
    fits_write_key(fptr, TDOUBLE, "FOCALLEN", &focal,   nullptr, &status);
    fits_write_key(fptr, TDOUBLE, "XPIXSZ",   &pxsz,    nullptr, &status);
    fits_write_key(fptr, TDOUBLE, "YPIXSZ",   &pxsz,    nullptr, &status);

    long fpixel = 1;
    fits_write_img(fptr, TUSHORT, fpixel,
                   static_cast<long>(params.width) * params.height,
                   const_cast<uint16_t *>(fb.pixels.data()), &status);

    fits_close_file(fptr, &status);
    return (status == 0);
}

// ---------------------------------------------------------------------------
// Main entry point
// ---------------------------------------------------------------------------

RenderedField renderStarField(const RenderParams &params, const QString &outputPath)
{
    RenderedField result;
    result.ra       = params.ra;
    result.dec      = params.dec;
    result.rotation = params.rotationDeg;

    // 206264.806 = arcsec per radian; divide by 1000 because focal is mm and
    // pixel is microns -> pixel/focal is already dimensionless * 1e-3.
    float const imageScaleX = static_cast<float>((params.pixelSizeUM / params.focalLengthMM) * 206.2648);
    float const imageScaleY = imageScaleX;
    result.pixscale = imageScaleX;

    double const theta = params.rotationDeg * (M_PI / 180.0);
    double const pprx  = params.focalLengthMM / params.pixelSizeUM * 1000.0;
    double const ppry  = pprx;

    double const pa =  pprx * std::cos(theta);
    double const pb =  ppry * std::sin(theta);
    double const pd = -pprx * std::sin(theta);
    double const pe =  ppry * std::cos(theta);
    double const pc =  params.width / 2.0;
    double const pf =  params.height / 2.0;
    double const ccdW = params.width;

    double const rar  = params.ra  * (M_PI / 180.0);
    double const decr = params.dec * (M_PI / 180.0);

    double radius = std::sqrt(
                        imageScaleX * imageScaleX * params.width / 2.0 * params.width / 2.0 +
                        imageScaleY * imageScaleY * params.height / 2.0 * params.height / 2.0) / 60.0;

    double lookuplimit = params.limitingMag;
    if (radius > 60)
        lookuplimit = 11;

    FrameBuffer fb;
    fb.resize(params.width, params.height);

    int drawn = 0;

    QString gscBin;
    {
        QStringList gscCandidates;
        QString kdeRoot = QString::fromLocal8Bit(qgetenv("KDEROOT"));
        if (!kdeRoot.isEmpty())
        {
            gscCandidates << kdeRoot + "/bin/KStars.app/Contents/MacOS/gsc";
            gscCandidates << kdeRoot + "/Applications/KDE/kstars.app/Contents/MacOS/gsc";
        }
        gscCandidates << QCoreApplication::applicationDirPath() + "/gsc";
        for (const auto &appDir : QStandardPaths::standardLocations(QStandardPaths::ApplicationsLocation))
        {
            gscCandidates << appDir + "/KDE/KStars.app/Contents/MacOS/gsc";
            gscCandidates << appDir + "/KStars.app/Contents/MacOS/gsc";
        }

        for (const auto &candidate : gscCandidates)
        {
            if (!candidate.isEmpty() && QFile::exists(candidate))
            {
                gscBin = candidate;
                break;
            }
        }
    }
    if (gscBin.isEmpty())
    {
        result.starsDrawn = -1;
        return result;
    }

    if (qEnvironmentVariableIsEmpty("GSCDAT"))
    {
        const QStringList candidates =
        {
            QDir::homePath() + "/Library/Application Support/kstars/gsc",
            QDir::homePath() + "/.local/share/kstars/gsc",
        };
        for (const auto &d : candidates)
            if (QDir(d).exists())
            {
                qputenv("GSCDAT", d.toLocal8Bit());
                break;
            }
    }

    {
        char *savedLocale = setlocale(LC_NUMERIC, nullptr);
        std::string localeBackup(savedLocale ? savedLocale : "");
        setlocale(LC_NUMERIC, "C");

        char gsccmd[256];
        snprintf(gsccmd, sizeof(gsccmd),
                 "%s -c %8.6f %+8.6f -r %4.1f -m 0 %4.2f -n 3000",
                 gscBin.toLocal8Bit().constData(),
                 range360(params.ra),
                 rangeDec(params.dec),
                 radius,
                 lookuplimit);

        FILE *pp = popen(gsccmd, "r");
        if (pp != nullptr)
        {
            char line[256];
            while (fgets(line, sizeof(line), pp) != nullptr)
            {
                char  id[20], plate[6], ob[6];
                float ra, dec, pose, mag, mage, dist;
                int   band, c, dir;

                int rc = sscanf(line, "%10s %f %f %f %f %f %d %d %4s %2s %f %d",
                                id, &ra, &dec, &pose, &mag, &mage,
                                &band, &c, plate, ob, &dist, &dir);
                if (rc != 12)
                    continue;

                double const srar  = ra  * (M_PI / 180.0);
                double const sdecr = dec * (M_PI / 180.0);

                double const denom = cos(decr) * cos(sdecr) * cos(srar - rar)
                                     + sin(decr) * sin(sdecr);
                if (denom <= 0)
                    continue;

                // Standard tangent-plane projection:
                //   xi  = cos(d)*sin(a-a0) / denom            (east-positive)
                //   eta = (sin(d)*cos(d0) - cos(d)*sin(d0)*cos(a-a0)) / denom
                //         (north-positive)
                double const sx = cos(sdecr) * sin(srar - rar) / denom;
                double const sy = (sin(sdecr) * cos(decr)
                                   - cos(sdecr) * sin(decr) * cos(srar - rar)) / denom;

                double const ccdx = ccdW - (pa * sx + pb * sy + pc);
                double const ccdy =          pd * sx + pe * sy + pf;

                drawn += drawStar(fb, mag,
                                  static_cast<float>(ccdx),
                                  static_cast<float>(ccdy),
                                  params.exposureSecs,
                                  imageScaleX, imageScaleY, params);
            }
            pclose(pp);
        }
        else
        {
            result.starsDrawn = -1;
            setlocale(LC_NUMERIC, localeBackup.c_str());
            return result;
        }

        setlocale(LC_NUMERIC, localeBackup.c_str());
    }

    for (int bsi = 0; bsi < s_BrightStarsCount; ++bsi)
    {
        float const bsRaDeg  = s_BrightStars[bsi].ra;
        float const bsDecDeg = s_BrightStars[bsi].dec;
        float const bsMag    = s_BrightStars[bsi].mag;

        double const srar  = bsRaDeg  * (M_PI / 180.0);
        double const sdecr = bsDecDeg * (M_PI / 180.0);

        double const denom = cos(decr) * cos(sdecr) * cos(srar - rar)
                             + sin(decr) * sin(sdecr);
        if (denom <= 0)
            continue;

        double const sx = cos(sdecr) * sin(srar - rar) / denom;
        double const sy = (sin(sdecr) * cos(decr)
                           - cos(sdecr) * sin(decr) * cos(srar - rar)) / denom;

        double const ccdx = ccdW - (pa * sx + pb * sy + pc);
        double const ccdy =          pd * sx + pe * sy + pf;

        drawn += drawStar(fb, bsMag,
                          static_cast<float>(ccdx),
                          static_cast<float>(ccdy),
                          params.exposureSecs,
                          imageScaleX, imageScaleY, params);
    }

    drawSkyGlow(fb, params.exposureSecs, imageScaleX, imageScaleY, params);
    applyReadNoise(fb, params.bias, params.readNoise, params.maxVal);

    result.starsDrawn = drawn;

    if (!writeFITS(fb, params, imageScaleX, outputPath))
    {
        result.starsDrawn = 0;
        result.filepath.clear();
        return result;
    }

    result.filepath = outputPath;
    return result;
}
