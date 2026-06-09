/*
    SPDX-FileCopyrightText: 2025 Christian Kemper <ckemper@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

#include "starfieldrenderer.h"

#include <QCoreApplication>
#include <QCommandLineParser>
#include <QFileInfo>
#include <QTextStream>

#include <cstdlib>

int main(int argc, char *argv[])
{
    QCoreApplication app(argc, argv);
    QCoreApplication::setApplicationName("render_starfield");
    QCoreApplication::setApplicationVersion("1.0");

    QCommandLineParser parser;
    parser.setApplicationDescription(
        "Render a synthetic FITS star field with WCS headers from GSC + bright-star catalogs.");
    parser.addHelpOption();
    parser.addVersionOption();

    const QList<QCommandLineOption> options =
    {
        { "out",          "Output FITS file path (required).",          "path"                       },
        { "ra",           "Field center RA in degrees.",                "deg",   "0.0"               },
        { "dec",          "Field center Dec in degrees.",               "deg",   "0.0"               },
        { "rotation",     "Position angle in degrees (east of north).", "deg",   "0.0"               },
        { "focal",        "Focal length in mm.",                        "mm",    "500.0"             },
        { "pixel",        "Pixel size in microns.",                     "um",    "5.0"               },
        { "width",        "Image width in pixels.",                     "px",    "1280"              },
        { "height",       "Image height in pixels.",                    "px",    "1024"              },
        { "exposure",     "Exposure time in seconds.",                  "s",     "300.0"             },
        { "seeing",       "Seeing FWHM in arcsec.",                     "as",    "8.0"               },
        { "sky-glow",     "Sky background magnitude.",                  "mag",   "19.5"              },
        { "read-noise",   "Read noise in ADU.",                         "adu",   "10"                },
        { "bias",         "Bias offset in ADU.",                        "adu",   "500"               },
        { "limit-mag",    "Faintest catalog magnitude rendered.",       "mag",   "14.0"              },
        { "saturation",   "Magnitude that saturates a pixel.",          "mag",   "2.0"               },
        { "max-val",      "Pixel value clamp (max).",                   "adu",   "65000"             },
    };
    parser.addOptions(options);
    parser.process(app);

    QTextStream err(stderr);
    QTextStream out(stdout);

    if (!parser.isSet("out"))
    {
        err << "error: --out is required\n";
        parser.showHelp(2);
    }

    RenderParams params;
    params.ra            = parser.value("ra").toDouble();
    params.dec           = parser.value("dec").toDouble();
    params.rotationDeg   = parser.value("rotation").toDouble();
    params.focalLengthMM = parser.value("focal").toDouble();
    params.pixelSizeUM   = parser.value("pixel").toDouble();
    params.width         = parser.value("width").toInt();
    params.height        = parser.value("height").toInt();
    params.exposureSecs  = parser.value("exposure").toFloat();
    params.seeingArcsec  = parser.value("seeing").toFloat();
    params.skyGlowMag    = parser.value("sky-glow").toFloat();
    params.readNoise     = parser.value("read-noise").toInt();
    params.bias          = parser.value("bias").toInt();
    params.limitingMag   = parser.value("limit-mag").toFloat();
    params.saturationMag = parser.value("saturation").toFloat();
    params.maxVal        = parser.value("max-val").toInt();

    const QString outPath = parser.value("out");

    RenderedField field = renderStarField(params, outPath);

    if (field.starsDrawn < 0)
    {
        err << "error: renderer failed (gsc binary or GSCDAT not found)\n";
        return 1;
    }
    if (field.filepath.isEmpty())
    {
        err << "error: failed to write FITS file: " << outPath << "\n";
        return 1;
    }

    out << "wrote=" << QFileInfo(field.filepath).absoluteFilePath() << "\n"
        << "ra="          << field.ra         << "\n"
        << "dec="         << field.dec        << "\n"
        << "rotation="    << field.rotation   << "\n"
        << "pixscale="    << field.pixscale   << "\n"
        << "stars_drawn=" << field.starsDrawn << "\n";

    return 0;
}
