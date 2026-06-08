/*
    SPDX-FileCopyrightText: 2025 Christian Kemper <ckemper@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

#pragma once

#include <QString>

struct RenderParams
{
    double ra              = 0.0;
    double dec             = 0.0;
    double focalLengthMM   = 500.0;
    double pixelSizeUM     = 5.0;
    int    width           = 1280;
    int    height          = 1024;
    float  exposureSecs    = 300.0f;
    float  seeingArcsec    = 8.0f;
    float  skyGlowMag      = 19.5f;
    int    readNoise       = 10;
    int    bias            = 500;
    double rotationDeg     = 0.0;
    float  limitingMag     = 14.0f;
    float  saturationMag   = 2.0f;
    int    maxVal          = 65000;
};

struct RenderedField
{
    QString filepath;
    double  ra         = 0.0;
    double  dec        = 0.0;
    double  pixscale   = 0.0;
    double  rotation   = 0.0;
    int     starsDrawn = 0;
};

RenderedField renderStarField(const RenderParams &params, const QString &outputPath);
