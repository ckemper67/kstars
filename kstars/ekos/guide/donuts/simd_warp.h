/*
    SPDX-FileCopyrightText: 2026 Christian Kemper <ckemper@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

#pragma once
#include "donuts.h"
#include <vector>

namespace Donuts
{

// Single-channel Catmull-Rom bicubic warp. Processes dst rows [ylo, yhi).
// m maps each destination pixel (xd, yd) to the source pixel (xs, ys):
//   xs = m.a*xd + m.b*yd + m.tx
//   ys = m.c*xd + m.d*yd + m.ty
// Out-of-bounds source pixels write 0.0f in dst.
void warpRows1(const std::vector<float> &src, int w, int h,
               const AffineMatrix &m,
               std::vector<float> &dst,
               int ylo, int yhi);

} // namespace Donuts
