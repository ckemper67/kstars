/*  Ekos Mount Model - Align Order Optimization
    SPDX-FileCopyrightText: 2026 Christian Kemper <ckemper@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

#pragma once

#include <vector>

namespace Ekos
{

// Default pier-flip penalty for GEM mounts: ~90 s flip at 3 deg/s = 270 deg.
// Used as the default for alignOrderOptimization's flipPenaltyDeg parameter.
constexpr double ALIGN_FLIP_PENALTY_DEG = 270.0;

/**
 * @brief Pre-parsed alignment point for the mount model sort.
 *
 * Equatorial coordinates (ra_h, dec_deg) are always populated.
 * Horizontal coordinates (az_deg, alt_deg) are populated only for alt-az
 * mounts; they are pre-computed once before the sort loop.
 */
struct AlignOrderPoint
{
    double ra_h    = 0;  // right ascension in hours
    double dec_deg = 0;  // declination in degrees
    double az_deg  = 0;  // azimuth in degrees  (alt-az mounts only)
    double alt_deg = 0;  // altitude in degrees (alt-az mounts only)
};

/**
 * @brief Sort mount model alignment points to minimize total slew cost.
 *
 * Returns a permutation of 0..pts.size()-1. Both mount axes slew simultaneously
 * so the dominant axis determines travel cost (Chebyshev distance):
 *   Alt-az  (isAltAz=true):  max(|delta_Az|, |delta_Alt|)
 *   Equatorial (isAltAz=false): max(|delta_HA|*15, |delta_Dec|) + flipPenaltyDeg
 *                               when the move crosses the meridian (HA sign changes).
 *
 * The flip penalty makes meridian-crossing moves expensive so the optimizer
 * naturally groups same-pier-side points and avoids unnecessary pier flips.
 * A typical GEM pier flip takes 60-120 s; at 3 deg/s that is 180-360 deg
 * equivalent. Set flipPenaltyDeg=0 to disable flip avoidance.
 *
 * Algorithm: greedy nearest-neighbor. For mount model sizes (N <= ~100) this
 * runs in microseconds.
 *
 * This function has no dependencies beyond the C++ standard library and can be
 * compiled and tested standalone (no KStars, Qt, or INDI headers required).
 *
 * @param pts            Pre-parsed alignment points.
 * @param start          Starting position in the same format as pts.
 * @param isAltAz        True for alt-az mounts; false for equatorial (GEM or fork).
 * @param lst_h          Local sidereal time in hours (used for HA = LST - RA).
 * @param flipPenaltyDeg Extra cost added to any equatorial move that crosses the
 *                       meridian. Ignored for alt-az mounts. Default 0.
 * @return               Permutation: result[k] is the index into pts[] to visit at step k.
 */
std::vector<int> alignOrderOptimization(
    const std::vector<AlignOrderPoint> &pts,
    const AlignOrderPoint &start,
    bool isAltAz,
    double lst_h,
    double flipPenaltyDeg = ALIGN_FLIP_PENALTY_DEG);

} // namespace Ekos
