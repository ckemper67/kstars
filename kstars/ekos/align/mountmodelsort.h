/*  Ekos Mount Model - Nearest-Neighbor Sort
    SPDX-FileCopyrightText: 2026 Christian Kemper <ckemper@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

#pragma once

#include <vector>

namespace Ekos
{

/**
 * @brief Pre-parsed alignment point for the mount model sort.
 *
 * Equatorial coordinates (ra_h, dec_deg) are always populated.
 * Horizontal coordinates (az_deg, alt_deg) are populated only for alt-az
 * mounts; they are pre-computed once before the O(N^2) sort loop.
 */
struct MountSortPoint
{
    double ra_h   = 0;   // right ascension in hours
    double dec_deg = 0;  // declination in degrees
    double az_deg  = 0;  // azimuth in degrees  (alt-az mounts only)
    double alt_deg = 0;  // altitude in degrees  (alt-az mounts only)
};

/**
 * @brief Greedy nearest-neighbor sort order for mount model alignment points.
 *
 * Returns a permutation of 0..pts.size()-1. Both mount axes slew simultaneously
 * so the dominant axis determines travel cost (Chebyshev distance):
 *   Alt-az  (isAltAz=true):  max(|delta_Az|, |delta_Alt|)
 *   Equatorial (isAltAz=false): max(|delta_HA|*15, |delta_Dec|)
 *
 * For GEMs, the HA-based metric makes cross-meridian moves expensive, naturally
 * grouping same-pier-side points and reducing unnecessary pier flips.
 *
 * This function has no dependencies beyond the C++ standard library and can be
 * compiled and tested standalone (no KStars, Qt, or INDI headers required).
 *
 * @param pts     Pre-parsed alignment points.
 * @param start   Starting position in the same format as pts.
 * @param isAltAz True for alt-az mounts; false for equatorial (GEM or fork).
 * @param lst_h   Local sidereal time in hours (used for HA = LST - RA).
 * @return        Permutation: result[k] is the index into pts[] to visit at step k.
 */
std::vector<int> mountModelNearestNeighborOrder(
    const std::vector<MountSortPoint> &pts,
    const MountSortPoint &start,
    bool isAltAz,
    double lst_h);

} // namespace Ekos
