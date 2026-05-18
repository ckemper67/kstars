/*  Ekos Mount Model - Nearest-Neighbor Sort
    SPDX-FileCopyrightText: 2026 Christian Kemper <ckemper@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

#include "mountmodelsort.h"

#include <algorithm>
#include <cmath>

namespace Ekos
{

// Normalize hour angle to [-12, 12).
static double rangeHours(double h)
{
    while (h >  12.0) h -= 24.0;
    while (h < -12.0) h += 24.0;
    return h;
}

// Shortest wrap-aware azimuth difference in [0, 180].
static double wrapAz(double d)
{
    return d > 180.0 ? 360.0 - d : d;
}

std::vector<int> mountModelNearestNeighborOrder(
    const std::vector<MountSortPoint> &pts,
    const MountSortPoint &start,
    bool isAltAz,
    double lst_h)
{
    const int n = static_cast<int>(pts.size());
    if (n == 0)
        return {};

    auto cost = [&](const MountSortPoint &a, const MountSortPoint &b) -> double
    {
        if (isAltAz)
            return std::max(wrapAz(std::abs(a.az_deg - b.az_deg)),
                            std::abs(a.alt_deg - b.alt_deg));
        const double haA = rangeHours(lst_h - a.ra_h);
        const double haB = rangeHours(lst_h - b.ra_h);
        return std::max(std::abs(haA - haB) * 15.0,
                        std::abs(a.dec_deg - b.dec_deg));
    };

    std::vector<bool> visited(n, false);
    std::vector<int>  order;
    order.reserve(n);

    // First point: closest to start.
    double best = 1e9;
    int    idx  = 0;
    for (int i = 0; i < n; i++)
    {
        double d = cost(start, pts[i]);
        if (d < best) { best = d; idx = i; }
    }
    order.push_back(idx);
    visited[idx] = true;

    // Greedy nearest-neighbor for the rest.
    for (int step = 1; step < n; step++)
    {
        const MountSortPoint &cur = pts[order.back()];
        best = 1e9;
        idx  = -1;
        for (int i = 0; i < n; i++)
        {
            if (visited[i]) continue;
            double d = cost(cur, pts[i]);
            if (d < best) { best = d; idx = i; }
        }
        order.push_back(idx);
        visited[idx] = true;
    }

    return order;
}

} // namespace Ekos
