/*  Ekos Mount Model - Align Order Optimization
    SPDX-FileCopyrightText: 2026 Christian Kemper <ckemper@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

#include "alignorderoptimization.h"

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

// Shortest wrap-aware azimuth delta in [0, 180].
static double wrapAz(double d)
{
    return d > 180.0 ? 360.0 - d : d;
}

std::vector<int> alignOrderOptimization(
    const std::vector<AlignOrderPoint> &pts,
    const AlignOrderPoint              &start,
    bool                               isAltAz,
    double                             lst_h,
    double                             flipPenaltyDeg)
{
    const int n = static_cast<int>(pts.size());
    if (n == 0)
        return {};

    auto cost = [&](const AlignOrderPoint &a, const AlignOrderPoint &b) -> double
    {
        if (isAltAz)
            return std::max(wrapAz(std::abs(a.az_deg - b.az_deg)),
                            std::abs(a.alt_deg - b.alt_deg));
        // HA is treated as a linear axis coordinate, not circular.
        // Equatorial mounts have mechanical limits (~+/-6h) that prevent
        // wrapping through HA=12, so the raw difference is the true travel.
        const double haA = rangeHours(lst_h - a.ra_h);
        const double haB = rangeHours(lst_h - b.ra_h);
        double c = std::max(std::abs(haA - haB) * 15.0,
                            std::abs(a.dec_deg - b.dec_deg));
        if (flipPenaltyDeg > 0.0 && haA * haB < 0.0)
            c += flipPenaltyDeg;
        return c;
    };

    // Greedy nearest-neighbor: at each step pick the unvisited point closest
    // to the current position under the axis-travel cost function.
    std::vector<bool> visited(n, false);
    std::vector<int>  order;
    order.reserve(n);

    double best = 1e9;
    int    idx  = 0;
    for (int i = 0; i < n; ++i)
    {
        double d = cost(start, pts[i]);
        if (d < best) { best = d; idx = i; }
    }
    order.push_back(idx);
    visited[idx] = true;

    for (int step = 1; step < n; ++step)
    {
        const AlignOrderPoint &cur = pts[order.back()];
        best = 1e9;
        idx  = -1;
        for (int i = 0; i < n; ++i)
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
