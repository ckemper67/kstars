/*  Ekos Mount Model - Sort
    SPDX-FileCopyrightText: 2026 Christian Kemper <ckemper@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

#include "alignorderoptimization.h"

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

// Shortest wrap-aware azimuth delta in [0, 180].
static double wrapAz(double d)
{
    return d > 180.0 ? 360.0 - d : d;
}

// 2-opt improvement on an open path: start -> pts[order[0]] -> ... -> pts[order[n-1]].
// Tries all (i, j) edge-swap pairs; on the first improvement reverses order[i..j]
// and restarts. Returns true if any improvement was found.
// The inner while loop converges to a 2-opt local optimum before returning.
template<typename CostFn>
static bool twoOpt(
    std::vector<int>              &order,
    const std::vector<AlignOrderPoint> &pts,
    const AlignOrderPoint          &start,
    CostFn                         cost)
{
    const int n   = static_cast<int>(order.size());
    bool      any = false;
    bool      found = true;

    while (found)
    {
        found = false;
        for (int i = 0; i < n - 1 && !found; ++i)
        {
            const AlignOrderPoint &prev = (i == 0) ? start : pts[order[i - 1]];
            for (int j = i + 1; j < n && !found; ++j)
            {
                // Reversing order[i..j] changes two edges:
                //   (prev, order[i]) and (order[j], order[j+1])
                //   -> (prev, order[j]) and (order[i], order[j+1])
                double cur = cost(prev, pts[order[i]]);
                double swp = cost(prev, pts[order[j]]);
                if (j < n - 1)
                {
                    cur += cost(pts[order[j]],     pts[order[j + 1]]);
                    swp += cost(pts[order[i]], pts[order[j + 1]]);
                }
                if (cur - swp > 1e-9)
                {
                    std::reverse(order.begin() + i, order.begin() + j + 1);
                    found = true;
                    any   = true;
                }
            }
        }
    }
    return any;
}

// Or-opt-1 improvement: for each point, find the globally best relocation
// (insert it at a different position in the tour). Restarts from the beginning
// after each successful move. Returns true if any improvement was found.
// Handles 2-opt blind spots such as single-point detours.
template<typename CostFn>
static bool orOpt1(
    std::vector<int>              &order,
    const std::vector<AlignOrderPoint> &pts,
    const AlignOrderPoint          &start,
    CostFn                         cost)
{
    const int n = static_cast<int>(order.size());
    if (n < 3)
        return false;

    bool any = false;

    for (int s = 0; s < n; ++s)
    {
        const AlignOrderPoint &ps = pts[order[s]];

        // reducedPt(r): point at position r in the (n-1)-element tour with s removed.
        // r = -1 always maps to start.
        auto reducedPt = [&](int r) -> const AlignOrderPoint &
        {
            if (r < 0) return start;
            return pts[order[r < s ? r : r + 1]];
        };

        // Net gain from removing ps from its current position:
        //   cost of edges removed - cost of bridge added.
        const AlignOrderPoint &Ls = (s > 0) ? pts[order[s - 1]] : start;
        double removal_gain;
        if (s < n - 1)
        {
            const AlignOrderPoint &Rs = pts[order[s + 1]];
            removal_gain = cost(Ls, ps) + cost(ps, Rs) - cost(Ls, Rs);
        }
        else
        {
            removal_gain = cost(Ls, ps);  // last point: remove the single trailing edge
        }

        // Try inserting ps at each gap g in the reduced tour.
        // g = -1         : before reduced[0]  (predecessor is start)
        // g = 0..n-3     : between reduced[g] and reduced[g+1]
        // g = n-2        : after  reduced[n-2]  (append at end)
        // g = s-1        : no-op (puts ps back where it was)
        double best_benefit = 1e-9;  // minimum improvement threshold
        int    best_g       = -2;    // sentinel: no improving move found

        for (int g = -1; g <= n - 2; ++g)
        {
            if (g == s - 1)
                continue;  // no-op: recreates original position

            const AlignOrderPoint &A = reducedPt(g);
            double ins;
            if (g < n - 2)
            {
                const AlignOrderPoint &B = reducedPt(g + 1);
                ins = cost(A, ps) + cost(ps, B) - cost(A, B);
            }
            else
            {
                ins = cost(A, ps);  // append: only one new edge
            }

            double benefit = removal_gain - ins;
            if (benefit > best_benefit)
            {
                best_benefit = benefit;
                best_g       = g;
            }
        }

        if (best_g != -2)
        {
            int saved = order[s];
            order.erase(order.begin() + s);
            order.insert(order.begin() + (best_g + 1), saved);
            any = true;
            s   = -1;  // restart scan from the beginning after any move
        }
    }
    return any;
}

// ---- Public API --------------------------------------------------------------

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

    // Phase 1: greedy nearest-neighbor initialization.
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

    // Phase 2: local search -- alternate 2-opt and Or-opt-1 until neither fires.
    // Each helper converges to its own local-optimum before returning; the outer
    // loop is needed because a move by one phase can expose new improvements for
    // the other (typically 1-2 outer iterations).
    // Bitwise OR: both phases must run before checking convergence.
    while (twoOpt(order, pts, start, cost) | orOpt1(order, pts, start, cost))
        ;

    return order;
}

} // namespace Ekos
