/*
 * Standalone unit tests for mountModelNearestNeighborOrder.
 *
 * Compiles and links mountmodelsort.cpp directly -- no KStars, Qt, or INDI
 * headers required.
 *
 * Build:
 *   g++ -std=c++17 -O2 \
 *       -I../../../kstars/ekos/align \
 *       standalone_mount_model_sort_algo_test.cpp \
 *       ../../../kstars/ekos/align/mountmodelsort.cpp \
 *       -o standalone_mount_model_sort_algo_test \
 *   && ./standalone_mount_model_sort_algo_test
 */

#include "mountmodelsort.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <vector>

using Ekos::MountSortPoint;

// ---- Helpers ----------------------------------------------------------------

static constexpr double DEG = M_PI / 180.0;

static double halton(int index, int base)
{
    double result = 0.0, f = 1.0;
    while (index > 0) { f /= base; result += (index % base) * f; index /= base; }
    return result;
}

static double rangeHA(double ha)
{
    while (ha >  12.0) ha -= 24.0;
    while (ha < -12.0) ha += 24.0;
    return ha;
}

static double rangeRA(double ra)
{
    while (ra <   0.0) ra += 24.0;
    while (ra >= 24.0) ra -= 24.0;
    return ra;
}

static double wrapAz(double d) { return d > 180.0 ? 360.0 - d : d; }

// Az/Alt -> RA/Dec (Az from North toward East, KStars convention).
static void horizontalToEquatorial(double az_deg, double alt_deg,
                                   double lst_h,  double lat_deg,
                                   double &ra_h,  double &dec_deg)
{
    const double az  = az_deg  * DEG;
    const double alt = alt_deg * DEG;
    const double lat = lat_deg * DEG;

    const double sin_dec = std::sin(lat) * std::sin(alt)
                         + std::cos(lat) * std::cos(alt) * std::cos(az);
    dec_deg = std::asin(std::max(-1.0, std::min(1.0, sin_dec))) / DEG;

    const double ha_rad = std::atan2(-std::sin(az) * std::cos(alt),
                                      std::sin(alt) * std::cos(lat)
                                      - std::cos(alt) * std::sin(lat) * std::cos(az));
    ra_h = rangeRA(lst_h - ha_rad / DEG / 15.0);
}

// Generate N Halton points in AltAz space and convert to MountSortPoint.
static std::vector<MountSortPoint> generateHaltonPoints(
    int n, double lst_h, double lat_deg, double minAlt, double maxAlt, double maxAbsDec)
{
    const double sinMin = std::sin(minAlt * DEG);
    const double sinMax = std::sin(maxAlt * DEG);

    std::vector<MountSortPoint> pts;
    pts.reserve(n);

    for (int i = 1; static_cast<int>(pts.size()) < n && i <= n * 10; i++)
    {
        const double az  = halton(i, 2) * 360.0;
        const double alt = std::asin(sinMin + halton(i, 3) * (sinMax - sinMin)) / DEG;

        double ra_h, dec_deg;
        horizontalToEquatorial(az, alt, lst_h, lat_deg, ra_h, dec_deg);
        if (std::abs(dec_deg) > maxAbsDec)
            continue;

        MountSortPoint p;
        p.ra_h    = ra_h;
        p.dec_deg = dec_deg;
        p.az_deg  = az;
        p.alt_deg = alt;
        pts.push_back(p);
    }
    return pts;
}

// Total cost of a visit order under the equatorial HA metric.
static double totalCostEquatorial(const std::vector<MountSortPoint> &pts,
                                  const std::vector<int> &order,
                                  const MountSortPoint &start,
                                  double lst_h)
{
    double total = 0;
    const MountSortPoint *prev = &start;
    for (int idx : order)
    {
        double haA = rangeHA(lst_h - prev->ra_h);
        double haB = rangeHA(lst_h - pts[idx].ra_h);
        total += std::max(std::abs(haA - haB) * 15.0,
                          std::abs(prev->dec_deg - pts[idx].dec_deg));
        prev = &pts[idx];
    }
    return total;
}

// Total cost of a visit order under the alt-az axis metric.
static double totalCostAltAz(const std::vector<MountSortPoint> &pts,
                              const std::vector<int> &order,
                              const MountSortPoint &start)
{
    double total = 0;
    const MountSortPoint *prev = &start;
    for (int idx : order)
    {
        total += std::max(wrapAz(std::abs(prev->az_deg - pts[idx].az_deg)),
                          std::abs(prev->alt_deg - pts[idx].alt_deg));
        prev = &pts[idx];
    }
    return total;
}

// Check that order is a valid permutation of 0..n-1.
static bool isValidPermutation(const std::vector<int> &order, int n)
{
    if (static_cast<int>(order.size()) != n) return false;
    std::vector<bool> seen(n, false);
    for (int i : order)
    {
        if (i < 0 || i >= n || seen[i]) return false;
        seen[i] = true;
    }
    return true;
}

// ---- Test infrastructure ----------------------------------------------------

static int g_passed = 0;
static int g_failed = 0;

#define EXPECT(cond, msg) \
    do { \
        if (cond) { g_passed++; } \
        else { g_failed++; std::printf("FAIL [%s:%d] %s\n", __func__, __LINE__, msg); } \
    } while(0)

// ---- Tests ------------------------------------------------------------------

static void testEmptyInput()
{
    std::vector<MountSortPoint> pts;
    MountSortPoint start;
    auto order = Ekos::mountModelNearestNeighborOrder(pts, start, false, 6.0);
    EXPECT(order.empty(), "empty input -> empty output");
}

static void testSinglePoint()
{
    MountSortPoint p; p.ra_h = 3.0; p.dec_deg = 20.0;
    MountSortPoint start;
    auto order = Ekos::mountModelNearestNeighborOrder({p}, start, false, 6.0);
    EXPECT(order.size() == 1 && order[0] == 0, "single point -> {0}");
}

static void testTwoPointsEquatorial_closerFirst()
{
    // start at HA=0, Dec=34. Point A is nearby (HA=1), point B is far (HA=5).
    const double lst = 6.0;
    MountSortPoint start; start.ra_h = rangeRA(lst - 0.0); start.dec_deg = 34.0;
    MountSortPoint a;     a.ra_h     = rangeRA(lst - 1.0); a.dec_deg     = 34.0;
    MountSortPoint b;     b.ra_h     = rangeRA(lst - 5.0); b.dec_deg     = 34.0;

    auto order = Ekos::mountModelNearestNeighborOrder({a, b}, start, false, lst);
    EXPECT(isValidPermutation(order, 2), "two points: valid permutation");
    EXPECT(order[0] == 0, "two points: closer point (A) visited first");
}

static void testTwoPointsAltAz_closerFirst()
{
    MountSortPoint start; start.az_deg = 90.0; start.alt_deg = 45.0;
    MountSortPoint a;     a.az_deg     = 100.0; a.alt_deg    = 48.0;  // nearby
    MountSortPoint b;     b.az_deg     = 270.0; b.alt_deg    = 20.0;  // far

    auto order = Ekos::mountModelNearestNeighborOrder({a, b}, start, true, 6.0);
    EXPECT(isValidPermutation(order, 2), "altaz two points: valid permutation");
    EXPECT(order[0] == 0, "altaz two points: closer point (A) visited first");
}

static void testPermutationIsValid_equatorial()
{
    const double lst = 6.0;
    auto pts = generateHaltonPoints(50, lst, 34.0, 20.0, 85.0, 80.0);
    MountSortPoint start; start.ra_h = rangeRA(lst); start.dec_deg = 34.0;

    auto order = Ekos::mountModelNearestNeighborOrder(pts, start, false, lst);
    EXPECT(isValidPermutation(order, static_cast<int>(pts.size())),
           "equatorial N=50: output is valid permutation");
}

static void testPermutationIsValid_altaz()
{
    const double lst = 6.0;
    auto pts = generateHaltonPoints(50, lst, 34.0, 20.0, 85.0, 80.0);
    MountSortPoint start; start.az_deg = 0.0; start.alt_deg = 56.0;

    auto order = Ekos::mountModelNearestNeighborOrder(pts, start, true, lst);
    EXPECT(isValidPermutation(order, static_cast<int>(pts.size())),
           "altaz N=50: output is valid permutation");
}

static void testFirstPointIsNearestEquatorial()
{
    // The first element of the returned order must be the point with minimum
    // cost from start -- that is the invariant of greedy nearest-neighbor.
    const double lst = 6.0;
    auto pts = generateHaltonPoints(50, lst, 34.0, 20.0, 85.0, 80.0);
    MountSortPoint start; start.ra_h = rangeRA(lst); start.dec_deg = 34.0;

    auto order = Ekos::mountModelNearestNeighborOrder(pts, start, false, lst);

    // Find the true nearest point to start.
    int nearestIdx = 0;
    double minCost = 1e9;
    for (int i = 0; i < static_cast<int>(pts.size()); i++)
    {
        double haS = rangeHA(lst - start.ra_h);
        double haI = rangeHA(lst - pts[i].ra_h);
        double d = std::max(std::abs(haS - haI) * 15.0,
                            std::abs(start.dec_deg - pts[i].dec_deg));
        if (d < minCost) { minCost = d; nearestIdx = i; }
    }
    EXPECT(order[0] == nearestIdx, "equatorial: first visited point is nearest to start");
}

static void testEquatorialCostBelowBaseline()
{
    // Greedy sort total cost must be <= cost of the identity order (unsorted).
    const double lst = 6.0;
    auto pts = generateHaltonPoints(100, lst, 34.0, 20.0, 85.0, 80.0);
    MountSortPoint start; start.ra_h = rangeRA(lst); start.dec_deg = 34.0;

    std::vector<int> identity(pts.size());
    std::iota(identity.begin(), identity.end(), 0);

    auto order = Ekos::mountModelNearestNeighborOrder(pts, start, false, lst);

    double sortedCost  = totalCostEquatorial(pts, order,    start, lst);
    double identityCost = totalCostEquatorial(pts, identity, start, lst);

    EXPECT(sortedCost <= identityCost,
           "equatorial N=100: sorted cost <= identity order cost");

    std::printf("  equatorial: sorted=%.1fdeg identity=%.1fdeg (%.0f%% improvement)\n",
                sortedCost, identityCost,
                100.0 * (identityCost - sortedCost) / identityCost);
}

static void testAltAzCostBelowBaseline()
{
    const double lst = 6.0;
    auto pts = generateHaltonPoints(100, lst, 34.0, 20.0, 85.0, 80.0);
    MountSortPoint start; start.az_deg = 0.0; start.alt_deg = 56.0;

    std::vector<int> identity(pts.size());
    std::iota(identity.begin(), identity.end(), 0);

    auto order = Ekos::mountModelNearestNeighborOrder(pts, start, true, lst);

    double sortedCost   = totalCostAltAz(pts, order,    start);
    double identityCost = totalCostAltAz(pts, identity, start);

    EXPECT(sortedCost <= identityCost,
           "altaz N=100: sorted cost <= identity order cost");

    std::printf("  altaz:       sorted=%.1fdeg identity=%.1fdeg (%.0f%% improvement)\n",
                sortedCost, identityCost,
                100.0 * (identityCost - sortedCost) / identityCost);
}

static void testEquatorialGroupsByPierSide()
{
    // For a GEM, the HA metric should prefer same-side-of-meridian moves.
    // Generate 10 points cleanly on each side and verify the sort does not
    // alternate sides more than once (greedy may still do one cross, but
    // should not interleave excessively).
    const double lst = 6.0;
    std::vector<MountSortPoint> pts;

    // 10 points on the east side (HA < 0, ra > lst)
    for (int i = 0; i < 10; i++)
    {
        MountSortPoint p;
        p.ra_h    = rangeRA(lst + 1.0 + i * 0.3);  // HA = -(1+i*0.3)
        p.dec_deg = 20.0 + i * 3.0;
        pts.push_back(p);
    }
    // 10 points on the west side (HA > 0, ra < lst)
    for (int i = 0; i < 10; i++)
    {
        MountSortPoint p;
        p.ra_h    = rangeRA(lst - 1.0 - i * 0.3);  // HA = +(1+i*0.3)
        p.dec_deg = 20.0 + i * 3.0;
        pts.push_back(p);
    }

    MountSortPoint start; start.ra_h = rangeRA(lst); start.dec_deg = 34.0;
    auto order = Ekos::mountModelNearestNeighborOrder(pts, start, false, lst);

    // Count meridian crossings.
    int crossings = 0;
    double prevHA = rangeHA(lst - start.ra_h);
    for (int idx : order)
    {
        double curHA = rangeHA(lst - pts[idx].ra_h);
        if (prevHA * curHA < 0.0) crossings++;
        prevHA = curHA;
    }

    // Expect at most 1 crossing: the sort should keep east/west grouped.
    EXPECT(crossings <= 1, "equatorial pier grouping: at most 1 meridian crossing");
    std::printf("  pier grouping: %d meridian crossing(s)\n", crossings);
}

static void testDeterminism()
{
    const double lst = 6.0;
    auto pts = generateHaltonPoints(50, lst, 34.0, 20.0, 85.0, 80.0);
    MountSortPoint start; start.ra_h = rangeRA(lst); start.dec_deg = 34.0;

    auto order1 = Ekos::mountModelNearestNeighborOrder(pts, start, false, lst);
    auto order2 = Ekos::mountModelNearestNeighborOrder(pts, start, false, lst);
    EXPECT(order1 == order2, "deterministic: same input -> same output");
}

// ---- Main -------------------------------------------------------------------

int main()
{
    std::printf("=== mountModelNearestNeighborOrder unit tests ===\n\n");

    testEmptyInput();
    testSinglePoint();
    testTwoPointsEquatorial_closerFirst();
    testTwoPointsAltAz_closerFirst();
    testPermutationIsValid_equatorial();
    testPermutationIsValid_altaz();
    testFirstPointIsNearestEquatorial();
    testEquatorialCostBelowBaseline();
    testAltAzCostBelowBaseline();
    testEquatorialGroupsByPierSide();
    testDeterminism();

    std::printf("\n%d passed, %d failed\n", g_passed, g_failed);
    return g_failed == 0 ? 0 : 1;
}
