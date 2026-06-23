/*
 * Standalone unit tests for alignOrderOptimization.
 *
 * Compiles and links alignorderoptimization.cpp directly -- no KStars, Qt, or INDI
 * headers required.
 *
 * Build:
 *   g++ -std=c++17 -O2 \
 *       -I../../../kstars/ekos/align \
 *       align_order_optimization_test.cpp \
 *       ../../../kstars/ekos/align/alignorderoptimization.cpp \
 *       -o align_order_optimization_test \
 *   && ./align_order_optimization_test
 */

#include "alignorderoptimization.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <vector>

using Ekos::AlignOrderPoint;

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

// Generate N Halton points in AltAz space and convert to AlignOrderPoint.
static std::vector<AlignOrderPoint> generateHaltonPoints(
    int n, double lst_h, double lat_deg, double minAlt, double maxAlt, double maxAbsDec)
{
    const double sinMin = std::sin(minAlt * DEG);
    const double sinMax = std::sin(maxAlt * DEG);

    std::vector<AlignOrderPoint> pts;
    pts.reserve(n);

    for (int i = 1; static_cast<int>(pts.size()) < n && i <= n * 10; i++)
    {
        const double az  = halton(i, 2) * 360.0;
        const double alt = std::asin(sinMin + halton(i, 3) * (sinMax - sinMin)) / DEG;

        double ra_h, dec_deg;
        horizontalToEquatorial(az, alt, lst_h, lat_deg, ra_h, dec_deg);
        if (std::abs(dec_deg) > maxAbsDec)
            continue;

        AlignOrderPoint p;
        p.ra_h    = ra_h;
        p.dec_deg = dec_deg;
        p.az_deg  = az;
        p.alt_deg = alt;
        pts.push_back(p);
    }
    return pts;
}

// Total cost of a visit order under the equatorial HA metric.
static double totalCostEquatorial(const std::vector<AlignOrderPoint> &pts,
                                  const std::vector<int> &order,
                                  const AlignOrderPoint &start,
                                  double lst_h)
{
    double total = 0;
    const AlignOrderPoint *prev = &start;
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
static double totalCostAltAz(const std::vector<AlignOrderPoint> &pts,
                              const std::vector<int> &order,
                              const AlignOrderPoint &start)
{
    double total = 0;
    const AlignOrderPoint *prev = &start;
    for (int idx : order)
    {
        total += std::max(wrapAz(std::abs(prev->az_deg - pts[idx].az_deg)),
                          std::abs(prev->alt_deg - pts[idx].alt_deg));
        prev = &pts[idx];
    }
    return total;
}

// Reference greedy nearest-neighbor: used only in tests to compare against
// the full local-search result.
static std::vector<int> greedyNNEquatorial(
    const std::vector<AlignOrderPoint> &pts,
    const AlignOrderPoint &start,
    double lst_h)
{
    const int n = static_cast<int>(pts.size());
    std::vector<bool> visited(n, false);
    std::vector<int>  order;
    order.reserve(n);

    auto cost = [&](const AlignOrderPoint &a, const AlignOrderPoint &b) -> double
    {
        double haA = rangeHA(lst_h - a.ra_h);
        double haB = rangeHA(lst_h - b.ra_h);
        return std::max(std::abs(haA - haB) * 15.0,
                        std::abs(a.dec_deg - b.dec_deg));
    };

    double best = 1e9; int idx = 0;
    for (int i = 0; i < n; ++i) {
        double d = cost(start, pts[i]);
        if (d < best) { best = d; idx = i; }
    }
    order.push_back(idx); visited[idx] = true;

    for (int step = 1; step < n; ++step) {
        const AlignOrderPoint &cur = pts[order.back()];
        best = 1e9; idx = -1;
        for (int i = 0; i < n; ++i) {
            if (visited[i]) continue;
            double d = cost(cur, pts[i]);
            if (d < best) { best = d; idx = i; }
        }
        order.push_back(idx); visited[idx] = true;
    }
    return order;
}

static std::vector<int> greedyNNAltAz(
    const std::vector<AlignOrderPoint> &pts,
    const AlignOrderPoint &start)
{
    const int n = static_cast<int>(pts.size());
    std::vector<bool> visited(n, false);
    std::vector<int>  order;
    order.reserve(n);

    auto cost = [](const AlignOrderPoint &a, const AlignOrderPoint &b) -> double
    {
        return std::max(wrapAz(std::abs(a.az_deg - b.az_deg)),
                        std::abs(a.alt_deg - b.alt_deg));
    };

    double best = 1e9; int idx = 0;
    for (int i = 0; i < n; ++i) {
        double d = cost(start, pts[i]);
        if (d < best) { best = d; idx = i; }
    }
    order.push_back(idx); visited[idx] = true;

    for (int step = 1; step < n; ++step) {
        const AlignOrderPoint &cur = pts[order.back()];
        best = 1e9; idx = -1;
        for (int i = 0; i < n; ++i) {
            if (visited[i]) continue;
            double d = cost(cur, pts[i]);
            if (d < best) { best = d; idx = i; }
        }
        order.push_back(idx); visited[idx] = true;
    }
    return order;
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
    std::vector<AlignOrderPoint> pts;
    AlignOrderPoint start;
    auto order = Ekos::alignOrderOptimization(pts, start, false, 6.0);
    EXPECT(order.empty(), "empty input -> empty output");
}

static void testSinglePoint()
{
    AlignOrderPoint p; p.ra_h = 3.0; p.dec_deg = 20.0;
    AlignOrderPoint start;
    auto order = Ekos::alignOrderOptimization({p}, start, false, 6.0);
    EXPECT(order.size() == 1 && order[0] == 0, "single point -> {0}");
}

static void testTwoPointsEquatorial_closerFirst()
{
    // start at HA=0, Dec=34. Point A is nearby (HA=1), point B is far (HA=5).
    const double lst = 6.0;
    AlignOrderPoint start; start.ra_h = rangeRA(lst - 0.0); start.dec_deg = 34.0;
    AlignOrderPoint a;     a.ra_h     = rangeRA(lst - 1.0); a.dec_deg     = 34.0;
    AlignOrderPoint b;     b.ra_h     = rangeRA(lst - 5.0); b.dec_deg     = 34.0;

    auto order = Ekos::alignOrderOptimization({a, b}, start, false, lst);
    EXPECT(isValidPermutation(order, 2), "two points: valid permutation");
    EXPECT(order[0] == 0, "two points: closer point (A) visited first");
}

static void testTwoPointsAltAz_closerFirst()
{
    AlignOrderPoint start; start.az_deg = 90.0; start.alt_deg = 45.0;
    AlignOrderPoint a;     a.az_deg     = 100.0; a.alt_deg    = 48.0;  // nearby
    AlignOrderPoint b;     b.az_deg     = 270.0; b.alt_deg    = 20.0;  // far

    auto order = Ekos::alignOrderOptimization({a, b}, start, true, 6.0);
    EXPECT(isValidPermutation(order, 2), "altaz two points: valid permutation");
    EXPECT(order[0] == 0, "altaz two points: closer point (A) visited first");
}

static void testPermutationIsValid_equatorial()
{
    const double lst = 6.0;
    auto pts = generateHaltonPoints(50, lst, 34.0, 20.0, 85.0, 80.0);
    AlignOrderPoint start; start.ra_h = rangeRA(lst); start.dec_deg = 34.0;

    auto order = Ekos::alignOrderOptimization(pts, start, false, lst);
    EXPECT(isValidPermutation(order, static_cast<int>(pts.size())),
           "equatorial N=50: output is valid permutation");
}

static void testPermutationIsValid_altaz()
{
    const double lst = 6.0;
    auto pts = generateHaltonPoints(50, lst, 34.0, 20.0, 85.0, 80.0);
    AlignOrderPoint start; start.az_deg = 0.0; start.alt_deg = 56.0;

    auto order = Ekos::alignOrderOptimization(pts, start, true, lst);
    EXPECT(isValidPermutation(order, static_cast<int>(pts.size())),
           "altaz N=50: output is valid permutation");
}

static void testEquatorialCostBelowBaseline()
{
    // Sort total cost must be <= cost of the identity order (unsorted).
    const double lst = 6.0;
    auto pts = generateHaltonPoints(100, lst, 34.0, 20.0, 85.0, 80.0);
    AlignOrderPoint start; start.ra_h = rangeRA(lst); start.dec_deg = 34.0;

    std::vector<int> identity(pts.size());
    std::iota(identity.begin(), identity.end(), 0);

    auto order = Ekos::alignOrderOptimization(pts, start, false, lst);

    double sortedCost   = totalCostEquatorial(pts, order,    start, lst);
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
    AlignOrderPoint start; start.az_deg = 0.0; start.alt_deg = 56.0;

    std::vector<int> identity(pts.size());
    std::iota(identity.begin(), identity.end(), 0);

    auto order = Ekos::alignOrderOptimization(pts, start, true, lst);

    double sortedCost   = totalCostAltAz(pts, order,    start);
    double identityCost = totalCostAltAz(pts, identity, start);

    EXPECT(sortedCost <= identityCost,
           "altaz N=100: sorted cost <= identity order cost");

    std::printf("  altaz:       sorted=%.1fdeg identity=%.1fdeg (%.0f%% improvement)\n",
                sortedCost, identityCost,
                100.0 * (identityCost - sortedCost) / identityCost);
}

static void testLocalSearchImprovesOverGreedyNN()
{
    // The 2-opt+Or-opt-1 local search must not be worse than plain greedy NN,
    // and should improve total cost on a non-trivial random-ish point set.
    const double lst = 6.0;
    auto pts = generateHaltonPoints(50, lst, 34.0, 20.0, 85.0, 80.0);
    AlignOrderPoint start; start.ra_h = rangeRA(lst); start.dec_deg = 34.0;

    auto orderFull = Ekos::alignOrderOptimization(pts, start, false, lst);
    auto orderNN   = greedyNNEquatorial(pts, start, lst);

    EXPECT(isValidPermutation(orderFull, static_cast<int>(pts.size())),
           "local-search N=50: valid permutation");

    double costFull = totalCostEquatorial(pts, orderFull, start, lst);
    double costNN   = totalCostEquatorial(pts, orderNN,   start, lst);

    EXPECT(costFull <= costNN + 1e-6,
           "local-search N=50: cost <= plain greedy-NN cost");

    std::printf("  local-search vs greedy-NN: %.1fdeg vs %.1fdeg (%.0f%% gain)\n",
                costFull, costNN,
                100.0 * (costNN - costFull) / (costNN > 0 ? costNN : 1.0));
}

static void testLocalSearchImprovesOverGreedyNN_altaz()
{
    const double lst = 6.0;
    auto pts = generateHaltonPoints(50, lst, 34.0, 20.0, 85.0, 80.0);
    AlignOrderPoint start;
    start.ra_h = rangeRA(lst); start.dec_deg = 34.0;
    start.az_deg = 180.0; start.alt_deg = 56.0;

    auto orderFull = Ekos::alignOrderOptimization(pts, start, true, lst);
    auto orderNN   = greedyNNAltAz(pts, start);

    EXPECT(isValidPermutation(orderFull, static_cast<int>(pts.size())),
           "local-search altaz N=50: valid permutation");

    double costFull = totalCostAltAz(pts, orderFull, start);
    double costNN   = totalCostAltAz(pts, orderNN,   start);

    EXPECT(costFull <= costNN + 1e-6,
           "local-search altaz N=50: cost <= plain greedy-NN cost");

    std::printf("  local-search altaz vs greedy-NN: %.1fdeg vs %.1fdeg (%.0f%% gain)\n",
                costFull, costNN,
                100.0 * (costNN - costFull) / (costNN > 0 ? costNN : 1.0));
}

static void testEquatorialGroupsByPierSide()
{
    // For a GEM, the HA metric should prefer same-side-of-meridian moves.
    // Generate 10 points cleanly on each side and verify the sort does not
    // alternate sides excessively.
    const double lst = 6.0;
    std::vector<AlignOrderPoint> pts;

    // 10 points on the east side (HA < 0, ra > lst)
    for (int i = 0; i < 10; i++)
    {
        AlignOrderPoint p;
        p.ra_h    = rangeRA(lst + 1.0 + i * 0.3);  // HA = -(1+i*0.3)
        p.dec_deg = 20.0 + i * 3.0;
        pts.push_back(p);
    }
    // 10 points on the west side (HA > 0, ra < lst)
    for (int i = 0; i < 10; i++)
    {
        AlignOrderPoint p;
        p.ra_h    = rangeRA(lst - 1.0 - i * 0.3);  // HA = +(1+i*0.3)
        p.dec_deg = 20.0 + i * 3.0;
        pts.push_back(p);
    }

    AlignOrderPoint start; start.ra_h = rangeRA(lst); start.dec_deg = 34.0;
    auto order = Ekos::alignOrderOptimization(pts, start, false, lst);

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

static void testAllIdenticalPoints()
{
    const double lst = 6.0;
    std::vector<AlignOrderPoint> pts(10);
    for (auto &p : pts) { p.ra_h = 3.0; p.dec_deg = 45.0; p.az_deg = 90.0; p.alt_deg = 50.0; }
    AlignOrderPoint start; start.ra_h = 6.0; start.dec_deg = 34.0;

    auto orderEq = Ekos::alignOrderOptimization(pts, start, false, lst);
    EXPECT(isValidPermutation(orderEq, 10), "all-identical equatorial: valid permutation");
    EXPECT(totalCostEquatorial(pts, orderEq, start, lst) >= 0.0,
           "all-identical equatorial: non-negative cost");

    start.az_deg = 180.0; start.alt_deg = 56.0;
    auto orderAz = Ekos::alignOrderOptimization(pts, start, true, lst);
    EXPECT(isValidPermutation(orderAz, 10), "all-identical altaz: valid permutation");
}

static void testThreePoints()
{
    const double lst = 6.0;
    std::vector<AlignOrderPoint> pts(3);
    pts[0].ra_h = 5.0; pts[0].dec_deg = 30.0;
    pts[1].ra_h = 7.0; pts[1].dec_deg = 30.0;
    pts[2].ra_h = 6.0; pts[2].dec_deg = 50.0;
    AlignOrderPoint start; start.ra_h = 6.0; start.dec_deg = 34.0;

    auto order = Ekos::alignOrderOptimization(pts, start, false, lst);
    EXPECT(isValidPermutation(order, 3), "N=3 equatorial: valid permutation");

    double costOpt = totalCostEquatorial(pts, order, start, lst);
    auto orderNN = greedyNNEquatorial(pts, start, lst);
    double costNN = totalCostEquatorial(pts, orderNN, start, lst);
    EXPECT(costOpt <= costNN + 1e-6, "N=3 equatorial: local search <= greedy NN");
}

static void testHAWrappingBoundary()
{
    // Points near HA +/-12h. The HA metric intentionally does not wrap,
    // treating HA as a linear axis with mechanical limits.
    const double lst = 6.0;
    std::vector<AlignOrderPoint> pts(4);
    pts[0].ra_h = rangeRA(lst - 11.5); pts[0].dec_deg = 30.0;  // HA = +11.5
    pts[1].ra_h = rangeRA(lst + 11.5); pts[1].dec_deg = 30.0;  // HA = -11.5
    pts[2].ra_h = rangeRA(lst - 2.0);  pts[2].dec_deg = 30.0;  // HA = +2
    pts[3].ra_h = rangeRA(lst + 2.0);  pts[3].dec_deg = 30.0;  // HA = -2
    AlignOrderPoint start; start.ra_h = lst; start.dec_deg = 30.0;

    auto order = Ekos::alignOrderOptimization(pts, start, false, lst);
    EXPECT(isValidPermutation(order, 4), "HA-boundary: valid permutation");

    double costOpt = totalCostEquatorial(pts, order, start, lst);
    auto orderNN = greedyNNEquatorial(pts, start, lst);
    double costNN = totalCostEquatorial(pts, orderNN, start, lst);
    EXPECT(costOpt <= costNN + 1e-6, "HA-boundary: local search <= greedy NN");
}

static void testDeterminism()
{
    const double lst = 6.0;
    auto pts = generateHaltonPoints(50, lst, 34.0, 20.0, 85.0, 80.0);
    AlignOrderPoint start; start.ra_h = rangeRA(lst); start.dec_deg = 34.0;

    auto order1 = Ekos::alignOrderOptimization(pts, start, false, lst);
    auto order2 = Ekos::alignOrderOptimization(pts, start, false, lst);
    EXPECT(order1 == order2, "deterministic: same input -> same output");
}

// Count strict HA sign changes (haPrev * haCur < 0) in an order sequence.
static int countMeridianCrossings(const std::vector<AlignOrderPoint> &pts,
                                  const std::vector<int> &order,
                                  const AlignOrderPoint &start,
                                  double lst_h)
{
    int crossings = 0;
    double prevHA = rangeHA(lst_h - start.ra_h);
    for (int idx : order)
    {
        double curHA = rangeHA(lst_h - pts[idx].ra_h);
        if (prevHA * curHA < 0.0)
            crossings++;
        prevHA = curHA;
    }
    return crossings;
}

static void testTwoPhaseValidPermutation()
{
    const double lst = 6.0;
    auto pts = generateHaltonPoints(60, lst, 34.0, 20.0, 85.0, 80.0);
    AlignOrderPoint start; start.ra_h = rangeRA(lst); start.dec_deg = 34.0;

    auto order = Ekos::alignOrderOptimizationTwoPhase(pts, start, lst);
    EXPECT(isValidPermutation(order, static_cast<int>(pts.size())),
           "two-phase N=60: output is valid permutation");
}

static void testTwoPhaseBothSidesOneFlip()
{
    // Points explicitly placed on both pier sides; two-phase must cross the
    // meridian exactly once.
    const double lst = 6.0;
    std::vector<AlignOrderPoint> pts;

    // 10 east points (HA = -1..-4)
    for (int i = 0; i < 10; i++)
    {
        AlignOrderPoint p;
        p.ra_h    = rangeRA(lst + 1.0 + i * 0.3);
        p.dec_deg = 20.0 + i * 3.0;
        pts.push_back(p);
    }
    // 10 west points (HA = +1..+4)
    for (int i = 0; i < 10; i++)
    {
        AlignOrderPoint p;
        p.ra_h    = rangeRA(lst - 1.0 - i * 0.3);
        p.dec_deg = 20.0 + i * 3.0;
        pts.push_back(p);
    }

    AlignOrderPoint start; start.ra_h = rangeRA(lst); start.dec_deg = 34.0;
    auto order = Ekos::alignOrderOptimizationTwoPhase(pts, start, lst);

    EXPECT(isValidPermutation(order, 20),
           "two-phase both-sides: valid permutation");

    int crossings = countMeridianCrossings(pts, order, start, lst);
    EXPECT(crossings == 1,
           "two-phase both-sides: exactly 1 meridian crossing");
}

static void testTwoPhaseSingleSideFallback()
{
    // All points are east of the meridian (HA < 0).  Two-phase falls back to
    // single-side optimization and must produce a valid permutation with 0
    // meridian crossings.
    const double lst = 6.0;
    std::vector<AlignOrderPoint> pts;
    for (int i = 0; i < 15; i++)
    {
        AlignOrderPoint p;
        p.ra_h    = rangeRA(lst + 1.0 + i * 0.4);  // HA = -(1+i*0.4), all east
        p.dec_deg = 15.0 + i * 4.0;
        pts.push_back(p);
    }

    AlignOrderPoint start; start.ra_h = rangeRA(lst); start.dec_deg = 34.0;
    auto order = Ekos::alignOrderOptimizationTwoPhase(pts, start, lst);

    EXPECT(isValidPermutation(order, 15),
           "two-phase single-side fallback: valid permutation");
    EXPECT(countMeridianCrossings(pts, order, start, lst) == 0,
           "two-phase single-side fallback: 0 meridian crossings");
}

// ---- Main -------------------------------------------------------------------

int main()
{
    std::printf("=== alignOrderOptimization unit tests ===\n\n");

    testEmptyInput();
    testSinglePoint();
    testTwoPointsEquatorial_closerFirst();
    testTwoPointsAltAz_closerFirst();
    testPermutationIsValid_equatorial();
    testPermutationIsValid_altaz();
    testEquatorialCostBelowBaseline();
    testAltAzCostBelowBaseline();
    testLocalSearchImprovesOverGreedyNN();
    testLocalSearchImprovesOverGreedyNN_altaz();
    testEquatorialGroupsByPierSide();
    testAllIdenticalPoints();
    testThreePoints();
    testHAWrappingBoundary();
    testDeterminism();
    testTwoPhaseValidPermutation();
    testTwoPhaseBothSidesOneFlip();
    testTwoPhaseSingleSideFallback();

    std::printf("\n%d passed, %d failed\n", g_passed, g_failed);
    return g_failed == 0 ? 0 : 1;
}
