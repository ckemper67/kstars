/*
    SPDX-FileCopyrightText: 2026 Christian Kemper <ckemper@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

#include "test_mountmodel_halton.h"
#include "auxiliary/dms.h"
#include "skyobjects/skypoint.h"
#include "../../kstars/skycomponents/linelist.h"
#include "../../kstars/skycomponents/artificialhorizoncomponent.h"
#include "../../kstars/ekos/align/haltonsequence.h"

#include <cmath>
#include <QtTest>

// Number of Halton points to generate per test case.
static constexpr int kNumPoints = 50;

// Matches the constants in mountmodel.cpp.
static constexpr double kMaxAlt    = 85.0;
static constexpr double kMaxAbsDec = 80.0;

TestMountModelHalton::TestMountModelHalton() : QObject() {}

// Replicate the halton() member from MountModel (3-line pure function).
double TestMountModelHalton::halton(int index, int base)
{
    double result = 0.0;
    double f      = 1.0;
    while (index > 0)
    {
        f /= static_cast<double>(base);
        result += (index % base) * f;
        index /= base;
    }
    return result;
}

// ---------------------------------------------------------------------------
// testPointsAboveHorizon
//
// For each (lat, LST, minAlt) combination, generate kNumPoints via the same
// AltAz-space Halton logic used in MountModel::slotWizardAlignmentPoints(),
// then convert the stored (ra, dec) back to AltAz and verify alt >= minAlt.
// ---------------------------------------------------------------------------

void TestMountModelHalton::testPointsAboveHorizon_data()
{
    QTest::addColumn<double>("lat");
    QTest::addColumn<double>("lst");   // hours
    QTest::addColumn<double>("minAlt");

    // Mid-northern latitudes
    QTest::newRow("lat=34 lst=0  minAlt=10") << 34.0 << 0.0  << 10.0;
    QTest::newRow("lat=34 lst=6  minAlt=20") << 34.0 << 6.0  << 20.0;
    QTest::newRow("lat=34 lst=12 minAlt=30") << 34.0 << 12.0 << 30.0;
    QTest::newRow("lat=34 lst=18 minAlt=30") << 34.0 << 18.0 << 30.0;

    // Mid-southern latitudes
    QTest::newRow("lat=-34 lst=0  minAlt=20") << -34.0 << 0.0  << 20.0;
    QTest::newRow("lat=-34 lst=12 minAlt=30") << -34.0 << 12.0 << 30.0;

    // Near-equatorial
    QTest::newRow("lat=0 lst=6 minAlt=20")  << 0.0  << 6.0  << 20.0;
    QTest::newRow("lat=5 lst=12 minAlt=20") << 5.0  << 12.0 << 20.0;

    // High northern latitude (NCP nearly overhead)
    QTest::newRow("lat=60 lst=0  minAlt=10") << 60.0 << 0.0  << 10.0;
    QTest::newRow("lat=60 lst=12 minAlt=20") << 60.0 << 12.0 << 20.0;

    // High southern latitude
    QTest::newRow("lat=-60 lst=6 minAlt=10") << -60.0 << 6.0 << 10.0;
}

void TestMountModelHalton::testPointsAboveHorizon()
{
    QFETCH(double, lat);
    QFETCH(double, lst);
    QFETCH(double, minAlt);

    dms latDms(lat);
    dms lstDms;
    lstDms.setH(lst);

    double sinMin = std::sin(minAlt   * dms::DegToRad);
    double sinMax = std::sin(kMaxAlt  * dms::DegToRad);

    for (int i = 1; i <= kNumPoints; i++)
    {
        double az  = halton(i, 2) * 360.0;
        double alt = std::asin(sinMin + halton(i, 3) * (sinMax - sinMin)) / dms::DegToRad;

        // Convert AltAz -> equatorial (same as the implementation).
        SkyPoint sp;
        sp.setAlt(alt);
        sp.setAz(az);
        sp.HorizontalToEquatorial(&lstDms, &latDms);

        double ra  = sp.ra().Hours();
        double dec = std::copysign(qMin(std::abs(sp.dec().Degrees()), kMaxAbsDec),
                                   sp.dec().Degrees());

        // Convert (ra, dec) back to AltAz to verify the point is above the horizon.
        SkyPoint check;
        check.setRA(ra);
        check.setDec(dec);
        check.EquatorialToHorizontal(&lstDms, &latDms);

        double checkAlt = check.alt().Degrees();
        QVERIFY2(checkAlt >= minAlt - 0.001,
                 qPrintable(QString("Point %1: az=%2 alt=%3 -> ra=%4 dec=%5 -> checkAlt=%6 < minAlt=%7")
                            .arg(i).arg(az, 0, 'f', 2).arg(alt, 0, 'f', 2)
                            .arg(ra, 0, 'f', 4).arg(dec, 0, 'f', 4)
                            .arg(checkAlt, 0, 'f', 4).arg(minAlt)));
    }
}

// ---------------------------------------------------------------------------
// testPointsAwayFromPole
//
// Verify that after Dec clamping, |Dec| <= kMaxAbsDec for all points and
// all observer locations.
// ---------------------------------------------------------------------------

void TestMountModelHalton::testPointsAwayFromPole_data()
{
    QTest::addColumn<double>("lat");
    QTest::addColumn<double>("lst");
    QTest::addColumn<double>("minAlt");

    QTest::newRow("lat=34  lst=0  minAlt=30") << 34.0  << 0.0  << 30.0;
    QTest::newRow("lat=-34 lst=12 minAlt=30") << -34.0 << 12.0 << 30.0;
    QTest::newRow("lat=60  lst=6  minAlt=10") << 60.0  << 6.0  << 10.0;
    QTest::newRow("lat=-60 lst=18 minAlt=10") << -60.0 << 18.0 << 10.0;
    QTest::newRow("lat=0   lst=0  minAlt=20") << 0.0   << 0.0  << 20.0;
}

void TestMountModelHalton::testPointsAwayFromPole()
{
    QFETCH(double, lat);
    QFETCH(double, lst);
    QFETCH(double, minAlt);

    dms latDms(lat);
    dms lstDms;
    lstDms.setH(lst);

    double sinMin = std::sin(minAlt  * dms::DegToRad);
    double sinMax = std::sin(kMaxAlt * dms::DegToRad);

    for (int i = 1; i <= kNumPoints; i++)
    {
        double az  = halton(i, 2) * 360.0;
        double alt = std::asin(sinMin + halton(i, 3) * (sinMax - sinMin)) / dms::DegToRad;

        SkyPoint sp;
        sp.setAlt(alt);
        sp.setAz(az);
        sp.HorizontalToEquatorial(&lstDms, &latDms);

        double dec = std::copysign(qMin(std::abs(sp.dec().Degrees()), kMaxAbsDec),
                                   sp.dec().Degrees());

        QVERIFY2(std::abs(dec) <= kMaxAbsDec + 0.001,
                 qPrintable(QString("Point %1: |dec|=%2 exceeds kMaxAbsDec=%3")
                            .arg(i).arg(std::abs(dec), 0, 'f', 4).arg(kMaxAbsDec)));
    }
}

void TestMountModelHalton::testStatefulHaltonSequence()
{
    Ekos::HaltonSequence hs2(2);
    Ekos::HaltonSequence hs3(3);
    for (int i = 1; i <= 2000; ++i)
    {
        double val2_stateful = hs2.next();
        double val2_stateless = halton(i, 2);
        QCOMPARE(hs2.index(), i);
        QVERIFY2(std::abs(val2_stateful - val2_stateless) < 1e-9,
                 qPrintable(QString("Base 2 mismatch at index %1: stateful=%2, stateless=%3")
                            .arg(i).arg(val2_stateful).arg(val2_stateless)));

        double val3_stateful = hs3.next();
        double val3_stateless = halton(i, 3);
        QCOMPARE(hs3.index(), i);
        QVERIFY2(std::abs(val3_stateful - val3_stateless) < 1e-9,
                 qPrintable(QString("Base 3 mismatch at index %1: stateful=%2, stateless=%3")
                            .arg(i).arg(val3_stateful).arg(val3_stateless)));
    }
}

void TestMountModelHalton::testHorizonRejection()
{
    ArtificialHorizon horizon;
    horizon.setTesting();

    // Setup a blocked region: Azimuth 45 to 135 degrees, altitude up to 50 degrees
    std::shared_ptr<LineList> list(new LineList());
    
    auto p1 = std::make_shared<SkyPoint>();
    p1->setAz(dms(45.0));
    p1->setAlt(dms(50.0));
    list->append(p1);

    auto p2 = std::make_shared<SkyPoint>();
    p2->setAz(dms(135.0));
    p2->setAlt(dms(50.0));
    list->append(p2);

    horizon.addRegion("BlockedRegion", true, list, false);
    QVERIFY(horizon.altitudeConstraintsExist());

    // Point generation simulation parameters
    constexpr double minAlt = 15.0;
    constexpr double maxAlt = 85.0;
    double sinMin = std::sin(minAlt * dms::DegToRad);
    double sinMax = std::sin(maxAlt * dms::DegToRad);
    
    Ekos::HaltonSequence haltonAz(2);
    Ekos::HaltonSequence haltonAlt(3);
    
    int targetPoints = 100;
    int generatedCount = 0;
    int rejectedCount = 0;
    
    constexpr int MAX_CANDIDATES = 10000;
    
    while (generatedCount < targetPoints && haltonAz.index() <= MAX_CANDIDATES)
    {
        double az = haltonAz.next() * 360.0;
        double alt = std::asin(sinMin + haltonAlt.next() * (sinMax - sinMin)) / dms::DegToRad;
        
        if (!horizon.isAltitudeOK(az, alt, nullptr))
        {
            rejectedCount++;
            if (az >= 45.0 && az <= 135.0)
            {
                double constraint = horizon.altitudeConstraint(az);
                QVERIFY(alt < constraint);
            }
            continue;
        }
        
        if (az >= 45.0 && az <= 135.0)
        {
            QVERIFY2(alt >= 50.0, qPrintable(QString("Accepted point in blocked zone! az=%1, alt=%2").arg(az).arg(alt)));
        }
        
        generatedCount++;
    }
    
    QCOMPARE(generatedCount, targetPoints);
    QVERIFY(rejectedCount > 0); // Ensure points were actually filtered and rejected
}

QTEST_GUILESS_MAIN(TestMountModelHalton)
