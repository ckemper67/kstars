/*
    SPDX-FileCopyrightText: 2026 Christian Kemper <ckemper@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

#include <QtGlobal>
#if QT_VERSION >= QT_VERSION_CHECK(6, 0, 0)
#include <QtTest/QTest>
#else
#include <QTest>
#endif

#include <QDir>
#include <QElapsedTimer>
#include <QTemporaryDir>

#include "fitsviewer/fitsdata.h"
#include "ekos/auxiliary/solverutils.h"
#include "ekos/auxiliary/stellarsolverprofile.h"
#include "Options.h"

#include "starfieldrenderer.h"

class TestStarFieldRenderer : public QObject
{
        Q_OBJECT

    private:
        static bool ensureIndexFiles();

    private Q_SLOTS:
        void initTestCase();
        void testRenderProducesFile();
        void testRoundTrip();
};

bool TestStarFieldRenderer::ensureIndexFiles()
{
    if (Options::astrometryIndexFolderList().isEmpty())
    {
        const QStringList candidates =
        {
            QDir::homePath() + "/Library/Application Support/kstars/astrometry",
            QDir::homePath() + "/.local/share/kstars/astrometry",
        };
        for (const auto &d : candidates)
            if (QDir(d).exists())
            {
                Options::setAstrometryIndexFolderList(QStringList() << d);
                break;
            }
    }
    return !Options::astrometryIndexFolderList().isEmpty();
}

void TestStarFieldRenderer::initTestCase()
{
    Options::setAutoDebayer(false);
}

// Verifies that renderStarField() produces a non-empty FITS file with stars.
void TestStarFieldRenderer::testRenderProducesFile()
{
    QTemporaryDir dir;
    QVERIFY(dir.isValid());

    RenderParams params;
    params.ra           = 83.8;  // Orion belt region
    params.dec          = -5.4;
    params.focalLengthMM = 500.0;
    params.pixelSizeUM  = 5.0;
    params.width        = 1280;
    params.height       = 1024;

    RenderedField field = renderStarField(params, dir.filePath("orion.fits"));

    QVERIFY2(field.starsDrawn >= 0, "renderer reported error (GSC/GSCDAT missing?)");
    QVERIFY2(!field.filepath.isEmpty(), "renderer returned empty filepath");
    QVERIFY2(QFileInfo::exists(field.filepath), "output FITS file was not created");
    QVERIFY2(field.starsDrawn > 5, qPrintable(
                 QString("too few stars drawn: %1").arg(field.starsDrawn)));

    // Reported metadata must round-trip through the struct.
    QCOMPARE_EQ(field.ra,  params.ra);
    QCOMPARE_EQ(field.dec, params.dec);
    QVERIFY2(field.pixscale > 0.0, "pixscale must be positive");
}

// Renders a star field and then plate-solves it; asserts that the solution
// agrees with the rendered field parameters within tolerances.
void TestStarFieldRenderer::testRoundTrip()
{
    if (!ensureIndexFiles())
        QSKIP("No astrometry index files found -- skipping round-trip solve");

    QTemporaryDir dir;
    QVERIFY(dir.isValid());

    // Use the Orion belt region: dense enough for easy solving at any algorithm.
    RenderParams params;
    params.ra            = 83.8;
    params.dec           = -5.4;
    params.focalLengthMM = 500.0;
    params.pixelSizeUM   = 5.0;
    params.width         = 1280;
    params.height        = 1024;
    params.exposureSecs  = 300.0f;
    params.seeingArcsec  = 6.0f;
    params.limitingMag   = 14.0f;

    RenderedField field = renderStarField(params, dir.filePath("orion_rt.fits"));
    if (field.starsDrawn < 0)
        QSKIP("Renderer failed -- GSC catalog data not available");
    QVERIFY(field.starsDrawn > 5);

    auto image = QSharedPointer<FITSData>(new FITSData());
    QVERIFY(image->loadFromFile(field.filepath).result());

    auto profiles = Ekos::getDefaultAlignOptionsProfiles();
    auto solverParams = profiles.at(3);
    solverParams.keepNum = 50;
    solverParams.minwidth = 0.1;
    solverParams.maxwidth = 10.0;

    QSharedPointer<SolverUtils> solver(new SolverUtils(solverParams, /*timeoutSecs=*/60),
                                       &QObject::deleteLater);
    solver->useScale(true, field.pixscale * 0.9, field.pixscale * 1.1, ARCSEC_PER_PIX);
    solver->usePosition(true, field.ra, field.dec);

    bool finished = false;
    bool success  = false;
    FITSImage::Solution solution;

    connect(solver.get(), &SolverUtils::done, this,
            [&](bool timedOut, bool ok, const FITSImage::Solution & sol, double)
    {
        finished = true;
        success  = !timedOut && ok;
        solution = sol;
    });

    solver->runSolver(image);

    QElapsedTimer waitTimer;
    waitTimer.start();
    while (!finished && waitTimer.elapsed() < 65000)
        QTest::qWait(200);

    QVERIFY2(finished, "solver did not finish within timeout");
    QVERIFY2(success,  "solver failed to find a solution");

    // RA tolerance: allow wrap at 360.
    double raDiff = std::abs(solution.ra - field.ra);
    if (raDiff > 180.0) raDiff = 360.0 - raDiff;
    QVERIFY2(raDiff < 0.5,
             qPrintable(QString("RA mismatch: expected %1, got %2 (diff %3 deg)")
                        .arg(field.ra).arg(solution.ra).arg(raDiff)));

    double decDiff = std::abs(solution.dec - field.dec);
    QVERIFY2(decDiff < 0.5,
             qPrintable(QString("Dec mismatch: expected %1, got %2 (diff %3 deg)")
                        .arg(field.dec).arg(solution.dec).arg(decDiff)));

    double scaleDiff = std::abs(solution.pixscale - field.pixscale) / field.pixscale;
    QVERIFY2(scaleDiff < 0.02,
             qPrintable(QString("pixscale mismatch: expected %1, got %2 (rel diff %3)")
                        .arg(field.pixscale).arg(solution.pixscale).arg(scaleDiff)));
}

QTEST_GUILESS_MAIN(TestStarFieldRenderer)

#include "teststarfieldrenderer.moc"
