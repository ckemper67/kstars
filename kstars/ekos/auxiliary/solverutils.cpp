/*
    SPDX-FileCopyrightText: 2022 Hy Murveit <hy@murveit.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

#include "solverutils.h"

#include "fitsviewer/fitsdata.h"
#include "Options.h"
#include <QRegularExpression>
#include <QUuid>

int SolverUtils::s_MultiAlgorithmOverride = -1;

SolverUtils::SolverUtils(const SSolver::Parameters &parameters, double timeoutSeconds,
                         SSolver::ProcessType type) :
    m_Parameters(parameters), m_TimeoutMilliseconds(timeoutSeconds * 1000.0), m_Type(type)
{
    connect(&m_Watcher, &QFutureWatcher<bool>::finished, this, &SolverUtils::executeSolver, Qt::UniqueConnection);
    connect(&m_SolverTimer, &QTimer::timeout, this, &SolverUtils::solverTimeout, Qt::UniqueConnection);

    m_StellarSolver.reset(new StellarSolver());
    m_ActiveSolver = m_StellarSolver.get();
}

SolverUtils::~SolverUtils()
{
    disconnect(&m_Watcher, &QFutureWatcher<bool>::finished, this, &SolverUtils::executeSolver);
    disconnect(&m_SolverTimer, &QTimer::timeout, this, &SolverUtils::solverTimeout);
    if (m_StellarSolver.get())
        disconnect(m_StellarSolver.get(), &StellarSolver::finished, this, &SolverUtils::solverDone);
    if (m_HedgeSolver.get())
        disconnect(m_HedgeSolver.get(), &StellarSolver::finished, this, &SolverUtils::solverDone);
}

void SolverUtils::executeSolver()
{
    runSolver(m_ImageData);
}

void SolverUtils::runSolver(const QString &filename)
{
    m_ImageData.reset(new FITSData(), &QObject::deleteLater);
    QFuture<bool> response = m_ImageData->loadFromFile(filename);
    m_Watcher.setFuture(response);
}

void SolverUtils::setHealpix(int indexToUse, int healpixToUse)
{
    m_IndexToUse = indexToUse;
    m_HealpixToUse = healpixToUse;
}

void SolverUtils::abort(bool wait)
{
    if (m_StellarSolver.get())
    {
        if (wait)
            m_StellarSolver->abortAndWait();
        else
            m_StellarSolver->abort();
    }
    if (m_HedgeSolver.get())
    {
        if (wait)
            m_HedgeSolver->abortAndWait();
        else
            m_HedgeSolver->abort();
    }
}

bool SolverUtils::isRunning() const
{
    if (m_StellarSolver && m_StellarSolver->isRunning()) return true;
    if (m_HedgeSolver && m_HedgeSolver->isRunning()) return true;
    return false;
}

void SolverUtils::getSolutionHealpix(int *indexUsed, int *healpixUsed) const
{
    *indexUsed = m_ActiveSolver->getSolutionIndexNumber();
    *healpixUsed = m_ActiveSolver->getSolutionHealpix();
}

void SolverUtils::configureSolver(StellarSolver *solver, const bool stack)
{
    if (solver->isRunning())
        solver->abort();
    solver->setProperty("ProcessType", m_Type);
    if (stack)
        solver->loadNewImageBuffer(m_ImageData->getStackStatistics(), m_ImageData->getStackImageBuffer());
    else
        solver->loadNewImageBuffer(m_ImageData->getStatistics(), m_ImageData->getImageBuffer());
    solver->setProperty("ExtractorType", Options::solveSextractorType());
    solver->setProperty("SolverType", Options::solverType());

    if (m_IndexToUse >= 0)
    {
        QStringList indexFiles = StellarSolver::getIndexFiles(
                                     Options::astrometryIndexFolderList(), m_IndexToUse, m_HealpixToUse);
        solver->setIndexFilePaths(indexFiles);
    }
    else
        solver->setIndexFolderPaths(Options::astrometryIndexFolderList());

    ExternalProgramPaths externalPaths;
    externalPaths.sextractorBinaryPath = Options::sextractorBinary();
    externalPaths.solverPath = Options::astrometrySolverBinary();
    externalPaths.astapBinaryPath = Options::aSTAPExecutable();
    externalPaths.watneyBinaryPath = Options::watneyBinary();
    externalPaths.wcsPath = Options::astrometryWCSInfo();
    solver->setExternalFilePaths(externalPaths);

    solver->setProperty("AutoGenerateAstroConfig", true);

    auto params = m_Parameters;
    params.partition = Options::stellarSolverPartition();
    solver->setParameters(params);

    if (m_UseScale)
        solver->setSearchScale(m_ScaleLow * 0.8, m_ScaleHigh * 1.2, m_ScaleUnits);
    else
        solver->setProperty("UseScale", false);

    if (m_UsePosition)
        solver->setSearchPositionInDegrees(m_raDegrees, m_decDegrees);
    else
        solver->setProperty("UsePosition", false);

    solver->setLogLevel(SSolver::LOG_NONE);
    solver->setSSLogLevel(SSolver::LOG_NORMAL);
    connect(solver, &StellarSolver::logOutput, this,
            [](const QString &msg) { qInfo("StellarSolver: %s", qPrintable(msg)); });
}

void SolverUtils::prepareSolver(const bool stack)
{
    m_HedgeActive = false;
    m_ActiveSolver = m_StellarSolver.get();

    configureSolver(m_StellarSolver.get(), stack);
    connect(m_StellarSolver.get(), &StellarSolver::finished, this, &SolverUtils::solverDone, Qt::UniqueConnection);

    m_TemporaryFilename.clear();

    const SSolver::SolverType type = static_cast<SSolver::SolverType>(m_StellarSolver->property("SolverType").toInt());
    if(type == SSolver::SOLVER_LOCALASTROMETRY || type == SSolver::SOLVER_ASTAP || type == SSolver::SOLVER_WATNEYASTROMETRY)
    {
        m_TemporaryFilename = QDir::tempPath() + QString("/solver%1.fits").arg(QUuid::createUuid().toString().remove(
                                  QRegularExpression("[-{}]")));
        m_ImageData->saveImage(m_TemporaryFilename);
        m_StellarSolver->setProperty("FileToProcess", m_TemporaryFilename);
    }
    else if (type == SSolver::SOLVER_ONLINEASTROMETRY )
    {
        m_TemporaryFilename = QDir::tempPath() + QString("/solver%1.fits").arg(QUuid::createUuid().toString().remove(
                                  QRegularExpression("[-{}]")));
        m_ImageData->saveImage(m_TemporaryFilename);
        m_StellarSolver->setProperty("FileToProcess", m_TemporaryFilename);
        m_StellarSolver->setProperty("AstrometryAPIKey", Options::astrometryAPIKey());
        m_StellarSolver->setProperty("AstrometryAPIURL", Options::astrometryAPIURL());
    }

    // Determine which multi-algorithm to use.
    int algo = resolveMultiAlgorithm(m_StellarSolver.get());

    const bool canHedge = (algo == MULTI_HEDGE)
                          && (type == SSolver::SOLVER_STELLARSOLVER)
                          && (m_Type == SSolver::SOLVE);

    if (canHedge)
    {
        auto p1 = m_StellarSolver->getCurrentParameters();
        p1.multiAlgorithm = MULTI_SCALES;
        m_StellarSolver->setParameters(p1);

        if (!m_HedgeSolver)
            m_HedgeSolver.reset(new StellarSolver());
        configureSolver(m_HedgeSolver.get(), stack);
        connect(m_HedgeSolver.get(), &StellarSolver::finished, this, &SolverUtils::solverDone, Qt::UniqueConnection);
        auto p2 = m_HedgeSolver->getCurrentParameters();
        p2.multiAlgorithm = MULTI_DEPTHS;
        m_HedgeSolver->setParameters(p2);

        m_HedgeActive = true;
    }
    else
    {
        auto params = m_StellarSolver->getCurrentParameters();
        if (algo != MULTI_HEDGE)
            params.multiAlgorithm = static_cast<SSolver::MultiAlgo>(algo);
        else
            params.multiAlgorithm = MULTI_SCALES;
        m_StellarSolver->setParameters(params);
        m_HedgeSolver.reset();
    }
}

void SolverUtils::runSolver(const QSharedPointer<FITSData> &data, const bool stack)
{
    m_SolverTimer.setSingleShot(true);
    m_SolverTimer.setInterval(m_TimeoutMilliseconds);
    m_SolverTimer.start();
    m_StartTime = QDateTime::currentMSecsSinceEpoch();

    m_ImageData = data;
    prepareSolver(stack);
    m_StellarSolver->start();
    if (m_HedgeActive)
        m_HedgeSolver->start();
}

SolverUtils &SolverUtils::useScale(bool useIt, double scaleLow, double scaleHigh, SSolver::ScaleUnits units)
{
    m_UseScale = useIt;
    m_ScaleLow = scaleLow;
    m_ScaleHigh = scaleHigh;
    m_ScaleUnits = units;
    return *this;
}

SolverUtils &SolverUtils::usePosition(bool useIt, double raDegrees, double decDegrees)
{
    m_UsePosition = useIt;
    m_raDegrees = raDegrees;
    m_decDegrees = decDegrees;
    return *this;
}

void SolverUtils::solverDone()
{
    StellarSolver *finished = qobject_cast<StellarSolver*>(sender());
    if (!finished)
        finished = m_StellarSolver.get();

    if (m_HedgeActive)
    {
        bool success;
        FITSImage::Solution solution;

        if (m_Type == SSolver::SOLVE)
        {
            success = finished->solvingDone() && !finished->failed();
            if (success)
                solution = finished->getSolution();
        }
        else
        {
            success = finished->extractionDone() && !finished->failed();
        }

        if (success)
        {
            const double elapsed = (QDateTime::currentMSecsSinceEpoch() - m_StartTime) / 1000.0;
            m_SolverTimer.stop();
            m_HedgeActive = false;
            m_ActiveSolver = finished;

            StellarSolver *other = (finished == m_StellarSolver.get())
                ? m_HedgeSolver.get() : m_StellarSolver.get();
            disconnect(other, &StellarSolver::finished, this, &SolverUtils::solverDone);
            other->abort();

            const char *result = (finished == m_StellarSolver.get()) ? "MULTI_SCALES" : "MULTI_DEPTHS";
            qInfo("Hedge result: %s in %.1fs", result, elapsed);

            Q_EMIT done(false, true, solution, elapsed);
            if (!m_TemporaryFilename.isEmpty())
                QFile::remove(m_TemporaryFilename);
            m_TemporaryFilename.clear();
            return;
        }

        // This solver failed. If the other is still running, wait for it.
        StellarSolver *other = (finished == m_StellarSolver.get())
            ? m_HedgeSolver.get() : m_StellarSolver.get();
        if (other && other->isRunning())
            return;

        // Both failed.
        const double elapsed = (QDateTime::currentMSecsSinceEpoch() - m_StartTime) / 1000.0;
        m_SolverTimer.stop();
        m_HedgeActive = false;
        Q_EMIT done(false, false, FITSImage::Solution(), elapsed);
        if (!m_TemporaryFilename.isEmpty())
            QFile::remove(m_TemporaryFilename);
        m_TemporaryFilename.clear();
        return;
    }

    // Non-hedge path (original logic).
    const double elapsed = (QDateTime::currentMSecsSinceEpoch() - m_StartTime) / 1000.0;
    m_SolverTimer.stop();

    if (m_Type == SSolver::SOLVE)
    {
        FITSImage::Solution solution;
        const bool success = m_StellarSolver->solvingDone() && !m_StellarSolver->failed();
        if (success)
            solution = m_StellarSolver->getSolution();
        Q_EMIT done(false, success, solution, elapsed);
    }
    else
    {
        const bool success = m_StellarSolver->extractionDone() && !m_StellarSolver->failed();
        Q_EMIT done(false, success, FITSImage::Solution(), elapsed);
    }
    if (!m_TemporaryFilename.isEmpty())
        QFile::remove(m_TemporaryFilename);
    m_TemporaryFilename.clear();
}

void SolverUtils::solverTimeout()
{
    m_SolverTimer.stop();

    disconnect(m_StellarSolver.get(), &StellarSolver::finished, this, &SolverUtils::solverDone);
    if (m_HedgeSolver)
        disconnect(m_HedgeSolver.get(), &StellarSolver::finished, this, &SolverUtils::solverDone);

    m_HedgeActive = false;
    abort();

    FITSImage::Solution empty;
    Q_EMIT done(true, false, empty, m_TimeoutMilliseconds / 1000.0);
    if (!m_TemporaryFilename.isEmpty())
        QFile::remove(m_TemporaryFilename);
    m_TemporaryFilename.clear();
}

int SolverUtils::resolveMultiAlgorithm(StellarSolver *solver)
{
    if (s_MultiAlgorithmOverride >= 0)
        return s_MultiAlgorithmOverride;

    if (!solver)
        return MULTI_HEDGE;

    auto algo = solver->getCurrentParameters().multiAlgorithm;
    if (algo == MULTI_AUTO || algo == MULTI_DEPTHS)
        return MULTI_HEDGE;

    return algo;
}

// Legacy entry point -- kept for callers outside SolverUtils.
void SolverUtils::patchMultiAlgorithm(StellarSolver *solver)
{
    if (!solver)
        return;

    int algo = resolveMultiAlgorithm(solver);
    if (algo == MULTI_HEDGE)
        algo = MULTI_SCALES;

    auto params = solver->getCurrentParameters();
    params.multiAlgorithm = static_cast<SSolver::MultiAlgo>(algo);
    solver->setParameters(params);
}
