/*
    SPDX-FileCopyrightText: 2022 Hy Murveit <hy@murveit.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

#pragma once

#include <stellarsolver.h>
#undef Const

#include <QObject>
#include <QString>
#include <QTimer>
#include <QFutureWatcher>
#include <mutex>
#include <memory>
#include <QSharedPointer>

#ifdef _WIN32
#undef Unused
#endif

class FITSData;

// This is a wrapper to make calling the StellarSolver solver a bit simpler.
// Must supply the imagedata and stellar solver parameters
// and connect to the signals. Remote solving not supported.
//
// When using the internal StellarSolver for plate solving, two solver
// instances hedge in parallel: one with MULTI_SCALES and one with
// MULTI_DEPTHS. The first successful result wins and the loser is
// cancelled. This is transparent to callers.
class SolverUtils : public QObject
{
        Q_OBJECT

    public:
        SolverUtils(const SSolver::Parameters &parameters, double timeoutSeconds = 15,
                    SSolver::ProcessType type = SSolver::SOLVE);
        ~SolverUtils();

        void runSolver(const QSharedPointer<FITSData> &data, const bool stack = false);
        void runSolver(const QString &filename);
        SolverUtils &useScale(bool useIt, double scaleLow, double scaleHigh, SSolver::ScaleUnits units = ARCSEC_PER_PIX);
        SolverUtils &usePosition(bool useIt, double raDegrees, double decDegrees);
        bool isRunning() const;
        void abort(bool wait = false);

        void setHealpix(int indexToUse = -1, int healpixToUse = -1);
        void getSolutionHealpix(int *indexUsed, int *healpixUsed) const;

        const FITSImage::Background &getBackground() const
        {
            if (!m_ActiveSolver) return *new FITSImage::Background();
            return m_ActiveSolver->getBackground();
        }
        const QList<FITSImage::Star> &getStarList() const
        {
            if (!m_ActiveSolver) return *new QList<FITSImage::Star>();
            return m_ActiveSolver->getStarList();
        }
        int getNumStarsFound() const
        {
            if (!m_ActiveSolver) return 0;
            return m_ActiveSolver->getNumStarsFound();
        };

        // Hedge MULTI_SCALES against MULTI_DEPTHS, take the first to succeed.
        // Defined outside StellarSolver's MultiAlgo enum range.
        static constexpr int MULTI_HEDGE = 100;

        static void patchMultiAlgorithm(StellarSolver *solver);

        // Override the multi-algorithm selection. Accepts any MultiAlgo value
        // or MULTI_HEDGE to force racing.
        static void setMultiAlgorithmOverride(int algo) { s_MultiAlgorithmOverride = algo; }
        static void clearMultiAlgorithmOverride() { s_MultiAlgorithmOverride = -1; }

    Q_SIGNALS:
        void done(bool timedOut, bool success, const FITSImage::Solution &solution, double elapsedSeconds);
        void newLog(const QString &logText);

    private:
        void solverDone();
        void solverTimeout();
        void executeSolver();
        void prepareSolver(const bool stack = false);
        void configureSolver(StellarSolver *solver, const bool stack);
        static int resolveMultiAlgorithm(StellarSolver *solver);

        std::unique_ptr<StellarSolver> m_StellarSolver;
        std::unique_ptr<StellarSolver> m_HedgeSolver;
        StellarSolver *m_ActiveSolver { nullptr };
        bool m_HedgeActive { false };

        qint64 m_StartTime;
        QTimer m_SolverTimer;
        SSolver::Parameters m_Parameters;
        const uint32_t m_TimeoutMilliseconds {0};
        QString m_TemporaryFilename;
        QFutureWatcher<bool> m_Watcher;
        double m_ScaleLow {0}, m_ScaleHigh {0};
        SSolver::ScaleUnits m_ScaleUnits { ARCSEC_PER_PIX };

        QSharedPointer<FITSData> m_ImageData;

        int m_IndexToUse { -1 };
        int m_HealpixToUse { -1 };

        bool m_UseScale { false };
        bool m_UsePosition { false };
        double m_raDegrees { 0.0 };
        double m_decDegrees { 0.0 };

        SSolver::ProcessType m_Type = SSolver::SOLVE;
        std::mutex deleteSolverMutex;

        static int s_MultiAlgorithmOverride;
};
