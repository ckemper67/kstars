/*
    SPDX-FileCopyrightText: 2026 Christian Kemper <ckemper@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

#pragma once

#include <QDateTime>
#include <QString>
#include <memory>
#include <Eigen/Dense>
#include "fftperiodestimator.h"

class MPCSolver;
class TelescopePlant;
class LaguerreNetwork;

class MPCGuider
{
    public:
        MPCGuider(const QString &id);
        ~MPCGuider();

        void setParameters(double Q, double R, double J = 1e-6, double Bf = 1.0, double Kt = 1.0, int N = 10, double alpha = 0.5);
        void setMinMove(double minMove) { m_MinMove = minMove; }
        // Configure gear backlash (arcsec). When >0, the solver's punch-through
        // path injects a one-shot bias on commanded direction reversals to
        // cross the deadband. Set to the known/measured mount backlash.
        void setBacklash(double b) { m_Backlash = b; m_Initialized = false; }

        // Time is implicitly computed. Returns correction in arcseconds.
        double guide(double offset);

        void reset();

        bool isHarmonicDetected() const { return m_IsHarmonicDetected; }

    private:
        QString m_ID;
        double m_MinMove { 0.0 };

        // Controller / Plant parameters
        double m_Q { 10.0 };
        double m_R { 0.1 };
        double m_J { 1e-6 };
        double m_Bf { 1.0 };
        double m_Kt { 1.0 };
        int m_N { 10 };
        double m_alpha { 0.5 };
        double m_Backlash { 0.0 };

        std::unique_ptr<TelescopePlant> m_Plant;
        std::unique_ptr<LaguerreNetwork> m_Network;
        std::unique_ptr<MPCSolver> m_Solver;

        // FFT active frequency learning
        std::unique_ptr<FFTPeriodEstimator> m_FFTEstimator;
        double m_Omega1 { 0.0 };
        double m_Omega2 { 0.0 };
        double m_LearnedT1 { 0.0 };
        double m_LearnedT2 { 0.0 };
        bool m_PeriodsLearned { false };

        double m_LastActiveR { -1.0 };

        // Luenberger Observer vectors
        Eigen::VectorXd m_Xhat;
        Eigen::VectorXd m_XhatPred;
        Eigen::VectorXd m_XhatPrev;

        bool m_Initialized { false };
        double m_LastDt { 0.0 };

        // State reconstruction history
        double m_PrevOffset { 0.0 };
        double m_PrevVelocity { 0.0 };
        bool m_HasPrevOffset { false };
        bool m_HasPrevVelocity { false };

        // Harmonic drive auto-detection passive learning
        std::vector<double> m_OffsetHistory;
        std::vector<double> m_VelocityHistory;
        bool m_IsHarmonicDetected { false };
        int m_HarmonicDetectionCounter { 0 };

        QDateTime m_LastGuideTime;
        int m_GuiderIteration { 0 };
};
