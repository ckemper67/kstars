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

        // Declare a compliant mount (motor-axis 2-mass coupling). Ks is spring
        // stiffness (Nm/rad), Bs is spring damping (Nm*s/rad), J_axis is axis
        // inertia (kg*m^2). When set such that Ks < 1e7 and J_axis > 1e-6, the
        // plant builds a 5-state flexible model instead of the default rigid
        // pure-integrator. This is a user-declared configuration, not online
        // detection. Setting Ks=0 (or leaving unset) restores rigid behavior.
        // NOTE: compliance and the FFT-driven IMP path are currently mutually
        // exclusive -- declaring compliance suppresses the 7-state IMP rebuild
        // so the 5-state flexible plant is used. Full flexible-IMP fusion
        // (9 states) is a follow-up.
        void setMechanicalParams(double Ks, double Bs, double J_axis)
        { m_Ks = Ks; m_Bs = Bs; m_J_axis = J_axis; m_Initialized = false; }

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
        // Mechanical compliance parameters (0 = rigid, default).
        double m_Ks { 0.0 };
        double m_Bs { 0.0 };
        double m_J_axis { 0.0 };

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
