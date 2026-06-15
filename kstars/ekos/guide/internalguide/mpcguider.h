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
        // Primary user-facing configuration: pick the mount class. Auto runs
        // detection across the servo-lag and harmonic-drive paths in parallel;
        // the first to converge sets the effective type. Explicit non-Auto
        // values pre-arm the appropriate detector and suppress the others.
        // Numeric setters (setServoLag, setMechanicalParams, setBacklash) are
        // advanced overrides; most users should pick a mount type and let
        // detection handle the rest.
        enum class MountType {
            Auto,         // default; detect class from observed behavior
            WormGear,     // rigid plant; FFT IMP only
            DirectDrive,  // enable servo-lag detection (or honor manual hook)
            StrainWave,   // pre-arm harmonic-drive R bump
            Belt,         // single-stage compliance; same gating as WormGear
        };

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
        {
            m_Ks = Ks; m_Bs = Bs; m_J_axis = J_axis;
            m_Initialized = false;
            // 5-state flexible plant has different topology than the rigid
            // 3-state where servo-lag detection is defined; turn it off.
            m_ServoLagDetectorActive = false;
        }

        // Declare a closed-loop velocity servo with first-order tracking lag
        // (direct-drive mounts: 10Micron, ASA, ZWO TC40, etc.). The rigid
        // 3-state plant already parameterizes the motor as a first-order rate
        // loop with mechanical time constant tau = J/Bf; setting J=tau, Bf=1,
        // Kt=1 makes the plant match a DD servo with first-order lag tau.
        // Supersedes setParameters' J/Bf/Kt for the rigid path; call order
        // does not matter. Locks the effective mount type to DirectDrive and
        // disables servo-lag auto-detection for the lifetime of this guider.
        // Most users: prefer setMountType(DirectDrive) and let auto-detect
        // handle tau; setServoLag is an advanced override (tests, power
        // users who measured tau on a bench).
        void setServoLag(double tau)
        {
            if (tau < 1e-6) tau = 1e-6;
            m_J = tau; m_Bf = 1.0; m_Kt = 1.0;
            m_Initialized = false;
            m_ServoLagManuallySet = true;
            m_ServoLagDetectorActive = false;
            m_EffectiveMountType = MountType::DirectDrive;
        }

        // Primary configuration: pick mount class. Defaults to Auto. See the
        // MountType enum for semantics. Calling this resets detector gating
        // to the appropriate defaults for the chosen class, but does not
        // clear previously-set numeric overrides (setServoLag, setMechanical-
        // Params, setBacklash).
        void setMountType(MountType type);

        // The class the user asked for (defaults to Auto).
        MountType getMountType() const { return m_MountType; }

        // The class currently in effect: same as getMountType() unless Auto
        // is set and a detector has converged, in which case this returns
        // the detected class.
        MountType getEffectiveMountType() const { return m_EffectiveMountType; }

        // Auto-detected servo lag in seconds. Returns 0 if detection has
        // not converged or is disabled.
        double getDetectedServoLag() const { return m_AutoServoLag; }

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

        // Mount type framework
        MountType m_MountType { MountType::Auto };
        MountType m_EffectiveMountType { MountType::Auto };
        bool m_ServoLagDetectorActive { true };
        bool m_HarmonicDetectorActive { true };

        // Servo-lag auto-detection state
        bool m_ServoLagManuallySet { false };
        double m_AutoServoLag { 0.0 };           // detected tau (0 = not yet)
        std::vector<double> m_ServoLagSamples;   // ring buffer of tau_est
        double m_PrevU { 0.0 };
        bool m_HasPrevU { false };

        QDateTime m_LastGuideTime;
        int m_GuiderIteration { 0 };
};
