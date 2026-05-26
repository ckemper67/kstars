/*
    SPDX-FileCopyrightText: 2026 Christian Kemper <ckemper@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

#include "mpcguider.h"
#include "MPCSolver.h"
#include "TelescopePlant.h"
#include "LaguerreNetwork.h"

#include <cmath>
#include "ekos_guide_debug.h"

MPCGuider::MPCGuider(const QString &id) : m_ID(id)
{
    reset();
}

MPCGuider::~MPCGuider()
{
}

void MPCGuider::reset()
{
    m_Initialized = false;
    m_LastDt = 0.0;
    m_PrevOffset = 0.0;
    m_PrevVelocity = 0.0;
    m_HasPrevOffset = false;
    m_HasPrevVelocity = false;
    m_LastGuideTime = QDateTime();
    m_GuiderIteration = 0;

    m_OffsetHistory.clear();
    m_VelocityHistory.clear();
    m_IsHarmonicDetected = false;
    m_HarmonicDetectionCounter = 0;

    m_Omega1 = 0.0;
    m_Omega2 = 0.0;
    m_LearnedT1 = 0.0;
    m_LearnedT2 = 0.0;
    m_PeriodsLearned = false;
    if (m_FFTEstimator) {
        m_FFTEstimator->reset();
    }
    m_Xhat = Eigen::VectorXd();
    m_XhatPred = Eigen::VectorXd();
    m_LastActiveR = -1.0;

    m_Solver.reset();
    m_Plant.reset();
    m_Network.reset();
}

void MPCGuider::setParameters(double Q, double R, double J, double Bf, double Kt, int N, double alpha)
{
    m_Q = Q;
    m_R = R;
    m_J = J;
    m_Bf = Bf;
    m_Kt = Kt;
    m_N = N;
    m_alpha = alpha;

    // Force re-initialization with new parameters on next guide tick
    m_Initialized = false;
}

double MPCGuider::guide(double offset)
{
    const QDateTime now = QDateTime::currentDateTime();
    double dt = 2.0; // Standard fallback (seconds)

    if (m_LastGuideTime.isValid())
    {
        constexpr int MAX_GUIDE_LAG = 30; // seconds
        const double seconds = m_LastGuideTime.msecsTo(now) / 1000.0;
        if (seconds > 0.1 && seconds < MAX_GUIDE_LAG)
        {
            dt = seconds;
        }
        else if (seconds > MAX_GUIDE_LAG || seconds < 0.0)
        {
            reset();
        }
    }
    m_LastGuideTime = now;

    // State reconstruction: delta_theta & velocity
    double delta_theta = m_HasPrevOffset ? (offset - m_PrevOffset) : 0.0;
    double velocity = delta_theta / dt;
    double delta_omega = m_HasPrevVelocity ? (velocity - m_PrevVelocity) : 0.0;

    // Save states for next iteration
    m_PrevOffset = offset;
    m_PrevVelocity = velocity;
    m_HasPrevOffset = true;
    m_HasPrevVelocity = true;

    // Store history for passive harmonic drive auto-detection
    m_OffsetHistory.push_back(offset);
    m_VelocityHistory.push_back(velocity);
    if (m_OffsetHistory.size() > 50)
    {
        m_OffsetHistory.erase(m_OffsetHistory.begin());
        m_VelocityHistory.erase(m_VelocityHistory.begin());
    }

    // Passive learning classification: ratio of rate of error (velocity) RMS to absolute error RMS
    if (m_OffsetHistory.size() >= 40 && !m_IsHarmonicDetected)
    {
        double sum_sq_offset = 0.0;
        double sum_sq_velocity = 0.0;
        for (size_t i = 0; i < m_OffsetHistory.size(); ++i)
        {
            sum_sq_offset += m_OffsetHistory[i] * m_OffsetHistory[i];
            sum_sq_velocity += m_VelocityHistory[i] * m_VelocityHistory[i];
        }
        double rms_offset = std::sqrt(sum_sq_offset / m_OffsetHistory.size());
        double rms_velocity = std::sqrt(sum_sq_velocity / m_OffsetHistory.size());
        double ratio = rms_velocity / (rms_offset + 1e-3);

        // ratio > 0.08 indicates dominant periodic error < 80 seconds (e.g. harmonic drives)
        // rms_offset > 0.15" ensures we are correcting substantial PE rather than pure small seeing noise
        if (ratio > 0.08 && rms_offset > 0.15)
        {
            m_HarmonicDetectionCounter++;
            if (m_HarmonicDetectionCounter >= 10) // Must persist for 10 frames
            {
                m_IsHarmonicDetected = true;
                qCDebug(KSTARS_EKOS_GUIDE) << QString("[MPCGuider %1] Passive Auto-Detection: Harmonic Drive mount signature identified (ratio=%2, RMS=%3\"). Dynamically raising control penalty R from %4 to 1.0 to prevent resonance.")
                                           .arg(m_ID).arg(ratio, 0, 'f', 3).arg(rms_offset, 0, 'f', 2).arg(m_R);
            }
        }
        else
        {
            m_HarmonicDetectionCounter = std::max(0, m_HarmonicDetectionCounter - 1);
        }
    }

    double activeR = m_R;
    if (activeR != m_LastActiveR)
    {
        m_Initialized = false; // Force dynamic rebuild with the new active R
        m_LastActiveR = activeR;
    }

    // Active learning frequency estimator update
    if (!m_FFTEstimator)
    {
        m_FFTEstimator = std::make_unique<FFTPeriodEstimator>();
    }
    double u_cum = m_Solver ? m_Solver->getCurrentU() : 0.0;
    double openLoopError = offset - u_cum;
    double elapsedSec = m_GuiderIteration * dt;
    m_FFTEstimator->addDataPoint(elapsedSec, openLoopError);

    if (m_GuiderIteration >= 40 && m_GuiderIteration % 10 == 0)
    {
        double t1 = 0.0, t2 = 0.0;
        // Gate: require totalTime > 200s before trusting FFT output.
        // This prevents spurious period detection on short pure-drift runs.
        bool ok = m_FFTEstimator->estimatePeriods(100.0, t1, t2);
        qDebug() << "estimatePeriods returned" << ok << "t1=" << t1 << "t2=" << t2;
        if (ok)
        {
            double alpha_freq = 0.1;
            if (m_LearnedT1 < 1.0)
            {
                m_LearnedT1 = t1;
                if (t2 > 1.0) m_LearnedT2 = t2;
            }
            else
            {
                m_LearnedT1 = alpha_freq * t1 + (1.0 - alpha_freq) * m_LearnedT1;
                if (t2 > 1.0) m_LearnedT2 = alpha_freq * t2 + (1.0 - alpha_freq) * m_LearnedT2;
            }

            double new_omega1 = 2.0 * M_PI / m_LearnedT1;
            double new_omega2 = (m_LearnedT2 > 1.0) ? 2.0 * M_PI / m_LearnedT2 : 0.0;

            bool shouldRebuild = (m_Omega1 < 1e-6)
                               || (std::abs(new_omega1 - m_Omega1) > 0.02 * m_Omega1)
                               || (new_omega2 > 1e-6 && m_Omega2 < 1e-6)
                               || (new_omega2 < 1e-6 && m_Omega2 > 1e-6)
                               || (new_omega2 > 1e-6 && std::abs(new_omega2 - m_Omega2) > 0.02 * m_Omega2);
            if (shouldRebuild)
            {
                m_Omega1 = new_omega1;
                m_Omega2 = new_omega2;
                m_PeriodsLearned = true;
                m_Initialized = false; // Force dynamic rebuild with new frequencies
                qCDebug(KSTARS_EKOS_GUIDE) << QString("[MPCGuider %1] Active Learning: Learned dominant periods T1=%2s, T2=%3s. Rebuilding solver.")
                                           .arg(m_ID).arg(m_LearnedT1, 0, 'f', 1).arg(m_LearnedT2, 0, 'f', 1);
            }
        }
    }

    // Dynamic rebuild of MPC matrices if dt, parameters, active R or frequencies changed
    if (!m_Initialized || std::abs(dt - m_LastDt) > 0.05 * m_LastDt)
    {
        double oldU = m_Solver ? m_Solver->getCurrentU() : 0.0;
        m_Kt = m_Bf / dt; // Enforce unit step gain for pure integrator plant in guiding
        m_Plant = std::make_unique<TelescopePlant>(m_J, m_Bf, m_Kt, dt);
        if (m_Backlash > 0.0)
        {
            // MPCSolver picks this up via plant.getBacklash() in rebuildMatrices
            // and activates its punch-through path on commanded direction changes.
            m_Plant->setBacklash(m_Backlash);
        }

        if (m_PeriodsLearned && m_Omega1 > 0.0)
        {
            m_Plant->setDisturbanceFrequencies(m_Omega1, m_Omega2);
        }

        m_Network = std::make_unique<LaguerreNetwork>(m_N, m_alpha);
        m_Solver = std::make_unique<MPCSolver>(*m_Plant, *m_Network, m_Q, activeR);
        m_Solver->setCurrentU(oldU);
        m_LastDt = dt;
        m_Initialized = true;

        int nx = m_Plant->getOrder();
        if (nx == 7)
        {
            const Eigen::RowVectorXd& kx = m_Solver->getKx();
            qDebug() << QString("[MPCGuider %1] IMP Kx=[%2 %3 %4 %5 %6 %7 %8] Kr=%9 omega1=%10 oldU=%11")
                        .arg(m_ID).arg(kx(0),6,'f',4).arg(kx(1),6,'f',4).arg(kx(2),6,'f',4)
                        .arg(kx(3),6,'f',4).arg(kx(4),6,'f',4).arg(kx(5),6,'f',4).arg(kx(6),6,'f',4)
                        .arg(m_Solver->getKr(),6,'f',4).arg(m_Omega1,6,'f',4).arg(oldU,6,'f',3);
        }
        if (nx == 7 && m_Xhat.size() != 6)
        {
            // Fresh transition into IMP mode.
            // theta_motor = oldU; d1 = offset - oldU (zero initial observer error).
            // d1_dot seeded from velocity so the Kx[3] velocity feedforward is
            // correct on the first IMP frame instead of starting cold at zero.
            m_Xhat = Eigen::VectorXd::Zero(6);
            m_Xhat(0) = oldU;
            m_Xhat(2) = offset - oldU;
            m_Xhat(3) = velocity;
            m_XhatPred = m_Plant->getAd6() * m_Xhat;
        }
    }

    m_GuiderIteration++;

    Eigen::VectorXd x_aug;
    int nx = m_Plant->getOrder();

    if (nx == 7)
    {
        // --- Luenberger Observer Update ---
        const Eigen::RowVectorXd& Cd = m_Plant->getCd6();

        double est_y = Cd.dot(m_XhatPred);
        double err = offset - est_y;

        Eigen::VectorXd L = Eigen::VectorXd::Zero(6);
        L(0) = 0.02;                                                  // motor position (extremely slow to prevent PE absorption)
        L(1) = 0.002 / dt;                                            // motor velocity
        L(2) = 0.60;                                                  // 1st harmonic displacement
        L(3) = 0.60 * m_Omega1;                                       // 1st harmonic velocity
        L(4) = 0.30;                                                  // 2nd harmonic displacement
        L(5) = 0.30 * m_Omega2;                                       // 2nd harmonic velocity

        m_Xhat = m_XhatPred + L * err;

        // incremental state update
        Eigen::VectorXd delta_xhat;
        if (m_XhatPrev.size() != 6)
        {
            m_XhatPrev = m_Xhat;
            delta_xhat = Eigen::VectorXd::Zero(6);
        }
        else
        {
            delta_xhat = m_Xhat - m_XhatPrev;
            m_XhatPrev = m_Xhat;
        }

        x_aug = Eigen::VectorXd::Zero(7);
        x_aug.segment<6>(0) = delta_xhat;
        x_aug(0) = 0.0; // Clear motor position increment to break the unstable observer-controller loop
        x_aug(1) = 0.0; // Clear the noisy velocity state increment to ensure high-frequency stability
        x_aug(6) = offset;
    }
    else if (nx == 5)
    {
        x_aug = Eigen::VectorXd::Zero(5);
        x_aug << 0.0, 0.0, 0.0, 0.0, offset;
    }
    else
    {
        x_aug = Eigen::VectorXd::Zero(3);
        x_aug << 0.0, 0.0, offset;
    }

    // Setpoint is 0.0 (regulate to zero drift)
    double prev_u = m_Solver->getCurrentU();
    double current_u = m_Solver->computeDeltaU(x_aug, 0.0);
    double guideVal = -(current_u - prev_u);

    // Predict next state for the observer
    if (nx == 7)
    {
        double delta_u = current_u - prev_u;
        m_XhatPred = m_Plant->getAd6() * m_Xhat + m_Plant->getBd6() * delta_u;
    }

    // Apply minMove threshold
    QString comment;
    if (std::abs(guideVal) > 0 && std::abs(guideVal) < m_MinMove)
    {
        comment = QString("Value %1 < minMove %2").arg(guideVal, 0, 'f', 2).arg(m_MinMove, 0, 'f', 2);
        guideVal = 0.0;
    }

    if (m_IsHarmonicDetected && comment.isEmpty())
    {
        comment = "Harmonic auto-tuned (R=1.0)";
    }

    qDebug() << QString("MPCGuide(%1,%2) Q=%3 R=%4 dt=%5s dx=%6 d_omega=%7 offset=%8 --> %9: %10")
                               .arg(m_ID, 3).arg(m_GuiderIteration, 3).arg(m_Q, 4, 'f', 1).arg(activeR, 4, 'f', 2)
                               .arg(dt, 4, 'f', 2).arg(delta_theta, 5, 'f', 2).arg(delta_omega, 5, 'f', 2)
                               .arg(offset, 5, 'f', 2).arg(guideVal, 5, 'f', 2).arg(comment);

    return guideVal;
}
