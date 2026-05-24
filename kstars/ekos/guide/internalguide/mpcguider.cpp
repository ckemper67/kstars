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

    double activeR = m_IsHarmonicDetected ? std::max(m_R, 1.0) : m_R;
    static double lastActiveR = -1.0;
    if (activeR != lastActiveR)
    {
        m_Initialized = false; // Force dynamic rebuild with the new active R
        lastActiveR = activeR;
    }

    // Dynamic rebuild of MPC matrices if dt, parameters or active R changed
    if (!m_Initialized || std::abs(dt - m_LastDt) > 0.05 * m_LastDt)
    {
        m_Kt = m_Bf / dt; // Enforce unit step gain for pure integrator plant in guiding
        m_Plant = std::make_unique<TelescopePlant>(m_J, m_Bf, m_Kt, dt);
        m_Network = std::make_unique<LaguerreNetwork>(m_N, m_alpha);
        m_Solver = std::make_unique<MPCSolver>(*m_Plant, *m_Network, m_Q, activeR);
        m_LastDt = dt;
        m_Initialized = true;
    }

    m_GuiderIteration++;

    Eigen::VectorXd x_aug = Eigen::VectorXd::Zero(3);
    x_aug << 0.0, 0.0, offset;

    // Setpoint is 0.0 (regulate to zero drift)
    double prev_u = m_Solver->getCurrentU();
    double current_u = m_Solver->computeDeltaU(x_aug, 0.0);
    double guideVal = -(current_u - prev_u);

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

    qCDebug(KSTARS_EKOS_GUIDE) << QString("MPCGuide(%1,%2) Q=%3 R=%4 dt=%5s dx=%6 d_omega=%7 offset=%8 --> %9: %10")
                               .arg(m_ID, 3).arg(m_GuiderIteration, 3).arg(m_Q, 4, 'f', 1).arg(activeR, 4, 'f', 2)
                               .arg(dt, 4, 'f', 2).arg(delta_theta, 5, 'f', 2).arg(delta_omega, 5, 'f', 2)
                               .arg(offset, 5, 'f', 2).arg(guideVal, 5, 'f', 2).arg(comment);

    return guideVal;
}
