/*
    SPDX-FileCopyrightText: 2026 Christian Kemper <ckemper@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

#pragma once

#include <QObject>

class TestGuideAlgorithms : public QObject
{
    Q_OBJECT

public:
    TestGuideAlgorithms();

private slots:
    void initTestCase();
    void cleanupTestCase();

    // LinearGuider
    void testLinearGuiderConstantDrift();
    void testLinearGuiderStepResponse();
    void testLinearGuiderPE();

    // HysteresisGuider
    void testHysteresisGuiderConstantDrift();
    void testHysteresisGuiderDampening();
    void testHysteresisGuiderPE();

    // GPG (GaussianProcessGuider direct, no Options dependency)
    void testGPGConstantDrift();
    void testGPGHarmonicDrivePE();
    void testGPGPeriodDetection();
    void testGPGMultiHarmonic();
    void testGPGSeeingNoise();

    // MPCGuider
    void testMPCGuiderConstantDrift();
    void testMPCGuiderPE();
};
