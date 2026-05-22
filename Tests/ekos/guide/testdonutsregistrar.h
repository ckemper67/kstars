/*
    SPDX-FileCopyrightText: 2024 Christian Kemper <ckemper@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

#pragma once

#include <QObject>
#include <vector>

class TestDonutsRegistrar : public QObject
{
    Q_OBJECT

public:
    TestDonutsRegistrar();

private Q_SLOTS:
    void initTestCase();
    void cleanupTestCase();
    void testMathRecovery();
    void testOffCenterRotation();
    void testTranslationMagnitudeSweep();
    void testTranslationDirections();
    void testStarDensity();

private:
    using Buffer = std::vector<double>;

    static Buffer makeStarField(int w, int h);
    static Buffer transformBuffer(const Buffer &src, int w, int h, double dx, double dy, double dthetaDeg);
};
