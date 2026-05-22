/*
    SPDX-FileCopyrightText: 2026 Christian Kemper <ckemper@gmail.com>

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

private slots:
    void initTestCase();
    void cleanupTestCase();

    void testMathRecovery();
    void testTranslationMagnitudeSweep();
    void testTranslationDirections();
    void testStarDensity();

private:
    using Buffer = std::vector<double>;

    static Buffer makeStarField(int w, int h, uint32_t seed = 42);
    static Buffer transformBuffer(const Buffer &src, int w, int h,
                                  double dx, double dy, double dthetaDeg = 0.0);
};
