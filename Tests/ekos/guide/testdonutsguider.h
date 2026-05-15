/*
    SPDX-FileCopyrightText: 2024 Christian Kemper <ckemper@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

#pragma once

#include <QObject>
#include <QSharedPointer>
#include <vector>

class FITSData;

class TestDonutsGuider : public QObject
{
    Q_OBJECT

public:
    TestDonutsGuider();

private Q_SLOTS:
    void initTestCase();
    void cleanupTestCase();
    void testMathRecovery();
    void testTranslationMagnitudeSweep();
    void testTranslationDirections();
    void testStarDensity();

private:
    // FITSData-based helper used by testMathRecovery
    QSharedPointer<FITSData> transformImage(const QSharedPointer<FITSData> &source, double dx, double dy, double dtheta);

    // Synthetic double-buffer helpers used by the translation tests
    using Buffer = std::vector<double>;
    static Buffer makeStarField(int w, int h);
    static Buffer transformBuffer(const Buffer &src, int w, int h, double dx, double dy, double dthetaDeg);
};
