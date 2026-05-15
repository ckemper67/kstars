/*
    SPDX-FileCopyrightText: 2024 Christian Kemper <ckemper@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

#pragma once

#include <QObject>
#include <QSharedPointer>

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

private:
    QSharedPointer<FITSData> transformImage(const QSharedPointer<FITSData> &source, double dx, double dy, double dtheta);
};
