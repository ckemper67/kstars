/*
    SPDX-FileCopyrightText: 2024 Christian Kemper <ckemper@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

#include "testdonutsguider.h"
#include "ekos/guide/internalguide/donutsguider.h"
#include "fitsviewer/fitsdata.h"
#include <QTest>
#include <QSignalSpy>
#include <cmath>

TestDonutsGuider::TestDonutsGuider() : QObject()
{
}

void TestDonutsGuider::initTestCase()
{
}

void TestDonutsGuider::cleanupTestCase()
{
}

// Simple bilinear interpolation to rotate/translate image
QSharedPointer<FITSData> TestDonutsGuider::transformImage(const QSharedPointer<FITSData> &source, double dx, double dy, double dtheta)
{
    int w = source->width();
    int h = source->height();
    QSharedPointer<FITSData> target(new FITSData(FITS_NORMAL));
    // Allocate buffer manually for simplicity in test
    // FITSData usually allocates on load, but we can set properties
    // Actually, FITSData is a bit complex to allocate from scratch.
    // Let's just clone and overwrite the buffer.
    target.reset(new FITSData(source));
    
    uint16_t *srcBuf = const_cast<uint16_t *>(reinterpret_cast<const uint16_t *>(source->getImageBuffer()));
    uint16_t *dstBuf = const_cast<uint16_t *>(reinterpret_cast<const uint16_t *>(target->getImageBuffer()));
    
    double cosT = cos(dtheta * M_PI / 180.0);
    double sinT = sin(dtheta * M_PI / 180.0);
    double cx = w / 2.0;
    double cy = h / 2.0;
    
    for (int y = 0; y < h; ++y)
    {
        for (int x = 0; x < w; ++x)
        {
            // Reverse transform to find source pixel
            // 1. Translation
            double x1 = x - dx;
            double y1 = y - dy;
            
            // 2. Rotation around center
            double x2 = (x1 - cx) * cosT + (y1 - cy) * sinT + cx;
            double y2 = -(x1 - cx) * sinT + (y1 - cy) * cosT + cy;
            
            if (x2 >= 0 && x2 < w - 1 && y2 >= 0 && y2 < h - 1)
            {
                int ix = floor(x2);
                int iy = floor(y2);
                double fx = x2 - ix;
                double fy = y2 - iy;
                
                double v00 = srcBuf[iy * w + ix];
                double v10 = srcBuf[iy * w + ix + 1];
                double v01 = srcBuf[(iy + 1) * w + ix];
                double v11 = srcBuf[(iy + 1) * w + ix + 1];
                
                double val = (1 - fx) * (1 - fy) * v00 + fx * (1 - fy) * v10 + (1 - fx) * fy * v01 + fx * fy * v11;
                dstBuf[y * w + x] = static_cast<uint16_t>(val);
            }
            else
            {
                dstBuf[y * w + x] = 0;
            }
        }
    }
    
    return target;
}

void TestDonutsGuider::testMathRecovery()
{
    QString fitsPath = QFINDTESTDATA("../../fitsviewer/m47_sim_stars.fits");
    if (fitsPath.isEmpty())
        QSKIP("M47 sim stars fits not found");
        
    QSharedPointer<FITSData> master(new FITSData());
    master->loadFromFile(fitsPath).waitForFinished();
    
    Ekos::DonutsGuider guider;
    guider.setReference(master);
    
    // Test 1: Pure translation
    double tx = 5.4, ty = -3.2;
    auto shifted = transformImage(master, tx, ty, 0);
    Ekos::DonutsGuider::Transform res1 = guider.calculateTransform(shifted);
    
    qDebug() << "Test 1 (Translation): Expected (" << tx << "," << ty << ") Got (" << res1.dx << "," << res1.dy << ")";
    QVERIFY(std::abs(res1.dx - tx) < 0.1);
    QVERIFY(std::abs(res1.dy - ty) < 0.1);
    
    // Test 2: Pure rotation
    double rot = 0.5; // degrees
    auto rotated = transformImage(master, 0, 0, rot);
    Ekos::DonutsGuider::Transform res2 = guider.calculateTransform(rotated);
    
    qDebug() << "Test 2 (Rotation): Expected" << rot << "deg Got" << res2.dtheta * 180.0 / M_PI << "deg";
    // For rotation, error might be higher due to simplified solver and bilinear artifacts
    QVERIFY(std::abs(res2.dtheta * 180.0 / M_PI - rot) < 0.1);

    // Test 3: Combined
    auto combined = transformImage(master, 2.0, 2.0, 0.2);
    Ekos::DonutsGuider::Transform res3 = guider.calculateTransform(combined);
    qDebug() << "Test 3 (Combined): Expected (2.0, 2.0, 0.2) Got (" << res3.dx << "," << res3.dy << "," << res3.dtheta * 180.0 / M_PI << ")";
    QVERIFY(std::abs(res3.dx - 2.0) < 0.1);
    QVERIFY(std::abs(res3.dy - 2.0) < 0.1);
    QVERIFY(std::abs(res3.dtheta * 180.0 / M_PI - 0.2) < 0.05);
}

QTEST_GUILESS_MAIN(TestDonutsGuider)
