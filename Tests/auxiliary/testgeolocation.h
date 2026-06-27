/*  KStars class tests
    SPDX-FileCopyrightText: 2020 Eric Dejouhanet <eric.dejouhanet@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

#ifndef TESTGEOLOCATION_H
#define TESTGEOLOCATION_H

#include <QtGlobal>
#if QT_VERSION >= QT_VERSION_CHECK(6, 0, 0)
#include <QtTest/QTest>
#else
#include <QTest>
#endif

#include <QObject>

class TestGeolocation : public QObject
{
        Q_OBJECT
    public:
        explicit TestGeolocation(QObject *parent = nullptr);

    private Q_SLOTS:
        void initTestCase();
        void cleanupTestCase();

        void init();
        void cleanup();

        // GeoLocation unit tests -- no database required
        void testTranslatedNameFromDb();
        void testTranslatedNameFallback();

        // End-to-end: real citydb.sqlite, real query, real GeoLocation
        void testCityI18nSchema();
        void testLoadCitiesFromDb();
        void testLoadCitiesFromDbEmptyReturnsFalse();
        void testI18nExactLocaleBeatsBase();
        void testPlainQueryNoTranslation();
        void testCologneGermanLocale();
        void testTokyoJapaneseLocale();
        void testDubaiArabicLocale();
};

#endif // TESTGEOLOCATION_H
