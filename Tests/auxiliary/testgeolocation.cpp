/*  KStars class tests
    SPDX-FileCopyrightText: 2020 Eric Dejouhanet <eric.dejouhanet@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

#include "testgeolocation.h"
#include "auxiliary/geolocation.h"
#include "citydb_helpers.h"
#include "dms.h"

#include <QFile>
#include <QSqlDatabase>
#include <QSqlError>

TestGeolocation::TestGeolocation(QObject *parent) : QObject(parent)
{
}

void TestGeolocation::initTestCase()
{
}

void TestGeolocation::cleanupTestCase()
{
}

void TestGeolocation::init()
{
}

void TestGeolocation::cleanup()
{
}

// GeoLocation stores a db-sourced translation and returns it from translatedName().
// Uses distinct name ("Cologne") and translation ("Koeln") to prove translatedName()
// returns the stored translation, not the city name.
// Fails if: GeonameId/TranslatedName fields are removed, or translatedName() stops
// checking them and always falls through to i18nc.
void TestGeolocation::testTranslatedNameFromDb()
{
    dms lat(48.0), lng(8.0);
    GeoLocation loc(lng, lat, "Cologne", "North Rhine-Westphalia", "Germany",
                    1.0, nullptr, 44.0, true, 4, 2886242, "Koeln");

    QCOMPARE(loc.geonameId(), 2886242);
    QCOMPARE(loc.dbTranslatedName(), QString("Koeln"));
    // translatedName() must return the stored translation, not the city name "Cologne"
    QCOMPARE(loc.translatedName(), QString("Koeln"));
    QVERIFY(loc.translatedName() != loc.name());
}

// When GeonameId is 0, translatedName() must fall back to i18nc, not return
// the stored translation string.
// Fails if: the geonameId > 0 guard is removed from translatedName().
void TestGeolocation::testTranslatedNameFallback()
{
    dms lat(48.0), lng(8.0);
    GeoLocation withId(lng, lat, "TestCity", "TestProv", "TestCountry",
                       1.0, nullptr, 0.0, true, 4, 42, "TranslatedName");
    GeoLocation noId(lng, lat, "TestCity", "TestProv", "TestCountry",
                     1.0, nullptr, 0.0, true, 4, 0, "TranslatedName");

    // geonameId > 0: stored translation wins
    QCOMPARE(withId.translatedName(), QString("TranslatedName"));

    // geonameId == 0: must NOT return the stored string -- falls back to i18nc
    QVERIFY(noId.translatedName() != QString("TranslatedName"));
    QVERIFY(!noId.translatedName().isEmpty());
}

// Helper: open citydb, run the i18n JOIN query for one city, return a GeoLocation.
// Caller owns the returned pointer.
static GeoLocation *queryOneCity(const QString &name, const QString &country,
                                 const QString &lang, const QString &baseLang,
                                 const QString &connectionName)
{
    GeoLocation *loc = nullptr;
    {
        QSqlDatabase db = QSqlDatabase::addDatabase("QSQLITE", connectionName);
        db.setDatabaseName(QString(KSTARS_CITYDB_PATH));
        if (db.open())
        {
            QSqlQuery q(db);
            q.prepare(QString(CITY_I18N_SQL) +
                      " WHERE c.Name = :name AND c.Country = :country AND c.GeonameId > 0");
            q.bindValue(":lang",      lang);
            q.bindValue(":base_lang", baseLang);
            q.bindValue(":name",      name);
            q.bindValue(":country",   country);
            if (q.exec() && q.next())
            {
                QMap<QString, TimeZoneRule> empty;
                loc = geoLocationFromCityRow(q, empty);
            }
            db.close();
        }
    }
    QSqlDatabase::removeDatabase(connectionName);
    return loc;
}

// Verify the city_i18n table is present in the built database.
// Fails if generate-citydb.py stops creating city_i18n, or if
// cities15000_i18n.tsv is not passed at build time.
void TestGeolocation::testCityI18nSchema()
{
    const QString dbPath = QString(KSTARS_CITYDB_PATH);
    QVERIFY2(QFile::exists(dbPath), qPrintable("citydb.sqlite not found at: " + dbPath));

    {
        QSqlDatabase db = QSqlDatabase::addDatabase("QSQLITE", "testgeolocation_schema");
        db.setDatabaseName(dbPath);
        QVERIFY2(db.open(), qPrintable("Failed to open: " + dbPath));
        QVERIFY2(db.tables().contains("city_i18n", Qt::CaseInsensitive),
                 "city_i18n table missing -- rebuild citydb.sqlite with --i18n-file");
        db.close();
    }
    QSqlDatabase::removeDatabase("testgeolocation_schema");
}

// End-to-end: loadCitiesFromDb() populates a list from the real citydb.sqlite.
// Verifies the integration loop in kstarsdata.cpp: table detection, query
// execution, row iteration, and GeoLocation construction all work together.
// Fails if: loadCitiesFromDb() is removed or broken, or Cologne disappears
// from citydb.tsv.
void TestGeolocation::testLoadCitiesFromDb()
{
    QList<GeoLocation *> cities;
    QMap<QString, TimeZoneRule> emptyRulebook;

    {
        QSqlDatabase db = QSqlDatabase::addDatabase("QSQLITE", "testgeolocation_load");
        db.setDatabaseName(QString(KSTARS_CITYDB_PATH));
        QVERIFY2(db.open(), "Failed to open citydb.sqlite");
        bool loaded = loadCitiesFromDb(db, "de_de", "de", cities, emptyRulebook);
        db.close();
        QVERIFY2(loaded, "loadCitiesFromDb returned false -- no cities loaded");
    }
    QSqlDatabase::removeDatabase("testgeolocation_load");

    QVERIFY2(!cities.isEmpty(), "City list is empty after loadCitiesFromDb");

    GeoLocation *cologne = nullptr;
    for (GeoLocation *loc : cities)
    {
        if (loc->name() == "Cologne" && loc->country() == "Germany")
        {
            cologne = loc;
            break;
        }
    }

    QVERIFY2(cologne != nullptr, "Cologne not found in loaded city list");
    QCOMPARE(cologne->geonameId(), 2886242);
    QCOMPARE(cologne->translatedName(), QString("Köln")); // Köln

    qDeleteAll(cities);
}

// loadCitiesFromDb must report failure (return false) when the query executes
// but yields no rows -- the caller treats false as "no city data" and aborts.
// Uses an in-memory db with an empty city table and no city_i18n, so the plain
// SELECT path runs and returns zero rows.
// Fails if: loadCitiesFromDb returns true on an empty result.
void TestGeolocation::testLoadCitiesFromDbEmptyReturnsFalse()
{
    QList<GeoLocation *> cities;
    QMap<QString, TimeZoneRule> emptyRulebook;
    bool loaded = true;

    {
        QSqlDatabase db = QSqlDatabase::addDatabase("QSQLITE", "testgeolocation_empty");
        db.setDatabaseName(":memory:");
        QVERIFY2(db.open(), "Failed to open in-memory db");
        QSqlQuery create(db);
        QVERIFY2(create.exec("CREATE TABLE city (id INTEGER, Name TEXT, Province TEXT, "
                             "Country TEXT, Latitude TEXT, Longitude TEXT, TZ REAL, "
                             "TZRule TEXT, Elevation REAL, GeonameId INTEGER)"),
                 qPrintable(create.lastError().text()));
        loaded = loadCitiesFromDb(db, "de_de", "de", cities, emptyRulebook);
        db.close();
    }
    QSqlDatabase::removeDatabase("testgeolocation_empty");

    QVERIFY2(!loaded, "loadCitiesFromDb returned true for an empty city table");
    QVERIFY2(cities.isEmpty(), "No GeoLocation should be appended for an empty table");
}

// Exercises the exact-locale LEFT JOIN in CITY_I18N_SQL: a city with both an
// exact-locale (de_de) and a base-language (de) translation must resolve to the
// exact one, while a city with only a base translation falls back to it. Uses an
// in-memory db because the shipped city_i18n carries only base-language codes, so
// the exact-locale join is otherwise never covered.
// Fails if: the ie/ib joins are swapped, the exact-over-base preference is
// inverted, or the translation columns move out of positions 10/11.
void TestGeolocation::testI18nExactLocaleBeatsBase()
{
    QList<GeoLocation *> cities;
    QMap<QString, TimeZoneRule> rulebook;

    {
        QSqlDatabase db = QSqlDatabase::addDatabase("QSQLITE", "testgeolocation_i18n");
        db.setDatabaseName(":memory:");
        QVERIFY2(db.open(), "Failed to open in-memory db");

        QSqlQuery q(db);
        QVERIFY2(q.exec("CREATE TABLE city (id INTEGER PRIMARY KEY, Name TEXT, Province TEXT, "
                        "Country TEXT, Latitude TEXT, Longitude TEXT, TZ REAL, TZRule TEXT, "
                        "Elevation REAL, GeonameId INTEGER)"),
                 qPrintable(q.lastError().text()));
        QVERIFY2(q.exec("CREATE TABLE city_i18n (geoname_id INTEGER, lang TEXT, name TEXT, "
                        "PRIMARY KEY (geoname_id, lang))"),
                 qPrintable(q.lastError().text()));
        QVERIFY2(q.exec("INSERT INTO city (Name, Province, Country, Latitude, Longitude, TZ, "
                        "TZRule, Elevation, GeonameId) VALUES "
                        "('Exactville','','Testland','0','0',1.0,'--',100.0,4242),"
                        "('Baseville','','Testland','0','0',1.0,'--',100.0,4243)"),
                 qPrintable(q.lastError().text()));
        QVERIFY2(q.exec("INSERT INTO city_i18n (geoname_id, lang, name) VALUES "
                        "(4242,'de_de','ExactDE'),(4242,'de','BaseDE'),"
                        "(4243,'de','OnlyBaseDE')"),
                 qPrintable(q.lastError().text()));

        QVERIFY2(loadCitiesFromDb(db, "de_de", "de", cities, rulebook),
                 "loadCitiesFromDb returned false");
        db.close();
    }
    QSqlDatabase::removeDatabase("testgeolocation_i18n");

    QString exact, base;
    for (GeoLocation *loc : cities)
    {
        if (loc->name() == "Exactville")
            exact = loc->dbTranslatedName();
        else if (loc->name() == "Baseville")
            base = loc->dbTranslatedName();
    }
    QCOMPARE(exact, QString("ExactDE"));   // exact-locale join wins over base
    QCOMPARE(base, QString("OnlyBaseDE")); // base-language join used when no exact

    qDeleteAll(cities);
}

// Exercises the plain-SELECT path: geoLocationFromCityRow on a real 10-column
// "SELECT * FROM city" result must read geonameId from column 9 and leave the
// translation empty (no city_i18n columns are present).
// Fails if: geonameId is read from the wrong column, or geoLocationFromCityRow
// populates a db translation when no translation columns exist.
// Note: this does NOT guard the cols >= 12 branch itself -- Qt returns an empty
// QVariant for out-of-range value(10)/value(11), so removing that guard would
// only emit a warning, not change any result asserted here.
void TestGeolocation::testPlainQueryNoTranslation()
{
    QSqlQuery q;
    {
        QSqlDatabase db = QSqlDatabase::addDatabase("QSQLITE", "testgeolocation_plain");
        db.setDatabaseName(QString(KSTARS_CITYDB_PATH));
        QVERIFY2(db.open(), "Failed to open citydb.sqlite");

        q = QSqlQuery(db);
        QVERIFY2(q.exec("SELECT * FROM city WHERE Name = 'Cologne' AND Country = 'Germany' "
                         "AND GeonameId > 0"),
                 qPrintable(q.lastError().text()));
        QVERIFY2(q.next(), "Cologne not found via plain SELECT");

        // 10 columns: cols >= 10 so geonameId is set; cols < 12 so no translation
        QCOMPARE(q.record().count(), 10);

        QMap<QString, TimeZoneRule> emptyRulebook;
        GeoLocation *cologne = geoLocationFromCityRow(q, emptyRulebook);
        QVERIFY(cologne != nullptr);

        // geonameId is populated from col 9
        QCOMPARE(cologne->geonameId(), 2886242);

        // no translation columns -- dbTranslatedName must be empty
        QVERIFY(cologne->dbTranslatedName().isEmpty());

        // translatedName() must NOT return a db value; falls back to i18nc
        QVERIFY(cologne->translatedName() != QString("Köln")); // Köln

        delete cologne;
        db.close();
    }
    QSqlDatabase::removeDatabase("testgeolocation_plain");
}

// Cologne with a German locale -- base-lang fallback ("de" matches, "de_de" does not).
// Fails if: citydb.tsv loses GeonameId 2886242, or translatedName() stops
// returning the db-sourced value.
void TestGeolocation::testCologneGermanLocale()
{
    GeoLocation *cologne = queryOneCity("Cologne", "Germany", "de_de", "de",
                                       "testgeolocation_cologne");
    QVERIFY2(cologne != nullptr, "Cologne/Germany not found in citydb.sqlite");
    QCOMPARE(cologne->geonameId(), 2886242);
    QCOMPARE(cologne->translatedName(), QString("Köln")); // Köln
    delete cologne;
}

// Tokyo queried with a Japanese locale should return the kanji name.
// Fails if: city_i18n loses the ja row for geoname_id 1850147, or
// geoLocationFromCityRow stops picking up the base-lang fallback.
void TestGeolocation::testTokyoJapaneseLocale()
{
    GeoLocation *tokyo = queryOneCity("Tokyo", "Japan", "ja_jp", "ja",
                                     "testgeolocation_tokyo");
    QVERIFY2(tokyo != nullptr, "Tokyo/Japan not found in citydb.sqlite");
    QCOMPARE(tokyo->geonameId(), 1850147);
    QCOMPARE(tokyo->translatedName(), QString("東京都")); // 東京都
    delete tokyo;
}

// Dubai queried with an Arabic locale should return the Arabic name.
// Fails if: city_i18n loses the ar row for geoname_id 292223, or
// geoLocationFromCityRow stops picking up the base-lang fallback.
void TestGeolocation::testDubaiArabicLocale()
{
    GeoLocation *dubai = queryOneCity("Dubai", "United Arab Emirates", "ar_ae", "ar",
                                     "testgeolocation_dubai");
    QVERIFY2(dubai != nullptr, "Dubai/UAE not found in citydb.sqlite");
    QCOMPARE(dubai->geonameId(), 292223);
    QCOMPARE(dubai->translatedName(), QString("دبي")); // دبي
    delete dubai;
}

QTEST_GUILESS_MAIN(TestGeolocation)
