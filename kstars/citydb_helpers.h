/*
    SPDX-FileCopyrightText: 2026 Christian Kemper <ckemper@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

#pragma once

// Internal helpers for city database loading shared between kstarsdata.cpp
// and the geolocation test.  Do not include from other production code.

#include "auxiliary/geolocation.h"
#include "timezonerule.h"

#include <QList>
#include <QMap>
#include <QSqlDatabase>
#include <QSqlQuery>
#include <QSqlRecord>
#include <QString>

// SQL query that joins city with city_i18n to fetch translated names.
// Bind :lang to the full locale (e.g. "de_de") and :base_lang to the
// language code (e.g. "de") before executing.
// Each LEFT JOIN matches at most one row -- city_i18n has PRIMARY KEY
// (geoname_id, lang) -- so there is no fan-out and a missing translation
// yields NULL, exactly as the previous correlated subqueries did. The same
// (geoname_id, lang) primary-key index serves both joins. Callers may append
// a WHERE clause on the c.* columns after this string.
static const char CITY_I18N_SQL[] =
    "SELECT c.id, c.Name, c.Province, c.Country, c.Latitude, c.Longitude, "
    "c.TZ, c.TZRule, c.Elevation, c.GeonameId, "
    "ie.name AS TranslationExact, "
    "ib.name AS TranslationBase "
    "FROM city c "
    "LEFT JOIN city_i18n ie ON ie.geoname_id = c.GeonameId AND ie.lang = :lang "
    "LEFT JOIN city_i18n ib ON ib.geoname_id = c.GeonameId AND ib.lang = :base_lang";

// Construct a GeoLocation from the current row of a city query.
// Works for both the i18n query (12 columns) and the plain fallback
// query (10 columns, no translation).  The caller owns the returned pointer.
// Pass a Rulebook so TZRule strings can be resolved; an unknown key auto-inserts
// a default empty rule, so the resulting TZrule is never null.
// readOnly is forwarded to GeoLocation: true for the shipped system database,
// false for the user's editable local database (mycitydb).
inline GeoLocation *geoLocationFromCityRow(const QSqlQuery &q,
                                           QMap<QString, TimeZoneRule> &rulebook,
                                           bool readOnly = true)
{
    const int cols = q.record().count();

    const QString name     = q.value(1).toString();
    const QString province = q.value(2).toString();
    const QString country  = q.value(3).toString();
    const dms     lat      = dms(q.value(4).toString());
    const dms     lng      = dms(q.value(5).toString());
    const double  TZ       = q.value(6).toDouble();
    // operator[] auto-inserts a default (empty, no-DST) TimeZoneRule for an
    // unknown key, so TZrule is never null. This matches the legacy
    // readCityData() behavior; downstream code dereferences tzrule() without
    // null checks, so returning nullptr here would crash on city selection.
    TimeZoneRule *TZrule   = &rulebook[q.value(7).toString()];
    const double  elevation = q.value(8).toDouble();

    int     geonameId     = 0;
    QString translatedName;

    if (cols >= 10)
        geonameId = q.value(9).toInt();

    if (cols >= 12)
    {
        const QString exact = q.value(10).toString();
        const QString base  = q.value(11).toString();
        translatedName = !exact.isEmpty() ? exact : base;
    }

    return new GeoLocation(lng, lat, name, province, country,
                           TZ, TZrule, elevation, readOnly, 4,
                           geonameId, translatedName);
}

// Load cities from an open database into out, using i18n translations for the
// given lang/baseLang if the city_i18n table is present, otherwise falling back
// to a plain SELECT.  Returns true if at least one city was loaded.
// The caller owns the GeoLocation objects appended to out.
inline bool loadCitiesFromDb(QSqlDatabase &db,
                             const QString &lang,
                             const QString &baseLang,
                             QList<GeoLocation *> &out,
                             QMap<QString, TimeZoneRule> &rulebook)
{
    QSqlQuery q(db);
    bool success = false;

    if (db.tables().contains("city_i18n", Qt::CaseInsensitive))
    {
        q.prepare(CITY_I18N_SQL);
        q.bindValue(":lang",      lang);
        q.bindValue(":base_lang", baseLang);
        success = q.exec();
    }

    if (!success)
        success = q.exec("SELECT * FROM city");

    if (!success)
        return false;

    bool anyLoaded = false;
    while (q.next())
    {
        out.append(geoLocationFromCityRow(q, rulebook));
        anyLoaded = true;
    }

    // Contract: true only if at least one city was loaded. A query that executed
    // against an empty city table is still a load failure for the caller, which
    // treats a false return as "no city data" -- matching the legacy readCityData.
    return anyLoaded;
}
