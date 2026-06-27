#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Christian Kemper <ckemper@gmail.com>
# SPDX-License-Identifier: GPL-2.0-or-later
"""Generate citydb.sqlite from KStars TSV city files.

Usage (build-time):
    python3 generate-citydb.py --countries countryInfo.txt --strict \
        --sqlite citydb.sqlite --min-population 50000 \
        --i18n-file cities15000_i18n.tsv \
        citydb.tsv cities15000.tsv

Usage (l10n extraction):
    python3 generate-citydb.py --countries countryInfo.txt --extract-l10n \
        citydb.tsv cities15000.tsv
"""
import sys
import os
import argparse
import sqlite3

from kstars_geo import DMS_RE, load_countries, DISPLAY_NAME_OVERRIDES

_errors = 0
_strict = False


def _warn(msg):
    global _errors
    if _strict:
        _errors += 1
        print("ERROR: " + msg, file=sys.stderr)
    else:
        print("WARNING: " + msg, file=sys.stderr)


def main():
    global _strict

    parser = argparse.ArgumentParser(description="Generate city database from TSV files")
    parser.add_argument("inputs", nargs="+", metavar="cities.tsv",
                        help="one or more KStars 8/9/10-column TSV files")
    parser.add_argument("-o", metavar="output.sql", dest="output_file", default=None)
    parser.add_argument("--sqlite", metavar="output.sqlite", dest="sqlite_file", default=None)
    parser.add_argument("--countries", metavar="countryInfo.txt", dest="countries_file", default=None,
                        help="GeoNames countryInfo.txt for country-name validation")
    parser.add_argument("--min-population", type=int, default=0, metavar="N",
                        help="exclude TSV rows with population below N (default: 0, no filter)")
    parser.add_argument("--strict", action="store_true",
                        help="treat validation warnings as errors (non-zero exit)")
    parser.add_argument("--extract-l10n", action="store_true",
                        help="print xi18nc() translation string declarations")
    parser.add_argument("--i18n-file", metavar="cities_i18n.tsv", dest="i18n_file", default=None,
                        help="TSV file with GeoNames translations (geonameid, lang, name)")
    args = parser.parse_args()

    _strict = args.strict

    if not args.output_file and not args.sqlite_file and not args.extract_l10n:
        parser.error("at least one of -o, --sqlite, or --extract-l10n is required")

    iso_countries = {}
    known_countries = set()
    if args.countries_file:
        iso_countries, known_countries = load_countries(args.countries_file)

    DDL = [
        "DROP TABLE IF EXISTS city",
        "CREATE TABLE city ("
        " id INTEGER DEFAULT NULL PRIMARY KEY AUTOINCREMENT,"
        " Name TEXT DEFAULT NULL,"
        " Province TEXT DEFAULT NULL,"
        " Country TEXT DEFAULT NULL,"
        " Latitude TEXT DEFAULT NULL,"
        " Longitude TEXT DEFAULT NULL,"
        " TZ REAL DEFAULT NULL,"
        " TZRule TEXT DEFAULT NULL,"
        " Elevation REAL NOT NULL DEFAULT -10,"
        " GeonameId INTEGER DEFAULT NULL"
        " )",
        "CREATE INDEX IF NOT EXISTS idx_name_country ON city (Name, Country)",
        "DROP TABLE IF EXISTS city_i18n",
        "CREATE TABLE city_i18n ("
        " geoname_id INTEGER,"
        " lang TEXT,"
        " name TEXT,"
        " PRIMARY KEY (geoname_id, lang)"
        " )",
        # No separate index on geoname_id: the (geoname_id, lang) primary key
        # already indexes geoname_id as its leftmost column, which serves both
        # the geoname_id+lang join in CITY_I18N_SQL and any geoname_id-only lookup.
    ]
    INSERT = ("INSERT INTO city"
              " (Name, Province, Country, Latitude, Longitude, TZ, TZRule, Elevation, GeonameId)"
              " VALUES")

    sql_out = None
    conn = None
    cur = None
    try:
        if args.output_file:
            sql_out = open(args.output_file, 'w', encoding='utf-8')
            sql_out.write("BEGIN TRANSACTION;\n")
            for stmt in DDL:
                sql_out.write(stmt + ";\n")

        if args.sqlite_file:
            if os.path.exists(args.sqlite_file):
                os.remove(args.sqlite_file)
            conn = sqlite3.connect(args.sqlite_file)
            cur = conn.cursor()
            for stmt in DDL:
                cur.execute(stmt)
            conn.commit()
            cur.execute("BEGIN")

        seen_cities = {}
        l10n_cities = []
        l10n_regions = set()
        l10n_countries = set()
        count = 0

        for cities_file in args.inputs:
            with open(cities_file, 'r', encoding='utf-8') as f:
                for lineno, line in enumerate(f, 1):
                    fields = line.strip().split('\t')
                    if len(fields) not in (8, 9, 10):
                        if len(fields) > 1:
                            _warn(f"{cities_file}:{lineno}: skipping row with"
                                  f" unexpected field count {len(fields)}")
                        continue

                    name, province, country = fields[0], fields[1], fields[2]
                    country = iso_countries.get(country, country)
                    lat_s, lon_s, tz_s, rule, elev_s = fields[3], fields[4], fields[5], fields[6], fields[7]

                    if len(fields) >= 9 and fields[8]:
                        try:
                            pop = int(fields[8])
                        except ValueError:
                            pop = 0
                        if pop < args.min_population:
                            continue

                    if known_countries and country and country.lower() not in known_countries:
                        _warn(f"{cities_file}:{lineno}: {name!r}: unknown country"
                              f" name {country!r}")
                    if not DMS_RE.match(lat_s) or not DMS_RE.match(lon_s):
                        _warn(f"{cities_file}:{lineno}: skipping {name!r} ({country}):"
                              f" malformed lat/lon lat={lat_s!r} lon={lon_s!r}")
                        continue
                    try:
                        tz_val = float(tz_s)
                        elev_val = float(elev_s) if elev_s else -10.0
                    except ValueError as e:
                        _warn(f"{cities_file}:{lineno}: skipping {name!r} ({country}): {e}")
                        continue
                    if not (-12 <= tz_val <= 14):
                        _warn(f"{cities_file}:{lineno}: skipping {name!r} ({country}):"
                              f" TZ {tz_val} out of range [-12, 14]")
                        continue
                    if not (-500 <= elev_val <= 9000):
                        _warn(f"{cities_file}:{lineno}: skipping {name!r} ({country}):"
                              f" elevation {elev_val} out of range [-500, 9000]")
                        continue

                    geoname_id = int(fields[9]) if len(fields) == 10 and fields[9] else None

                    city_key = (name.lower(), province.lower(), country.lower())
                    if city_key in seen_cities:
                        continue
                    seen_cities[city_key] = geoname_id
                    row = (name, province, country, lat_s, lon_s, tz_val, rule, elev_val, geoname_id)

                    if args.extract_l10n:
                        l10n_cities.append((name, province, country))
                        if province:
                            l10n_regions.add((province, country))
                        if country:
                            l10n_countries.add(country)

                    if sql_out:
                        nm, pv, ct, la, lo, tz_v, ru, el, gid = row
                        nm = nm.replace("'", "''")
                        pv = pv.replace("'", "''")
                        ct = ct.replace("'", "''")
                        la = la.replace("'", "''")
                        lo = lo.replace("'", "''")
                        gid_s = str(gid) if gid is not None else "NULL"
                        sql_out.write(f"{INSERT} ('{nm}', '{pv}', '{ct}', '{la}', '{lo}',"
                                      f" {tz_v}, '{ru}', {el}, {gid_s});\n")
                    if cur:
                        cur.execute(f"{INSERT} (?,?,?,?,?,?,?,?,?)", row)

                    count += 1
                    if count % 1000 == 0:
                        if sql_out:
                            sql_out.write("COMMIT; BEGIN TRANSACTION;\n")
                        if conn:
                            conn.commit()
                            cur.execute("BEGIN")

        if args.i18n_file:
            with open(args.i18n_file, 'r', encoding='utf-8') as f:
                for line in f:
                    fields = line.strip('\r\n').split('\t')
                    if len(fields) >= 3:
                        try:
                            geoname_id = int(fields[0])
                        except ValueError:
                            continue
                        lang = fields[1]
                        translated_name = fields[2]
                        if sql_out:
                            lang_esc = lang.replace("'", "''")
                            name_esc = translated_name.replace("'", "''")
                            sql_out.write(f"INSERT OR REPLACE INTO city_i18n"
                                          f" (geoname_id, lang, name) VALUES"
                                          f" ({geoname_id}, '{lang_esc}', '{name_esc}');\n")
                        if cur:
                            cur.execute("INSERT OR REPLACE INTO city_i18n"
                                        " (geoname_id, lang, name) VALUES (?, ?, ?)",
                                        (geoname_id, lang, translated_name))

        if sql_out:
            sql_out.write("COMMIT;\n")
        if conn:
            conn.commit()
    finally:
        if sql_out:
            sql_out.close()
        if conn:
            conn.close()

    if args.extract_l10n:
        for nm, pv, ct in sorted(l10n_cities):
            place = "{0} {1}".format(pv, ct) if pv else ct
            print('xi18nc("City in {0}", "{1}")'.format(place, nm))
        for pv, ct in sorted(l10n_regions):
            print('xi18nc("Region/state in {0}", "{1}")'.format(ct, pv))
        for ct in sorted(l10n_countries):
            print('xi18nc("Country name", "{0}")'.format(ct))

    if count == 0:
        print("ERROR: no valid rows were produced -- check input files and warnings above",
              file=sys.stderr)
        sys.exit(1)

    if _errors > 0:
        print(f"ERROR: {_errors} validation error(s) in strict mode", file=sys.stderr)
        sys.exit(1)

    print(f"Generated {count} cities", file=sys.stderr)


if __name__ == "__main__":
    main()
