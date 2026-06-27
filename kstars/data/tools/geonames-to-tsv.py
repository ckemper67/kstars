#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Christian Kemper <ckemper@gmail.com>
# SPDX-License-Identifier: GPL-2.0-or-later
"""Convert a raw GeoNames cities file to KStars 10-column TSV format.

Column 3 is the ISO 3166-1 alpha-2 country code (taken straight from the
GeoNames country-code field). generate-citydb.py resolves it to a display
name at build time, so both citydb.tsv and cities15000.tsv share one
unambiguous country representation.

Usage:
    python3 geonames-to-tsv.py admin1CodesASCII.txt timeZones.txt \
        /path/to/cities15000.txt > cities15000.tsv
"""
import argparse

from kstars_geo import (deg_to_dms, get_tz_rule, load_admin1,
                        load_dst_observed, load_timezones)


def main():
    parser = argparse.ArgumentParser(
        description="Convert raw GeoNames cities file to KStars 10-column TSV")
    parser.add_argument("admin1_file", metavar="admin1CodesASCII.txt")
    parser.add_argument("tz_file", metavar="timeZones.txt")
    parser.add_argument("cities_file", metavar="cities15000.txt")
    args = parser.parse_args()

    admin1 = load_admin1(args.admin1_file)
    timezones = load_timezones(args.tz_file)
    dst_observed = load_dst_observed(args.tz_file)

    rows = []
    with open(args.cities_file, 'r', encoding='utf-8') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 19:
                continue
            try:
                population = int(fields[14])
            except ValueError:
                population = 0

            country_code = fields[8]
            admin1_code = fields[10]
            lat = float(fields[4])
            lng = float(fields[5])
            elevation = fields[15]
            iana_tz = fields[17]

            if not elevation or elevation == "":
                elevation = "-10"
            province_key = f"{country_code}.{admin1_code}"
            province = admin1.get(province_key, "")
            lat_dms = deg_to_dms(lat, True)
            lng_dms = deg_to_dms(lng, False)
            tz_offset = timezones.get(iana_tz, "0.0")
            tz_rule = get_tz_rule(iana_tz, country_code,
                                  dst_observed.get(iana_tz))

            city_name = fields[2]  # asciiname (always ASCII)
            geoname_id = fields[0]

            rows.append(f"{city_name}\t{province}\t{country_code}\t{lat_dms}\t{lng_dms}\t{tz_offset}\t{tz_rule}\t{elevation}\t{population}\t{geoname_id}")

    rows.sort(key=lambda r: tuple(r.split('\t')[i].lower() for i in (2, 1, 0)))
    for row in rows:
        print(row)


if __name__ == "__main__":
    main()
