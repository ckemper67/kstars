#!/usr/bin/env python3
# Tool to export citydb.sqlite to citydb.tsv and optionally resolve GeonameIds
# for legacy entries using GeoNames alternate-name data.
#
# Provenance of citydb.tsv:
#   1. Run this script against the old citydb.sqlite to produce an initial TSV.
#   2. Manually fix 10 malformed lat/lon entries that existed in the original binary DB:
#      - Altenstadt (Germany): lat was decimal degrees (47.8339), converted to DMS  47 50' 2"
#      - Adak (USA):           lat missing leading space sign character
#      - Amami Island (Japan): lat missing leading space; lon had trailing space
#      - Amchitka (USA):       lon missing leading space sign character
#      - Bastia (France):      lat missing leading space sign character
#      - Bosscha (Indonesia):  lon missing leading space sign character
#      - Cape Canaveral (USA): lat missing leading space sign character
#      - Cape May (USA):       lat missing leading space sign character
#      - Sewerqia (Saudi Arabia): lat and lon missing leading space sign character
#      - Preston (United Kingdom): lat and lon had no spaces in DMS tokens and no sign;
#                                  lon corrected to negative (west longitude)
#   3. Elevation values were rounded to 2 decimal places to remove float32
#      artifacts from the original SQLite database (e.g. 1788.29004 -> 1788.29).
#   4. Three elevations were corrected from feet to meters (audited against GeoNames):
#      - Bryce Canyon (USA): 9100 -> 2774
#      - Hollywood, CA (USA): 276 -> 84
#      - Olympia, WA (USA): 119 -> 36
#
# The valid KStars DMS format is: leading ' ' (positive) or '-' (negative) followed by
# DD deg MM' SS". generate-citydb.py validates this at build time and skips any entry
# that does not conform.
import sqlite3
import argparse
from collections import defaultdict


def build_en_name_lookups(geonames_tsv, i18n_file):
    """Build (en_name, country) -> geoname_id lookups from GeoNames data.

    Returns two dicts:
      en_name_lookup:          (en_name_lower, country_lower)                -> geoname_id
      en_name_lookup_province: (en_name_lower, province_lower, country_lower) -> geoname_id
    """
    # Step 1: geoname_id -> (country, province) from 10-column TSV
    geonameid_to_loc = {}
    with open(geonames_tsv, 'r', encoding='utf-8') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) == 10 and fields[9]:
                try:
                    geonameid_to_loc[int(fields[9])] = (fields[2], fields[1])
                except ValueError:
                    pass

    # Step 2: collect English alternate names grouped by (name, country)
    candidates = defaultdict(set)  # (name, country) -> set of (geoname_id, province)
    with open(i18n_file, 'r', encoding='utf-8') as f:
        for line in f:
            fields = line.strip('\r\n').split('\t')
            if len(fields) >= 3 and fields[1] == 'en':
                try:
                    geoname_id = int(fields[0])
                except ValueError:
                    continue
                loc = geonameid_to_loc.get(geoname_id)
                if not loc:
                    continue
                country, province = loc
                key = (fields[2].lower(), country.lower())
                candidates[key].add((geoname_id, province))

    # Step 3: unambiguous country-level matches; province lookup for the rest
    en_name_lookup = {}
    en_name_lookup_province = {}
    ambiguous = 0
    for key, cands in candidates.items():
        if len(cands) == 1:
            en_name_lookup[key] = next(iter(cands))[0]
        else:
            for geoname_id, province in cands:
                en_name_lookup_province[(key[0], province.lower(), key[1])] = geoname_id
            ambiguous += 1
    if ambiguous:
        print(f"WARNING: {ambiguous} English names ambiguous at country level "
              f"(will resolve via province when available)")

    return en_name_lookup, en_name_lookup_province


def resolve_geoname_id(name, province, country, en_name_lookup, en_name_lookup_province):
    nl, pl, cl = name.lower(), province.lower(), country.lower()
    return (en_name_lookup_province.get((nl, pl, cl))
            or en_name_lookup.get((nl, cl)))


def main():
    parser = argparse.ArgumentParser(
        description="Export citydb.sqlite to citydb.tsv, optionally resolving GeonameIds")
    parser.add_argument("input", metavar="citydb.sqlite")
    parser.add_argument("output", metavar="citydb.tsv")
    parser.add_argument("--geonames-tsv", metavar="cities15000.tsv",
                        help="10-column GeoNames TSV used to build the GeonameId lookup")
    parser.add_argument("--i18n-file", metavar="cities15000_i18n.tsv",
                        help="GeoNames alternate-names TSV used to build the GeonameId lookup")
    args = parser.parse_args()

    if bool(args.geonames_tsv) != bool(args.i18n_file):
        parser.error("--geonames-tsv and --i18n-file must be provided together")

    en_name_lookup = {}
    en_name_lookup_province = {}
    if args.geonames_tsv:
        en_name_lookup, en_name_lookup_province = build_en_name_lookups(
            args.geonames_tsv, args.i18n_file)

    conn = sqlite3.connect(args.input)
    cur = conn.cursor()
    cur.execute("SELECT Name, Province, Country, Latitude, Longitude, TZ, TZRule, Elevation, GeonameId FROM city ORDER BY Country, Name")
    rows = cur.fetchall()
    conn.close()

    matched = missing = mismatched = 0
    with open(args.output, 'w', encoding='utf-8') as f:
        for name, province, country, lat, lon, tz, rule, elev, geoname_id_db in rows:
            elev_s = f"{elev:.2f}" if elev is not None else ""
            fields = [
                name or "", province or "", country or "",
                lat or "", lon or "", str(tz) if tz is not None else "",
                rule or "", elev_s,
            ]

            inferred = None
            if en_name_lookup or en_name_lookup_province:
                inferred = resolve_geoname_id(
                    name or "", province or "", country or "",
                    en_name_lookup, en_name_lookup_province)

            geoname_id_out = geoname_id_db
            if inferred:
                if geoname_id_db is None:
                    print(f"WARNING: {name!r} ({province}, {country})"
                          f" missing GeonameId -- expected {inferred}")
                    geoname_id_out = inferred
                    matched += 1
                elif geoname_id_db != inferred:
                    print(f"ERROR: {name!r} ({province}, {country})"
                          f" GeonameId mismatch: db={geoname_id_db} expected={inferred}")
                    mismatched += 1
                else:
                    matched += 1
            elif geoname_id_db is None:
                missing += 1

            # Column 9 is population (empty for legacy entries), column 10 is GeonameId.
            # This matches the cities15000.tsv 10-column layout so generate-citydb.py
            # reads both files uniformly.
            fields.append("")  # population: not available for legacy entries
            fields.append(str(geoname_id_out) if geoname_id_out is not None else "")
            f.write("\t".join(fields) + "\n")

    total = len(rows)
    print(f"Exported {total} cities to {args.output}")
    if en_name_lookup or en_name_lookup_province:
        print(f"  GeonameId: {matched} resolved, {missing} unmatched in GeoNames, "
              f"{mismatched} mismatches (see ERRORs above)")


if __name__ == "__main__":
    main()
