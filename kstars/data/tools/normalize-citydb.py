#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Christian Kemper <ckemper@gmail.com>
# SPDX-License-Identifier: GPL-2.0-or-later
"""Backfill GeonameIds into a KStars city TSV using GeoNames data.

Column 3 of the input TSV is an ISO 3166-1 alpha-2 country code -- the same
code GeoNames uses (field 9 of the raw cities file). Matching is done on
(name, ISO code): the ISO code is unambiguous, so no country-name resolution
is needed. City names are matched in any language or script via
alternateNames.txt, with the GeoNames primary and ASCII names as fallbacks
and an ASCII-folded retry for diacritics.

Several cities in one country can share a name (e.g. "Ambler" in Alaska and
Pennsylvania, or the generic Turkish district name "Merkez"), so a name+ISO
match alone is not enough. The candidate whose GeoNames coordinates are
closest to the legacy entry's coordinates is chosen, and only accepted if it
is within MAX_COORD_DELTA_DEG. Coordinates are a robust disambiguator: they
need no province-name string matching and reject the wrong same-name city even
when the intended one is absent from the GeoNames input. An unmatched entry is
simply left without a GeonameId (its English name still displays); a wrong
GeonameId would attach wrong translations, so the matcher errs toward leaving
entries unresolved.

For each legacy entry (empty population field) without a GeonameId, the
matching record's id is written to column 10. City names, provinces, and
country codes are left unchanged.

Usage:
    python3 normalize-citydb.py \
        --cities-raw /path/to/cities1000.txt \
        --alt-names  /path/to/alternateNames.txt \
        --input      citydb.tsv \
        --output     citydb.tsv

    # Dry-run (stats + unresolved entries, no write):
    python3 normalize-citydb.py ... --input citydb.tsv --dry-run

    # Restrict to one ISO country code:
    python3 normalize-citydb.py ... --country TR
"""
import re
import sys
import argparse
import math
import unicodedata
from collections import defaultdict

# A name+ISO match is accepted only if the candidate's GeoNames coordinates are
# close enough to the legacy entry's coordinates. Two tolerances apply:
#
#   AMBIGUOUS_MAX_DEG - when several GeoNames cities in the country share the
#       name, the closest must be within this (strict) bound. This disambiguates
#       reliably, e.g. picking the right "Springfield".
#   UNIQUE_MAX_DEG    - when the name is unique in the country, a looser bound
#       tolerates imprecise legacy coordinates (legacy entries are often only
#       arcminute-accurate, e.g. Rio de Janeiro is ~0.7 deg off) while still
#       rejecting a different same-name place when the intended city is absent
#       from the GeoNames input (e.g. legacy "Ambler, Alaska" must not match the
#       only GeoNames "Ambler", which is in Pennsylvania, ~87 deg away).
AMBIGUOUS_MAX_DEG = 0.5
UNIQUE_MAX_DEG = 2.0

_DMS_RE = re.compile(r'^\s*([ -]?)(\d+)\xb0\s*(\d+)\'\s*(\d+)"')


def dms_to_deg(s):
    """Parse a KStars DMS string to signed decimal degrees, or None."""
    m = _DMS_RE.match(s)
    if not m:
        return None
    sign = -1 if m.group(1) == '-' else 1
    return sign * (int(m.group(2)) + int(m.group(3)) / 60.0 + int(m.group(4)) / 3600.0)


def coord_delta(lat1, lon1, lat2, lon2):
    """Approximate angular distance in degrees, handling longitude wraparound."""
    dlon = abs(lon1 - lon2)
    dlon = min(dlon, 360.0 - dlon)
    return math.hypot(lat1 - lat2, dlon)


def to_ascii(s):
    """Strip diacritics via NFD, with explicit Turkish dotless-i handling."""
    s = s.replace('ı', 'i').replace('İ', 'I')
    nfd = unicodedata.normalize('NFD', s)
    return ''.join(c for c in nfd if unicodedata.category(c) != 'Mn')


def load_cities_raw(cities_raw_file, iso_filter=None):
    """Load a raw GeoNames cities file (19-column format).

    Returns gid_to_loc: geoname_id_str -> (iso_code, lat, lon).
    """
    gid_to_loc = {}
    with open(cities_raw_file, encoding='utf-8') as f:
        for line in f:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 19:
                continue
            iso = fields[8]
            if iso_filter and iso != iso_filter:
                continue
            try:
                lat = float(fields[4])
                lon = float(fields[5])
            except ValueError:
                continue
            gid_to_loc[fields[0]] = (iso, lat, lon)
    return gid_to_loc


def build_name_index(gid_to_loc, alt_names_file, cities_raw_file):
    """Build (name_lower, iso_lower) -> set of geoname_id from all name sources.

    Sources: every alternateNames.txt entry plus the GeoNames primary (UTF-8)
    and ASCII names from the raw cities file. No uniqueness pruning happens
    here; ambiguity is resolved later by coordinates.
    """
    city_ids = set(gid_to_loc)
    name_to_gids = defaultdict(set)

    print("Scanning alternateNames.txt...", file=sys.stderr)
    with open(alt_names_file, encoding='utf-8') as f:
        for lineno, line in enumerate(f, 1):
            if lineno % 2000000 == 0:
                print(f"  {lineno // 1000000}M lines...", file=sys.stderr)
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 4:
                continue
            gid = fields[1]
            if gid not in city_ids:
                continue
            alt_name = fields[3]
            if alt_name:
                iso_l = gid_to_loc[gid][0].lower()
                name_to_gids[(alt_name.lower(), iso_l)].add(gid)

    print("Indexing primary and ASCII names from raw cities file...", file=sys.stderr)
    with open(cities_raw_file, encoding='utf-8') as f:
        for line in f:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 3:
                continue
            gid = fields[0]
            if gid not in city_ids:
                continue
            iso_l = gid_to_loc[gid][0].lower()
            for nm in (fields[1], fields[2]):  # primary UTF-8 name, ASCII name
                if nm:
                    name_to_gids[(nm.lower(), iso_l)].add(gid)

    return name_to_gids


def resolve_city(name, lat, lon, iso_l, name_to_gids, gid_to_loc):
    """Return (geoname_id, coord_delta) for the best match, or (None, None).

    Candidates are cities sharing the name (exact, then ASCII-folded) within the
    country. The closest by coordinate is chosen, accepted only within tolerance:
    a strict bound when the name is ambiguous (several candidates), a looser one
    when it is unique (see the module constants).
    """
    candidates = (name_to_gids.get((name.lower(), iso_l))
                  or name_to_gids.get((to_ascii(name).lower(), iso_l)))
    if not candidates:
        return None, None
    best_d = UNIQUE_MAX_DEG if len(candidates) == 1 else AMBIGUOUS_MAX_DEG
    best_gid = None
    for gid in candidates:
        _, glat, glon = gid_to_loc[gid]
        d = coord_delta(lat, lon, glat, glon)
        if d <= best_d:
            best_d = d
            best_gid = gid
    return best_gid, (best_d if best_gid else None)


def main():
    parser = argparse.ArgumentParser(
        description="Backfill GeonameIds into a KStars city TSV (matches on ISO + coordinates)")
    parser.add_argument("--cities-raw", required=True, metavar="cities1000.txt",
                        help="Raw GeoNames cities file (19-column format)")
    parser.add_argument("--alt-names", required=True, metavar="alternateNames.txt",
                        help="GeoNames alternateNames.txt (for name matching)")
    parser.add_argument("--input", required=True, metavar="citydb.tsv")
    parser.add_argument("--output", metavar="citydb.tsv",
                        help="Output path (omit for dry-run)")
    parser.add_argument("--country", metavar="ISO",
                        help="Restrict to this ISO alpha-2 country code (e.g. TR)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print stats and unresolved entries without writing output")
    args = parser.parse_args()

    if not args.output and not args.dry_run:
        parser.error("specify --output or --dry-run")

    iso_filter = args.country.upper() if args.country else None

    print("Loading cities from raw GeoNames file...", file=sys.stderr)
    gid_to_loc = load_cities_raw(args.cities_raw, iso_filter)

    name_to_gids = build_name_index(gid_to_loc, args.alt_names, args.cities_raw)

    # First pass: resolve each row into [fields, status, gid, delta], where
    # status is 'skip', 'already', 'match' or 'nomatch'.
    rows = []
    with open(args.input, encoding='utf-8') as f:
        for line in f:
            fields = line.rstrip('\n').split('\t')

            if len(fields) < 8:
                rows.append([fields, 'skip', None, None])
                continue
            iso = fields[2]
            if iso_filter and iso.upper() != iso_filter:
                rows.append([fields, 'skip', None, None])
                continue
            is_legacy = (len(fields) < 9 or fields[8] == '')
            gid_existing = fields[9] if len(fields) >= 10 else ''
            if not is_legacy or gid_existing:
                rows.append([fields, 'already', None, None])
                continue

            lat = dms_to_deg(fields[3])
            lon = dms_to_deg(fields[4])
            gid = delta = None
            if lat is not None and lon is not None:
                gid, delta = resolve_city(fields[0], lat, lon, iso.lower(),
                                          name_to_gids, gid_to_loc)
            rows.append([fields, 'match' if gid else 'nomatch', gid, delta])

    # Collision pass: a geoname id identifies one city. Entries that match it
    # within the strict bound are confidently that city (co-located duplicates
    # such as "Jakarta"/"Djakarta") and are all kept. Entries accepted only by
    # the looser unique-name bound are redundant-or-wrong if a confident match
    # already claims the id, and at most one survives otherwise -- this drops
    # generic names like several Turkish "Merkez" district centers collapsing
    # onto one GeoNames "Merkez".
    by_gid = defaultdict(list)
    for idx, row in enumerate(rows):
        if row[1] == 'match':
            by_gid[row[2]].append(idx)
    for gid, idxs in by_gid.items():
        if len(idxs) == 1:
            continue
        has_close = any(rows[i][3] <= AMBIGUOUS_MAX_DEG for i in idxs)
        keep = None if has_close else min(idxs, key=lambda i: rows[i][3])
        for i in idxs:
            close = rows[i][3] <= AMBIGUOUS_MAX_DEG
            if not (close if has_close else i == keep):
                rows[i][1] = 'collision'

    # Emit pass: write the GeonameId for rows that kept their match.
    counts = defaultdict(int)
    lines_out = []
    for fields, status, gid, _ in rows:
        counts[status] += 1
        if status == 'match':
            while len(fields) < 9:
                fields.append('')
            if len(fields) == 9:
                fields.append('')
            fields[9] = gid
            lines_out.append('\t'.join(fields))
        else:
            lines_out.append('\t'.join(fields))
            if args.dry_run and status in ('nomatch', 'collision'):
                print(f"  {status}: {fields[0]!r} / {fields[1]!r} / {fields[2]!r}")

    print(f"Total rows: {len(lines_out)}  (skipped: {counts['skip']})")
    print(f"  already had GeonameId:          {counts['already']}")
    print(f"  resolved (GeonameId filled):    {counts['match']}")
    print(f"  unresolved (no GeoNames match): {counts['nomatch']}")
    print(f"  dropped (id claimed by closer entry): {counts['collision']}")

    if not args.dry_run and args.output:
        with open(args.output, 'w', encoding='utf-8') as f:
            for line in lines_out:
                f.write(line + '\n')
        print(f"Wrote {len(lines_out)} rows to {args.output}")


if __name__ == '__main__':
    main()
