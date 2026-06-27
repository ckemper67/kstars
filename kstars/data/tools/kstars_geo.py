#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Christian Kemper <ckemper@gmail.com>
# SPDX-License-Identifier: GPL-2.0-or-later
"""Shared helpers for KStars location-database tools."""

import re

# Valid KStars DMS format: optional leading ' ' (positive) or '-' (negative),
# then digits, degree symbol, space, minutes, quote, seconds, double-quote.
DMS_RE = re.compile(r'^[ -]?\d+\xb0 \d+\' \d+"$')

# GeoNames canonical country names that KStars displays differently.
# These are uncontroversial short-form preferences, not political overrides.
# All other countries use the GeoNames canonical English name from countryInfo.txt.
DISPLAY_NAME_OVERRIDES = {
    "United States": "USA",
    "Czech Republic": "Czechia",
    "Turkey": "Türkiye",
}

# Reverse lookup: KStars display name -> GeoNames canonical name.
_DISPLAY_TO_CANONICAL = {v: k for k, v in DISPLAY_NAME_OVERRIDES.items()}


def deg_to_dms(deg, is_lat):
    """Convert decimal degrees to KStars DMS string."""
    abs_deg = abs(deg)
    d = int(abs_deg)
    m = int((abs_deg - d) * 60)
    s = int(round((abs_deg - d - m / 60.0) * 3600))
    if s == 60:
        m += 1
        s = 0
    if m == 60:
        d += 1
        m = 0
    sign = " "
    if deg < 0:
        sign = "-"
    return f"{sign}{d}\xb0 {m:02d}' {s:02d}\""


def get_tz_rule(iana_id, country_code, observes_dst=None):
    """Map an IANA timezone ID to a KStars DST rule code.

    observes_dst, when not None, is the authoritative DST flag for the zone
    derived from GeoNames timeZones.txt (its January and July offsets differ).
    A False value forces the no-DST rule "--" regardless of the prefix heuristic
    below, so zones in countries that abolished DST (Russia, Turkey, Belarus,
    and any future case) are handled from data instead of a hardcoded list.
    When None (zone not found in timeZones.txt) the heuristic alone is used.

    The heuristic only ever picks which DST *schedule* applies (EU, US, ...);
    GeoNames does not carry the transition dates, so it cannot replace it.
    """
    if not iana_id:
        return "--"
    if observes_dst is False:
        return "--"
    rule = "--"
    if iana_id.startswith("Europe/"):
        if iana_id in [
            "Europe/Helsinki", "Europe/Kyiv", "Europe/Riga",
            "Europe/Tallinn", "Europe/Vilnius", "Europe/Sofia",
            "Europe/Bucharest", "Europe/Athens",
        ]:
            rule = "EE"
        else:
            rule = "EU"
    elif iana_id.startswith("America/"):
        if country_code in ("US", "CA"):
            rule = "US"
        elif country_code == "MX":
            rule = "MX"
        elif country_code == "BR":
            rule = "BZ"
        elif country_code == "CL":
            rule = "CL"
        elif country_code == "PY":
            rule = "PY"
    elif iana_id.startswith("Australia/"):
        rule = "AU"
    elif iana_id.startswith("Pacific/"):
        if "Auckland" in iana_id:
            rule = "NZ"
        elif "Chatham" in iana_id:
            rule = "CZ"
        elif "Tongatapu" in iana_id:
            rule = "TG"
    elif iana_id == "Asia/Tehran":
        rule = "IR"
    elif iana_id == "Asia/Amman":
        rule = "JD"
    elif iana_id == "Asia/Beirut":
        rule = "LB"
    elif iana_id == "Asia/Damascus":
        rule = "SY"
    elif iana_id in ("Asia/Jerusalem", "Asia/Gaza", "Asia/Hebron"):
        rule = "ZN"
    elif iana_id.startswith("Africa/Windhoek"):
        rule = "NB"
    elif iana_id.startswith("Atlantic/Stanley"):
        rule = "FK"
    return rule


def country_name(iso_code, canonical_name):
    """Return the KStars display name for a country, applying overrides."""
    return DISPLAY_NAME_OVERRIDES.get(canonical_name, canonical_name)


def canonical_country_name(display_name):
    """Return the GeoNames canonical name for a KStars display name."""
    return _DISPLAY_TO_CANONICAL.get(display_name, display_name)


def load_countries(countries_file):
    """Load countryInfo.txt -> dict mapping ISO code to KStars display name.

    Also returns the set of valid display names (lowercased) for validation.
    """
    countries = {}
    with open(countries_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 5:
                iso = fields[0]
                countries[iso] = country_name(iso, fields[4])
    known = {name.lower() for name in countries.values()}
    return countries, known


def load_admin1(admin1_file):
    """Load admin1CodesASCII.txt -> dict mapping 'CC.ADM1' to province name.

    The KStars copy is a 5-column reformat of the GeoNames file:
    country code, admin1 code, name, ASCII name, geonameid. The ASCII name
    is used so province values stay ASCII-clean.
    """
    admin1 = {}
    with open(admin1_file, 'r', encoding='utf-8') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 4:
                admin1[f"{fields[0]}.{fields[1]}"] = fields[3]
    return admin1


def load_timezones(tz_file):
    """Load timeZones.txt -> dict mapping IANA timezone ID to UTC offset string."""
    timezones = {}
    with open(tz_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('CountryCode'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 5:
                timezones[fields[1]] = fields[4]
    return timezones


def load_dst_observed(tz_file):
    """Load timeZones.txt -> dict mapping IANA timezone ID to a DST-observed bool.

    GeoNames records each zone's offset on 1 January and 1 July. The two differ
    only when the zone shifts its clocks, so unequal offsets mean DST is observed
    (true in either hemisphere -- only the inequality matters, not its sign).
    Equal offsets mean permanent standard time. This is the authoritative source
    for whether a zone has DST; see get_tz_rule's observes_dst argument.
    """
    dst = {}
    with open(tz_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('CountryCode'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 4:
                try:
                    jan = float(fields[2])
                    jul = float(fields[3])
                except ValueError:
                    continue
                dst[fields[1]] = (jan != jul)
    return dst
