# KStars Location Database

This directory contains the source data and tools used to generate `citydb.sqlite`, the geographic location database for KStars.

## Data Files

*   `citydb.tsv`: Contains the original ~3,400 cities from the legacy KStars database, augmented with GeonameIds where available.
*   `cities15000.tsv`: Contains all ~34,000 entries from the GeoNames `cities15000` dataset of all cities with a population >15,000 worldwide.
*   `cities15000_i18n.tsv`: Localized city names derived from the GeoNames `alternateNames` dataset, used to populate the `city_i18n` table for KStars UI translation.
*   `countryInfo.txt`, `admin1CodesASCII.txt`, `timeZones.txt`: Support files used for mapping GeoNames codes to full names and KStars DST rules.

## Country names

Column 3 of both `citydb.tsv` and `cities15000.tsv` is an ISO 3166-1 alpha-2
country code (`US`, `DE`, `TR`, ...) -- the same code GeoNames uses. ISO codes
have exactly one correct value per country, so the source data cannot drift into
inconsistencies like "Turkiye" vs "Turkey" or "Ivory coast" vs "Ivory Coast".

`generate-citydb.py` resolves each ISO code to a display name at build time using
the GeoNames canonical English name from `countryInfo.txt`, with two display-name
overrides defined in `tools/kstars_geo.py`:

| GeoNames canonical | KStars display name | Rationale |
|---|---|---|
| United States | USA | Universally recognized short form |
| Czech Republic | Czechia | Official short form (adopted 2016) |
| Turkey | Türkiye | Official English short form (adopted 2022) |

All other countries use the GeoNames canonical name as-is. The C++ code continues
to see human-readable display names; only the source TSVs carry ISO codes.

## Build Instructions

The binary `citydb.sqlite` is automatically generated at build time by CMake using strict validation. You can also generate it manually from the `kstars/data/` directory:

```bash
python3 tools/generate-citydb.py \
    --countries countryInfo.txt --strict \
    --sqlite citydb.sqlite --min-population 50000 \
    --i18n-file cities15000_i18n.tsv \
    citydb.tsv cities15000.tsv
```

In `--strict` mode, validation warnings (unknown country names, malformed
coordinates, etc.) become fatal errors and the build fails. This is the
default in the CMake build to prevent bad data from reaching users.

## How to Contribute (Fixing or Adding Cities)

Contributing to the location database no longer requires knowledge of SQL. You can now edit the data using a standard text editor or spreadsheet application.

### To add or fix a city:

1.  **Edit the source**: Open `citydb.tsv` in your favorite text editor or spreadsheet program (use Tab as the delimiter).
2.  **Make your changes**: Find the row you want to fix, or add a new row at the bottom.
    The file uses a 10-column format:

    | # | Field | Example |
    |---|-------|---------|
    | 1 | Name | Ankara |
    | 2 | Province | |
    | 3 | Country (ISO 3166-1 alpha-2) | TR |
    | 4 | Latitude (DMS) | 39 deg 55' 48" |
    | 5 | Longitude (DMS) | 32 deg 51' 00" |
    | 6 | TZ offset | 3.0 |
    | 7 | TZ rule | -- |
    | 8 | Elevation (m) | 938.00 |
    | 9 | Population | *(leave empty for legacy entries)* |
    | 10 | GeonameId | 323786 |

    Column 3 must be an ISO 3166-1 alpha-2 country code (e.g. `TR`, not `Turkey`).
    `generate-citydb.py` resolves it to a display name at build time, and `--strict`
    fails the build on an unknown code. Run a build after editing to catch errors.

3.  **Prioritization**: Changes made in `citydb.tsv` **automatically override** any data from the larger GeoNames dataset.
4.  **Rebuild and Verify**: CMake rebuilds `citydb.sqlite` automatically on the next build,
    or run the generator manually (see [Build Instructions](#build-instructions)).

## Population cutoff

The GeoNames `cities15000` dataset contains every city with population > 15,000 worldwide
(~34,000 entries). Including all of them would flood the KStars location picker with obscure
towns. The table below shows how database size varies with the cutoff:

| Min population | GeoNames cities | Total (+ legacy) |
|---------------:|----------------:|-----------------:|
|         15,000 |          33,521 |           36,950 |
|         25,000 |          22,436 |           25,865 |
|         50,000 |          12,128 |           15,557 |
|         75,000 |           8,118 |           11,547 |
|        100,000 |           6,063 |            9,492 |
|        150,000 |           3,953 |            7,382 |
|        200,000 |           2,989 |            6,418 |
|        500,000 |           1,163 |            4,592 |
|      1,000,000 |             561 |            3,990 |

**50,000 was chosen as the default.** It adds roughly 3.5x the number of legacy cities while
keeping the list navigable. Cities below 50,000 that matter (observatories, polar stations,
etc.) are already curated in `citydb.tsv`.

To change the population cutoff, update the `--min-population` value in `kstars/data/CMakeLists.txt`.

## How to update the GeoNames data

If you want to update `cities15000.tsv` from a new release of GeoNames:

1.  Download `cities15000.zip` from [GeoNames.org](http://download.geonames.org/export/dump/) and unzip it to get `cities15000.txt`.
2.  Run the converter from the `kstars/data/` directory:
    ```bash
    python3 tools/geonames-to-tsv.py admin1CodesASCII.txt timeZones.txt \
        /path/to/cities15000.txt > cities15000.tsv
    ```
3.  Regenerate the i18n file:
    ```bash
    python3 tools/extract-geonames-i18n.py \
        /path/to/alternateNames.txt cities15000.tsv > cities15000_i18n.tsv
    ```
4.  If the GeoNames support files have also changed, download and commit the updated
    `countryInfo.txt`, `admin1CodesASCII.txt`, and `timeZones.txt` from the same page.
5.  Commit the updated files.

## How to backfill GeonameIds in citydb.tsv

Linking a legacy entry to its GeoNames record (via the GeonameId in column 10)
lets the build attach localized names from `cities15000_i18n.tsv`. Use
`tools/normalize-citydb.py` to backfill the GeonameId for legacy entries that
lack one. The script matches each entry on `(name, ISO code)` -- names are
matched in any language or script via `alternateNames.txt` -- and confirms the
match by coordinates, so it never attaches a wrong id (an unmatched entry is
simply left without one; its English name still displays). City names,
provinces, and country codes are left unchanged.

It requires the raw GeoNames data files (not included in the repository; download from
[GeoNames.org](http://download.geonames.org/export/dump/)):

```bash
python3 tools/normalize-citydb.py \
    --cities-raw /path/to/cities1000.txt \
    --alt-names  /path/to/alternateNames.txt \
    --input      citydb.tsv \
    --output     citydb.tsv
```

To restrict processing to one country, pass its ISO code (useful when adding a
batch of cities for a single country):

```bash
python3 tools/normalize-citydb.py ... --country TR
```

Run with `--dry-run` instead of `--output` to preview what would change without writing.

## Testing the database on macOS (dev builds)

`citydb.sqlite` is read-only data regenerated at build time, so KStars reads it
directly from the app bundle's `Contents/Resources/kstars/`. On macOS
`KSPaths::locate` searches that directory explicitly (see
`kstars/auxiliary/kspaths.cpp`), and startup removes any stale writable copy
from `~/Library/Application Support/kstars/` (see `setupMacKStarsIfNeeded` in
`kstars/auxiliary/ksutils.cpp`). A rebuilt database is therefore picked up
automatically -- no manual copy into Application Support is needed.

## Provenance of citydb.tsv

`citydb.tsv` was created from the original KStars binary `citydb.sqlite` using
`tools/export-citydb.py`. Ten entries in the original binary database had malformed
lat/lon strings and were corrected manually during the export:

| City | Country | Problem | Fix |
|------|---------|---------|-----|
| Altenstadt | Germany | lat was decimal degrees (47.8339) | Converted to DMS: 47 deg 50' 02" |
| Adak | USA | lat missing leading sign character | Added leading space |
| Amami Island | Japan | lat missing leading space; lon had trailing space | Trimmed and padded |
| Amchitka | USA | lon missing leading sign character | Added leading space |
| Bastia | France | lat missing leading sign character | Added leading space |
| Bosscha | Indonesia | lon missing leading sign character | Added leading space |
| Cape Canaveral | USA | lat missing leading sign character | Added leading space |
| Cape May | USA | lat missing leading sign character | Added leading space |
| Sewerqia | Saudi Arabia | lat and lon missing leading sign character | Added leading space to both |
| Preston | United Kingdom | lat and lon had no spaces in DMS tokens and no sign | Reformatted; lon corrected to negative (west) |

The valid KStars DMS format is: optional leading ' ' (positive) or '-' (negative),
followed by `DD deg MM' SS"`. The leading sign character is optional; the C++ parser
accepts bare `DD deg MM' SS"` as positive. `generate-citydb.py` validates this at build
time and skips any entry that does not conform.

## Tools

*   `tools/kstars_geo.py`: Shared module with DMS conversion, DST rule mapping,
    country-name overrides, and GeoNames support-file loaders. Used by all other tools.
*   `tools/generate-citydb.py`: Generates `citydb.sqlite` from one or more TSV files. Validates
    DMS format, TZ range, elevation, and country names at import time. Use `--strict` to make
    validation errors fatal (this is the default in the CMake build). Run with `--help` for usage.
*   `tools/test-generate-citydb.py`: Unit tests for `generate-citydb.py`. Run with `python3`.
*   `tools/geonames-to-tsv.py`: Converts raw GeoNames cities files into the KStars 10-column TSV format.
*   `tools/export-citydb.py`: The one-time bootstrap tool that created `citydb.tsv`
    from the original `citydb.sqlite` binary (see [Provenance](#provenance-of-citydbtsv)).
    Kept for historical reference; it predates the ISO-code format and writes display
    names in column 3, so do not use it to regenerate the current `citydb.tsv`.
*   `tools/extract-geonames-i18n.py`: Filters GeoNames `alternateNames.txt` down to the localized
    city names needed for KStars UI translation, producing `cities15000_i18n.tsv`.
*   `tools/normalize-citydb.py`: Backfills GeonameIds for legacy `citydb.tsv` entries,
    matching on ISO code and name (any script) and confirming the match by coordinates.
