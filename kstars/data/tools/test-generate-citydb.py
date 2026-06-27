#!/usr/bin/env python3
"""
Tests for generate-citydb.py.

Run with:  python3 kstars/data/tools/test-generate-citydb.py
"""
import os
import sys
import sqlite3
import subprocess
import tempfile
import unittest

SCRIPT = os.path.join(os.path.dirname(__file__), "generate-citydb.py")
EXTRACT_SCRIPT = os.path.join(os.path.dirname(__file__), "extract-geonames-i18n.py")

# Generic valid 8-column row used by tests that need a passing row as a control.
GOOD_ROW = ("TestCity", "TestProvince", "TestCountry", " 37\xb0 20' 26\"", " 10\xb0 00' 00\"", "1.0", "EU", "100.0")

# One representative input per detectable validation error class.
# error: the invalid row; valid: the corrected row that must be accepted.
VALIDATION_CASES = [
    {
        "label": "malformed latlon (Altenstadt: lat is decimal degrees not DMS)",
        "error": ("Altenstadt", "", "Germany", "47.8339", " 10\xb0 52' 04\"", "1.0", "EU", "739.00"),
        "valid": ("Altenstadt", "", "Germany", " 47\xb0 50' 02\"", " 10\xb0 52' 04\"", "1.0", "EU", "739.00"),
        "name":  "'Altenstadt'",
    },
    {
        "label": "elevation out of range (Bryce Canyon: 9100 m exceeds 9000 m limit)",
        "error": ("Bryce Canyon National Park (Tropic) IDS", "Utah", "USA", " 37\xb0 59' 30\"", "-112\xb0 19' 11\"", "-6.0", "US", "9100.00"),
        "valid": ("Bryce Canyon National Park (Tropic) IDS", "Utah", "USA", " 37\xb0 59' 30\"", "-112\xb0 19' 11\"", "-6.0", "US", "2774.00"),
        "name":  "'Bryce Canyon National Park (Tropic) IDS'",
    },
    {
        "label": "non-numeric TZ (Adak: TZ field is not a number)",
        "error": ("Adak", "Alaska", "USA", " 51\xb0 52' 46\"", "-176\xb0 43' 48\"", "notanumber", "US", "1.37"),
        "valid": ("Adak", "Alaska", "USA", " 51\xb0 52' 46\"", "-176\xb0 43' 48\"", "-9.0",       "US", "1.37"),
        "name":  "'Adak'",
    },
    {
        "label": "non-numeric elevation (Bryce Canyon: elevation field is not a number)",
        "error": ("Bryce Canyon National Park (Tropic) IDS", "Utah", "USA", " 37\xb0 59' 30\"", "-112\xb0 19' 11\"", "-6.0", "US", "notanumber"),
        "valid": ("Bryce Canyon National Park (Tropic) IDS", "Utah", "USA", " 37\xb0 59' 30\"", "-112\xb0 19' 11\"", "-6.0", "US", "2774.00"),
        "name":  "'Bryce Canyon National Park (Tropic) IDS'",
    },
    {
        "label": "TZ out of range (Adak: TZ 99.0 outside [-12, 14])",
        "error": ("Adak", "Alaska", "USA", " 51\xb0 52' 46\"", "-176\xb0 43' 48\"", "99.0", "US", "1.37"),
        "valid": ("Adak", "Alaska", "USA", " 51\xb0 52' 46\"", "-176\xb0 43' 48\"", "-9.0", "US", "1.37"),
        "name":  "'Adak'",
    },
]


def run_script(args):
    r = subprocess.run([sys.executable, SCRIPT] + args, capture_output=True, text=True)
    return r


def write_tsv(rows, f):
    for row in rows:
        f.write("\t".join(row) + "\n")
    f.flush()


class TestValidationErrors(unittest.TestCase):

    def test_validation_errors(self):
        """Each detectable error class must warn, skip the bad row, and accept the corrected row."""
        for case in VALIDATION_CASES:
            with self.subTest(case["label"]):
                with tempfile.TemporaryDirectory() as tmpdir:
                    tsv = os.path.join(tmpdir, "cities.tsv")
                    db  = os.path.join(tmpdir, "out.sqlite")
                    with open(tsv, "w", encoding="utf-8") as f:
                        write_tsv([case["error"], case["valid"]], f)
                    r = run_script([tsv, "--sqlite", db])
                    self.assertEqual(r.returncode, 0, r.stderr)
                    self.assertIn("cities.tsv:1:", r.stderr)
                    self.assertIn(case["name"], r.stderr)
                    count = sqlite3.connect(db).execute("SELECT COUNT(*) FROM city").fetchone()[0]
                    self.assertEqual(count, 1)

    def test_strict_mode_exits_nonzero_on_unknown_country(self):
        """In strict mode, an unknown country name must cause a non-zero exit."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tsv = os.path.join(tmpdir, "cities.tsv")
            db  = os.path.join(tmpdir, "out.sqlite")
            countries = os.path.join(tmpdir, "countryInfo.txt")
            with open(countries, "w", encoding="utf-8") as f:
                f.write("US\tUSA\tUSA\t840\tUnited States\t...\n")
            with open(tsv, "w", encoding="utf-8") as f:
                write_tsv([
                    ("Ankara", "", "Turkiye", " 39\xb0 55' 48\"", " 32\xb0 51' 00\"", "3.0", "--", "938.0"),
                ], f)
            r = run_script([tsv, "--sqlite", db, "--countries", countries, "--strict"])
            self.assertNotEqual(r.returncode, 0)
            self.assertIn("ERROR", r.stderr)
            self.assertIn("Turkiye", r.stderr)

    def test_zero_rows_exits_nonzero(self):
        """An input with no valid rows must exit with a non-zero return code.
        Altenstadt (Germany): lat stored as decimal in the original binary db."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tsv = os.path.join(tmpdir, "cities.tsv")
            db  = os.path.join(tmpdir, "out.sqlite")
            with open(tsv, "w", encoding="utf-8") as f:
                write_tsv([
                    ("Altenstadt", "", "Germany", "47.8339", " 10\xb0 52' 04\"", "1.0", "EU", "739.00"),
                ], f)
            r = run_script([tsv, "--sqlite", db])
            self.assertNotEqual(r.returncode, 0)
            self.assertIn("ERROR", r.stderr)
            self.assertIn("cities.tsv:1:", r.stderr)
            self.assertIn("'Altenstadt'", r.stderr)


class TestSQLOutput(unittest.TestCase):

    def test_sql_output_apostrophe_escaping(self):
        """City names with apostrophes must be properly escaped in SQL output."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tsv = os.path.join(tmpdir, "cities.tsv")
            sql = os.path.join(tmpdir, "out.sql")
            db  = os.path.join(tmpdir, "out.sqlite")
            with open(tsv, "w", encoding="utf-8") as f:
                write_tsv([
                    ("Coeur d'Alene", "Idaho", "USA", " 47\xb0 40' 37\"", "-116\xb0 46' 45\"", "-7.0", "US", "665.0"),
                ], f)
            r = run_script([tsv, "-o", sql])
            self.assertEqual(r.returncode, 0, r.stderr)
            conn = sqlite3.connect(db)
            with open(sql, encoding="utf-8") as f:
                conn.executescript(f.read())
            row = conn.execute("SELECT Name FROM city").fetchone()
            conn.close()
            self.assertIsNotNone(row)
            self.assertEqual(row[0], "Coeur d'Alene")

    def test_sql_output_row_count_matches_sqlite(self):
        """SQL output and --sqlite output must produce the same row count."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tsv = os.path.join(tmpdir, "cities.tsv")
            sql = os.path.join(tmpdir, "out.sql")
            db_sql    = os.path.join(tmpdir, "from_sql.sqlite")
            db_direct = os.path.join(tmpdir, "direct.sqlite")
            with open(tsv, "w", encoding="utf-8") as f:
                write_tsv([
                    ("Adak", "Alaska", "USA", " 51\xb0 52' 46\"", "-176\xb0 43' 48\"", "-9.0", "US", "1.37"),
                    ("Coeur d'Alene", "Idaho", "USA", " 47\xb0 40' 37\"", "-116\xb0 46' 45\"", "-7.0", "US", "665.0"),
                ], f)
            run_script([tsv, "-o", sql])
            conn = sqlite3.connect(db_sql)
            with open(sql, encoding="utf-8") as f:
                conn.executescript(f.read())
            count_sql = conn.execute("SELECT COUNT(*) FROM city").fetchone()[0]
            conn.close()

            run_script([tsv, "--sqlite", db_direct])
            count_direct = sqlite3.connect(db_direct).execute("SELECT COUNT(*) FROM city").fetchone()[0]

            self.assertEqual(count_sql, count_direct)


class TestProvinceDedup(unittest.TestCase):

    def test_same_name_different_province_both_kept(self):
        """Two cities with same name and country but different province must both be inserted."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tsv = os.path.join(tmpdir, "cities.tsv")
            db  = os.path.join(tmpdir, "out.sqlite")
            with open(tsv, "w", encoding="utf-8") as f:
                write_tsv([
                    ("Springfield", "Illinois", "USA", " 39\xb0 47' 58\"", "-89\xb0 39' 00\"", "-6.0", "US", "168.0"),
                    ("Springfield", "Missouri", "USA", " 37\xb0 13' 00\"", "-93\xb0 18' 00\"", "-6.0", "US", "417.0"),
                ], f)
            run_script([tsv, "--sqlite", db])
            rows = sqlite3.connect(db).execute(
                "SELECT Province FROM city WHERE Name='Springfield' ORDER BY Province"
            ).fetchall()
            self.assertEqual([r[0] for r in rows], ["Illinois", "Missouri"])

    def test_true_duplicate_dropped(self):
        """Exact duplicate (same name+province+country) must be inserted only once."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tsv = os.path.join(tmpdir, "cities.tsv")
            db  = os.path.join(tmpdir, "out.sqlite")
            entry = ("Springfield", "Illinois", "USA", " 39\xb0 47' 58\"", "-89\xb0 39' 00\"", "-6.0", "US", "168.0")
            with open(tsv, "w", encoding="utf-8") as f:
                write_tsv([entry, entry], f)
            run_script([tsv, "--sqlite", db])
            count = sqlite3.connect(db).execute("SELECT COUNT(*) FROM city").fetchone()[0]
            self.assertEqual(count, 1)


class TestMinPopulation(unittest.TestCase):

    def test_9col_below_cutoff_excluded(self):
        """A 9-column row with population below --min-population must be skipped."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tsv = os.path.join(tmpdir, "cities.tsv")
            db  = os.path.join(tmpdir, "out.sqlite")
            low  = GOOD_ROW + ("30000",)
            high = ("BigCity", "P", "C", " 10\xb0 00' 00\"", " 20\xb0 00' 00\"", "1.0", "EU", "0.0", "80000")
            with open(tsv, "w", encoding="utf-8") as f:
                write_tsv([low, high], f)
            r = run_script([tsv, "--sqlite", db, "--min-population", "50000"])
            self.assertEqual(r.returncode, 0, r.stderr)
            count = sqlite3.connect(db).execute("SELECT COUNT(*) FROM city").fetchone()[0]
            self.assertEqual(count, 1)
            name = sqlite3.connect(db).execute("SELECT Name FROM city").fetchone()[0]
            self.assertEqual(name, "BigCity")

    def test_8col_unaffected_by_min_population(self):
        """An 8-column row must pass through regardless of --min-population."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tsv = os.path.join(tmpdir, "cities.tsv")
            db  = os.path.join(tmpdir, "out.sqlite")
            with open(tsv, "w", encoding="utf-8") as f:
                write_tsv([GOOD_ROW], f)
            r = run_script([tsv, "--sqlite", db, "--min-population", "999999999"])
            self.assertEqual(r.returncode, 0, r.stderr)
            count = sqlite3.connect(db).execute("SELECT COUNT(*) FROM city").fetchone()[0]
            self.assertEqual(count, 1)

    def test_9col_no_filter_when_default(self):
        """Without --min-population both low- and high-population 9-column rows are kept."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tsv = os.path.join(tmpdir, "cities.tsv")
            db  = os.path.join(tmpdir, "out.sqlite")
            low  = GOOD_ROW + ("1000",)
            high = ("BigCity", "P", "C", " 10\xb0 00' 00\"", " 20\xb0 00' 00\"", "1.0", "EU", "0.0", "80000")
            with open(tsv, "w", encoding="utf-8") as f:
                write_tsv([low, high], f)
            r = run_script([tsv, "--sqlite", db])
            self.assertEqual(r.returncode, 0, r.stderr)
            count = sqlite3.connect(db).execute("SELECT COUNT(*) FROM city").fetchone()[0]
            self.assertEqual(count, 2)


class TestDegToDms(unittest.TestCase):

    def test_single_digit_degree_passes_validation(self):
        """deg_to_dms output must match DMS_RE for all degree magnitudes.
        Regression: :2d padding bug produced '- 9 deg' instead of '-9 deg'."""
        from kstars_geo import deg_to_dms, DMS_RE
        cases = [
            (-9.66, "negative single-digit"),
            (5.25,  "positive single-digit"),
            (-0.5,  "negative sub-one-degree"),
            (51.5,  "normal two-digit positive"),
            (-33.9, "normal two-digit negative"),
            (179.9, "near 180 positive"),
            (-179.9, "near 180 negative"),
        ]
        for deg, label in cases:
            dms = deg_to_dms(deg, True)
            self.assertRegex(dms, DMS_RE, f"deg_to_dms({deg}) = {dms!r} [{label}]")


class TestGeonameIdAndI18n(unittest.TestCase):

    def test_geoname_id_and_i18n_table(self):
        """10-column rows must import GeonameId, and --i18n-file must populate city_i18n."""
        with tempfile.TemporaryDirectory() as tmpdir:
            cities_tsv = os.path.join(tmpdir, "cities.tsv")
            i18n_tsv   = os.path.join(tmpdir, "cities_i18n.tsv")
            db         = os.path.join(tmpdir, "out.sqlite")

            row1 = GOOD_ROW + ("30000", "123456")
            row2 = ("OtherCity", "OtherProvince", "OtherCountry", "-12\xb0 00' 00\"", " 45\xb0 00' 00\"", "2.0", "--", "50.0")

            with open(cities_tsv, "w", encoding="utf-8") as f:
                write_tsv([row1, row2], f)

            with open(i18n_tsv, "w", encoding="utf-8") as f:
                f.write("123456\tde\tTestStadt\n")
                f.write("123456\tes\tCiudadPrueba\n")

            r = run_script([cities_tsv, "--sqlite", db, "--i18n-file", i18n_tsv])
            self.assertEqual(r.returncode, 0, r.stderr)

            conn = sqlite3.connect(db)

            rows = conn.execute("SELECT Name, GeonameId FROM city ORDER BY Name").fetchall()
            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[0][0], "OtherCity")
            self.assertIsNone(rows[0][1])
            self.assertEqual(rows[1][0], "TestCity")
            self.assertEqual(rows[1][1], 123456)

            i18n_rows = conn.execute("SELECT geoname_id, lang, name FROM city_i18n ORDER BY lang").fetchall()
            self.assertEqual(len(i18n_rows), 2)
            self.assertEqual(i18n_rows[0], (123456, "de", "TestStadt"))
            self.assertEqual(i18n_rows[1], (123456, "es", "CiudadPrueba"))

            conn.close()


class TestTzRule(unittest.TestCase):

    @staticmethod
    def _valid_rule_codes():
        """Rule codes defined in data/TZrules.dat (first token of each data line)."""
        path = os.path.join(os.path.dirname(__file__), "..", "TZrules.dat")
        codes = set()
        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                codes.add(line.split()[0])
        return codes

    @staticmethod
    def _tz_file():
        return os.path.join(os.path.dirname(__file__), "..", "timeZones.txt")

    def test_get_tz_rule_only_emits_known_codes(self):
        """For every zone in timeZones.txt, get_tz_rule must return a code that
        exists in TZrules.dat. A missing code resolves to a null TimeZoneRule in
        C++ and crashes on city selection (regression: Asian-Russia once emitted
        the bogus "RU"). Driven by the real DST map, this also guards every zone,
        not just a hand-picked sample."""
        from kstars_geo import get_tz_rule, load_dst_observed
        valid = self._valid_rule_codes()
        dst = load_dst_observed(self._tz_file())
        self.assertTrue(dst, "load_dst_observed returned nothing")
        for iana_id, observes in dst.items():
            # country_code is only consulted by the America/ branch; "" exercises
            # the default path for every zone and never widens the emitted set.
            code = get_tz_rule(iana_id, "", observes)
            self.assertIn(code, valid,
                          f"get_tz_rule({iana_id!r}, observes_dst={observes}) "
                          f"= {code!r} not in TZrules.dat")

    def test_no_dst_zones_resolve_to_dashes(self):
        """Zones whose Jan/Jul offsets are equal in timeZones.txt must resolve to
        the no-DST rule '--', even when they sit under the Europe/ prefix that the
        heuristic would otherwise map to EU. Covers Turkey/Belarus/Russia, which
        abolished DST, while a genuine DST zone (Berlin) still gets EU."""
        from kstars_geo import get_tz_rule, load_dst_observed
        dst = load_dst_observed(self._tz_file())
        no_dst = ["Europe/Istanbul", "Europe/Minsk", "Europe/Moscow",
                  "Asia/Novosibirsk", "Asia/Vladivostok"]
        for iana_id in no_dst:
            self.assertIn(iana_id, dst, f"{iana_id} missing from timeZones.txt")
            self.assertFalse(dst[iana_id], f"{iana_id} should read as no-DST")
            self.assertEqual(get_tz_rule(iana_id, "", dst[iana_id]), "--",
                             f"{iana_id} must resolve to '--'")
        self.assertTrue(dst.get("Europe/Berlin"), "Berlin should read as DST")
        self.assertEqual(get_tz_rule("Europe/Berlin", "DE", dst["Europe/Berlin"]),
                         "EU")

    def test_observes_dst_false_overrides_heuristic(self):
        """The observes_dst gate, not the prefix heuristic, decides no-DST. Proves
        the gate is wired in: the same zone yields EU when DST is observed and '--'
        when it is not."""
        from kstars_geo import get_tz_rule
        self.assertEqual(get_tz_rule("Europe/Berlin", "DE", True), "EU")
        self.assertEqual(get_tz_rule("Europe/Berlin", "DE", False), "--")


class TestExtractI18nLangNormalization(unittest.TestCase):

    def test_regional_locale_separator_normalized_to_underscore(self):
        """extract-geonames-i18n.py must store regional locales with an underscore
        (zh-CN -> zh_cn) so the lang codes match QLocale::name() (e.g. "zh_cn"),
        which CITY_I18N_SQL binds for the exact-locale join. A hyphen would never
        match and the exact translation would be silently unreachable."""
        with tempfile.TemporaryDirectory() as tmpdir:
            cities = os.path.join(tmpdir, "cities.tsv")
            altnames = os.path.join(tmpdir, "alt.txt")

            # 10-column KStars TSV: geonameid is the 10th column.
            with open(cities, "w", encoding="utf-8") as f:
                f.write("Testville\t\tXX\t0\t0\t1.0\t--\t100\t30000\t12345\n")

            # GeoNames alternateNames columns: id, geonameid, isolanguage, name, ...
            with open(altnames, "w", encoding="utf-8") as f:
                f.write("1\t12345\tzh-CN\tExactName\t0\t0\t0\t0\n")
                f.write("2\t12345\tde\tBaseName\t0\t0\t0\t0\n")

            r = subprocess.run([sys.executable, EXTRACT_SCRIPT, altnames, cities],
                               capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)

            langs = {line.split("\t")[1] for line in r.stdout.splitlines() if line.strip()}
            self.assertIn("zh_cn", langs, f"expected normalized 'zh_cn' in {langs}")
            self.assertNotIn("zh-cn", langs, "hyphenated locale leaked into output")
            self.assertIn("de", langs, "base-language row missing")


if __name__ == "__main__":
    unittest.main(verbosity=2)
