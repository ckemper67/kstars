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

# Entries that were malformed in the original binary database and have been corrected.
# Each dict has: error_kind, name, province, country, bad (lat, lon, tz, rule, elev),
# good (lat, lon, tz, rule, elev). The test assembles full rows as (name, province, country) + bad/good.
FORMERLY_BROKEN = [
    {"error_kind": "latlon", "name": "Altenstadt", "province": "", "country": "Germany",
     "bad":  ("47.8339",            " 10\xb0 52' 04\"",  "1.0",  "EU", "739.00"),
     "good": (" 47\xb0 50' 2\"",    " 10\xb0 52' 04\"",  "1.0",  "EU", "739.00")},
    {"error_kind": "latlon", "name": "Adak", "province": "Alaska", "country": "USA",
     "bad":  ("51\xb0 52' 46\"",    "-176\xb0 43' 48\"", "-9.0", "US", "1.37"),
     "good": (" 51\xb0 52' 46\"",   "-176\xb0 43' 48\"", "-9.0", "US", "1.37")},
    {"error_kind": "latlon", "name": "Amami Island", "province": "", "country": "Japan",
     "bad":  ("28\xb0 19' 15\"",    "129\xb0 8' 55\" ",  "9.0",  "--", "219.57"),
     "good": (" 28\xb0 19' 15\"",   " 129\xb0 8' 55\"",  "9.0",  "--", "219.57")},
    {"error_kind": "latlon", "name": "Amchitka", "province": "Alaska", "country": "USA",
     "bad":  (" 51\xb0 22' 48\"",   "179\xb0 15' 35\"",  "-9.0", "US", "59.99"),
     "good": (" 51\xb0 22' 48\"",   " 179\xb0 15' 35\"", "-9.0", "US", "59.99")},
    {"error_kind": "latlon", "name": "Bastia", "province": "Haute-Corse", "country": "France",
     "bad":  ("42\xb0 42' 8\"",     " 9\xb0 30' 00\"",   "1.0",  "EU", "296.44"),
     "good": (" 42\xb0 42' 8\"",    " 9\xb0 30' 00\"",   "1.0",  "EU", "296.44")},
    {"error_kind": "latlon", "name": "Bosscha", "province": "", "country": "Indonesia",
     "bad":  ("-6\xb0 53' 21\"",    "107\xb0 35' 56\"",  "7.5",  "--", "961.34"),
     "good": ("-6\xb0 53' 21\"",    " 107\xb0 35' 56\"", "7.5",  "--", "961.34")},
    {"error_kind": "latlon", "name": "Cape Canaveral", "province": "Florida", "country": "USA",
     "bad":  ("28\xb0 30' 27\"",    "-80\xb0 39' 9\"",   "-5.0", "US", "1.00"),
     "good": (" 28\xb0 30' 27\"",   "-80\xb0 39' 9\"",   "-5.0", "US", "1.00")},
    {"error_kind": "latlon", "name": "Cape May", "province": "New Jersey", "country": "USA",
     "bad":  ("38\xb0 56' 25\"",    "-74\xb0 56' 20\"",  "-5.0", "US", "3.00"),
     "good": (" 38\xb0 56' 25\"",   "-74\xb0 56' 20\"",  "-5.0", "US", "3.00")},
    {"error_kind": "latlon", "name": "Sewerqia", "province": "Madina", "country": "Saudi Arabia",
     "bad":  ("23\xb0 20' 36\"",    "40\xb0 19' 30\"",   "3.0",  "--", "905.07"),
     "good": (" 23\xb0 20' 36\"",   " 40\xb0 19' 30\"",  "3.0",  "--", "905.07")},
    {"error_kind": "latlon", "name": "Preston", "province": "Lancashire", "country": "United Kingdom",
     "bad":  ("53\xb045'35\"",      "2\xb042'25\"",      "0.0",  "EU", "-10.00"),
     "good": (" 53\xb0 45' 35\"",   "-2\xb0 42' 25\"",   "0.0",  "EU", "-10.00")},
    {"error_kind": "elev", "name": "Bryce Canyon National Park (Tropic) IDS", "province": "Utah", "country": "USA",
     "bad":  (" 37\xb0 59' 30\"",   "-112\xb0 19' 11\"", "-6.0", "US", "9100.00"),
     "good": (" 37\xb0 59' 30\"",   "-112\xb0 19' 11\"", "-6.0", "US", "2774.00")},
    {"error_kind": "elev", "name": "Hollywood", "province": "California", "country": "USA",
     "bad":  (" 34\xb0 07' 30\"",   "-118\xb0 20' 13\"", "-8.0", "US", "276.46"),
     "good": (" 34\xb0 07' 30\"",   "-118\xb0 20' 13\"", "-8.0", "US", "84.00")},
    {"error_kind": "elev", "name": "Olympia", "province": "Washington", "country": "USA",
     "bad":  (" 47\xb0 01' 59\"",   "-122\xb0 58' 01\"", "-8.0", "US", "118.83"),
     "good": (" 47\xb0 01' 59\"",   "-122\xb0 58' 01\"", "-8.0", "US", "36.20")},
]

GOOD_ROW = ("TestCity", "TestProvince", "TestCountry", " 37\xb0 20' 26\"", " 10\xb0 00' 00\"", "1.0", "EU", "100.0")


def run_script(args):
    r = subprocess.run([sys.executable, SCRIPT] + args, capture_output=True, text=True)
    return r


def write_tsv(rows, f):
    for row in rows:
        f.write("\t".join(row) + "\n")
    f.flush()


class TestDMSValidation(unittest.TestCase):

    def test_broken_originals_are_rejected(self):
        """Entries with detectable errors must be rejected with a WARNING."""
        detectable = [e for e in FORMERLY_BROKEN
                      if e["error_kind"] == "latlon"
                      or float(e["bad"][4]) < -500 or float(e["bad"][4]) > 9000]
        with tempfile.TemporaryDirectory() as tmpdir:
            tsv = os.path.join(tmpdir, "cities.tsv")
            db  = os.path.join(tmpdir, "out.sqlite")
            with open(tsv, "w", encoding="utf-8") as f:
                write_tsv(
                    [(e["name"], e["province"], e["country"]) + e["bad"]
                     for e in detectable],
                    f
                )
            r = run_script([tsv, "--sqlite", db])
            self.assertNotEqual(r.returncode, 0)
            self.assertEqual(r.stderr.count("WARNING"), len(detectable))
            for e in detectable:
                self.assertIn(e["name"], r.stderr)
            self.assertIn("cities.tsv:", r.stderr)

    def test_fixed_entries_are_accepted(self):
        """The corrected versions of all entries must be accepted."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tsv = os.path.join(tmpdir, "cities.tsv")
            db  = os.path.join(tmpdir, "out.sqlite")
            with open(tsv, "w", encoding="utf-8") as f:
                write_tsv(
                    [(e["name"], e["province"], e["country"]) + e["good"]
                     for e in FORMERLY_BROKEN],
                    f
                )
            r = run_script([tsv, "--sqlite", db])
            self.assertEqual(r.returncode, 0, r.stderr)
            count = sqlite3.connect(db).execute("SELECT COUNT(*) FROM city").fetchone()[0]
            self.assertEqual(count, len(FORMERLY_BROKEN))
            self.assertEqual(r.stderr.count("WARNING"), 0, r.stderr)

    def test_bad_tz_warns_and_skips(self):
        """A row with non-numeric TZ must warn and be skipped, not crash."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tsv = os.path.join(tmpdir, "cities.tsv")
            db  = os.path.join(tmpdir, "out.sqlite")
            with open(tsv, "w", encoding="utf-8") as f:
                write_tsv([
                    ("BadTZ", "P", "C", " 10\xb0 00' 00\"", " 20\xb0 00' 00\"", "notanumber", "EU", "0.0"),
                    GOOD_ROW,
                ], f)
            r = run_script([tsv, "--sqlite", db])
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn("cities.tsv:1:", r.stderr)
            self.assertIn("'BadTZ'", r.stderr)
            count = sqlite3.connect(db).execute("SELECT COUNT(*) FROM city").fetchone()[0]
            self.assertEqual(count, 1)

    def test_bad_elevation_warns_and_skips(self):
        """A row with non-numeric elevation must warn and be skipped, not crash."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tsv = os.path.join(tmpdir, "cities.tsv")
            db  = os.path.join(tmpdir, "out.sqlite")
            with open(tsv, "w", encoding="utf-8") as f:
                write_tsv([
                    ("BadElev", "P", "C", " 10\xb0 00' 00\"", " 20\xb0 00' 00\"", "1.0", "EU", "notanumber"),
                    GOOD_ROW,
                ], f)
            r = run_script([tsv, "--sqlite", db])
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn("cities.tsv:1:", r.stderr)
            self.assertIn("'BadElev'", r.stderr)
            count = sqlite3.connect(db).execute("SELECT COUNT(*) FROM city").fetchone()[0]
            self.assertEqual(count, 1)

    def test_tz_out_of_range_warns_and_skips(self):
        """A TZ outside [-12, 14] must warn and be skipped."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tsv = os.path.join(tmpdir, "cities.tsv")
            db  = os.path.join(tmpdir, "out.sqlite")
            with open(tsv, "w", encoding="utf-8") as f:
                write_tsv([
                    ("BadTZRange", "P", "C", " 10\xb0 00' 00\"", " 20\xb0 00' 00\"", "99.0", "EU", "0.0"),
                    GOOD_ROW,
                ], f)
            r = run_script([tsv, "--sqlite", db])
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn("cities.tsv:1:", r.stderr)
            self.assertIn("'BadTZRange'", r.stderr)
            count = sqlite3.connect(db).execute("SELECT COUNT(*) FROM city").fetchone()[0]
            self.assertEqual(count, 1)

    def test_zero_rows_exits_nonzero(self):
        """An input with no valid rows must exit with a non-zero return code."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tsv = os.path.join(tmpdir, "cities.tsv")
            db  = os.path.join(tmpdir, "out.sqlite")
            with open(tsv, "w", encoding="utf-8") as f:
                write_tsv([
                    ("BadLat", "P", "C", "51.5085", " 20\xb0 00' 00\"", "1.0", "EU", "0.0"),
                ], f)
            r = run_script([tsv, "--sqlite", db])
            self.assertNotEqual(r.returncode, 0)
            self.assertIn("ERROR", r.stderr)
            self.assertIn("cities.tsv:1:", r.stderr)
            self.assertIn("'BadLat'", r.stderr)


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
            # Apply the generated SQL to a fresh database and verify the row survives
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
                    GOOD_ROW,
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
                    ("Springfield", "Illinois",  "USA", " 39\xb0 47' 58\"", "-89\xb0 39' 00\"", "-6.0", "US", "168.0"),
                    ("Springfield", "Missouri",  "USA", " 37\xb0 13' 00\"", "-93\xb0 18' 00\"", "-6.0", "US", "417.0"),
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
            row = ("Springfield", "Illinois", "USA", " 39\xb0 47' 58\"", "-89\xb0 39' 00\"", "-6.0", "US", "168.0")
            with open(tsv, "w", encoding="utf-8") as f:
                write_tsv([row, row], f)
            run_script([tsv, "--sqlite", db])
            count = sqlite3.connect(db).execute("SELECT COUNT(*) FROM city").fetchone()[0]
            self.assertEqual(count, 1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
