#!/usr/bin/env python3
import sys
import argparse
import re

# Match valid language codes: 2-3 characters (e.g., de, fr) or regional locales (e.g., zh-cn, pt-br)
# Excludes post, link, wkdt, iata, icao, etc.
LANG_RE = re.compile(r'^[a-z]{2,3}([-_][a-zA-Z0-9]+)?$')

def main():
    parser = argparse.ArgumentParser(description="Extract and filter alternate names from GeoNames alternateNames.txt")
    parser.add_argument("alternate_names_file", help="Path to raw alternateNames.txt")
    parser.add_argument("cities_file", help="Path to cities15000.tsv (KStars TSV) or raw cities15000.txt (GeoNames)")
    args = parser.parse_args()

    # 1. Load active geoname IDs from the cities file
    active_ids = set()
    print(f"Reading active city IDs from {args.cities_file}...", file=sys.stderr)
    with open(args.cities_file, 'r', encoding='utf-8') as f:
        for line in f:
            fields = line.strip('\r\n').split('\t')
            if not fields or fields[0].startswith('#'):
                continue
            if len(fields) >= 19:
                # Raw GeoNames file: geonameid is the 1st column
                active_ids.add(fields[0])
            elif len(fields) >= 10:
                # KStars 10-column TSV: geonameid is the 10th column
                active_ids.add(fields[9])
            elif len(fields) == 9:
                # Fallback check if it's cities15000.tsv but 9-column
                # (if we haven't updated it yet)
                pass

    print(f"Loaded {len(active_ids)} active city IDs.", file=sys.stderr)
    if not active_ids:
        print("Error: No city IDs found. Check if the cities file format is correct.", file=sys.stderr)
        sys.exit(1)

    # 2. Process alternateNames.txt line-by-line to avoid loading the huge file into memory
    print(f"Filtering alternate names from {args.alternate_names_file}...", file=sys.stderr)
    
    # Store: (geonameid, lang) -> (alternate_name, score)
    translations = {}

    with open(args.alternate_names_file, 'r', encoding='utf-8') as f:
        for lineno, line in enumerate(f, 1):
            if lineno % 1000000 == 0:
                print(f"Processed {lineno} lines...", file=sys.stderr)

            fields = line.strip('\r\n').split('\t')
            if len(fields) < 4:
                continue

            geonameid = fields[1]
            if geonameid not in active_ids:
                continue

            lang = fields[2].lower()
            if not LANG_RE.match(lang):
                continue

            # Convert to standard format with underscores, e.g. zh-cn -> zh_cn
            lang = lang.replace('-', '_')

            name = fields[3]
            if not name:
                continue

            # Extract flags for scoring
            is_preferred = fields[4] if len(fields) > 4 else '0'
            is_short = fields[5] if len(fields) > 5 else '0'
            is_colloquial = fields[6] if len(fields) > 6 else '0'
            is_historic = fields[7] if len(fields) > 7 else '0'

            # Calculate a score to determine the best translation for this language
            score = 0
            if is_preferred == '1':
                score += 10
            if is_short == '1':
                score += 5
            if is_colloquial == '1':
                score -= 5
            if is_historic == '1':
                score -= 20  # Strongly avoid historical names

            # Keep only the translation with the highest score
            key = (geonameid, lang)
            if key in translations:
                existing_name, existing_score = translations[key]
                if score > existing_score:
                    translations[key] = (name, score)
                elif score == existing_score:
                    # Tie-breaker: prefer shorter names to keep database clean
                    if len(name) < len(existing_name):
                        translations[key] = (name, score)
            else:
                translations[key] = (name, score)

    # 3. Output the filtered alternate names
    print(f"Writing {len(translations)} translation rows...", file=sys.stderr)
    for (geonameid, lang), (name, score) in sorted(translations.items()):
        # Output format: geonameid \t lang \t alternate_name
        # Keep historic names out if the final score is very low
        if score > -10:
            print(f"{geonameid}\t{lang}\t{translations[(geonameid, lang)][0]}")

    print("Done!", file=sys.stderr)

if __name__ == "__main__":
    main()
