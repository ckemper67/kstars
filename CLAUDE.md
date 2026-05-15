# KStars Development Guide (macOS, Apple Silicon)

## Overview

KStars and all its dependencies (Qt6, KDE Frameworks, INDI, StellarSolver, etc.) are built
and installed into a self-contained craft prefix at `~/work/craft-root`. The craft meta-build
system manages all dependencies. Local source overrides are configured in
`~/work/craft-root/etc/BlueprintSettings.ini`.

## Source layout

    ~/work/kstars             - kstars (this repo, active development branches)
    ~/work/kstars-staging     - git worktree on branch kstars-staging (craft builds from here)
    ~/work/indi               - INDI core (local changes)
    ~/work/indi-3rdparty      - INDI 3rd party drivers (local changes)
    ~/work/craft-root         - craft prefix (all deps install here)

## Staging worktree

Craft builds kstars from `~/work/kstars-staging`, not from `~/work/kstars` directly.
`kstars-staging` is a git worktree that shares the same object store as `~/work/kstars`.

To integrate branches into a build:

    cd ~/work/kstars-staging
    git merge donuts-guiding   # merge any feature branch
    # then rebuild (see below)

`~/work/kstars` is for active development. `~/work/kstars-staging` is the integration point
for craft builds. They share history -- no duplicate git objects.

## Building

**IMPORTANT**: Source `craftenv.sh` only once per shell session. Sourcing it multiple times
prepends `CRAFT: ` to PS1 repeatedly. Open a new terminal tab instead of re-sourcing.

### Full build via craft (first time or after dependency changes)

    source ~/work/craft-root/craft/craftenv.sh
    craft -i kde/applications/kstars

This builds from `~/work/kstars-staging` (srcDir override in BlueprintSettings.ini).
Takes ~10-30 min with cached packages; much longer on first run.

### Incremental rebuild (daily workflow after changing kstars source)

    cd ~/work/craft-root/build/kde/applications/kstars/work/build
    ninja -j$(sysctl -n hw.ncpu)

Run ninja directly in the craft build directory -- skips craft's fetch/configure overhead.
When ninja finishes cleanly, run `craft -i kde/applications/kstars` once to let craft do
the install step into the app bundle.

### After changing INDI core

    source ~/work/craft-root/craft/craftenv.sh
    cd ~/work/indi/build-mac
    ninja -j$(sysctl -n hw.ncpu) && ninja install

### After changing INDI 3rd party

    source ~/work/craft-root/craft/craftenv.sh
    # Libraries first:
    cd ~/work/indi-3rdparty/build-mac-libs
    ninja -j$(sysctl -n hw.ncpu) && ninja install
    # Then drivers:
    cd ~/work/indi-3rdparty/build-mac
    ninja -j$(sysctl -n hw.ncpu) && ninja install

### Clean the craft kstars build dir

    rm -rf ~/work/craft-root/build/kde/applications/kstars

## Running

    source ~/work/craft-root/craft/craftenv.sh
    ~/work/craft-root/Applications/KDE/KStars.app/Contents/MacOS/kstars

Sourcing craftenv.sh is required at runtime -- it sets DYLD_LIBRARY_PATH and Qt plugin
paths so the binary finds craft's libraries.

## Running tests

Tests require a separate cmake build with `BUILD_TESTING=ON`. The craft build uses
`BUILD_TESTING=OFF`, so tests are not built there.

### One-time test build setup

    source ~/work/craft-root/craft/craftenv.sh
    mkdir -p ~/work/kstars/build-mac-tests
    cd ~/work/kstars/build-mac-tests
    cmake \
      -G Ninja \
      -DCMAKE_PREFIX_PATH="$KDEROOT" \
      -DCMAKE_INSTALL_PREFIX="$KDEROOT" \
      -DBUILD_WITH_QT6=ON \
      -DBUILD_TESTING=ON \
      -DCMAKE_BUILD_TYPE=Debug \
      ~/work/kstars
    ninja -j$(sysctl -n hw.ncpu)

### Running tests

    # Full unit test suite (excludes UI tests that need a display):
    ctest --test-dir ~/work/kstars/build-mac-tests/Tests -LE ui --output-on-failure

    # Stable tests only:
    ctest --test-dir ~/work/kstars/build-mac-tests/Tests -L stable --output-on-failure

    # Run a specific test binary directly (fastest for development):
    ~/work/kstars/build-mac-tests/Tests/ekos/guide/testdonutsguider.app/Contents/MacOS/testdonutsguider -v2

    # Run a single test function:
    ~/work/kstars/build-mac-tests/Tests/ekos/guide/testdonutsguider.app/Contents/MacOS/testdonutsguider testDonutsBasic -v2

    # Python test runner with formatted output:
    cd ~/work/kstars
    python3 Tests/run_tests.py --no-build -k guide

See `Tests/README.md` for the full test taxonomy, CTest labels, and available test binaries.

### Standalone tests (no cmake required)

Some ekos tests can be built and run standalone:

    cd ~/work/kstars/Tests/ekos/guide
    # standalone_donuts_test binary is committed here for quick iteration

## craft configuration files

    ~/work/craft-root/etc/CraftSettings.ini      - craft global settings (Python path, cache, etc.)
    ~/work/craft-root/etc/BlueprintSettings.ini  - per-package overrides (srcDir, buildTests, etc.)

Key BlueprintSettings.ini overrides:
- `[kde/applications/kstars] srcDir` -> `~/work/kstars-staging`
- `[libs/indilib/indi] srcDir` -> `~/work/indi`
- `[libs/indilib/indi-3rdparty] srcDir` / `[libs/indilib/indi-3rdparty-libs] srcDir` -> `~/work/indi-3rdparty`
- `[libs/libftdi] buildTests = False` (boost_unit_test_framework not available in craft)

See `~/work/kstars-on-osx-craft/BUILD_NOTES.md` for full history of pitfalls and fixes.
