#!/usr/bin/env bash
# compare-solvers.sh — compare arcsec vs ASTAP on a FITS image.
# Usage: compare-solvers.sh <fits_file> [--ra <hours>] [--spd <degrees>] [--fov <degrees>]
#                           [--db <path>] [--db-name <name>] [--astap <path>] [--arcsec <path>]
#
# --ra   right ascension in hours (0-24)   — passed to both solvers
# --spd  south-pole-distance in degrees    — passed to both solvers
# --fov  field-of-view in degrees          — if omitted, each solver auto-detects

set -uo pipefail

# ── Defaults ────────────────────────────────────────────────────────────────
FITS_FILE=""
RA=""
SPD=""
FOV=""
DB_PATH=""
DB_NAME="d80"
ASTAP_BIN="${ASTAP:-astap}"
ARCSEC_BIN="${ARCSEC:-$(dirname "$0")/../target/release/arcsec}"
if [[ ! -x "$ARCSEC_BIN" ]]; then
    ARCSEC_BIN="$(dirname "$0")/../target/debug/arcsec"
fi

# ── Argument parsing ─────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --ra)      RA="$2";       shift 2 ;;
        --spd)     SPD="$2";      shift 2 ;;
        --fov)     FOV="$2";      shift 2 ;;
        --db)      DB_PATH="$2";  shift 2 ;;
        --db-name) DB_NAME="$2";  shift 2 ;;
        --astap)   ASTAP_BIN="$2"; shift 2 ;;
        --arcsec)  ARCSEC_BIN="$2"; shift 2 ;;
        -*)        echo "Unknown option: $1" >&2; exit 1 ;;
        *)         FITS_FILE="$1"; shift ;;
    esac
done

if [[ -z "$FITS_FILE" ]]; then
    echo "Usage: compare-solvers.sh <fits_file> [--ra <hours>] [--spd <deg>] [--fov <deg>] [--db <path>]" >&2
    exit 1
fi

if [[ ! -f "$FITS_FILE" ]]; then
    echo "Error: FITS file not found: $FITS_FILE" >&2
    exit 1
fi

STEM="${FITS_FILE%.*}"
DIR="$(dirname "$FITS_FILE")"
BASE="$(basename "$STEM")"

# ── Helper: parse .ini file (always exits 0 even if field is absent) ─────────
parse_ini() {
    local file="$1" key="$2"
    { grep -i "^${key}=" "$file" 2>/dev/null || true; } | head -1 | cut -d= -f2- | tr -d ' '
}

# ── Run ASTAP ────────────────────────────────────────────────────────────────
ASTAP_ELAPSED="N/A"
ASTAP_EXIT=127
A_CRVAL1="" A_CRVAL2="" A_CDELT1="" A_CDELT2="" A_CROTA2="" A_RMS=""
ASTAP_INI="${DIR}/${BASE}.ini"
ASTAP_WCS="${DIR}/${BASE}.wcs"

if command -v "$ASTAP_BIN" &>/dev/null || [[ -x "$ASTAP_BIN" ]]; then
    ASTAP_ARGS=("-f" "$FITS_FILE")
    [[ -n "$RA"      ]] && ASTAP_ARGS+=("-ra"  "$RA")
    [[ -n "$SPD"     ]] && ASTAP_ARGS+=("-spd" "$SPD")
    [[ -n "$FOV"     ]] && ASTAP_ARGS+=("-fov" "$FOV")
    [[ -n "$DB_PATH" ]] && ASTAP_ARGS+=("-d"   "$DB_PATH")
    ASTAP_ARGS+=("-update")

    echo "── Running ASTAP ──────────────────────────────────────────────"
    ASTAP_TIME_FILE=$(mktemp)
    set +e
    /usr/bin/time -f "%e" "$ASTAP_BIN" "${ASTAP_ARGS[@]}" 2>"$ASTAP_TIME_FILE"
    ASTAP_EXIT=$?
    set -e
    ASTAP_ELAPSED=$(cat "$ASTAP_TIME_FILE")
    rm -f "$ASTAP_TIME_FILE"

    # Extract ASTAP values NOW before arcsec overwrites the .ini
    if [[ $ASTAP_EXIT -eq 0 && -f "$ASTAP_INI" ]]; then
        A_CRVAL1=$(parse_ini "$ASTAP_INI" CRVAL1)
        A_CRVAL2=$(parse_ini "$ASTAP_INI" CRVAL2)
        A_CDELT1=$(parse_ini "$ASTAP_INI" CDELT1)
        A_CDELT2=$(parse_ini "$ASTAP_INI" CDELT2)
        A_CROTA2=$(parse_ini "$ASTAP_INI" CROTA2)
        A_RMS=$(parse_ini    "$ASTAP_INI" RMS)
    fi
else
    echo "── ASTAP not found (set ASTAP env var or pass --astap) ─────────"
fi

# ── Run arcsec ───────────────────────────────────────────────────────────────
if [[ ! -x "$ARCSEC_BIN" ]]; then
    echo "Error: arcsec binary not found at $ARCSEC_BIN" >&2
    echo "Run: cargo build --release -p arcsec"
    exit 1
fi

# Write arcsec output to a temp dir to avoid overwriting ASTAP's files.
ARCSEC_TMPDIR=$(mktemp -d)
ARCSEC_OUT_BASE="$ARCSEC_TMPDIR/result"

ARCSEC_ARGS=("-f" "$FITS_FILE" "-o" "$ARCSEC_OUT_BASE")
[[ -n "$RA"      ]] && ARCSEC_ARGS+=("--ra"  "$RA")
[[ -n "$SPD"     ]] && ARCSEC_ARGS+=("--spd" "$SPD")
[[ -n "$FOV"     ]] && ARCSEC_ARGS+=("--fov" "$FOV")
[[ -n "$DB_PATH" ]] && ARCSEC_ARGS+=("-d"    "$DB_PATH")
ARCSEC_ARGS+=("-D" "$DB_NAME")

echo "── Running arcsec ─────────────────────────────────────────────"
ARCSEC_TIME_FILE=$(mktemp)
set +e
/usr/bin/time -f "%e" "$ARCSEC_BIN" "${ARCSEC_ARGS[@]}" 2>"$ARCSEC_TIME_FILE"
ARCSEC_EXIT=$?
set -e
ARCSEC_ELAPSED=$(cat "$ARCSEC_TIME_FILE")
rm -f "$ARCSEC_TIME_FILE"

ARCSEC_INI="$ARCSEC_OUT_BASE.ini"
ARCSEC_WCS="$ARCSEC_OUT_BASE.wcs"
P_CRVAL1="" P_CRVAL2="" P_CDELT1="" P_CDELT2="" P_CROTA2="" P_RMS=""
if [[ $ARCSEC_EXIT -eq 0 && -f "$ARCSEC_INI" ]]; then
    P_CRVAL1=$(parse_ini "$ARCSEC_INI" CRVAL1)
    P_CRVAL2=$(parse_ini "$ARCSEC_INI" CRVAL2)
    P_CDELT1=$(parse_ini "$ARCSEC_INI" CDELT1)
    P_CDELT2=$(parse_ini "$ARCSEC_INI" CDELT2)
    P_CROTA2=$(parse_ini "$ARCSEC_INI" CROTA2)
    P_RMS=$(parse_ini    "$ARCSEC_INI" RMS)
fi

# ── Report ────────────────────────────────────────────────────────────────────
echo ""
echo "════════════════════════════════════════════════════════════════"
echo "  COMPARISON REPORT: $(basename "$FITS_FILE")"
echo "════════════════════════════════════════════════════════════════"
echo ""

echo "── Solution status ─────────────────────────────────────────────"
if [[ $ASTAP_EXIT -eq 0 ]]; then
    ASTAP_SOLVED="YES"
elif [[ $ASTAP_EXIT -eq 127 ]]; then
    ASTAP_SOLVED="N/A (not installed)"
else
    ASTAP_SOLVED="NO (exit $ASTAP_EXIT)"
fi
if [[ $ARCSEC_EXIT -eq 0 ]]; then ARCSEC_SOLVED="YES"; else ARCSEC_SOLVED="NO (exit $ARCSEC_EXIT)"; fi
printf "  ASTAP:  %-25s  arcsec: %s\n" "$ASTAP_SOLVED" "$ARCSEC_SOLVED"
echo ""

echo "── WCS solution ────────────────────────────────────────────────"
[[ -n "$A_CRVAL1" ]] && printf "  ASTAP   RA=%-16s Dec=%-16s CDELT2=%-14s CROTA2=%-10s RMS=%s\"\n" \
    "$A_CRVAL1" "$A_CRVAL2" "$A_CDELT2" "$A_CROTA2" "${A_RMS:-?}"
[[ -n "$P_CRVAL1" ]] && printf "  arcsec  RA=%-16s Dec=%-16s CDELT2=%-14s CROTA2=%-10s RMS=%s\"\n" \
    "$P_CRVAL1" "$P_CRVAL2" "$P_CDELT2" "$P_CROTA2" "${P_RMS:-?}"

if [[ -n "$A_CRVAL1" && -n "$P_CRVAL1" && $ASTAP_EXIT -eq 0 && $ARCSEC_EXIT -eq 0 ]]; then
    echo ""
    echo "── Deltas (arcsec − ASTAP) ─────────────────────────────────────"
    awk -v ac1="$A_CRVAL1" -v pc1="$P_CRVAL1" \
        -v ac2="$A_CRVAL2" -v pc2="$P_CRVAL2" \
        -v acd2="$A_CDELT2" -v pcd2="$P_CDELT2" \
    'BEGIN {
        dra_as  = ((pc1+0) - (ac1+0)) * 3600
        ddec_as = ((pc2+0) - (ac2+0)) * 3600
        # Use abs(CDELT) to avoid sign-convention differences
        abs_acd2 = acd2 < 0 ? -acd2 : acd2
        abs_pcd2 = pcd2 < 0 ? -pcd2 : pcd2
        dcd2_pct = abs_acd2 != 0 ? (abs_pcd2 - abs_acd2) / abs_acd2 * 100 : 0
        flag_ra   = (dra_as  >  1 || dra_as  < -1) ? "  <<< >1 arcsec" : ""
        flag_dec  = (ddec_as >  1 || ddec_as < -1) ? "  <<< >1 arcsec" : ""
        flag_cd2  = (dcd2_pct > 0.1 || dcd2_pct < -0.1) ? "  <<< >0.1%" : ""
        printf "  ΔRA   = %+.4f arcsec%s\n", dra_as,   flag_ra
        printf "  ΔDec  = %+.4f arcsec%s\n", ddec_as,  flag_dec
        printf "  |ΔCDELT2| = %+.4f%%%s\n",  dcd2_pct, flag_cd2
    }'
fi
echo ""

echo "── Star counts ─────────────────────────────────────────────────"
[[ -n "$A_CRVAL1" ]] && printf "  ASTAP   NSTARS=%s  NQUADS=%s\n" \
    "$(parse_ini "$ASTAP_INI" NSTARS)" "$(parse_ini "$ASTAP_INI" NQUADS)"
if [[ -f "$ARCSEC_INI" ]]; then
    printf "  arcsec  NSTARS=%s  NQUADS=%s\n" \
        "$(parse_ini "$ARCSEC_INI" NSTARS)" "$(parse_ini "$ARCSEC_INI" NQUADS)"
fi
echo ""

echo "── WCS header diff (ASTAP vs arcsec — numeric values only) ─────"
if [[ -f "$ASTAP_WCS" && -f "$ARCSEC_WCS" ]]; then
    diff --label ASTAP --label arcsec \
        <(grep -E "^(CRPIX|CRVAL|CD[12]_[12]|CDELT|CROTA|CTYPE|CUNIT)" "$ASTAP_WCS") \
        <(grep -E "^(CRPIX|CRVAL|CD[12]_[12]|CDELT|CROTA|CTYPE|CUNIT)" "$ARCSEC_WCS") || true
elif [[ -f "$ARCSEC_WCS" ]]; then
    echo "  (ASTAP .wcs not available; arcsec .wcs written to $ARCSEC_WCS)"
fi
echo ""

echo "── Performance ─────────────────────────────────────────────────"
printf "  ASTAP   elapsed: %s s\n" "$ASTAP_ELAPSED"
printf "  arcsec  elapsed: %s s\n" "$ARCSEC_ELAPSED"
if [[ "$ASTAP_ELAPSED" != "N/A" && "$ARCSEC_ELAPSED" != "N/A" ]]; then
    awk -v at="$ASTAP_ELAPSED" -v pt="$ARCSEC_ELAPSED" \
    'BEGIN { if (at+0 > 0) printf "  ratio   arcsec/ASTAP = %.2fx\n", (pt+0)/(at+0) }'
fi
echo ""
echo "════════════════════════════════════════════════════════════════"

# Cleanup temp dir
rm -rf "$ARCSEC_TMPDIR"
