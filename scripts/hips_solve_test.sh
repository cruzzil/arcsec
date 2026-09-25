#!/usr/bin/env bash
# hips_solve_test.sh — download 20 test fields (1000x1000 px) from HiPS2FITS, run the
# arcsec blind solver (-i) on each, and report solve rate and positional accuracy against
# the known ground truth. Run from the repository root.
#
# Needs the astrometry.net 4100 indexes (`arcsec catalog install anet-4100`) and an ASTAP
# star database in the catalogue directory, since the blind estimate is refined by the
# catalogue solve.
#
# Usage:
#   scripts/hips_solve_test.sh [--fov <deg>] [--index <path>] [--arcsec <path>]
#                              [--work-dir <dir>] [--survey <hips>] [--skip-download]
#
# Options:
#   --fov <deg>       image FOV in degrees [default: 1.5]
#   --index <path>    index dir or file [default: the output of `arcsec catalog path`]
#   --arcsec <path>   arcsec binary [default: ./target/release/arcsec]
#   --work-dir <dir>  where to store downloaded FITS and outputs [default: /tmp/hips_test]
#   --survey <hips>   HiPS survey [default: CDS/P/DSS2/blue]
#   --skip-download   reuse FITS files already in work-dir

set -uo pipefail

ARCSEC_BIN="./target/release/arcsec"
INDEX_PATH=""
FOV="1.5"
WORK_DIR="/tmp/hips_test"
SKIP_DOWNLOAD=0
HIPS="CDS/P/DSS2/blue"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --fov)       FOV="$2"; shift 2 ;;
        --index)     INDEX_PATH="$2"; shift 2 ;;
        --arcsec)    ARCSEC_BIN="$2"; shift 2 ;;
        --work-dir)  WORK_DIR="$2"; shift 2 ;;
        --survey)    HIPS="$2"; shift 2 ;;
        --skip-download) SKIP_DOWNLOAD=1; shift ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# ── Ground-truth positions: label ra_deg dec_deg ────────────────────────────
declare -a TARGETS=(
    "Orion-M42        83.82   -5.39"
    "Andromeda-M31    10.68   41.27"
    "Virgo-cluster   187.71   12.39"
    "Cygnus-X        305.55   40.73"
    "Auriga-M36       78.00   43.00"
    "Perseus-h-chi    51.40   49.85"
    "Coma-Berenices  195.00   28.00"
    "Centaurus-A     202.47  -43.02"
    "Ursa-Major      180.00   65.00"
    "Aquarius        350.00   -1.00"
    "Leo-triplet     169.73   13.39"
    "Taurus-Pleiades  56.75   24.11"
    "Gemini-M35      100.00   20.00"
    "Pegasus         345.00   15.00"
    "Sculptor        015.00  -33.00"
    "Pisces          020.00   15.00"
    "Eridanus        055.00  -22.00"
    "Corvus          185.00  -17.00"
    "Bootes          218.00   35.00"
    "Piscis-Aus      344.41  -29.62"
)

WIDTH=1000
HEIGHT=1000

# Default to the catalogue directory, where `arcsec catalog install anet-4100` puts the
# astrometry.net index files.
[[ -n "$INDEX_PATH" ]] || INDEX_PATH="$("$ARCSEC_BIN" catalog path)"

mkdir -p "$WORK_DIR"

pass=0; fail=0; total=0
sum_ra_err=0; sum_dec_err=0; n_accurate=0

echo "Arcsec:  $ARCSEC_BIN"
echo "Index:   $INDEX_PATH"
echo "FOV:     ${FOV}°   Image: ${WIDTH}×${HEIGHT}"
echo "Work:    $WORK_DIR"
echo ""
printf "%-22s  %8s  %8s  %8s  %8s  %s\n" "Target" "RA_true" "DEC_true" "ΔRA\"" "ΔDEC\"" "Result"
echo "$(printf '%0.s-' {1..80})"

for entry in "${TARGETS[@]}"; do
    label=$(echo "$entry" | awk '{print $1}')
    ra_true=$(echo "$entry" | awk '{print $2}')
    dec_true=$(echo "$entry" | awk '{print $3}')

    fits_file="$WORK_DIR/${label}.fits"
    wcs_file="${fits_file%.fits}.wcs"
    log_file="${fits_file%.fits}.log"

    # ── Download ─────────────────────────────────────────────────────────────
    if [[ $SKIP_DOWNLOAD -eq 0 || ! -f "$fits_file" ]]; then
        url="https://alasky.cds.unistra.fr/hips-image-services/hips2fits"
        url+="?hips=$(python3 -c "import urllib.parse; print(urllib.parse.quote('$HIPS'))")"
        url+="&ra=${ra_true}&dec=${dec_true}&fov=${FOV}"
        url+="&width=${WIDTH}&height=${HEIGHT}&projection=TAN&coordsys=icrs&format=fits"

        if ! curl -sf --max-time 30 "$url" -o "$fits_file"; then
            printf "%-22s  %8s  %8s  %8s  %8s  %s\n" "$label" "$ra_true" "$dec_true" "-" "-" "DOWNLOAD_FAIL"
            (( fail++ )) || true
            (( total++ )) || true
            continue
        fi
    fi

    # ── Solve ─────────────────────────────────────────────────────────────────
    rm -f "$wcs_file" "$log_file"
    solve_start=$(date +%s%3N)
    blind_out=$("$ARCSEC_BIN" -f "$fits_file" -i "$INDEX_PATH" --fov "$FOV" -t 0.02 2>/dev/null)
    solve_end=$(date +%s%3N)
    elapsed=$(( solve_end - solve_start ))

    # ── Parse blind position estimate from stdout ──────────────────────────────
    # Line format: "Index position estimate: RA=56.750°, Dec=24.108°"
    ra_got=$(echo "$blind_out" | grep -o 'RA=[0-9.-]*' | head -1 | cut -d= -f2)
    dec_got=$(echo "$blind_out" | grep -o 'Dec=[0-9.-]*' | head -1 | cut -d= -f2)

    if [[ -z "$ra_got" || -z "$dec_got" ]]; then
        printf "%-22s  %8s  %8s  %8s  %8s  %s\n" "$label" "$ra_true" "$dec_true" "-" "-" "NO_SOLVE  (${elapsed}ms)"
        (( fail++ )) || true
        (( total++ )) || true
        continue
    fi

    # Angular errors in arcseconds (cos-DEC correction on RA)
    delta_ra=$(python3 -c "
import math
ra1, ra2 = math.radians($ra_true), math.radians($ra_got)
dec = math.radians($dec_true)
dra = abs(ra1 - ra2)
if dra > math.pi: dra = 2*math.pi - dra
print(f'{dra * math.cos(dec) * 3600 * 180 / math.pi:.2f}')
")
    delta_dec=$(python3 -c "print(f'{abs($dec_true - $dec_got) * 3600:.2f}')")

    # Success: blind solve within 30 arcmin of true position on both axes
    ok=$(python3 -c "print('1' if float('$delta_ra') < 1800 and float('$delta_dec') < 1800 else '0')")
    if [[ "$ok" == "1" ]]; then
        (( pass++ )) || true
        sum_ra_err=$(python3 -c "print($sum_ra_err + $delta_ra)")
        sum_dec_err=$(python3 -c "print($sum_dec_err + $delta_dec)")
        (( n_accurate++ )) || true
        status="OK  (${elapsed}ms)"
    else
        (( fail++ )) || true
        status="WRONG (${elapsed}ms)"
    fi
    (( total++ )) || true

    printf "%-22s  %8s  %8s  %8s  %8s  %s\n" \
        "$label" "$ra_true" "$dec_true" "$delta_ra" "$delta_dec" "$status"
done

echo "$(printf '%0.s-' {1..80})"
echo ""
echo "Results: $pass/$total solved  ($(( pass * 100 / (total > 0 ? total : 1) ))% solve rate)"
if [[ $n_accurate -gt 0 ]]; then
    mean_ra=$(python3 -c "print(f'{$sum_ra_err / $n_accurate:.2f}')")
    mean_dec=$(python3 -c "print(f'{$sum_dec_err / $n_accurate:.2f}')")
    echo "Accuracy: mean ΔRA=${mean_ra}\"  mean ΔDEC=${mean_dec}\"  (n=$n_accurate)"
fi
