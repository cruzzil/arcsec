#!/usr/bin/env bash
# hips_fov_sweep.sh — test the blind solver across multiple FOV sizes.
# Uses a fixed set of 8 star fields; downloads a 1000x1000 image of each at
# 0.5, 1, 1.5, 2, 3 and 5 degrees. Same requirements as hips_solve_test.sh.
# Run from the repository root.
#
# Usage: scripts/hips_fov_sweep.sh [--index <path>] [--arcsec <path>]
#                                  [--work-dir <dir>] [--skip-download]
#   --index <path>    index dir or file [default: the output of `arcsec catalog path`]
#   --arcsec <path>   arcsec binary [default: ./target/release/arcsec]
#   --work-dir <dir>  download and output directory [default: /tmp/hips_fov_sweep]
#   --skip-download   reuse FITS files already in work-dir

set -uo pipefail

ARCSEC_BIN="./target/release/arcsec"
INDEX_PATH=""
WORK_DIR="/tmp/hips_fov_sweep"
SKIP_DOWNLOAD=0
HIPS="CDS/P/DSS2/blue"
WIDTH=1000
HEIGHT=1000

while [[ $# -gt 0 ]]; do
    case "$1" in
        --arcsec)    ARCSEC_BIN="$2"; shift 2 ;;
        --index)     INDEX_PATH="$2"; shift 2 ;;
        --work-dir)  WORK_DIR="$2"; shift 2 ;;
        --skip-download) SKIP_DOWNLOAD=1; shift ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# Fields known to solve reliably at 1.5° (pure star fields)
declare -a TARGETS=(
    "Andromeda   10.68  41.27"
    "Virgo      187.71  12.39"
    "Cygnus-X   305.55  40.73"
    "Perseus     51.40  49.85"
    "Ursa-Major 180.00  65.00"
    "Aquarius   350.00  -1.00"
    "Corvus     185.00 -17.00"
    "Bootes     218.00  35.00"
)

# FOVs to sweep; each column in the output
declare -a FOVS=("0.5" "1.0" "1.5" "2.0" "3.0" "5.0")

# Default to the catalogue directory, where `arcsec catalog install anet-4100` puts the
# astrometry.net index files.
[[ -n "$INDEX_PATH" ]] || INDEX_PATH="$("$ARCSEC_BIN" catalog path)"

mkdir -p "$WORK_DIR"

echo "Arcsec:  $ARCSEC_BIN"
echo "Index:   $INDEX_PATH"
echo "Survey:  $HIPS  (${WIDTH}×${HEIGHT})"
echo ""
printf "%-14s" "Target"
for fov in "${FOVS[@]}"; do printf "  %7s°" "$fov"; done
echo ""
echo "$(printf '%0.s-' {1..80})"

for entry in "${TARGETS[@]}"; do
    label=$(echo "$entry" | awk '{print $1}')
    ra_true=$(echo "$entry" | awk '{print $2}')
    dec_true=$(echo "$entry" | awk '{print $3}')

    printf "%-14s" "$label"

    for fov in "${FOVS[@]}"; do
        fits_file="$WORK_DIR/${label}_fov${fov}.fits"

        if [[ $SKIP_DOWNLOAD -eq 0 || ! -f "$fits_file" ]]; then
            url="https://alasky.cds.unistra.fr/hips-image-services/hips2fits"
            url+="?hips=$(python3 -c "import urllib.parse; print(urllib.parse.quote('$HIPS'))")"
            url+="&ra=${ra_true}&dec=${dec_true}&fov=${fov}"
            url+="&width=${WIDTH}&height=${HEIGHT}&projection=TAN&coordsys=icrs&format=fits"

            if ! curl -sf --max-time 30 "$url" -o "$fits_file"; then
                printf "  %8s" "DL_FAIL"
                continue
            fi
        fi

        solve_start=$(date +%s%3N)
        blind_out=$("$ARCSEC_BIN" -f "$fits_file" -i "$INDEX_PATH" --fov "$fov" -t 0.02 2>/dev/null)
        solve_end=$(date +%s%3N)
        elapsed=$(( solve_end - solve_start ))

        ra_got=$(echo "$blind_out" | grep -o 'RA=[0-9.-]*' | head -1 | cut -d= -f2)
        dec_got=$(echo "$blind_out" | grep -o 'Dec=[0-9.-]*' | head -1 | cut -d= -f2)

        if [[ -z "$ra_got" || -z "$dec_got" ]]; then
            printf "  %8s" "NO_SOLVE"
            continue
        fi

        ok=$(python3 -c "
import math
ra1,ra2 = math.radians($ra_true),math.radians($ra_got)
dec = math.radians($dec_true)
dra = abs(ra1-ra2); dra = min(dra, 2*math.pi-dra)
dra_as = dra * math.cos(dec) * 3600 * 180 / math.pi
ddec_as = abs($dec_true - $dec_got) * 3600
print('OK' if dra_as < 1800 and ddec_as < 1800 else 'WRONG')
")
        printf "  %5s(%ds)" "$ok" "$(( elapsed / 1000 ))"
    done
    echo ""
done

echo "$(printf '%0.s-' {1..80})"
echo "FOV columns: ${FOVS[*]} degrees"
