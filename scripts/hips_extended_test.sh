#!/usr/bin/env bash
# hips_extended_test.sh — 20 additional diverse sky fields at 1.5° FOV, blind solver.
# Complements hips_solve_test.sh (different positions, field types); same requirements.
# Run from the repository root.
#
# Usage: scripts/hips_extended_test.sh [--fov <deg>] [--index <path>] [--arcsec <path>]
#                                      [--work-dir <dir>] [--survey <hips>] [--skip-download]
#   --fov <deg>       image FOV in degrees [default: 1.5]
#   --index <path>    index dir or file [default: the output of `arcsec catalog path`]
#   --arcsec <path>   arcsec binary [default: ./target/release/arcsec]
#   --work-dir <dir>  download and output directory [default: /tmp/hips_ext_test]
#   --survey <hips>   HiPS survey string [default: CDS/P/DSS2/blue]
#   --skip-download   reuse FITS files already in work-dir

set -uo pipefail

ARCSEC_BIN="./target/release/arcsec"
INDEX_PATH=""
FOV="1.5"
WORK_DIR="/tmp/hips_ext_test"
SKIP_DOWNLOAD=0
HIPS="CDS/P/DSS2/blue"
WIDTH=1000
HEIGHT=1000

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

# ── Ground-truth positions ───────────────────────────────────────────────────
# Mix: pure star fields (expect OK), clusters/nebulae/galaxies (expect NO_SOLVE)
declare -a TARGETS=(
    # Pure star fields — various galactic latitudes and RA coverage
    "Lyra-Vega         279.23   38.78"   # b=+19°, Vega nearby, dense-ish
    "Cepheus           340.00   70.00"   # b=+16°, sparse polar field
    "Ophiuchus         261.00   -7.00"   # b=+13°, near galactic plane
    "Libra             228.00  -15.00"   # b=+35°, clean mid-latitude
    "Columba            84.00  -35.00"   # b=-26°, southern hemisphere
    "Volans            128.00  -70.00"   # b=-13°, southern sparse
    "Lupus             237.00  -35.00"   # b=+22°, relatively empty
    "Phoenix            15.00  -50.00"   # b=-67°, high southern lat, sparse
    "Ara               261.77  -56.60"   # b=-10°, galactic plane dense
    "Pyxis             134.00  -27.00"   # b=+25°, clean star field
    # Galactic plane / dense fields
    "Sagittarius-arm   274.00  -20.00"   # b=-5°, galactic arm, very dense
    "Cassiopeia        350.85   58.81"   # b=-2°, galactic plane
    "Monoceros          96.00    4.00"   # b=-4°, galactic plane
    # Cluster fields (expect NO_SOLVE — open clusters confuse star detection)
    "Cancer-M44        130.10   19.67"   # Beehive open cluster
    "Hercules-M13      250.42   36.46"   # Rich globular cluster
    "Omega-Cen         201.70  -47.48"   # Omega Centauri (huge globular)
    # Galaxy fields (expect NO_SOLVE — galaxy nuclei ≠ catalog stars)
    "Triangulum-M33     23.46   30.66"   # M33 spiral galaxy
    "Fornax-cluster     54.62  -35.45"   # Fornax galaxy cluster
    # Other challenging fields
    "M17-Nebula        275.09  -16.18"   # Omega nebula (emission)
    "LMC-center         80.89  -69.75"   # Large Magellanic Cloud
)

# Default to the catalogue directory, where `arcsec catalog install anet-4100` puts the
# astrometry.net index files.
[[ -n "$INDEX_PATH" ]] || INDEX_PATH="$("$ARCSEC_BIN" catalog path)"

mkdir -p "$WORK_DIR"

pass=0; fail=0; total=0
sum_ra_err=0; sum_dec_err=0; n_accurate=0

survey_label=$(echo "$HIPS" | tr '/' '-')
echo "Arcsec:  $ARCSEC_BIN"
echo "Index:   $INDEX_PATH"
echo "Survey:  $HIPS"
echo "FOV:     ${FOV}°   Image: ${WIDTH}×${HEIGHT}"
echo "Work:    $WORK_DIR"
echo ""
printf "%-22s  %8s  %8s  %8s  %8s  %s\n" "Target" "RA_true" "DEC_true" "ΔRA\"" "ΔDEC\"" "Result"
echo "$(printf '%0.s-' {1..80})"

for entry in "${TARGETS[@]}"; do
    label=$(echo "$entry" | awk '{print $1}')
    ra_true=$(echo "$entry" | awk '{print $2}')
    dec_true=$(echo "$entry" | awk '{print $3}')

    fits_file="$WORK_DIR/${label}_${survey_label}.fits"

    if [[ $SKIP_DOWNLOAD -eq 0 || ! -f "$fits_file" ]]; then
        url="https://alasky.cds.unistra.fr/hips-image-services/hips2fits"
        url+="?hips=$(python3 -c "import urllib.parse; print(urllib.parse.quote('$HIPS'))")"
        url+="&ra=${ra_true}&dec=${dec_true}&fov=${FOV}"
        url+="&width=${WIDTH}&height=${HEIGHT}&projection=TAN&coordsys=icrs&format=fits"

        if ! curl -sf --max-time 30 "$url" -o "$fits_file"; then
            printf "%-22s  %8s  %8s  %8s  %8s  %s\n" "$label" "$ra_true" "$dec_true" "-" "-" "DOWNLOAD_FAIL"
            (( fail++ )) || true; (( total++ )) || true; continue
        fi
    fi

    solve_start=$(date +%s%3N)
    blind_out=$("$ARCSEC_BIN" -f "$fits_file" -i "$INDEX_PATH" --fov "$FOV" -t 0.02 2>/dev/null)
    solve_end=$(date +%s%3N)
    elapsed=$(( solve_end - solve_start ))

    ra_got=$(echo "$blind_out" | grep -o 'RA=[0-9.-]*' | head -1 | cut -d= -f2)
    dec_got=$(echo "$blind_out" | grep -o 'Dec=[0-9.-]*' | head -1 | cut -d= -f2)

    if [[ -z "$ra_got" || -z "$dec_got" ]]; then
        printf "%-22s  %8s  %8s  %8s  %8s  %s\n" "$label" "$ra_true" "$dec_true" "-" "-" "NO_SOLVE  (${elapsed}ms)"
        (( fail++ )) || true; (( total++ )) || true; continue
    fi

    delta_ra=$(python3 -c "
import math
ra1, ra2 = math.radians($ra_true), math.radians($ra_got)
dec = math.radians($dec_true)
dra = abs(ra1 - ra2)
if dra > math.pi: dra = 2*math.pi - dra
print(f'{dra * math.cos(dec) * 3600 * 180 / math.pi:.2f}')
")
    delta_dec=$(python3 -c "print(f'{abs($dec_true - $dec_got) * 3600:.2f}')")

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
