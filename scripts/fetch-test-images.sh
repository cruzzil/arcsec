#!/usr/bin/env bash
# fetch-test-images.sh — materialise the benchmark corpus described in docs/test-images.md
# from the manifest in scripts/test-images.tsv.
#
# Usage:
#   scripts/fetch-test-images.sh                  # fetch everything
#   scripts/fetch-test-images.sh --tier A         # one tier (A, B or D); repeatable
#   scripts/fetch-test-images.sh --id fov_1p00    # one entry; repeatable
#   scripts/fetch-test-images.sh --list           # print the manifest, download nothing
#   scripts/fetch-test-images.sh --out <dir>      # default resources/testset
#   scripts/fetch-test-images.sh --force          # re-download files that already exist
#
# Downloads are skipped when the target already exists, so re-running is cheap.
# A truth.tsv is written alongside the images with the ground-truth centre and FOV of
# every tier-A entry (for tiers B and D the truth is the delivered header). Its `file`
# column is a bare filename, resolved against truth.tsv's own directory - an absolute
# path there would bake one machine's layout into the corpus.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MANIFEST="$SCRIPT_DIR/test-images.tsv"
OUT_DIR="$SCRIPT_DIR/../resources/testset"
LIST_ONLY=0
FORCE=0
declare -a WANT_TIERS=()
declare -a WANT_IDS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --tier)     WANT_TIERS+=("$2"); shift 2 ;;
        --id)       WANT_IDS+=("$2");   shift 2 ;;
        --out)      OUT_DIR="$2";       shift 2 ;;
        --manifest) MANIFEST="$2";      shift 2 ;;
        --list)     LIST_ONLY=1;        shift ;;
        --force)    FORCE=1;            shift ;;
        -h|--help)  sed -n '2,18p' "${BASH_SOURCE[0]}"; exit 0 ;;
        *)          echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

[[ -f "$MANIFEST" ]] || { echo "Manifest not found: $MANIFEST" >&2; exit 1; }

wanted() {
    local id="$1" tier="$2"
    if [[ ${#WANT_IDS[@]} -gt 0 ]]; then
        local w; for w in "${WANT_IDS[@]}"; do [[ "$w" == "$id" ]] && return 0; done
        return 1
    fi
    if [[ ${#WANT_TIERS[@]} -gt 0 ]]; then
        local w; for w in "${WANT_TIERS[@]}"; do [[ "$w" == "$tier" ]] && return 0; done
        return 1
    fi
    return 0
}

# curl with a generous timeout; prints the HTTP code on stdout.
fetch() {
    local url="$1" dest="$2"
    curl -sS -L --retry 2 --max-time 300 -o "$dest" -w '%{http_code}' "$url" 2>/dev/null
}

# A FITS file must start with "SIMPLE  =" and be a multiple of 2880 bytes.
valid_fits() {
    local f="$1"
    [[ -s "$f" ]] || return 1
    [[ "$(head -c 9 "$f")" == "SIMPLE  =" ]] || return 1
    local sz; sz=$(stat -c %s "$f")
    (( sz % 2880 == 0 ))
}

build_url() {
    local source="$1" ra="$2" dec="$3" fov="$4" w="$5" h="$6" extra="$7"
    case "$source" in
        hips2fits)
            local hips_enc="${extra//\//%2F}"
            printf '%s' "https://alasky.cds.unistra.fr/hips-image-services/hips2fits?hips=${hips_enc}&width=${w}&height=${h}&fov=${fov}&projection=TAN&coordsys=icrs&ra=${ra}&dec=${dec}&format=fits"
            ;;
        skyview)
            printf '%s' "https://skyview.gsfc.nasa.gov/current/cgi/pskcall?Survey=${extra}&Position=${ra},${dec}&Size=${fov}&Pixels=${w}&Return=FITS"
            ;;
        legacysurvey)
            # pixscale in arcsec/px derived from the requested fov and width
            local ps; ps=$(awk -v f="$fov" -v w="$w" 'BEGIN{printf "%.4f", f*3600.0/w}')
            printf '%s' "https://www.legacysurvey.org/viewer/cutout.fits?ra=${ra}&dec=${dec}&layer=ls-dr10&pixscale=${ps}&size=${w}&bands=${extra}"
            ;;
        sdss)
            # extra = rerun/run/camcol/band/field
            local rerun run camcol band field
            IFS='/' read -r rerun run camcol band field <<< "$extra"
            printf '%s' "https://data.sdss.org/sas/dr17/eboss/photoObj/frames/${rerun}/${run}/${camcol}/frame-${band}-$(printf '%06d' "$run")-${camcol}-$(printf '%04d' "$field").fits.bz2"
            ;;
        *) return 1 ;;
    esac
}

# Pan-STARRS needs a filename lookup first, so it gets its own path.
fetch_panstarrs() {
    local ra="$1" dec="$2" w="$3" band="$4" dest="$5"
    local list fname
    list=$(curl -sS -L --max-time 120 \
        "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py?ra=${ra}&dec=${dec}&size=${w}&format=fits&filters=${band}" 2>/dev/null)
    fname=$(printf '%s\n' "$list" | awk 'NR>1 {print $8; exit}')
    [[ -n "$fname" ]] || { echo "no PS1 skycell"; return 1; }
    fetch "https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?red=${fname}&format=fits&size=${w}&ra=${ra}&dec=${dec}" "$dest"
}

if [[ "$LIST_ONLY" -eq 0 ]]; then
    mkdir -p "$OUT_DIR"
    TRUTH="$OUT_DIR/truth.tsv"
    printf '# id\ttier\tra_deg\tdec_deg\tfov_deg\tfile\ttruth_source\n' > "$TRUTH"
fi

n_ok=0; n_skip=0; n_fail=0
declare -a FAILED=()

while IFS=$'\t' read -r id tier ra dec fov w h source extra; do
    [[ -z "${id:-}" || "$id" == \#* ]] && continue
    # tolerate space-aligned manifests as well as strict TSV
    if [[ -z "${extra:-}" ]]; then
        read -r id tier ra dec fov w h source extra <<< "$id $tier $ra $dec $fov $w $h $source ${extra:-}"
    fi
    wanted "$id" "$tier" || continue

    if [[ "$LIST_ONLY" -eq 1 ]]; then
        printf '%-16s %s  ra=%-10s dec=%-9s fov=%-6s %sx%s  %s %s\n' \
            "$id" "$tier" "$ra" "$dec" "$fov" "$w" "$h" "$source" "$extra"
        continue
    fi

    dest="$OUT_DIR/${id}.fits"
    if [[ -f "$dest" && "$FORCE" -eq 0 ]]; then
        echo "skip    $id (already present)"
        n_skip=$((n_skip + 1))
    else
        printf 'fetch   %-16s %-13s ' "$id" "$source"
        rc=""
        case "$source" in
            panstarrs)
                rc=$(fetch_panstarrs "$ra" "$dec" "$w" "$extra" "$dest") ;;
            sdss)
                url=$(build_url "$source" "$ra" "$dec" "$fov" "$w" "$h" "$extra")
                rc=$(fetch "$url" "${dest}.bz2")
                if [[ "$rc" == "200" ]]; then
                    bunzip2 -f "${dest}.bz2" 2>/dev/null && rc=200 || rc="bunzip2-failed"
                fi ;;
            *)
                url=$(build_url "$source" "$ra" "$dec" "$fov" "$w" "$h" "$extra")
                rc=$(fetch "$url" "$dest") ;;
        esac

        if [[ "$rc" == "200" ]] && valid_fits "$dest"; then
            echo "ok ($(du -h "$dest" | cut -f1))"
            n_ok=$((n_ok + 1))
        else
            echo "FAILED (http=$rc)"
            rm -f "$dest" "${dest}.bz2"
            FAILED+=("$id")
            n_fail=$((n_fail + 1))
            continue
        fi
    fi

    if [[ "$tier" == "A" || "$tier" == "D" ]]; then
        truth_src="requested"
    else
        truth_src="header"
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$id" "$tier" "$ra" "$dec" "$fov" "${dest##*/}" "$truth_src" >> "$TRUTH"
done < <(grep -v '^[[:space:]]*#' "$MANIFEST" | grep -v '^[[:space:]]*$' | tr -s ' \t' '\t')

[[ "$LIST_ONLY" -eq 1 ]] && exit 0

echo ""
echo "=============================================================="
echo " ${n_ok} downloaded, ${n_skip} already present, ${n_fail} failed"
echo " Images:  $OUT_DIR"
echo " Truth:   $OUT_DIR/truth.tsv"
if [[ ${#FAILED[@]} -gt 0 ]]; then
    echo " Failed:  ${FAILED[*]}"
    echo " (transient service errors are common — re-run to retry just those)"
fi
echo "=============================================================="
