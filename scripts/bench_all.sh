#!/usr/bin/env bash
# Parallel benchmark: run arcsec and astap-cli on all resources/*.fits
# Usage: bench_all.sh [concurrency] [method]
#   method: quads (default), tetra, quads+blind
# Results land in /tmp/arcsec_bench/<method>/

set -euo pipefail

ARCSEC="./target/release/arcsec"
ASTAP="~/astap-cli"
DB="~/star_database"
INDEX_DIR="resources/indexes/4100"
RESOURCES="resources"
METHOD="${2:-quads}"
# blind mode is CPU-intensive (sorts 580K-entry index tables per image);
# keep concurrency low to avoid CPU saturation and memory pressure.
DEFAULT_CONCURRENCY=10
[[ "$METHOD" == "quads+blind" ]] && DEFAULT_CONCURRENCY=4
CONCURRENCY="${1:-$DEFAULT_CONCURRENCY}"
BENCH="/tmp/arcsec_bench/${METHOD}"

mkdir -p "$BENCH"

mapfile -t ALL_FITS < <(ls "$RESOURCES"/*.fits | sort)
TOTAL="${#ALL_FITS[@]}"
echo "=== Arcsec batch benchmark ==="
echo "Method: $METHOD   Files: $TOTAL   Concurrency: $CONCURRENCY"
echo "DB: $DB   INDEX: $INDEX_DIR"
echo ""

# ── run_arcsec_one <fits_file> ───────────────────────────────────────────────
run_arcsec_one() {
    local file="$1"
    local name="${file##*/}"
    local stem="${name%.fits}"
    local out="$BENCH/p_${stem}"

    local extra_args=""
    case "$METHOD" in
        quads+blind) extra_args="-i $INDEX_DIR" ;;
        tetra)       extra_args="--method tetra" ;;
        quads)       extra_args="" ;;
    esac

    # blind method gets more time (blind step ~2s + catalog ~1s = well under 120s)
    local per_timeout=60
    [[ "$METHOD" == "quads+blind" ]] && per_timeout=30

    local t0; t0=$(date +%s%3N)
    timeout "$per_timeout" "$ARCSEC" -f "$file" --fov 2.61 -r 180 \
        -d "$DB" -o "$out" $extra_args \
        > "$out.stdout" 2>&1
    local ec=$?
    local t1; t1=$(date +%s%3N)
    local ms=$(( t1 - t0 ))

    local solved="N" crval1="" crval2=""
    if [[ $ec -eq 0 ]] && grep -q "PLTSOLVD=T" "$out.ini" 2>/dev/null; then
        solved="Y"
        crval1=$(grep "^CRVAL1=" "$out.ini" | cut -d= -f2)
        crval2=$(grep "^CRVAL2=" "$out.ini" | cut -d= -f2)
    fi

    echo "$stem,$solved,$crval1,$crval2,$ms" > "$BENCH/p_${stem}.result"
    printf '  arcsec(%-10s) %-50s %s  %5dms\n' "$METHOD" "$name" "$solved" "$ms" >&2
}

# ── run_astap_one <fits_file> ────────────────────────────────────────────────
run_astap_one() {
    local file="$1"
    local name="${file##*/}"
    local stem="${name%.fits}"
    local out="$BENCH/a_${stem}"

    local t0; t0=$(date +%s%3N)
    timeout 120 "$ASTAP" -f "$file" -fov 2.61 -r 180 \
        -d "$DB" -o "$out" \
        > "$out.stdout" 2>&1
    local ec=$?
    local t1; t1=$(date +%s%3N)
    local ms=$(( t1 - t0 ))

    local solved="N" crval1="" crval2=""
    if grep -q "PLTSOLVD=T" "$out.ini" 2>/dev/null; then
        solved="Y"
        crval1=$(grep "^CRVAL1=" "$out.ini" | cut -d= -f2)
        crval2=$(grep "^CRVAL2=" "$out.ini" | cut -d= -f2)
    fi

    echo "$stem,$solved,$crval1,$crval2,$ms" > "$BENCH/a_${stem}.result"
    printf '  astap  %-50s %s  %5dms\n' "$name" "$solved" "$ms" >&2
}

export -f run_arcsec_one run_astap_one
export ARCSEC ASTAP DB INDEX_DIR BENCH METHOD

# ── Phase 1: run arcsec on all files ────────────────────────────────────────
echo "--- Phase 1: arcsec($METHOD) on all $TOTAL files ---"
printf '%s\n' "${ALL_FITS[@]}" \
    | xargs -P "$CONCURRENCY" -I{} bash -c 'run_arcsec_one "$@"' _ {}
echo ""

# ── Phase 2: run astap on files without existing reference .ini ──────────────
ASTAP_BENCH="/tmp/arcsec_bench/astap"
mkdir -p "$ASTAP_BENCH"
mapfile -t NO_INI < <(
    for f in "${ALL_FITS[@]}"; do
        ini="${f%.fits}.ini"
        ares="$ASTAP_BENCH/a_${f##*/}"
        ares="${ares%.fits}.result"
        [[ ! -f "$ini" && ! -f "$ares" ]] && echo "$f"
    done
)
echo "--- Phase 2: astap on ${#NO_INI[@]} files without existing reference ---"
if [[ ${#NO_INI[@]} -gt 0 ]]; then
    printf '%s\n' "${NO_INI[@]}" \
        | xargs -P "$CONCURRENCY" -I{} bash -c \
            'run_astap_one "$@"' _ {} 2>/dev/null || true
fi
echo ""

# ── Phase 3: aggregate results ───────────────────────────────────────────────
echo "--- Phase 3: aggregating results ---"

p_solved=0; p_fail=0
p_time_total=0; p_time_min=99999999; p_time_max=0
p_succ_time=0; p_succ_n=0; p_fail_time=0; p_fail_n=0
p_err_ra_sum=0; p_err_dec_sum=0; p_err_n=0
p_err_max_ra=0; p_err_max_dec=0
p_acc1=0  # within 1 arcsec
p_acc5=0  # 1-5 arcsec
p_acc_bad=0  # >5 arcsec

a_solved=0; a_fail=0

SUMMARY_CSV="$BENCH/summary.csv"
echo "file,p_solved,p_crval1,p_crval2,p_ms,a_solved,a_crval1,a_crval2,a_ms,err_ra_arcsec,err_dec_arcsec" > "$SUMMARY_CSV"

for f in "${ALL_FITS[@]}"; do
    name="${f##*/}"
    stem="${name%.fits}"
    pres="$BENCH/p_${stem}.result"

    p_s="N"; p_r1=""; p_r2=""; p_ms=0
    if [[ -f "$pres" ]]; then
        IFS=',' read -r _ p_s p_r1 p_r2 p_ms < "$pres"
    fi

    a_s="N"; a_r1=""; a_r2=""; a_ms=0
    orig_ini="${f%.fits}.ini"
    ares="$ASTAP_BENCH/a_${stem}.result"
    if [[ -f "$orig_ini" ]] && grep -q "PLTSOLVD=T" "$orig_ini" 2>/dev/null; then
        a_s="Y"
        a_r1=$(grep "^CRVAL1=" "$orig_ini" | cut -d= -f2)
        a_r2=$(grep "^CRVAL2=" "$orig_ini" | cut -d= -f2)
    elif [[ -f "$ares" ]]; then
        IFS=',' read -r _ a_s a_r1 a_r2 a_ms < "$ares"
    fi

    if [[ "$p_s" == "Y" ]]; then (( p_solved++ )); else (( p_fail++ )); fi
    if [[ "$a_s" == "Y" ]]; then (( a_solved++ )); else (( a_fail++ )); fi

    if [[ -n "$p_ms" && "$p_ms" -gt 0 ]]; then
        (( p_time_total += p_ms ))
        (( p_ms < p_time_min )) && p_time_min=$p_ms
        (( p_ms > p_time_max )) && p_time_max=$p_ms
        if [[ "$p_s" == "Y" ]]; then
            (( p_succ_time += p_ms )); (( p_succ_n++ ))
        else
            (( p_fail_time += p_ms )); (( p_fail_n++ ))
        fi
    fi

    err_ra=""; err_dec=""
    if [[ "$p_s" == "Y" && "$a_s" == "Y" && -n "$p_r1" && -n "$a_r1" ]]; then
        result=$(python3 -c "
import math
pr,ar = $p_r1, $a_r1
pd,ad = $p_r2, $a_r2
dra = abs(pr-ar)*math.cos(math.radians((pd+ad)/2))*3600
ddec = abs(pd-ad)*3600
sep = math.sqrt(dra**2+ddec**2)
print(f'{dra:.2f},{ddec:.2f},{sep:.2f}')
" 2>/dev/null || echo ",,")
        IFS=',' read -r era edec esep <<< "$result"
        err_ra="$era"; err_dec="$edec"

        if [[ -n "$era" && -n "$esep" ]]; then
            (( p_err_n++ ))
            p_err_ra_sum=$(python3 -c "print($p_err_ra_sum + $era)")
            p_err_dec_sum=$(python3 -c "print($p_err_dec_sum + $edec)")
            max_ra_check=$(python3 -c "print(1 if $era > $p_err_max_ra else 0)")
            [[ "$max_ra_check" == "1" ]] && p_err_max_ra=$era
            max_dec_check=$(python3 -c "print(1 if $edec > $p_err_max_dec else 0)")
            [[ "$max_dec_check" == "1" ]] && p_err_max_dec=$edec
            acc=$(python3 -c "
sep=$esep
if sep<1: print('1')
elif sep<5: print('5')
else: print('bad')
")
            case "$acc" in
                1)   (( p_acc1++ )) ;;
                5)   (( p_acc5++ )) ;;
                bad) (( p_acc_bad++ )) ;;
            esac
        fi
    fi

    echo "$stem,$p_s,$p_r1,$p_r2,$p_ms,$a_s,$a_r1,$a_r2,$a_ms,$err_ra,$err_dec" >> "$SUMMARY_CSV"
done

# ── Print summary ─────────────────────────────────────────────────────────────
p_total=$(( p_solved + p_fail ))
a_total=$(( a_solved + a_fail ))
p_pct=$(python3 -c "print(f'{100*$p_solved/$p_total:.1f}')")
a_pct=$(python3 -c "print(f'{100*$a_solved/$a_total:.1f}')")
p_avg=$(python3 -c "print(f'{$p_time_total/$p_total/1000:.2f}')")
p_succ_avg=$(python3 -c "print(f'{$p_succ_time/max(1,$p_succ_n)/1000:.2f}')")
p_fail_avg=$(python3 -c "print(f'{$p_fail_time/max(1,$p_fail_n)/1000:.2f}')")
mean_ra=$(python3 -c "print(f'{$p_err_ra_sum/max(1,$p_err_n):.2f}')")
mean_dec=$(python3 -c "print(f'{$p_err_dec_sum/max(1,$p_err_n):.2f}')")

echo "================================================================"
echo " BENCHMARK RESULTS — method=${METHOD} — $(date)"
echo "================================================================"
printf " %-35s  %d / %d  (%s%%)\n" "Arcsec ($METHOD) solved" "$p_solved" "$p_total" "$p_pct"
printf " %-35s  %d / %d  (%s%%)\n" "ASTAP reference solved" "$a_solved" "$a_total" "$a_pct"
echo ""
printf " %-35s  %ss avg total\n" "Arcsec solve time" "$p_avg"
printf " %-35s  %ss avg (n=%d)\n" "  successes" "$p_succ_avg" "$p_succ_n"
printf " %-35s  %ss avg (n=%d)\n" "  failures/timeouts" "$p_fail_avg" "$p_fail_n"
echo ""
if [[ "$p_err_n" -gt 0 ]]; then
    printf " %-35s  ΔRA=%.2f\"  ΔDec=%.2f\"  (n=%d)\n" \
        "Position error vs ASTAP (mean)" "$mean_ra" "$mean_dec" "$p_err_n"
    printf " %-35s  ΔRA=%.2f\"  ΔDec=%.2f\"\n" \
        "Max error" "$p_err_max_ra" "$p_err_max_dec"
    printf " %-35s  %d excellent / %d good / %d poor\n" \
        "Accuracy (<1\"/1-5\"/>5\")" "$p_acc1" "$p_acc5" "$p_acc_bad"
fi
echo ""
both=$(( p_solved < a_solved ? p_solved : a_solved ))
arcsec_only=$(( p_solved > both ? p_solved - both : 0 ))
astap_only=$(( a_solved > both ? a_solved - both : 0 ))
neither=$(( p_total - both - arcsec_only - astap_only ))
printf " %-35s  %d\n" "Both solved" "$both"
printf " %-35s  %d\n" "Arcsec only" "$arcsec_only"
printf " %-35s  %d\n" "ASTAP only" "$astap_only"
printf " %-35s  %d\n" "Neither" "$neither"
echo "================================================================"
echo "Full CSV: $SUMMARY_CSV"
