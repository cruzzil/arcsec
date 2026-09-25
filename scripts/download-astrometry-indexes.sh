#!/usr/bin/env bash
# Superseded by `arcsec catalog`, which resumes interrupted downloads, reports
# progress, verifies the result and puts files where the solver already looks.
#
#   arcsec catalog install anet-4100     # Tycho-2, scales 07-19, ~160 MB
#   arcsec catalog install anet-5200     # Gaia LITE, 3 x 48 files, ~8.8 GB
#
# This shim forwards so existing muscle memory keeps working.
set -euo pipefail

ARCSEC="${ARCSEC:-$(dirname "$0")/../target/release/arcsec}"
[[ -x "$ARCSEC" ]] || ARCSEC="$(dirname "$0")/../target/debug/arcsec"
[[ -x "$ARCSEC" ]] || { echo "build arcsec first: cargo build --release" >&2; exit 1; }

echo "This script is superseded by 'arcsec catalog install'. Forwarding..." >&2
exec "$ARCSEC" catalog install "${@:-anet-4100}"
