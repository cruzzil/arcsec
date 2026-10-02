#!/usr/bin/env bash
# Run the C test (tests/c/smoke.c, via tests/c_api.rs) in every configuration CI
# checks on Linux: shared and static linking, AddressSanitizer + LeakSanitizer on
# the C side, and the whole process under valgrind memcheck.
#
#   libarcsec/tests/c/run-checks.sh            # all of them
#   libarcsec/tests/c/run-checks.sh asan       # just one: plain|static|asan|valgrind
#
# Needs gcc or clang; valgrind for the last. Run from anywhere in the repository.
set -euo pipefail

cd "$(dirname "$0")/../../.."
run() {
  echo "── $1 ──"
  shift
  env "$@" cargo test --locked -p arcsec-capi --test c_api -- --nocapture
}

want=${1:-all}
if [[ $want == all || $want == plain ]]; then
  run "shared library" ARCSEC_C_LINK=shared
fi
if [[ $want == all || $want == static ]]; then
  run "static library" ARCSEC_C_LINK=static
fi
if [[ $want == all || $want == asan ]]; then
  # The C program is instrumented; the Rust library is not, but its allocations
  # go through the same malloc, so LeakSanitizer sees a handle that is never
  # freed and ASan sees C overrunning anything the library hands out.
  run "AddressSanitizer + LeakSanitizer" \
    ARCSEC_C_CFLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer -fno-sanitize-recover=all" \
    ASAN_OPTIONS="detect_leaks=1:abort_on_error=1" \
    UBSAN_OPTIONS="print_stacktrace=1"
fi
if [[ $want == all || $want == valgrind ]]; then
  if command -v valgrind > /dev/null; then
    run "valgrind memcheck" \
      ARCSEC_C_RUNNER="valgrind --error-exitcode=99 --leak-check=full --errors-for-leak-kinds=definite,indirect --track-origins=yes"
  else
    echo "valgrind not installed; skipped" >&2
  fi
fi
