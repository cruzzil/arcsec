#!/usr/bin/env bash
# Build libarcsec and lay it out as an installable, relocatable tree:
#
#   DEST/include/arcsec.h
#   DEST/lib/libarcsec.so.1, libarcsec.so -> libarcsec.so.1, libarcsec.a     Linux and other ELF
#   DEST/lib/libarcsec.dylib, libarcsec.a                                      macOS
#   DEST/bin/arcsec.dll, DEST/lib/arcsec.lib (import), arcsec_static.lib       Windows (MSVC)
#   DEST/bin/arcsec.dll, DEST/lib/libarcsec.dll.a, libarcsec.a                 Windows (GNU)
#   DEST/lib/pkgconfig/arcsec.pc
#   DEST/lib/cmake/arcsec/arcsecConfig.cmake, arcsecConfigVersion.cmake
#   DEST/share/doc/arcsec/README.md, LICENSE
#
# Usage: libarcsec/dist.sh [--target TRIPLE] [--profile PROFILE] DEST
#
#   --target    a Rust target triple; default: the host
#   --profile   cargo profile; default: dist (release without debug info)
#
# To install system-wide: libarcsec/dist.sh /tmp/stage && sudo cp -a /tmp/stage/. /usr/local/
# (then `sudo ldconfig` on Linux). The pkg-config and CMake files find everything
# relative to themselves, so the tree works wherever it is copied.
set -euo pipefail

target=""
profile="dist"
dest=""
while [ $# -gt 0 ]; do
  case "$1" in
    --target) target="$2"; shift 2 ;;
    --profile) profile="$2"; shift 2 ;;
    -h|--help) sed -n '2,24p' "$0"; exit 0 ;;
    -*) echo "unknown option $1" >&2; exit 2 ;;
    *) dest="$1"; shift ;;
  esac
done
[ -n "$dest" ] || { echo "usage: $0 [--target TRIPLE] [--profile PROFILE] DEST" >&2; exit 2; }

here="$(cd "$(dirname "$0")" && pwd)"
root="$(cd "$here/.." && pwd)"
cd "$root"

host="$(rustc -vV | sed -n 's/^host: //p')"
[ -n "$target" ] || target="$host"
version="$(sed -n 's/^version = "\(.*\)"/\1/p' Cargo.toml | head -1)"
abi="$(sed -n 's/^pub const ARCSEC_ABI_VERSION: u32 = \([0-9]*\);/\1/p' libarcsec/src/lib.rs)"

# One build gives the shared and static libraries; the static library's system
# dependencies come from rustc's --print native-static-libs.
log="$(mktemp)"
cargo rustc --locked -p libarcsec --lib --profile "$profile" --target "$target" \
  -- --print native-static-libs 2> >(tee "$log" >&2)
native="$(sed -n 's/.*native-static-libs: //p' "$log" | tail -1)"
rm -f "$log"
out="target/$target/$profile"

mkdir -p "$dest/include" "$dest/lib/pkgconfig" "$dest/lib/cmake/arcsec" "$dest/share/doc/arcsec"
dest="$(cd "$dest" && pwd)"
cp libarcsec/include/arcsec.h "$dest/include/"
cp libarcsec/README.md "$dest/share/doc/arcsec/README.md"
cp LICENSE "$dest/share/doc/arcsec/LICENSE"

libs_extra=""
case "$target" in
  *-windows-msvc)
    mkdir -p "$dest/bin"
    cp "$out/arcsec.dll" "$dest/bin/"
    cp "$out/arcsec.dll.lib" "$dest/lib/arcsec.lib"
    cp "$out/arcsec.lib" "$dest/lib/arcsec_static.lib"
    shared=arcsec.dll; implib=arcsec.lib; static=arcsec_static.lib; soname="" ;;
  *-windows-gnu*)
    mkdir -p "$dest/bin"
    cp "$out/arcsec.dll" "$dest/bin/"
    cp "$out/libarcsec.dll.a" "$dest/lib/"
    cp "$out/libarcsec.a" "$dest/lib/"
    shared=arcsec.dll; implib=libarcsec.dll.a; static=libarcsec.a; soname="" ;;
  *-apple-*)
    cp "$out/libarcsec.dylib" "$out/libarcsec.a" "$dest/lib/"
    # The install name is @rpath/libarcsec.dylib, so a program finds the library
    # through its own rpath; the .pc file adds one for libdir.
    libs_extra=' -Wl,-rpath,${libdir}'
    shared=libarcsec.dylib; implib=""; static=libarcsec.a; soname="@rpath/libarcsec.dylib" ;;
  *)
    soname="libarcsec.so.$abi"
    cp "$out/libarcsec.so" "$dest/lib/$soname"
    ln -sf "$soname" "$dest/lib/libarcsec.so"
    cp "$out/libarcsec.a" "$dest/lib/"
    shared="$soname"; implib=""; static=libarcsec.a ;;
esac

case "$target" in
  i686-*|i586-*|armv7-*|arm-*|thumbv7*) ptr=4 ;;
  *) ptr=8 ;;
esac
native_cmake="$(printf '%s' "$native" | sed 's/ \+/;/g')"
fill() {
  sed -e "s|@VERSION@|$version|g" -e "s|@ABI@|$abi|g" \
      -e "s|@LIBS_EXTRA@|$libs_extra|g" -e "s|@NATIVE_STATIC_LIBS@|$native|g" \
      -e "s|@NATIVE_STATIC_LIBS_CMAKE@|$native_cmake|g" \
      -e "s|@SHARED_NAME@|$shared|g" -e "s|@IMPLIB_NAME@|$implib|g" \
      -e "s|@STATIC_NAME@|$static|g" -e "s|@SONAME@|$soname|g" \
      -e "s|@SIZEOF_VOID_P@|$ptr|g" "$1"
}
fill libarcsec/pkg/arcsec.pc.in > "$dest/lib/pkgconfig/arcsec.pc"
fill libarcsec/pkg/arcsecConfig.cmake.in > "$dest/lib/cmake/arcsec/arcsecConfig.cmake"
fill libarcsec/pkg/arcsecConfigVersion.cmake.in > "$dest/lib/cmake/arcsec/arcsecConfigVersion.cmake"

echo "libarcsec $version (ABI $abi) for $target in $dest:"
(cd "$dest" && find . -type f -o -type l | sort | sed 's|^\./|  |')
