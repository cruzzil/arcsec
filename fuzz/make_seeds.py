#!/usr/bin/env python3
"""Write the seed corpus for the fuzz targets into fuzz/seeds/<target>/.

Every seed is synthetic, so the corpus carries no third-party data: small valid
FITS images (every BITPIX, a cube, an image extension, gzip and bzip2 wrapping,
and tile-compressed variants when CFITSIO's `fpack` is on PATH or given with
--fpack), XISF and ASDF images, ASTAP .1476/.290/.001 tiles, catalogue archives
(zip, .deb, tar, xz), raw pixel arrays and solver parameters. The `anet_index` and
`arcsecix` seeds are index files written by arcsec-core's own test fixtures
(`RawIndex::fits_bytes`, `build_index`), checked in as they are.

Usage: fuzz/make_seeds.py [--fpack PATH]
"""
import argparse, bz2, gzip, io, lzma, math, os, random, shutil, struct, subprocess
import tarfile, tempfile, zipfile, zlib

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "seeds")
random.seed(7)


def put(target, name, data):
    d = os.path.join(ROOT, target)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, name), "wb") as f:
        f.write(data)


def star_image(w, h, n=6):
    img = [[1000 + random.randint(0, 30) for _ in range(w)] for _ in range(h)]
    for _ in range(n):
        cx, cy = random.uniform(3, w - 4), random.uniform(3, h - 4)
        peak = random.uniform(3000, 20000)
        for y in range(h):
            for x in range(w):
                img[y][x] += peak * math.exp(-((x - cx) ** 2 + (y - cy) ** 2) / 3.0)
    return img


# ── FITS ─────────────────────────────────────────────────────────────────────

def card(key, value=None):
    if value is None:
        return key.ljust(80)
    if isinstance(value, bool):
        v = ("T" if value else "F").rjust(20)
    elif isinstance(value, str):
        v = ("'" + value.ljust(8) + "'").ljust(20)
    else:
        v = str(value).rjust(20)
    return f"{key:<8}= {v}".ljust(80)


def header(cards):
    h = "".join(card(*c) for c in cards) + card("END")
    return (h + " " * ((-len(h)) % 2880)).encode("ascii")


def pad(data):
    return data + b"\0" * ((-len(data)) % 2880)


WCS = [("CTYPE1", "RA---TAN"), ("CTYPE2", "DEC--TAN"), ("CRVAL1", 150.1), ("CRVAL2", 2.2),
       ("CRPIX1", 16.5), ("CRPIX2", 16.5), ("CD1_1", -0.0003), ("CD1_2", 0.0), ("CD2_1", 0.0),
       ("CD2_2", 0.0003), ("RA", 150.0), ("DEC", 2.0), ("FOCALLEN", 500.0), ("XPIXSZ", 3.76),
       ("XBINNING", 1)]
FORMATS = {8: ">B", 16: ">h", 32: ">i", 64: ">q", -32: ">f", -64: ">d"}


def pixels(bitpix, img, bzero):
    out = bytearray()
    for row in img:
        for v in row:
            if bitpix == 8:
                v = min(255, int(v / 100))
            elif bitpix == 16:
                v = max(-32768, min(32767, int(v) - (32768 if bzero else 0)))
            elif bitpix in (32, 64):
                v = int(v)
            out += struct.pack(FORMATS[bitpix], v)
    return bytes(out)


def fits(bitpix, w, h, planes=1, extra=(), bzero=False):
    img = star_image(w, h)
    cards = [("SIMPLE", True), ("BITPIX", bitpix), ("NAXIS", 3 if planes > 1 else 2),
             ("NAXIS1", w), ("NAXIS2", h)] + ([("NAXIS3", planes)] if planes > 1 else [])
    if bzero:
        cards += [("BZERO", 32768), ("BSCALE", 1)]
    return header(cards + list(extra)) + pad(pixels(bitpix, img, bzero) * planes)


def fits_seeds(fpack):
    seeds = {
        "u16_wcs.fits": fits(16, 32, 32, extra=WCS, bzero=True),
        "u8.fits": fits(8, 24, 20),
        "i32.fits": fits(32, 16, 16, extra=WCS[:6]),
        "f32.fits": fits(-32, 32, 24, extra=WCS),
        "f64.fits": fits(-64, 12, 12),
        "i64.fits": fits(64, 8, 8),
        "rgb.fits": fits(16, 16, 16, planes=3, extra=WCS[10:]),
        "header_only.fits": header([("SIMPLE", True), ("BITPIX", 8), ("NAXIS", 2),
                                    ("NAXIS1", 4), ("NAXIS2", 4)]) + b"\0" * 2880,
        "ext_image.fits": header([("SIMPLE", True), ("BITPIX", 8), ("NAXIS", 0), ("EXTEND", True)])
        + header([("XTENSION", "IMAGE"), ("BITPIX", -32), ("NAXIS", 2), ("NAXIS1", 16),
                  ("NAXIS2", 16), ("PCOUNT", 0), ("GCOUNT", 1)] + WCS)
        + pad(pixels(-32, star_image(16, 16), False)),
    }
    for name, data in seeds.items():
        put("image_fits", name, data)
    put("image_fits", "u16_wcs.fits.gz", gzip.compress(seeds["u16_wcs.fits"], mtime=0))
    put("image_fits", "f32.fits.bz2", bz2.compress(seeds["f32.fits"]))
    if not fpack:
        print("fpack not found: no tile-compressed FITS seeds")
        return
    with tempfile.TemporaryDirectory() as tmp:
        for opt, src in [("-r", "u16_wcs.fits"), ("-g", "u16_wcs.fits"), ("-h", "u16_wcs.fits"),
                         ("-p", "u8.fits"), ("-r", "f32.fits"), ("-g2", "i32.fits")]:
            inp, out = os.path.join(tmp, src), os.path.join(tmp, f"{opt}{src}.fz")
            with open(inp, "wb") as f:
                f.write(seeds[src])
            subprocess.run([fpack, opt, "-O", out, inp], check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            with open(out, "rb") as f:
                put("image_fits", f"fpack{opt}_{src}.fz", f.read())


# ── XISF ─────────────────────────────────────────────────────────────────────

def xisf(fmt, w, h, channels, data, keywords=(), compress=False, offset=None):
    attrs = f'geometry="{w}:{h}:{channels}" sampleFormat="{fmt}" colorSpace="{"RGB" if channels == 3 else "Gray"}"'
    if offset is not None:
        attrs += f' offset="{offset}"'
    block = zlib.compress(data) if compress else data
    kws = "".join(f'<FITSKeyword name="{k}" value="{v}" comment=""/>' for k, v in keywords)

    def header_xml(pos):
        loc = f'location="attachment:{pos}:{len(block)}"'
        if compress:
            loc += f' compression="zlib:{len(data)}"'
        return (f'<?xml version="1.0" encoding="UTF-8"?><xisf version="1.0" '
                f'xmlns="http://www.pixinsight.com/xisf"><Image {attrs} {loc}>{kws}</Image></xisf>').encode()

    # The data follows the header, whose length depends on the position it names:
    # leave room for the digits and round up.
    pos = (16 + len(header_xml(0)) + 32 + 63) // 64 * 64
    head = header_xml(pos)
    out = b"XISF0100" + struct.pack("<II", len(head), 0) + head
    out += b"\0" * (pos - len(out))
    return out + block


def xisf_seeds():
    img = star_image(17, 11)
    u16 = b"".join(struct.pack("<H", min(65535, int(v))) for r in img for v in r)
    f32 = b"".join(struct.pack("<f", v / 70000.0) for r in img for v in r)
    kws = [("RA", "150.1"), ("DEC", "'2.2'"), ("FOCALLEN", "500"), ("XPIXSZ", "3.76"),
           ("XBINNING", "2"), ("CRVAL1", "150.0"), ("CRVAL2", "2.0"), ("CRPIX1", "8.5"),
           ("CRPIX2", "5.5"), ("CD1_1", "-0.0003"), ("CD2_2", "0.0003")]
    put("image_xisf", "u16.xisf", xisf("UInt16", 17, 11, 1, u16))
    put("image_xisf", "u16_keywords.xisf", xisf("UInt16", 17, 11, 1, u16, kws, offset=100))
    put("image_xisf", "u16_zlib.xisf", xisf("UInt16", 17, 11, 1, u16, compress=True))
    put("image_xisf", "f32.xisf", xisf("Float32", 17, 11, 1, f32, kws[:4]))
    put("image_xisf", "u8_rgb.xisf", xisf("UInt8", 17, 11, 3, bytes(min(255, int(v / 100)) for _ in range(3) for r in img for v in r)))
    put("image_xisf", "f64.xisf", xisf("Float64", 4, 3, 1, struct.pack("<12d", *range(12))))
    put("image_xisf", "u32.xisf", xisf("UInt32", 4, 3, 1, struct.pack("<12I", *range(12))))


# ── ASDF ─────────────────────────────────────────────────────────────────────

def asdf_block(data, compression=b"\0\0\0\0"):
    import hashlib
    payload = zlib.compress(data) if compression == b"zlib" else data
    return (b"\xd3BLK" + struct.pack(">HI4sQQQ16s", 48, 0, compression, len(payload),
                                     len(payload), len(data), hashlib.md5(data).digest()) + payload)


def asdf(tree, blocks):
    head = (b"#ASDF 1.0.0\n#ASDF_STANDARD 1.5.0\n%YAML 1.1\n%TAG ! tag:stsci.edu:asdf/\n"
            b"--- !core/asdf-1.1.0\n")
    return head + tree.encode() + b"\n...\n" + b"".join(blocks)


def nd(src, dtype, shape):
    return f"!core/ndarray-1.0.0 {{source: {src}, datatype: {dtype}, byteorder: little, shape: {shape}}}"


def asdf_seeds():
    w, h = 16, 12
    pix = b"".join(struct.pack("<f", v) for r in star_image(w, h) for v in r)
    meta = "meta:\n  wcsinfo: {ra_ref: 150.1, dec_ref: 2.2, pixel_scale: 0.11}\n"
    put("image_asdf", "plain.asdf", asdf(f"data: {nd(0, 'float32', [h, w])}\n" + meta, [asdf_block(pix)]))
    put("image_asdf", "zlib.asdf", asdf(f"data: {nd(0, 'float32', [h, w])}\n" + meta, [asdf_block(pix, b"zlib")]))
    put("image_asdf", "roman.asdf", asdf(
        f"roman:\n  data: {nd(0, 'float32', [h, w])}\n  meta:\n    wcsinfo: {{ra_ref: 10.0, dec_ref: -5.0}}\n",
        [asdf_block(pix)]))
    put("image_asdf", "cube_u16.asdf", asdf(f"sci:\n  data: {nd(0, 'uint16', [2, 4, 4])}\n",
                                           [asdf_block(struct.pack("<32H", *range(32)))]))
    put("image_asdf", "nested.asdf", asdf(
        f"a:\n  b:\n    - {nd(0, 'int32', [3, 5])}\n    - {nd(1, 'float64', [2, 2])}\nfocallen: 500\nxpixsz: 3.8\n",
        [asdf_block(struct.pack("<15i", *range(15))), asdf_block(struct.pack("<4d", 1, 2, 3, 4))]))
    put("image_asdf", "inline.asdf", asdf(
        "data: !core/ndarray-1.0.0 {data: [[1, 2, 3], [4, 5, 6]], datatype: int64, shape: [2, 3]}\nra: 1.0\ndec: 2.0\n", []))


# ── Catalogue tiles ──────────────────────────────────────────────────────────

def tile(record_size, dec_deg=None, groups=6, per=8):
    body = bytearray()
    for g in range(groups):
        if dec_deg is None:
            hi, base = random.choice([0, 1, 255, 254, 129]), None
        else:
            base = int(dec_deg / 90.0 * 8388607)
            hi = (base >> 16) & 0xFF
        body += b"\xff\xff\xff" + bytes([(hi + 128) & 0xFF, 76 + g * 5]) + b"\0" * (record_size - 5)
        for _ in range(per):
            ra = struct.pack("<I", random.randrange(0, 0xFFFFFE))[:3]
            if base is None:
                dec = bytes([random.randrange(256), random.randrange(256)])
            else:
                dec = struct.pack("<i", base + random.randrange(-20000, 20000))[:2]
            body += ra + dec + b"\0" * (record_size - 5)
    return b"arcsec fuzz seed tile".ljust(109, b" ") + bytes([record_size]) + bytes(body)


def catalog_seeds():
    put("catalog_tiles", "d50_5", b"\x00" + tile(5))
    put("catalog_tiles", "d50_6", b"\x00" + tile(6))
    put("catalog_tiles", "d50_south", b"\x00" + tile(5, -89.0))
    put("catalog_tiles", "g05_6", b"\x01" + tile(6))
    put("catalog_tiles", "g05_south", b"\x01" + tile(6, -85.0))
    w08 = struct.pack("<I", 20) + b"".join(
        struct.pack("<fff", -15 + i * 4, random.uniform(0, 6.28), random.uniform(-1.5, -0.8)) for i in range(20))
    put("catalog_tiles", "w08", b"\x02" + w08)


# ── Archives ─────────────────────────────────────────────────────────────────

def ar(name, payload):
    hdr = (name.encode().ljust(16) + b"0".ljust(12) + b"0".ljust(6) + b"0".ljust(6) + b"100644".ljust(8)
           + str(len(payload)).encode().ljust(10) + b"`\n")
    return hdr + payload + (b"\n" if len(payload) % 2 else b"")


def archive_seeds():
    bio = io.BytesIO()
    with zipfile.ZipFile(bio, "w", zipfile.ZIP_DEFLATED) as z:
        for name, data in [("d50/d50_0101.1476", b"tile data"), ("../../escape.1476", b"x"),
                           ("/abs/path.1476", b"y"), ("skip_readme.txt", b"no"), ("dir/", b"")]:
            z.writestr(name, data)
    put("archive", "zip", b"\x00" + bio.getvalue())
    bio = io.BytesIO()
    with tarfile.open(fileobj=bio, mode="w", format=tarfile.GNU_FORMAT) as t:
        for name, data in [("opt/astap/v05_0101.290", b"t" * 300), ("../../up.290", b"u"), ("skip_copyright", b"c")]:
            ti = tarfile.TarInfo(name)
            ti.size = len(data)
            t.addfile(ti, io.BytesIO(data))
        for name, kind, target in [("evil_link", tarfile.SYMTYPE, "/etc/passwd"),
                                   ("opt/astap/hard", tarfile.LNKTYPE, "../../x")]:
            ti = tarfile.TarInfo(name)
            ti.type, ti.linkname = kind, target
            t.addfile(ti)
    tar = bio.getvalue()
    xz = lzma.compress(tar, format=lzma.FORMAT_XZ, check=lzma.CHECK_CRC64)
    put("archive", "deb_raw", b"\x01" + b"!<arch>\n" + ar("debian-binary", b"2.0\n")
        + ar("control.tar.xz", b"xx") + ar("data.tar.xz", xz))
    put("archive", "deb_plain_tar", b"\x02" + tar)
    put("archive", "deb_xz", b"\x03" + xz)
    put("archive", "deb_xz_crc32", b"\x03" + lzma.compress(tar, format=lzma.FORMAT_XZ, check=lzma.CHECK_CRC32))


# ── Pixels and solver parameters ─────────────────────────────────────────────

def detect_seeds():
    img = star_image(40, 40)
    u16 = b"".join(struct.pack("<H", min(65535, int(v))) for r in img for v in r)
    put("detect", "u16_40", bytes([39, 39, 1]) + u16)
    put("detect", "tetra", bytes([39, 39, 0x81]) + u16)
    put("detect", "f32_40", bytes([39, 39, 0]) + b"".join(struct.pack("<f", v / 70000.0) for r in img for v in r))
    put("detect", "u8_40_bin2", bytes([39, 39, 5]) + bytes(min(255, int(v / 100)) for r in img for v in r))
    put("detect", "nan", bytes([20, 20, 0]) + struct.pack("<ffff", math.nan, math.inf, -math.inf, 1e38))


def solve_params_seeds():
    def p(name, values, flags):
        put("solve_params", name, b"".join(struct.pack("<d", v) for v in values) + bytes([flags]))
    p("plain", [1.0, 0.3, 0.01, 0.02, 0.007, 0.8, 0.5, 20.0], 0x4B)
    p("pole", [6.0, 1.5707, 0.05, 0.1, 0.007, 1.5, 1.0, 2.0], 0x0C)
    p("wide", [3.0, -0.5, 0.3, 0.0, 0.01, 0.8, 5.0, 50.0], 0x82)
    p("nan", [math.nan, math.inf, 1e-300, 1e300, -1.0, math.nan, 0.0, -1.0], 0xFF)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--fpack", default=shutil.which("fpack"), help="CFITSIO's fpack")
    args = ap.parse_args()
    fits_seeds(args.fpack)
    xisf_seeds()
    asdf_seeds()
    catalog_seeds()
    archive_seeds()
    detect_seeds()
    solve_params_seeds()


if __name__ == "__main__":
    main()
