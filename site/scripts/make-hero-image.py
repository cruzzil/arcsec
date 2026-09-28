#!/usr/bin/env python3
"""Render a solved FITS image as the website's hero picture, with its RA/Dec grid.

Pure Python 3 (standard library only), so it runs anywhere without numpy or Pillow.

    arcsec -f sdss_d.fits -d ~/star_database -fov 0.164 -r 3 -o /tmp/frame
    site/scripts/make-hero-image.py sdss_d.fits /tmp/frame.wcs --out site/src/assets/hero

writes two files into --out:

  sky.png        the image, 2x2 binned and asinh-stretched to 8-bit greyscale. Astro
                 re-encodes it as WebP at build time, so this is only the source.
  sky-grid.json  the coordinate grid computed from the solved WCS (TAN projection):
                 SVG path data for the lines of constant RA and Dec, in the PNG's pixel
                 coordinates, label positions as percentages, and a short summary of
                 the solution for the caption.

The grid is drawn by the page (src/components/SkyImage.astro), not baked into the
pixels, so it stays sharp and its colours follow the site theme.

The published image is an SDSS DR17 r-band frame (run 5115, camcol 4, field 150) from
the benchmark corpus (scripts/test-images.tsv, entry sdss_d), a field in Coma Berenices. SDSS imagery is free to
use with credit: "Image: Sloan Digital Sky Survey".
"""

import argparse
import array
import json
import math
import os
import struct
import sys
import zlib


# ── FITS ────────────────────────────────────────────────────────────────────────


def parse_cards(cards_iter):
    cards = {}
    for card in cards_iter:
        key = card[:8].strip()
        if key == "END":
            return cards
        if card[8:10] == "= " and "'" not in card[10:]:
            cards[key] = card[10:].split("/")[0].strip()
    raise ValueError("FITS header has no END card")


def read_header(f):
    """Read one FITS header (80-byte cards in 2880-byte blocks); keyword -> raw value."""

    def cards():
        while True:
            block = f.read(2880)
            if len(block) < 2880:
                raise ValueError("truncated FITS header")
            for i in range(0, 2880, 80):
                yield block[i : i + 80].decode("ascii", "replace")

    return parse_cards(cards())


def num(cards, key, default=None):
    if key not in cards:
        if default is None:
            raise KeyError(key)
        return default
    return float(cards[key].replace("D", "E"))


def read_image(path):
    """Return (width, height, pixels) of the primary HDU as a flat list, row 0 = FITS row 1."""
    with open(path, "rb") as f:
        cards = read_header(f)
        bitpix = int(num(cards, "BITPIX"))
        w = int(num(cards, "NAXIS1"))
        h = int(num(cards, "NAXIS2"))
        bzero = num(cards, "BZERO", 0.0)
        bscale = num(cards, "BSCALE", 1.0)
        code = {8: "B", 16: "h", 32: "i", -32: "f", -64: "d"}[bitpix]
        data = array.array(code)
        data.frombytes(f.read(w * h * abs(bitpix) // 8))
    if sys.byteorder == "little" and bitpix != 8:
        data.byteswap()
    if bscale != 1.0 or bzero != 0.0:
        data = [v * bscale + bzero for v in data]
    return w, h, data


def read_wcs(path):
    # arcsec's .wcs, like ASTAP's, is a FITS header with a newline after each card.
    with open(path, encoding="ascii", errors="replace") as f:
        cards = parse_cards(line.rstrip("\n") for line in f)
    return {k: num(cards, k) for k in ("CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CD1_1", "CD1_2", "CD2_1", "CD2_2")}


# ── TAN projection ──────────────────────────────────────────────────────────────


class Tan:
    """The gnomonic projection of a FITS WCS, in 1-based FITS pixel coordinates."""

    def __init__(self, wcs):
        self.w = wcs
        self.ra0 = math.radians(wcs["CRVAL1"])
        self.dec0 = math.radians(wcs["CRVAL2"])
        det = wcs["CD1_1"] * wcs["CD2_2"] - wcs["CD1_2"] * wcs["CD2_1"]
        self.inv = (wcs["CD2_2"] / det, -wcs["CD1_2"] / det, -wcs["CD2_1"] / det, wcs["CD1_1"] / det)

    def pix_to_sky(self, x, y):
        w = self.w
        dx, dy = x - w["CRPIX1"], y - w["CRPIX2"]
        xi = math.radians(w["CD1_1"] * dx + w["CD1_2"] * dy)
        eta = math.radians(w["CD2_1"] * dx + w["CD2_2"] * dy)
        d = math.cos(self.dec0) - eta * math.sin(self.dec0)
        ra = self.ra0 + math.atan2(xi, d)
        dec = math.atan2(math.sin(self.dec0) + eta * math.cos(self.dec0), math.hypot(xi, d))
        return math.degrees(ra) % 360.0, math.degrees(dec)

    def sky_to_pix(self, ra, dec):
        ra, dec = math.radians(ra), math.radians(dec)
        da = ra - self.ra0
        cosc = math.sin(self.dec0) * math.sin(dec) + math.cos(self.dec0) * math.cos(dec) * math.cos(da)
        if cosc <= 0:
            return None
        xi = math.degrees(math.cos(dec) * math.sin(da) / cosc)
        eta = math.degrees((math.cos(self.dec0) * math.sin(dec) - math.sin(self.dec0) * math.cos(dec) * math.cos(da)) / cosc)
        a, b, c, d = self.inv
        return self.w["CRPIX1"] + a * xi + b * eta, self.w["CRPIX2"] + c * xi + d * eta


# ── Image ───────────────────────────────────────────────────────────────────────


def binned(w, h, data, b):
    """Average b x b blocks; the result is flipped so row 0 is the top of the picture."""
    bw, bh = w // b, h // b
    out = [0.0] * (bw * bh)
    inv = 1.0 / (b * b)
    for by in range(bh):
        rows = [data[(by * b + j) * w : (by * b + j) * w + bw * b] for j in range(b)]
        dst = (bh - 1 - by) * bw
        for bx in range(bw):
            s = 0.0
            for r in rows:
                s += sum(r[bx * b : bx * b + b])
            out[dst + bx] = s * inv
    return bw, bh, out


def stretch(pixels):
    """Asinh stretch to 0-255, black point just above the sky so the background is quiet."""
    finite = [v for v in pixels[:: max(1, len(pixels) // 200_000)] if v == v]
    finite.sort()
    med = finite[len(finite) // 2]
    mad = sorted(abs(v - med) for v in finite)[len(finite) // 2]
    sigma = 1.4826 * mad or 1.0
    black = med + 0.8 * sigma
    white = finite[int(len(finite) * 0.9995)]
    soft = 6.0 * sigma  # where the curve turns from linear to logarithmic
    top = math.asinh((white - black) / soft)
    out = bytearray(len(pixels))
    for i, v in enumerate(pixels):
        if v != v or v <= black:
            continue
        t = math.asinh((v - black) / soft) / top
        out[i] = 255 if t >= 1 else int(t * 255 + 0.5)
    return out


def write_png(path, w, h, grey):
    raw = bytearray()
    prev = bytearray(w)
    for y in range(h):
        row = grey[y * w : (y + 1) * w]
        # Filter 2 (Up) suits a mostly black sky with sparse stars.
        raw.append(2)
        raw.extend((row[i] - prev[i]) & 0xFF for i in range(w))
        prev = row

    def chunk(tag, body):
        c = struct.pack(">I", len(body)) + tag + body
        return c + struct.pack(">I", zlib.crc32(tag + body) & 0xFFFFFFFF)

    with open(path, "wb") as f:
        f.write(b"\x89PNG\r\n\x1a\n")
        f.write(chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 0, 0, 0, 0)))
        f.write(chunk(b"IDAT", zlib.compress(bytes(raw), 9)))
        f.write(chunk(b"IEND", b""))


# ── Grid ────────────────────────────────────────────────────────────────────────

# Candidate spacings: Dec in arcseconds, RA in seconds of time.
DEC_STEPS = [10, 15, 20, 30, 60, 120, 300, 600, 900, 1200, 1800, 3600, 7200, 18000, 36000]
RA_STEPS = [1, 2, 5, 10, 15, 20, 30, 60, 120, 300, 600, 900, 1200, 1800, 3600, 7200]


def nice_step(span, steps, want=4):
    for s in steps:
        if span / s <= want:
            return s
    return steps[-1]


def fmt_ra(sec):
    sec = round(sec) % 86400
    h, m, s = sec // 3600, sec // 60 % 60, sec % 60
    return f"{h}h{m:02d}m{s:02d}s" if s else f"{h}h{m:02d}m"


def fmt_dec(arcsec):
    sign = "−" if arcsec < 0 else "+"
    a = round(abs(arcsec))
    d, m, s = a // 3600, a // 60 % 60, a % 60
    txt = f"{sign}{d}°{m:02d}′"
    return txt + f"{s:02d}″" if s else txt


def grid(tan, w, h, b):
    """Grid lines in binned, top-down pixel coordinates."""
    bw, bh = w / b, h / b

    def to_view(x, y):  # FITS 1-based pixel centre -> binned picture coordinates
        return (x - 0.5) / b, bh - (y - 0.5) / b

    # Sky extent from a ring of points around the border.
    ring = []
    for i in range(41):
        t = i / 40
        ring += [(0.5 + t * w, 0.5), (0.5 + t * w, h + 0.5), (0.5, 0.5 + t * h), (w + 0.5, 0.5 + t * h)]
    sky = [tan.pix_to_sky(x, y) for x, y in ring]
    ra_c, dec_c = tan.pix_to_sky(w / 2 + 0.5, h / 2 + 0.5)
    ras = [((r - ra_c + 180) % 360) - 180 for r, _ in sky]  # unwrap around the centre
    decs = [d for _, d in sky]
    dec_lo, dec_hi = min(decs), max(decs)
    ra_lo, ra_hi = ra_c + min(ras), ra_c + max(ras)

    dec_step = nice_step((dec_hi - dec_lo) * 3600, DEC_STEPS)
    ra_step = nice_step((ra_hi - ra_lo) * 240, RA_STEPS)

    lines = []
    samples = 64

    def trace(points):
        pts = [p for p in points if p is not None]
        return [to_view(x, y) for x, y in pts]

    d = math.ceil(dec_lo * 3600 / dec_step) * dec_step
    while d <= dec_hi * 3600:
        pts = trace(tan.sky_to_pix(ra_lo + (ra_hi - ra_lo) * i / samples, d / 3600) for i in range(samples + 1))
        lines.append(("dec", fmt_dec(d), pts))
        d += dec_step
    r = math.ceil(ra_lo * 240 / ra_step) * ra_step
    while r <= ra_hi * 240:
        pts = trace(tan.sky_to_pix(r / 240, dec_lo + (dec_hi - dec_lo) * i / samples) for i in range(samples + 1))
        lines.append(("ra", fmt_ra(r), pts))
        r += ra_step

    def inside(p):
        return 0 <= p[0] <= bw and 0 <= p[1] <= bh

    def edge_of(p):
        dists = {"left": p[0], "right": bw - p[0], "top": p[1], "bottom": bh - p[1]}
        return min(dists, key=dists.get)

    # Clip each line to the picture and find where it leaves it, for the label.
    clipped = []
    for kind, label, pts in lines:
        keep = [p for p in pts if inside(p)]
        if len(keep) < 2:
            continue
        clipped.append((kind, label, keep, [(edge_of(keep[0]), keep[0]), (edge_of(keep[-1]), keep[-1])]))

    # Each family is labelled along one edge: RA along the one it crosses most (bottom
    # first), Dec along the best of what is left.
    labels = []
    used = set()
    for kind in ("ra", "dec"):
        fam = [c for c in clipped if c[0] == kind]
        counts = {}
        for c in fam:
            for e, _ in c[3]:
                counts[e] = counts.get(e, 0) + 1
        order = ["bottom", "left", "top", "right"]
        edge = max((e for e in order if e not in used), key=lambda e: (counts.get(e, 0), -order.index(e)))
        used.add(edge)
        for _, label, _, ends in fam:
            for e, p in ends:
                if e == edge:
                    labels.append({"text": label, "edge": edge, "x": round(100 * p[0] / bw, 2), "y": round(100 * p[1] / bh, 2)})

    paths = [
        {"kind": kind, "d": "M" + " L".join(f"{x:.1f},{y:.1f}" for x, y in pts)}
        for kind, _, pts, _ in clipped
    ]
    return paths, labels, (ra_c, dec_c)


# ── Main ────────────────────────────────────────────────────────────────────────


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("fits")
    ap.add_argument("wcs", help="the .wcs file arcsec wrote for this image")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--bin", type=int, default=2, help="binning factor (default 2)")
    args = ap.parse_args()

    w, h, data = read_image(args.fits)
    wcs = read_wcs(args.wcs)
    tan = Tan(wcs)

    bw, bh, pix = binned(w, h, data, args.bin)
    grey = stretch(pix)
    os.makedirs(args.out, exist_ok=True)
    write_png(os.path.join(args.out, "sky.png"), bw, bh, grey)

    paths, labels, (ra_c, dec_c) = grid(tan, w, h, args.bin)
    scale = 3600 * math.sqrt(abs(wcs["CD1_1"] * wcs["CD2_2"] - wcs["CD1_2"] * wcs["CD2_1"]))
    meta = {
        "width": bw,
        "height": bh,
        "centre": {"ra": fmt_ra(ra_c * 240), "dec": fmt_dec(dec_c * 3600)},
        "scale": f"{scale:.3f}″/px",
        "field": f"{w * scale / 60:.1f}′ × {h * scale / 60:.1f}′",
        "paths": paths,
        "labels": labels,
    }
    with open(os.path.join(args.out, "sky-grid.json"), "w", encoding="utf-8") as f:
        json.dump(meta, f, ensure_ascii=False, indent=1)
        f.write("\n")
    print(f"{bw}x{bh}, centre {meta['centre']}, {meta['scale']}, {len(paths)} lines, {len(labels)} labels")


if __name__ == "__main__":
    main()
