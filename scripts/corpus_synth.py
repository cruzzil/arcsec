"""Synthetic tiers of the benchmark corpus: tier S (simulated camera artefacts) and the
generated tier-D negative controls. Standard library only; deterministic for a given
parent image and seed.

A recipe is a comma-separated list of operations, applied left to right to the parent
image (or to a generated one). Each operation keeps the ground truth exact: geometric
operations transform the WCS with the pixels, and `distort` encodes its own radial
distortion as SIP terms, so the derived truth is as exact as the parent's.

    geometry     crop:WxH+X+Y  bin:N  flipx  flipy  transpose  rot90:K
                 distort:K           radial distortion r' = r(1 + K r^2) about CRPIX
                                     (K > 0 barrel, K < 0 pincushion; px^-2)
    optics/sky   vignette:F          multiply by 1 - F (r / r_corner)^2
                 gradient:G:ANGLE    add a linear sky ramp of G x (image robust range)
                 sky:L               add a pedestal of L x (sky noise)
                 glow:A              amp glow: exponential brightening from one corner
                 clouds:A:SCALE      multiply by a smooth random transmission (1-A..1)
                 blur:SIGMA          Gaussian-ish blur (three box passes), px
                 trail:LEN:ANGLE     star trailing: mean of LEN+1 shifted copies
    detector     noise:S             add Gaussian noise, S x the measured sky noise
                 hot:N               N hot pixels at full scale
                 hotcol:N            N hot columns
                 satellite:N         N bright straight trails
                 saturate:P          clip everything above the P-th percentile
                 bayer:PAT           colour-camera raw: 2x2 gains (RGGB, GRBG, ...)
    output       u8 u16 i32 f32 f64  sample type (u16 writes BZERO = 32768, as cameras do)
                 gzip                write .fits.gz
                 xisf                write PixInsight XISF instead of FITS
                 keepwcs             leave the truth WCS in the image header (default:
                                     strip it, as a raw camera frame has none; the truth
                                     then lives in <id>.truth)

Generators (tier D, recipe starts with gen=...):
    gen=noise          Gaussian noise on a flat sky, no sources
    gen=fakestars:N    N Gaussian stars at random positions: star-like, no real asterisms
    gen=flat           a smooth gradient and nothing else
    shuffle:B          (an op) permute BxB blocks of a real field: real stars, no real
                       asterisms
"""

import array
import gzip
import math
import os
import random
import struct

import fitslite
from fitslite import Wcs


class Img:
    def __init__(self, w, h, data, wcs, cards):
        self.w, self.h, self.data, self.wcs, self.cards = w, h, data, wcs, cards
        self.bitpix = -32
        self.fmt = "fits"
        self.keepwcs = False
        self.bayer = None
        self.saturate = None


# ── statistics ──────────────────────────────────────────────────────────────


def robust_stats(img, samples=200000, rng=None):
    """(median, sigma, p99.9) from a strided sample."""
    n = img.w * img.h
    step = max(1, n // samples)
    v = sorted(img.data[i] for i in range(0, n, step))
    m = v[len(v) // 2]
    mad = sorted(abs(x - m) for x in v)[len(v) // 2]
    return m, max(1.4826 * mad, 1e-6), v[min(len(v) - 1, int(0.999 * len(v)))], v[-1], v[0]


# ── geometry ────────────────────────────────────────────────────────────────


def _remap_linear(img, m, crpix_new):
    """Replace the linear part: u_old = M u_new."""
    w = img.wcs
    cd = w.cd
    ncd = [[cd[0][0] * m[0][0] + cd[0][1] * m[1][0], cd[0][0] * m[0][1] + cd[0][1] * m[1][1]],
           [cd[1][0] * m[0][0] + cd[1][1] * m[1][0], cd[1][0] * m[0][1] + cd[1][1] * m[1][1]]]
    img.wcs = Wcs(crpix_new[0], crpix_new[1], w.crval1, w.crval2, ncd, proj=w.proj)


def _need_linear(img, op):
    if img.wcs.distorted:
        raise ValueError(f"{op}: apply orientation/binning before distort (parent has SIP/TPV)")


def op_crop(img, arg):
    size, x0, y0 = arg.split("+")
    cw, ch = (int(v) for v in size.split("x"))
    x0, y0 = int(x0), int(y0)
    if x0 + cw > img.w or y0 + ch > img.h:
        raise ValueError("crop outside the image")
    out = array.array("f", bytes(4 * cw * ch))
    for r in range(ch):
        s = (y0 + r) * img.w + x0
        out[r * cw:(r + 1) * cw] = img.data[s:s + cw]
    w = img.wcs
    img.wcs = Wcs(w.crpix1 - x0, w.crpix2 - y0, w.crval1, w.crval2, w.cd, proj=w.proj,
                  sip=w.sip, tpv=w.tpv)
    img.w, img.h, img.data = cw, ch, out


def op_bin(img, arg):
    _need_linear(img, "bin")
    n = int(arg)
    nw, nh = img.w // n, img.h // n
    out = array.array("f", bytes(4 * nw * nh))
    inv = 1.0 / (n * n)
    d, W = img.data, img.w
    for y in range(nh):
        for x in range(nw):
            s = 0.0
            for j in range(n):
                base = (y * n + j) * W + x * n
                for i in range(n):
                    s += d[base + i]
            out[y * nw + x] = s * inv
    w = img.wcs
    cd = [[c * n for c in row] for row in w.cd]
    img.wcs = Wcs((w.crpix1 - 0.5) / n + 0.5, (w.crpix2 - 0.5) / n + 0.5, w.crval1, w.crval2,
                  cd, proj=w.proj)
    img.w, img.h, img.data = nw, nh, out


def op_flipx(img, arg=None):
    _need_linear(img, "flipx")
    W = img.w
    out = array.array("f", img.data)
    for r in range(img.h):
        row = img.data[r * W:(r + 1) * W]
        row.reverse()
        out[r * W:(r + 1) * W] = row
    img.data = out
    _remap_linear(img, [[-1, 0], [0, 1]], (W + 1 - img.wcs.crpix1, img.wcs.crpix2))


def op_flipy(img, arg=None):
    _need_linear(img, "flipy")
    W, H = img.w, img.h
    out = array.array("f", bytes(4 * W * H))
    for r in range(H):
        out[r * W:(r + 1) * W] = img.data[(H - 1 - r) * W:(H - r) * W]
    img.data = out
    _remap_linear(img, [[1, 0], [0, -1]], (img.wcs.crpix1, H + 1 - img.wcs.crpix2))


def op_transpose(img, arg=None):
    _need_linear(img, "transpose")
    W, H = img.w, img.h
    out = array.array("f", bytes(4 * W * H))
    d = img.data
    for y in range(H):
        for x in range(W):
            out[x * H + y] = d[y * W + x]
    img.data = out
    img.w, img.h = H, W
    _remap_linear(img, [[0, 1], [1, 0]], (img.wcs.crpix2, img.wcs.crpix1))


def op_rot90(img, arg="1"):
    # counter-clockwise quarter turns, as transpose + flips
    for _ in range(int(arg) % 4):
        op_transpose(img)
        op_flipx(img)


def _bilinear(d, W, H, x, y):
    """0-based pixel coordinates; 0 outside."""
    if x < 0 or y < 0 or x > W - 1 or y > H - 1:
        return 0.0
    ix, iy = int(x), int(y)
    if ix >= W - 1:
        ix = W - 2
    if iy >= H - 1:
        iy = H - 2
    fx, fy = x - ix, y - iy
    i = iy * W + ix
    return ((d[i] * (1 - fx) + d[i + 1] * fx) * (1 - fy)
            + (d[i + W] * (1 - fx) + d[i + W + 1] * fx) * fy)


def op_distort(img, arg):
    """Radial distortion about CRPIX, encoded exactly as SIP."""
    _need_linear(img, "distort")
    k = float(arg)
    W, H = img.w, img.h
    c1, c2 = img.wcs.crpix1, img.wcs.crpix2
    d = img.data
    med = robust_stats(img)[0]
    out = array.array("f", bytes(4 * W * H))
    for y in range(H):
        v = (y + 1) - c2
        v2 = v * v
        base = y * W
        for x in range(W):
            u = (x + 1) - c1
            f = 1.0 + k * (u * u + v2)
            sx = c1 + u * f - 1.0
            sy = c2 + v * f - 1.0
            if sx < 0 or sy < 0 or sx > W - 1 or sy > H - 1:
                out[base + x] = med
            else:
                out[base + x] = _bilinear(d, W, H, sx, sy)
    img.data = out
    w = img.wcs
    a = {(3, 0): k, (1, 2): k}
    b = {(2, 1): k, (0, 3): k}
    img.wcs = Wcs(w.crpix1, w.crpix2, w.crval1, w.crval2, w.cd, proj=w.proj, sip=(a, b))


# ── optics / sky / detector ─────────────────────────────────────────────────


def op_vignette(img, arg):
    f = float(arg)
    W, H = img.w, img.h
    cx, cy = (W - 1) / 2.0, (H - 1) / 2.0
    rmax2 = cx * cx + cy * cy
    med = robust_stats(img)[0]
    d = img.data
    for y in range(H):
        dy2 = (y - cy) ** 2
        base = y * W
        for x in range(W):
            g = 1.0 - f * ((x - cx) ** 2 + dy2) / rmax2
            # vignetting dims sky and stars alike; keep the zero point at the sky
            d[base + x] = (d[base + x] - med) * g + med * g


def op_gradient(img, arg):
    parts = arg.split(":")
    g = float(parts[0])
    ang = math.radians(float(parts[1]) if len(parts) > 1 else 30.0)
    m, s, p999, mx, mn = robust_stats(img)
    amp = g * (p999 - m + 10 * s)
    W, H = img.w, img.h
    ca, sa = math.cos(ang), math.sin(ang)
    norm = 1.0 / max(W, H)
    d = img.data
    for y in range(H):
        base = y * W
        ty = y * sa
        for x in range(W):
            d[base + x] += amp * ((x * ca + ty) * norm)


def op_sky(img, arg):
    m, s, *_ = robust_stats(img)
    add = float(arg) * s
    d = img.data
    for i in range(len(d)):
        d[i] += add


def op_glow(img, arg):
    a = float(arg)
    m, s, p999, *_ = robust_stats(img)
    amp = a * (p999 - m)
    W, H = img.w, img.h
    L = 0.12 * max(W, H)
    d = img.data
    for y in range(H):
        base = y * W
        for x in range(W):
            r = math.hypot(x, y)
            if r < 6 * L:
                d[base + x] += amp * math.exp(-r / L)


def op_clouds(img, arg, rng):
    parts = arg.split(":")
    a = float(parts[0])
    scale = float(parts[1]) if len(parts) > 1 else 0.3
    W, H = img.w, img.h
    waves = []
    for _ in range(6):
        k = 2 * math.pi / (scale * max(W, H) * rng.uniform(0.6, 1.6))
        th = rng.uniform(0, 2 * math.pi)
        waves.append((k * math.cos(th), k * math.sin(th), rng.uniform(0, 2 * math.pi)))
    m = robust_stats(img)[0]
    d = img.data
    for y in range(H):
        base = y * W
        for x in range(W):
            t = 0.0
            for kx, ky, ph in waves:
                t += math.cos(kx * x + ky * y + ph)
            t = 0.5 + t / 12.0  # 0..1
            trans = 1.0 - a * t
            # scattered light: the cloud also glows a little
            d[base + x] = (d[base + x] - m) * trans + m * (1.0 + 0.3 * a * t)


def _box_rows(d, W, H, r):
    out = array.array("f", bytes(4 * W * H))
    n = 2 * r + 1
    for y in range(H):
        base = y * W
        s = 0.0
        for x in range(-r, r + 1):
            s += d[base + min(max(x, 0), W - 1)]
        for x in range(W):
            out[base + x] = s / n
            s += d[base + min(x + r + 1, W - 1)] - d[base + max(x - r, 0)]
    return out


def _box_cols(d, W, H, r):
    out = array.array("f", bytes(4 * W * H))
    n = 2 * r + 1
    for x in range(W):
        s = 0.0
        for y in range(-r, r + 1):
            s += d[min(max(y, 0), H - 1) * W + x]
        for y in range(H):
            out[y * W + x] = s / n
            s += d[min(y + r + 1, H - 1) * W + x] - d[max(y - r, 0) * W + x]
    return out


def op_blur(img, arg):
    sigma = float(arg)
    # three box passes approximate a Gaussian: box width from sigma
    r = max(1, int(round((math.sqrt(12 * sigma * sigma / 3 + 1) - 1) / 2)))
    d = img.data
    for _ in range(3):
        d = _box_rows(d, img.w, img.h, r)
        d = _box_cols(d, img.w, img.h, r)
    img.data = d


def op_trail(img, arg):
    parts = arg.split(":")
    ln = int(parts[0])
    ang = math.radians(float(parts[1]) if len(parts) > 1 else 20.0)
    W, H = img.w, img.h
    src = img.data
    acc = array.array("f", bytes(4 * W * H))
    shifts = [(int(round(t * math.cos(ang))), int(round(t * math.sin(ang)))) for t in range(ln + 1)]
    inv = 1.0 / len(shifts)
    for dx, dy in shifts:
        for y in range(H):
            sy = min(max(y - dy, 0), H - 1)
            sb = sy * W
            b = y * W
            if dx >= 0:
                for x in range(W):
                    acc[b + x] += src[sb + (x - dx if x >= dx else 0)] * inv
            else:
                for x in range(W):
                    xx = x - dx
                    acc[b + x] += src[sb + (xx if xx < W else W - 1)] * inv
    img.data = acc
    # the photocentre moves by the mean shift; move the truth with it
    mx = sum(s[0] for s in shifts) * inv
    my = sum(s[1] for s in shifts) * inv
    w = img.wcs
    img.wcs = Wcs(w.crpix1 + mx, w.crpix2 + my, w.crval1, w.crval2, w.cd, proj=w.proj,
                  sip=w.sip, tpv=w.tpv)


def op_noise(img, arg, rng):
    s = float(arg) * robust_stats(img)[1]
    d = img.data
    g = rng.gauss
    for i in range(len(d)):
        d[i] += g(0.0, s)


def op_hot(img, arg, rng):
    m, s, p999, mx, mn = robust_stats(img)
    top = max(mx, p999 * 4, m + 1000 * s)
    for _ in range(int(arg)):
        i = rng.randrange(len(img.data))
        img.data[i] = top * rng.uniform(0.3, 1.0)


def op_hotcol(img, arg, rng):
    m, s, *_ = robust_stats(img)
    for _ in range(int(arg)):
        x = rng.randrange(img.w)
        y0 = rng.randrange(img.h // 2)
        add = rng.uniform(5, 40) * s
        for y in range(y0, img.h):
            img.data[y * img.w + x] += add


def op_satellite(img, arg, rng):
    m, s, p999, *_ = robust_stats(img)
    W, H = img.w, img.h
    for _ in range(int(arg)):
        x0, y0 = rng.uniform(0, W), rng.uniform(0, H)
        ang = rng.uniform(0, math.pi)
        dx, dy = math.cos(ang), math.sin(ang)
        amp = rng.uniform(0.2, 1.0) * (p999 - m)
        L = int(2 * math.hypot(W, H))
        for t in range(-L // 2, L // 2):
            x, y = x0 + t * dx * 0.5, y0 + t * dy * 0.5
            ix, iy = int(x), int(y)
            for ox in (-1, 0, 1):
                xx, yy = ix + ox, iy
                if 0 <= xx < W and 0 <= yy < H:
                    img.data[yy * W + xx] = max(img.data[yy * W + xx], m + amp * (1 - 0.4 * abs(ox)))


def op_saturate(img, arg):
    p = float(arg)
    n = len(img.data)
    step = max(1, n // 400000)
    v = sorted(img.data[i] for i in range(0, n, step))
    lim = v[min(len(v) - 1, int(p / 100.0 * len(v)))]
    d = img.data
    for i in range(n):
        if d[i] > lim:
            d[i] = lim
    img.saturate = lim


BAYER_GAINS = {"R": 0.55, "G": 1.0, "B": 0.42}


def op_bayer(img, arg):
    pat = arg.upper()
    if len(pat) != 4 or set(pat) - set("RGB"):
        raise ValueError(f"bad bayer pattern {arg}")
    m = robust_stats(img)[0]
    W = img.w
    d = img.data
    # FITS row 1 is the bottom row; BAYERPAT is conventionally given top-down
    # (ROWORDER='TOP-DOWN' written below), so the pattern is laid on the rows as stored.
    for y in range(img.h):
        base = y * W
        for x in range(W):
            g = BAYER_GAINS[pat[(y % 2) * 2 + (x % 2)]]
            d[base + x] = (d[base + x] - m) * g + m * g
    img.bayer = pat


def op_shuffle(img, arg, rng):
    b = int(arg)
    nx, ny = img.w // b, img.h // b
    blocks = [(i, j) for j in range(ny) for i in range(nx)]
    perm = blocks[:]
    rng.shuffle(perm)
    src = img.data
    out = array.array("f", src)
    W = img.w
    for (di, dj), (si, sj) in zip(blocks, perm):
        for r in range(b):
            so = (sj * b + r) * W + si * b
            do = (dj * b + r) * W + di * b
            out[do:do + b] = src[so:so + b]
    img.data = out


# ── generators ──────────────────────────────────────────────────────────────


def gen_image(kind, w, h, ra, dec, fov, rng):
    """A generated frame with a nominal WCS (the hint the benchmark will pass)."""
    scale = fov / max(w, h)
    wcs = Wcs((w + 1) / 2.0, (h + 1) / 2.0, ra, dec, [[-scale, 0.0], [0.0, scale]])
    n = w * h
    sky, noise = 1000.0, 12.0
    d = array.array("f", bytes(4 * n))
    g = rng.gauss
    name, _, arg = kind.partition(":")
    if name == "flat":
        for y in range(h):
            for x in range(w):
                d[y * w + x] = sky + 300.0 * (x / w) + 150.0 * (y / h) + g(0, noise)
    else:
        for i in range(n):
            d[i] = sky + g(0, noise)
    if name == "fakestars":
        nstars = int(arg or 400)
        for _ in range(nstars):
            x0, y0 = rng.uniform(5, w - 6), rng.uniform(5, h - 6)
            flux = 10 ** rng.uniform(2.5, 5.5)
            sig = rng.uniform(1.2, 2.0)
            r = int(4 * sig) + 1
            norm = flux / (2 * math.pi * sig * sig)
            for yy in range(max(0, int(y0) - r), min(h, int(y0) + r + 1)):
                for xx in range(max(0, int(x0) - r), min(w, int(x0) + r + 1)):
                    d[yy * w + xx] += norm * math.exp(-((xx - x0) ** 2 + (yy - y0) ** 2) / (2 * sig * sig))
    elif name not in ("noise", "flat"):
        raise ValueError(f"unknown generator {kind}")
    return Img(w, h, d, wcs, [])


# ── recipe driver ───────────────────────────────────────────────────────────

_SIMPLE = {"crop": op_crop, "bin": op_bin, "flipx": op_flipx, "flipy": op_flipy,
           "transpose": op_transpose, "rot90": op_rot90, "distort": op_distort,
           "vignette": op_vignette, "gradient": op_gradient, "sky": op_sky,
           "glow": op_glow, "blur": op_blur, "trail": op_trail, "saturate": op_saturate,
           "bayer": op_bayer}
_RANDOM = {"clouds": op_clouds, "noise": op_noise, "hot": op_hot, "hotcol": op_hotcol,
           "satellite": op_satellite, "shuffle": op_shuffle}
_OUTPUT = {"u8": 8, "u16": 16, "i32": 32, "f32": -32, "f64": -64}


def apply_ops(img, ops, seed):
    rng = random.Random(seed)
    for op in [o for o in ops.split(",") if o]:
        name, _, arg = op.partition(":")
        if name in _SIMPLE:
            _SIMPLE[name](img, arg)
        elif name in _RANDOM:
            _RANDOM[name](img, arg, rng)
        elif name in _OUTPUT:
            img.bitpix = _OUTPUT[name]
        elif name == "gzip":
            img.fmt = "fits.gz"
        elif name == "xisf":
            img.fmt = "xisf"
        elif name == "keepwcs":
            img.keepwcs = True
        else:
            raise ValueError(f"unknown op {name}")
    return img


def load_parent(path):
    cards, hdr, w, h, data = fitslite.read_image(path)
    wcs = Wcs.from_header(hdr)
    if wcs is None:
        raise ValueError(f"{path}: parent has no usable WCS")
    if wcs.tpv:
        raise ValueError(f"{path}: TPV parents are not supported (use a TAN/SIP parent)")
    data = array.array("f", data)  # writable copy
    # NaNs (blank HiPS pixels) -> median, so arithmetic stays finite
    med = sorted(data[i] for i in range(0, len(data), max(1, len(data) // 100000)) if data[i] == data[i])
    fill = med[len(med) // 2] if med else 0.0
    for i in range(len(data)):
        if data[i] != data[i]:
            data[i] = fill
    return Img(w, h, data, wcs, fitslite.strip_wcs(cards))


def _scale_for_int(img):
    """Map physical values onto a camera-like ADU range for integer output."""
    n = len(img.data)
    step = max(1, n // 400000)
    v = sorted(img.data[i] for i in range(0, n, step))
    lo = v[int(0.001 * len(v))]
    hi = img.saturate if img.saturate is not None else v[-1]
    full = {8: 250.0, 16: 65000.0, 32: 4.0e6}[img.bitpix]
    ped = {8: 5.0, 16: 500.0, 32: 1000.0}[img.bitpix]
    k = (full - ped) / max(hi - lo, 1e-9)
    d = img.data
    for i in range(n):
        x = (d[i] - lo) * k + ped
        d[i] = 0.0 if x < 0 else (full if x > full else x)
    if img.saturate is not None:
        img.saturate = full


def write(img, dest_base, recipe_note=""):
    """Write the image (and a .truth sidecar). Returns the image path."""
    truth_cards = fitslite.wcs_cards(img.wcs)
    if img.bitpix > 0:
        _scale_for_int(img)
    extra = list(img.cards)
    extra.append(fitslite.card("ROWORDER", "BOTTOM-UP"))
    if img.bayer:
        extra.append(fitslite.card("BAYERPAT", img.bayer, "synthetic colour-camera raw"))
        extra.append(fitslite.card("XBAYROFF", 0))
        extra.append(fitslite.card("YBAYROFF", 0))
    if img.saturate is not None:
        extra.append(fitslite.card("SATURATE", float(img.saturate)))
    if recipe_note:
        extra.append(("HISTORY arcsec corpus: " + recipe_note)[:80].ljust(80))
    if img.keepwcs:
        extra += truth_cards
    # sidecar truth, same card syntax as arcsec's .wcs files
    with open(dest_base + ".truth", "w") as f:
        f.write(f"NAXIS1  = {img.w:>20d}\nNAXIS2  = {img.h:>20d}\n")
        for c in truth_cards:
            f.write(c.rstrip() + "\n")
        f.write("END\n")
    if img.fmt == "xisf":
        path = dest_base + ".xisf"
        write_xisf(path, img, truth_cards if img.keepwcs else [])
        return path
    path = dest_base + ".fits"
    bz = None if img.bitpix != 8 else 0.0
    fitslite.write_image(path, img.w, img.h, img.data, bitpix=img.bitpix, extra_cards=extra,
                         bzero=bz)
    if img.fmt == "fits.gz":
        with open(path, "rb") as fi, open(path + ".gz", "wb") as raw, \
                gzip.GzipFile(fileobj=raw, mode="wb", compresslevel=6, mtime=0) as fo:
            fo.write(fi.read())
        os.remove(path)
        path += ".gz"
    return path


def write_xisf(path, img, wcs_cards):
    """Monolithic XISF 1.0, one greyscale image, little-endian samples.

    Rows are stored in the same order as the FITS data (arcsec reads both formats
    without flipping), so the FITS-convention truth in <id>.truth applies as is."""
    samp = {8: ("UInt8", "B"), 16: ("UInt16", "H"), 32: ("UInt32", "I"),
            -32: ("Float32", "f"), -64: ("Float64", "d")}[img.bitpix]
    W, H = img.w, img.h
    if samp[1] in "fd":
        a = array.array(samp[1], img.data)
    else:
        a = array.array(samp[1], (int(round(v)) for v in img.data))
    if a.itemsize > 1 and struct.pack("=H", 1) != struct.pack("<H", 1):
        a.byteswap()
    blob = a.tobytes()
    bounds = ' bounds="0:1"' if samp[1] in "fd" else ""
    kw = []
    for c in wcs_cards:
        k = c[:8].strip()
        v = c[10:].split("/")[0].strip()
        kw.append(f'<FITSKeyword name="{k}" value="{v}" comment=""/>')
    header_tmpl = ('<?xml version="1.0" encoding="UTF-8"?>'
                   '<xisf version="1.0" xmlns="http://www.pixinsight.com/xisf">'
                   '<Image geometry="{W}:{H}:1" sampleFormat="{fmt}" colorSpace="Gray"{bounds}'
                   ' location="attachment:{pos}:{size}">{kw}</Image></xisf>')
    pos = 0
    for _ in range(3):
        hdr = header_tmpl.format(W=W, H=H, fmt=samp[0], bounds=bounds, pos=pos, size=len(blob),
                                 kw="".join(kw)).encode("utf-8")
        need = 16 + len(hdr)
        newpos = (need + 4095) // 4096 * 4096
        if newpos == pos:
            break
        pos = newpos
    with open(path, "wb") as f:
        f.write(b"XISF0100")
        f.write(struct.pack("<I", len(hdr)))
        f.write(b"\0\0\0\0")
        f.write(hdr)
        f.write(b"\0" * (pos - 16 - len(hdr)))
        f.write(blob)
