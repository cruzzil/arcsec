"""Minimal FITS and WCS helpers for the benchmark tooling. Standard library only.

Shared by scripts/benchmark.py (truth WCS, scoring) and scripts/fetch-corpus.py
(validation, TESS cropping, the synthetic tiers). It is deliberately small: enough
FITS to read and write a single 2-D image, and enough WCS to turn a pixel into a sky
position for every projection the corpus contains.

WCS support:
  * linear part: CDi_j, PCi_j + CDELTi, legacy PCiiijjj + CDELTi, or CROTA2 + CDELTi
  * projections: TAN and SIN (WISE L1b frames are SIN-SIP)
  * distortion: SIP (A_p_q / B_p_q, forward direction) and TPV (PVi_k, k <= 39)

Anything else (ZEA, ARC, ...) makes Wcs.from_header return None, so the caller sees
"no truth" rather than a silently wrong one.
"""

import array
import bz2
import gzip
import math
import sys

BLOCK = 2880

# ── low-level header reading ────────────────────────────────────────────────


def _open(path):
    with open(path, "rb") as f:
        magic = f.read(3)
    if magic[:2] == b"\x1f\x8b":
        return gzip.open(path, "rb")
    if magic == b"BZh":
        return bz2.open(path, "rb")
    return open(path, "rb")


def _parse_value(val):
    val = val.strip()
    if val.startswith("'"):
        # quoted string; '' is an escaped quote
        end = 1
        out = []
        while end < len(val):
            c = val[end]
            if c == "'":
                if end + 1 < len(val) and val[end + 1] == "'":
                    out.append("'")
                    end += 2
                    continue
                break
            out.append(c)
            end += 1
        return "".join(out).strip()
    v = val.split("/")[0].strip()
    if v in ("T", "F"):
        return v == "T"
    try:
        if any(c in v for c in ".EeDd"):
            return float(v.replace("D", "E").replace("d", "e"))
        return int(v)
    except ValueError:
        try:
            return float(v)
        except ValueError:
            return v


def parse_card(card):
    """(key, value) for a value card, (key, None) for commentary cards."""
    key = card[:8].strip()
    if card[8:10] != "= ":
        return key, None
    return key, _parse_value(card[10:])


def read_hdus(f, max_hdus=8, max_blocks=400):
    """Yield (cards, header_dict, data_bytes_len) for successive HDUs of an open file,
    skipping over the data. Stops quietly at EOF."""
    for _ in range(max_hdus):
        cards = []
        hdr = {}
        done = False
        for _ in range(max_blocks):
            block = f.read(BLOCK)
            if len(block) < BLOCK:
                return
            for i in range(0, BLOCK, 80):
                card = block[i:i + 80].decode("latin-1")
                key, val = parse_card(card)
                if key == "END":
                    done = True
                    break
                cards.append(card)
                if val is not None and key not in hdr:
                    hdr[key] = val
            if done:
                break
        if not done:
            return
        nbytes = data_size(hdr)
        yield cards, hdr, nbytes
        pad = (nbytes + BLOCK - 1) // BLOCK * BLOCK
        if pad:
            try:
                f.seek(pad, 1)
            except (OSError, ValueError):
                f.read(pad)


def data_size(h):
    naxis = int(h.get("NAXIS", 0))
    if naxis == 0:
        return 0
    n = 1
    for i in range(1, naxis + 1):
        n *= int(h.get(f"NAXIS{i}", 0))
    n = abs(int(h.get("BITPIX", 8))) // 8 * n
    # binary tables carry a heap after the main table
    n += int(h.get("PCOUNT", 0))
    return n


def is_image_hdu(h):
    if h.get("ZIMAGE") is True:
        return True
    return (h.get("XTENSION", "IMAGE") in ("IMAGE", "IMAGE   ") or "SIMPLE" in h) and \
        int(h.get("NAXIS", 0)) >= 2 and int(h.get("NAXIS1", 0)) > 0


def read_header(path, max_blocks=400):
    """Header of the first HDU holding an image, as {keyword: value}.

    Follows CFITSIO's fits_open_image, which is what arcsec uses: a primary HDU with
    NAXIS = 0 is skipped (TESS FFIs, fpacked LCO frames). For a tile-compressed image
    the Z-keywords are mapped back onto NAXISn / BITPIX. Keywords from an empty
    primary header are inherited when the image HDU does not repeat them.
    """
    first = None
    with _open(path) as f:
        for cards, hdr, _ in read_hdus(f, max_blocks=max_blocks):
            if first is None:
                first = hdr
            if is_image_hdu(hdr):
                if hdr.get("ZIMAGE") is True:
                    hdr = dict(hdr)
                    hdr["NAXIS"] = hdr.get("ZNAXIS", 2)
                    for i in (1, 2, 3):
                        if f"ZNAXIS{i}" in hdr:
                            hdr[f"NAXIS{i}"] = hdr[f"ZNAXIS{i}"]
                    hdr["BITPIX"] = hdr.get("ZBITPIX", hdr.get("BITPIX"))
                if first is not hdr:
                    merged = dict(first)
                    merged.update(hdr)
                    hdr = merged
                return hdr
    return first or {}


# ── image data ──────────────────────────────────────────────────────────────

_TYPECODES = {8: "B", 16: "h", 32: "i", -32: "f", -64: "d"}


def _typecode_for(bitpix):
    tc = _TYPECODES[bitpix]
    # 'i' / 'l' sizes vary by platform; pick the 4-byte one
    if bitpix == 32 and array.array("i").itemsize != 4:
        tc = "l"
    return tc


def read_image(path):
    """Read the first uncompressed image HDU.

    Returns (cards, header, width, height, data) where data is an array('f') of
    physical values (BZERO/BSCALE applied), row-major, first row = FITS row 1. A
    3-D cube is averaged over its planes (colour HiPS cutouts)."""
    with _open(path) as f:
        # read_hdus seeks past data; re-implement the walk so we can read it
        while True:
            cards = []
            hdr = {}
            done = False
            while not done:
                block = f.read(BLOCK)
                if len(block) < BLOCK:
                    raise ValueError(f"{path}: no image HDU")
                for i in range(0, BLOCK, 80):
                    card = block[i:i + 80].decode("latin-1")
                    key, val = parse_card(card)
                    if key == "END":
                        done = True
                        break
                    cards.append(card)
                    if val is not None and key not in hdr:
                        hdr[key] = val
            nbytes = data_size(hdr)
            if hdr.get("ZIMAGE") is True:
                raise ValueError(f"{path}: tile-compressed images are not supported here")
            if is_image_hdu(hdr):
                break
            pad = (nbytes + BLOCK - 1) // BLOCK * BLOCK
            f.read(pad)
        bitpix = int(hdr["BITPIX"])
        w, h = int(hdr["NAXIS1"]), int(hdr["NAXIS2"])
        planes = int(hdr.get("NAXIS3", 1)) if int(hdr.get("NAXIS", 2)) >= 3 else 1
        raw = f.read(abs(bitpix) // 8 * w * h * planes)
    a = array.array(_typecode_for(bitpix))
    a.frombytes(raw)
    if sys.byteorder == "little" and a.itemsize > 1:
        a.byteswap()
    bzero = float(hdr.get("BZERO", 0.0))
    bscale = float(hdr.get("BSCALE", 1.0))
    n = w * h
    out = array.array("f", bytes(4 * n))
    if planes == 1:
        if bzero == 0.0 and bscale == 1.0 and a.typecode == "f":
            out = a
        else:
            for i in range(n):
                out[i] = a[i] * bscale + bzero
    else:
        inv = 1.0 / planes
        for i in range(n):
            s = 0.0
            for p in range(planes):
                s += a[p * n + i]
            out[i] = (s * inv) * bscale + bzero
    return cards, hdr, w, h, out


def card(key, value, comment=""):
    """Format one 80-character header card."""
    if isinstance(value, bool):
        v = f"{'T' if value else 'F':>20}"
    elif isinstance(value, int):
        v = f"{value:>20d}"
    elif isinstance(value, float):
        s = f"{value:.15G}"
        if "." not in s and "E" not in s:
            s += "."
        v = f"{s:>20}"
    else:
        s = str(value).replace("'", "''")
        v = f"'{s:<8}'"
    c = f"{key:<8}= {v}"
    if comment:
        c += f" / {comment}"
    return c[:80].ljust(80)


# keywords write_image controls itself
_STRUCTURAL = {"SIMPLE", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "NAXIS3", "EXTEND",
               "BZERO", "BSCALE", "XTENSION", "PCOUNT", "GCOUNT", "END", "CHECKSUM",
               "DATASUM", "BLANK"}


def write_image(path, w, h, data, bitpix=-32, extra_cards=(), bzero=None, bscale=1.0):
    """Write a single-HDU 2-D image. `data` holds physical values, row-major.

    Integer types are stored with the conventional BZERO for unsigned data when
    `bzero` is None (32768 for 16-bit, 2^31 for 32-bit); values are rounded and
    clipped to the representable range."""
    cards = [card("SIMPLE", True, "conforms to FITS standard"),
             card("BITPIX", bitpix), card("NAXIS", 2),
             card("NAXIS1", w), card("NAXIS2", h)]
    if bitpix == 16 and bzero is None:
        bzero = 32768.0
    if bitpix == 32 and bzero is None:
        bzero = 2147483648.0
    if bitpix in (16, 32, 8) and bzero:
        cards.append(card("BZERO", float(bzero), "physical = stored * BSCALE + BZERO"))
        cards.append(card("BSCALE", float(bscale)))
    for c in extra_cards:
        key = c[:8].strip()
        if key in _STRUCTURAL:
            continue
        cards.append(c.ljust(80)[:80])
    cards.append("END".ljust(80))
    head = "".join(cards).encode("latin-1")
    head += b" " * ((-len(head)) % BLOCK)

    n = w * h
    tc = _typecode_for(bitpix)
    if bitpix in (-32, -64):
        out = array.array(tc, data) if not (isinstance(data, array.array) and data.typecode == tc) \
            else array.array(tc, data)
    else:
        lo, hi = {8: (0, 255), 16: (-32768, 32767), 32: (-2147483648, 2147483647)}[bitpix]
        z = float(bzero or 0.0)
        inv = 1.0 / float(bscale)
        out = array.array(tc, bytes(array.array(tc).itemsize * n))
        for i in range(n):
            v = int(round((data[i] - z) * inv))
            out[i] = lo if v < lo else hi if v > hi else v
    if sys.byteorder == "little" and out.itemsize > 1:
        out.byteswap()
    body = out.tobytes()
    body += b"\0" * ((-len(body)) % BLOCK)
    with open(path, "wb") as f:
        f.write(head)
        f.write(body)


# ── WCS ─────────────────────────────────────────────────────────────────────


def _tpv_terms(x, y):
    """The 40 TPV monomials, in PV index order (0..39)."""
    r = math.hypot(x, y)
    x2, y2 = x * x, y * y
    r2 = r * r
    t = [1.0, x, y, r,
         x2, x * y, y2,
         x2 * x, x2 * y, x * y2, y2 * y, r2 * r,
         x2 * x2, x2 * x * y, x2 * y2, x * y2 * y, y2 * y2,
         x2 * x2 * x, x2 * x2 * y, x2 * x * y2, x2 * y2 * y, x * y2 * y2, y2 * y2 * y, r2 * r2 * r,
         ]
    # 6th and 7th order (indices 24..39)
    for n in (6, 7):
        for k in range(n + 1):
            t.append(x ** (n - k) * y ** k)
        if n == 7:
            t.append(r ** 7)
    return t


class Wcs:
    """Celestial WCS for a 2-D image: linear CD matrix, TAN/SIN, optional SIP/TPV."""

    def __init__(self, crpix1, crpix2, crval1, crval2, cd, proj="TAN", sip=None, tpv=None):
        self.crpix1, self.crpix2 = crpix1, crpix2
        self.crval1, self.crval2 = crval1, crval2
        self.cd = cd  # [[cd1_1, cd1_2], [cd2_1, cd2_2]]
        self.proj = proj
        self.sip = sip  # (A dict {(p,q): c}, B dict)
        self.tpv = tpv  # (PV1 list[40], PV2 list[40])

    @property
    def distorted(self):
        return bool(self.sip or self.tpv)

    @classmethod
    def from_header(cls, h):
        try:
            crpix1, crpix2 = float(h["CRPIX1"]), float(h["CRPIX2"])
            crval1, crval2 = float(h["CRVAL1"]), float(h["CRVAL2"])
        except (KeyError, TypeError, ValueError):
            return None

        ctype1 = str(h.get("CTYPE1", "RA---TAN")).upper()
        ctype2 = str(h.get("CTYPE2", "DEC--TAN")).upper()
        if not ctype1.startswith("RA") or not ctype2.startswith("DEC"):
            # e.g. GLON/GLAT - the corpus always asks for equatorial frames
            if ctype1.startswith("DEC") and ctype2.startswith("RA"):
                return None
            if "CTYPE1" in h:
                return None
        code = ctype1[5:8] if len(ctype1) >= 8 else "TAN"
        if code not in ("TAN", "SIN", "TPV"):
            return None
        proj = "SIN" if code == "SIN" else "TAN"

        if "CD1_1" in h or "CD2_2" in h:
            cd = [[h.get("CD1_1", 0.0), h.get("CD1_2", 0.0)],
                  [h.get("CD2_1", 0.0), h.get("CD2_2", 0.0)]]
        else:
            cdelt1 = h.get("CDELT1")
            cdelt2 = h.get("CDELT2")
            if cdelt1 is None or cdelt2 is None:
                return None
            if "PC1_1" in h or "PC2_2" in h:
                pc = [[h.get("PC1_1", 1.0), h.get("PC1_2", 0.0)],
                      [h.get("PC2_1", 0.0), h.get("PC2_2", 1.0)]]
            elif "PC001001" in h:
                # Legacy PCiiijjj form, still used by Pan-STARRS skycell headers.
                # Missing it silently flips the RA axis sign and, with the large
                # off-image CRPIX those cutouts carry, throws the corners out by
                # over a degree.
                pc = [[h.get("PC001001", 1.0), h.get("PC001002", 0.0)],
                      [h.get("PC002001", 0.0), h.get("PC002002", 1.0)]]
            else:
                rot = math.radians(h.get("CROTA2", 0.0))
                pc = [[math.cos(rot), -math.sin(rot)],
                      [math.sin(rot), math.cos(rot)]]
            cd = [[cdelt1 * pc[0][0], cdelt1 * pc[0][1]],
                  [cdelt2 * pc[1][0], cdelt2 * pc[1][1]]]
        cd = [[float(v) for v in row] for row in cd]
        if abs(cd[0][0] * cd[1][1] - cd[0][1] * cd[1][0]) == 0.0:
            return None

        sip = None
        if ctype1.endswith("-SIP") or "A_ORDER" in h:
            a, b = {}, {}
            for p in range(0, 10):
                for q in range(0, 10 - p):
                    if f"A_{p}_{q}" in h:
                        a[(p, q)] = float(h[f"A_{p}_{q}"])
                    if f"B_{p}_{q}" in h:
                        b[(p, q)] = float(h[f"B_{p}_{q}"])
            if a or b:
                sip = (a, b)

        tpv = None
        if code == "TPV":
            pv1 = [0.0] * 40
            pv2 = [0.0] * 40
            pv1[1] = pv2[1] = 1.0
            for k in range(40):
                if f"PV1_{k}" in h:
                    pv1[k] = float(h[f"PV1_{k}"])
                if f"PV2_{k}" in h:
                    pv2[k] = float(h[f"PV2_{k}"])
            tpv = (pv1, pv2)
        return cls(crpix1, crpix2, crval1, crval2, cd, proj=proj, sip=sip, tpv=tpv)

    def intermediate(self, x, y):
        """1-based pixel -> intermediate world coordinates (xi, eta) in degrees."""
        u = x - self.crpix1
        v = y - self.crpix2
        if self.sip:
            a, b = self.sip
            du = dv = 0.0
            for (p, q), c in a.items():
                du += c * u ** p * v ** q
            for (p, q), c in b.items():
                dv += c * u ** p * v ** q
            u += du
            v += dv
        xi = self.cd[0][0] * u + self.cd[0][1] * v
        eta = self.cd[1][0] * u + self.cd[1][1] * v
        if self.tpv:
            pv1, pv2 = self.tpv
            t1 = _tpv_terms(xi, eta)
            t2 = _tpv_terms(eta, xi)
            xi, eta = (sum(c * t for c, t in zip(pv1, t1) if c),
                       sum(c * t for c, t in zip(pv2, t2) if c))
        return xi, eta

    def pix2sky(self, x, y):
        """1-based FITS pixel coords -> (ra_deg, dec_deg)."""
        xi_d, eta_d = self.intermediate(x, y)
        xi = math.radians(xi_d)
        eta = math.radians(eta_d)
        ra0 = math.radians(self.crval1)
        dec0 = math.radians(self.crval2)
        sd0, cd0 = math.sin(dec0), math.cos(dec0)
        if self.proj == "SIN":
            z = math.sqrt(max(0.0, 1.0 - xi * xi - eta * eta))
        else:
            z = 1.0
        denom = z * cd0 - eta * sd0
        ra = ra0 + math.atan2(xi, denom)
        dec = math.atan2(z * sd0 + eta * cd0, math.hypot(denom, xi))
        return math.degrees(ra) % 360.0, math.degrees(dec)

    def pixscale(self):
        """Degrees per pixel along the first image axis (linear part)."""
        return math.hypot(self.cd[0][0], self.cd[1][0])

    def rotation(self):
        return math.degrees(math.atan2(self.cd[1][0], self.cd[1][1]))

    def parity(self):
        """Sign of the CD determinant: negative for the usual sky orientation."""
        return 1 if (self.cd[0][0] * self.cd[1][1] - self.cd[0][1] * self.cd[1][0]) > 0 else -1


def angsep(ra1, dec1, ra2, dec2):
    """Angular separation in arcseconds."""
    r1, d1 = math.radians(ra1), math.radians(dec1)
    r2, d2 = math.radians(ra2), math.radians(dec2)
    c = math.sin(d1) * math.sin(d2) + math.cos(d1) * math.cos(d2) * math.cos(r1 - r2)
    return math.degrees(math.acos(max(-1.0, min(1.0, c)))) * 3600.0


def _project(ra, dec, ra0, dec0):
    """Gnomonic projection about (ra0, dec0); all degrees in, degrees out."""
    r, d = math.radians(ra), math.radians(dec)
    r0, d0 = math.radians(ra0), math.radians(dec0)
    cosc = math.sin(d0) * math.sin(d) + math.cos(d0) * math.cos(d) * math.cos(r - r0)
    xi = math.cos(d) * math.sin(r - r0) / cosc
    eta = (math.cos(d0) * math.sin(d) - math.sin(d0) * math.cos(d) * math.cos(r - r0)) / cosc
    return math.degrees(xi), math.degrees(eta)


def _solve3(m, v):
    """Solve a 3x3 linear system by Cramer's rule."""
    def det(a):
        return (a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
                - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
                + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]))
    d = det(m)
    out = []
    for k in range(3):
        mk = [row[:] for row in m]
        for i in range(3):
            mk[i][k] = v[i]
        out.append(det(mk) / d)
    return out


def best_linear(truth, naxis1, naxis2, n=9):
    """Best linear TAN approximation of `truth` over the image.

    A plate solver that fits only a linear plate model cannot do better than this.
    Returns (wcs, floor) where floor is the worst corner error, in arcsec, of the
    least-squares linear fit (tangent point at the image centre) against the truth
    sampled on an n x n grid. For an undistorted TAN truth the floor is ~0."""
    cx, cy = (naxis1 + 1) / 2.0, (naxis2 + 1) / 2.0
    ra0, dec0 = truth.pix2sky(cx, cy)
    rows = []
    for i in range(n):
        for j in range(n):
            x = 1.0 + (naxis1 - 1.0) * i / (n - 1)
            y = 1.0 + (naxis2 - 1.0) * j / (n - 1)
            ra, dec = truth.pix2sky(x, y)
            xi, eta = _project(ra, dec, ra0, dec0)
            rows.append((x - cx, y - cy, xi, eta))
    ata = [[0.0] * 3 for _ in range(3)]
    atx = [0.0] * 3
    aty = [0.0] * 3
    for u, v, xi, eta in rows:
        f = (1.0, u, v)
        for a in range(3):
            atx[a] += f[a] * xi
            aty[a] += f[a] * eta
            for b in range(3):
                ata[a][b] += f[a] * f[b]
    px = _solve3(ata, atx)
    py = _solve3(ata, aty)
    # fold the constant terms into CRPIX by keeping CRVAL at the centre: a small
    # offset in the tangent plane is absorbed by solving for the pixel it maps to
    cd = [[px[1], px[2]], [py[1], py[2]]]
    det = cd[0][0] * cd[1][1] - cd[0][1] * cd[1][0]
    du = (-px[0] * cd[1][1] + py[0] * cd[0][1]) / det
    dv = (-py[0] * cd[0][0] + px[0] * cd[1][0]) / det
    lin = Wcs(cx + du, cy + dv, ra0, dec0, cd)
    corners = [(1.0, 1.0), (naxis1, 1.0), (1.0, naxis2), (naxis1, naxis2)]
    floor = 0.0
    for (x, y) in corners:
        a = truth.pix2sky(x, y)
        b = lin.pix2sky(x, y)
        floor = max(floor, angsep(a[0], a[1], b[0], b[1]))
    return lin, floor


def wcs_cards(w, extra_sip=None):
    """Header cards for a TAN (optionally TAN-SIP) WCS."""
    sip = extra_sip if extra_sip is not None else w.sip
    suffix = "-SIP" if sip else ""
    proj = w.proj
    cards = [card("CTYPE1", f"RA---{proj}{suffix}"), card("CTYPE2", f"DEC--{proj}{suffix}"),
             card("CRPIX1", float(w.crpix1)), card("CRPIX2", float(w.crpix2)),
             card("CRVAL1", float(w.crval1)), card("CRVAL2", float(w.crval2)),
             card("CD1_1", float(w.cd[0][0])), card("CD1_2", float(w.cd[0][1])),
             card("CD2_1", float(w.cd[1][0])), card("CD2_2", float(w.cd[1][1])),
             card("RADESYS", "ICRS")]
    if sip:
        a, b = sip
        order = max([p + q for (p, q) in list(a) + list(b)] + [2])
        cards.append(card("A_ORDER", order))
        cards.append(card("B_ORDER", order))
        for (p, q), c in sorted(a.items()):
            cards.append(card(f"A_{p}_{q}", float(c)))
        for (p, q), c in sorted(b.items()):
            cards.append(card(f"B_{p}_{q}", float(c)))
    return cards


WCS_KEY_PREFIXES = ("CTYPE", "CRPIX", "CRVAL", "CDELT", "CROTA", "CD1_", "CD2_", "PC1_", "PC2_",
                    "PC00", "PV1_", "PV2_", "A_", "B_", "AP_", "BP_", "CUNIT", "WCSAXES",
                    "RADESYS", "RADECSYS", "EQUINOX", "LONPOLE", "LATPOLE", "WCSNAME",
                    "CRDER", "CSYER", "IMAGEW", "IMAGEH")


def strip_wcs(cards):
    """Drop every WCS-related card (for writing a new WCS in its place)."""
    out = []
    for c in cards:
        k = c[:8].strip()
        if any(k.startswith(p) for p in WCS_KEY_PREFIXES):
            continue
        out.append(c)
    return out
