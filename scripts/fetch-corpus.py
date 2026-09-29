#!/usr/bin/env python3
"""Materialise the expanded benchmark corpus (scripts/corpus.tsv) into resources/corpus.

    scripts/fetch-corpus.py                     # everything (resumable; re-run to retry)
    scripts/fetch-corpus.py --list              # show the manifest, download nothing
    scripts/fetch-corpus.py --stats             # coverage tables for the manifest
    scripts/fetch-corpus.py --tier B --source ztf --set core --id foo   (repeatable)
    scripts/fetch-corpus.py --force             # re-fetch files that already exist
    scripts/fetch-corpus.py --verify            # re-hash files against SHA256SUMS
    scripts/fetch-corpus.py --jobs 4            # parallel downloads (max 4)
    scripts/fetch-corpus.py --reuse-dir resources/testset   # hard-link v1 images already
                                                            # fetched by fetch-test-images.sh

Being a good citizen to the archives: at most --jobs (<= 4) requests in flight, at most
2 per host, a minimum interval between requests to one host, a descriptive User-Agent,
retries with exponential backoff (honouring Retry-After), and nothing re-fetched that is
already present and valid. Downloads land in `<file>.part` and are renamed only once
validated, so an interrupted run resumes cleanly.

Every file is validated before it is kept: FITS signature and size, a WCS the benchmark
can read, a centre within the requested field, and the checksum the archive publishes
when it publishes one (LCO md5, AWS ETag for TESS). SHA256 of every kept file is
recorded in <out>/SHA256SUMS; --verify re-checks it. Cutout services stamp dates into
headers, so a re-fetch is not byte-identical and the sums are local integrity checks,
not global ones.

Derived entries (source `synth`, tier S and generated tier-D controls) are built after
their parents, deterministically from the parent and the recipe's seed; see
scripts/corpus_synth.py.

Standard library only.
"""

import argparse
import bz2
import concurrent.futures
import gzip
import hashlib
import json
import os
import random
import re
import shutil
import sys
import threading
import time
import urllib.error
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fitslite  # noqa: E402

UA = "arcsec-corpus-fetch/1.0 (+https://github.com/cruzzil/arcsec; benchmark corpus, low rate)"
MAX_JOBS = 4
PER_HOST = 2
MIN_INTERVAL = {"default": 0.5, "archive-api.lco.global": 1.0, "skyview.gsfc.nasa.gov": 1.0,
                "api.skymapper.nci.org.au": 1.0, "ps1images.stsci.edu": 1.0}
IMAGE_EXTS = (".fits", ".fits.fz", ".fits.gz", ".xisf")

# ── manifest ────────────────────────────────────────────────────────────────


def parse_opts(extra):
    out = {}
    if "=" not in extra:
        return out
    for part in extra.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def load_manifest(path):
    out = []
    with open(path) as f:
        for ln, line in enumerate(f, 1):
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            p = re.split(r"[ \t]+", s)
            if len(p) < 9:
                raise SystemExit(f"{path}:{ln}: expected at least 9 columns")
            e = {"id": p[0], "tier": p[1], "ra": float(p[2]), "dec": float(p[3]),
                 "fov": float(p[4]), "w": int(p[5]), "h": int(p[6]), "source": p[7],
                 "extra": p[8], "dataset": p[9] if len(p) > 9 else "",
                 "sets": p[10].split(",") if len(p) > 10 else [],
                 "truth": p[11] if len(p) > 11 else ""}
            e["opts"] = parse_opts(e["extra"])
            out.append(e)
    ids = [e["id"] for e in out]
    dup = {i for i in ids if ids.count(i) > 1}
    if dup:
        raise SystemExit(f"duplicate ids in manifest: {sorted(dup)}")
    return out


# ── polite HTTP ─────────────────────────────────────────────────────────────


class Http:
    def __init__(self):
        self.lock = threading.Lock()
        self.host_sem = {}
        self.host_last = {}

    def _sem(self, host):
        with self.lock:
            if host not in self.host_sem:
                self.host_sem[host] = threading.Semaphore(PER_HOST)
            return self.host_sem[host]

    def _pace(self, host):
        gap = MIN_INTERVAL.get(host, MIN_INTERVAL["default"])
        while True:
            with self.lock:
                now = time.time()
                last = self.host_last.get(host, 0.0)
                if now - last >= gap:
                    self.host_last[host] = now
                    return
                wait = gap - (now - last)
            time.sleep(wait)

    def get(self, url, dest=None, tries=5, timeout=300, accept_codes=(200, 206)):
        """GET url. With dest, stream to dest (resuming a .part via Range when the
        server allows it) and return the path; otherwise return the body bytes."""
        host = urllib.parse.urlparse(url).hostname or "?"
        delay = 5.0
        last_err = None
        for attempt in range(tries):
            with self._sem(host):
                self._pace(host)
                try:
                    headers = {"User-Agent": UA}
                    have = 0
                    if dest and os.path.exists(dest):
                        have = os.path.getsize(dest)
                        if have:
                            headers["Range"] = f"bytes={have}-"
                    req = urllib.request.Request(url, headers=headers)
                    with urllib.request.urlopen(req, timeout=timeout) as r:
                        if r.status not in accept_codes:
                            raise urllib.error.HTTPError(url, r.status, "unexpected", r.headers, None)
                        if dest is None:
                            return r.read()
                        mode = "ab" if (have and r.status == 206) else "wb"
                        with open(dest, mode) as f:
                            shutil.copyfileobj(r, f, 1 << 20)
                        return dest
                except urllib.error.HTTPError as e:
                    last_err = e
                    if e.code == 416 and dest:  # range not satisfiable: start over
                        os.remove(dest)
                        continue
                    if e.code not in (408, 425, 429, 500, 502, 503, 504):
                        raise
                    ra = e.headers.get("Retry-After") if e.headers else None
                    if ra and ra.isdigit():
                        delay = max(delay, float(ra))
                except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as e:
                    last_err = e
            if attempt + 1 < tries:
                time.sleep(delay * random.uniform(0.8, 1.2))
                delay = min(delay * 3, 180.0)
        raise RuntimeError(f"giving up on {url}: {last_err}")


HTTP = Http()

# ── validation helpers ──────────────────────────────────────────────────────


def fits_ok(path):
    with open(path, "rb") as f:
        head = f.read(9)
    if head != b"SIMPLE  =":
        return False
    return os.path.getsize(path) % 2880 == 0


def check_wcs(path, e, tol_fields=0.75, need_centre=True):
    """Header WCS readable, and the image centre where the manifest says it is."""
    h = fitslite.read_header(path)
    w = fitslite.Wcs.from_header(h)
    if w is None:
        raise ValueError("no usable WCS in the delivered header")
    n1, n2 = int(h.get("NAXIS1", 0)), int(h.get("NAXIS2", 0))
    if n1 < 16 or n2 < 16:
        raise ValueError(f"image too small ({n1}x{n2})")
    if need_centre:
        ra, dec = w.pix2sky((n1 + 1) / 2.0, (n2 + 1) / 2.0)
        sep = fitslite.angsep(ra, dec, e["ra"], e["dec"]) / 3600.0
        fov = w.pixscale() * max(n1, n2)
        if sep > max(tol_fields * fov, 0.02):
            raise ValueError(f"centre {sep:.3f} deg from the requested position (fov {fov:.3f})")
    return h, w


def check_not_blank(path, max_blank=0.25):
    """Reject cutouts that fall (mostly) outside a survey's footprint."""
    _, _, w, h, data = fitslite.read_image(path)
    n = w * h
    step = max(1, n // 50000)
    tot = blank = 0
    for i in range(0, n, step):
        v = data[i]
        tot += 1
        if v != v or v == 0.0:
            blank += 1
    if tot and blank / tot > max_blank:
        raise ValueError(f"{100.0 * blank / tot:.0f}% blank pixels (outside the survey footprint?)")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def decompress(src, dest, kind):
    op = gzip.open if kind == "gz" else bz2.open
    with op(src, "rb") as fi, open(dest, "wb") as fo:
        shutil.copyfileobj(fi, fo, 1 << 20)
    os.remove(src)


def is_gzip(path):
    with open(path, "rb") as f:
        return f.read(2) == b"\x1f\x8b"


# ── fetchers: each writes <base>.<ext>.part and returns the final path ──────


def q(s):
    return urllib.parse.quote(str(s), safe="")


def f_hips2fits(e, base):
    hips = e["opts"].get("hips", e["extra"] if "=" not in e["extra"] else "")
    url = ("https://alasky.cds.unistra.fr/hips-image-services/hips2fits?"
           f"hips={q(hips)}&width={e['w']}&height={e['h']}&fov={e['fov']}"
           f"&projection=TAN&coordsys=icrs&ra={e['ra']}&dec={e['dec']}&format=fits")
    if "rot" in e["opts"]:
        url += f"&rotation_angle={e['opts']['rot']}"
    return simple(url, e, base)


def f_skyview(e, base):
    survey = urllib.parse.unquote(e["opts"].get("survey", e["extra"] if "=" not in e["extra"] else ""))
    pix = f"{e['w']},{e['h']}" if e["w"] != e["h"] else f"{e['w']}"
    size = e["fov"] if e["w"] == e["h"] else f"{e['fov']},{e['fov'] * e['h'] / e['w']:.6f}"
    url = ("https://skyview.gsfc.nasa.gov/current/cgi/pskcall?"
           f"Survey={q(survey)}&Position={e['ra']},{e['dec']}&Size={size}&Pixels={pix}"
           "&Projection=Tan&Coordinates=J2000&Return=FITS")
    if "rot" in e["opts"]:
        url += f"&Rotation={e['opts']['rot']}"
    return simple(url, e, base)


def f_legacysurvey(e, base):
    o = e["opts"]
    band = o.get("band", e["extra"] if "=" not in e["extra"] else "r")
    layer = o.get("layer", "ls-dr10")
    ps = e["fov"] * 3600.0 / e["w"]
    url = ("https://www.legacysurvey.org/viewer/cutout.fits?"
           f"ra={e['ra']}&dec={e['dec']}&layer={layer}&pixscale={ps:.4f}&size={e['w']}&bands={band}")
    return simple(url, e, base)


def f_panstarrs(e, base):
    o = e["opts"]
    band = o.get("band", e["extra"] if "=" not in e["extra"] else "r")
    size = int(o.get("size", e["w"]))
    lst = HTTP.get("https://ps1images.stsci.edu/cgi-bin/ps1filenames.py?"
                   f"ra={e['ra']}&dec={e['dec']}&size={size}&format=fits&filters={band}").decode()
    lines = [ln.split() for ln in lst.splitlines()[1:] if ln.strip()]
    if not lines:
        raise ValueError("no PS1 skycell")
    fname = lines[0][7]
    url = ("https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
           f"red={fname}&format=fits&size={size}&ra={e['ra']}&dec={e['dec']}")
    if "out" in o:
        url += f"&output_size={o['out']}"
    return simple(url, e, base)


def f_sdss(e, base):
    path = e["opts"].get("frame", e["extra"])
    rerun, run, camcol, band, field = path.split("/")
    url = (f"https://data.sdss.org/sas/dr17/eboss/photoObj/frames/{rerun}/{run}/{camcol}/"
           f"frame-{band}-{int(run):06d}-{camcol}-{int(field):04d}.fits.bz2")
    part = base + ".fits.bz2.part"
    HTTP.get(url, part)
    decompress(part, base + ".fits.part", "bz2")
    return finish_fits(base, e, need_centre=False)


def f_ztf(e, base):
    o = e["opts"]
    fn = o["file"]
    m = re.match(r"ztf_(\d{4})(\d{4})(\d{6})_(\d{6})_(z[gri])_c(\d\d)_([a-z])_q(\d)_sciimg\.fits", fn)
    if not m:
        raise ValueError(f"bad ZTF file name {fn}")
    y, md, frac = m.group(1), m.group(2), m.group(3)
    url = f"https://irsa.ipac.caltech.edu/ibe/data/ztf/products/sci/{y}/{md}/{frac}/{fn}"
    if "size" in o:
        url += f"?center={e['ra']},{e['dec']}&size={o['size']}pix"
    return simple(url, e, base)


def f_wise(e, base):
    o = e["opts"]
    fr = o["frame"]  # e.g. 05150a195
    scan, num = fr[:6], fr[6:]
    band = o.get("band", "1")
    url = (f"https://irsa.ipac.caltech.edu/ibe/data/wise/allsky/4band_p1bm_frm/{scan[-2:]}/"
           f"{scan}/{num}/{scan}{num}-w{band}-int-1b.fits")
    return simple(url, e, base, tol=1.0)


def f_skymapper(e, base):
    o = e["opts"]
    size = o.get("size", str(e["fov"]))
    url = ("https://api.skymapper.nci.org.au/public/siap/dr4/get_image?"
           f"image={o['image']}&format=fits&pos={e['ra']},{e['dec']}&size={size},{size}")
    return simple(url, e, base)


TESS_CACHE_LOCK = threading.Lock()
TESS_LOCKS = {}


def f_tess(e, base, cache_dir):
    """Crop a calibrated TESS FFI (SIP WCS) to a single-HDU float32 image."""
    o = e["opts"]
    key = o["key"]
    size = int(o.get("size", e["w"]))
    x0, y0 = int(o.get("x0", 44)), int(o.get("y0", 0))
    fn = os.path.join(cache_dir, "tess", key.replace("/", "_"))
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    with TESS_CACHE_LOCK:
        lk = TESS_LOCKS.setdefault(fn, threading.Lock())
    with lk:
        if not os.path.exists(fn):
            url = f"https://stpubdata.s3.amazonaws.com/tess/public/ffi/{key}"
            HTTP.get(url, fn + ".part")
            etag = o.get("md5")
            if etag and md5(fn + ".part") != etag:
                os.remove(fn + ".part")
                raise ValueError("TESS FFI md5 mismatch")
            os.replace(fn + ".part", fn)
    cards, hdr, w, h, data = fitslite.read_image(fn)
    if fitslite.Wcs.from_header(hdr) is None:
        raise ValueError("TESS FFI carries no celestial WCS")
    if x0 + size > w or y0 + size > h:
        raise ValueError("crop outside the FFI")
    import array
    out = array.array("f", bytes(4 * size * size))
    for r in range(size):
        s = (y0 + r) * w + x0
        out[r * size:(r + 1) * size] = data[s:s + size]
    keep = []
    for c in cards:
        k = c[:8].strip()
        if k in ("CRPIX1", "CRPIX2"):
            continue
        keep.append(c)
    keep.append(fitslite.card("CRPIX1", float(hdr["CRPIX1"]) - x0, "shifted for the corpus crop"))
    keep.append(fitslite.card("CRPIX2", float(hdr["CRPIX2"]) - y0, "shifted for the corpus crop"))
    keep.append(("HISTORY arcsec corpus: crop of " + key.split("/")[-1])[:80].ljust(80))
    keep.append((f"HISTORY   at x0={x0} y0={y0} size={size} (0-based FFI pixels)")[:80].ljust(80))
    part = base + ".fits.part"
    fitslite.write_image(part, size, size, out, bitpix=-32, extra_cards=keep)
    return finish_fits(base, e, need_centre=False)


def f_lco(e, base):
    o = e["opts"]
    meta = json.loads(HTTP.get(f"https://archive-api.lco.global/frames/{o['frame']}/"))
    if not meta.get("url"):
        raise ValueError("LCO frame has no download url (not public?)")
    part = base + ".fits.fz.part"
    HTTP.get(meta["url"], part)
    want = o.get("md5") or (meta.get("version_set") or [{}])[0].get("md5")
    if want and md5(part) != want:
        os.remove(part)
        raise ValueError("LCO md5 mismatch")
    h = fitslite.read_header(part)
    if int(h.get("WCSERR", 0) or 0) != 0:
        os.remove(part)
        raise ValueError("LCO pipeline WCS flagged bad (WCSERR != 0)")
    check_wcs(part, e)
    final = base + ".fits.fz"
    os.replace(part, final)
    return final


def simple(url, e, base, tol=0.75):
    part = base + ".fits.part"
    if os.path.exists(part):
        os.remove(part)  # cutout services do not support Range
    HTTP.get(url, part)
    if is_gzip(part):
        decompress(part, part + ".x", "gz")
        os.replace(part + ".x", part)
    return finish_fits(base, e, tol=tol)


def finish_fits(base, e, need_centre=True, tol=0.75):
    part = base + ".fits.part"
    if not fits_ok(part):
        with open(part, "rb") as f:
            snippet = f.read(200).decode("latin-1", "replace").replace("\n", " ")
        os.remove(part)
        raise ValueError(f"not a FITS file: {snippet[:120]!r}")
    try:
        check_wcs(part, e, tol_fields=tol, need_centre=need_centre)
        check_not_blank(part)
    except ValueError:
        os.remove(part)
        raise
    final = base + ".fits"
    os.replace(part, final)
    return final


def f_synth(e, base, out_dir):
    import corpus_synth as cs
    o = e["opts"]
    seed = int(o.get("seed", 1))
    ops = o.get("ops", "")
    if "gen" in o:
        img = cs.gen_image(o["gen"], e["w"], e["h"], e["ra"], e["dec"], e["fov"], random.Random(seed))
    else:
        ppath = existing(out_dir, o["parent"])
        if not ppath:
            raise ValueError(f"parent {o['parent']} not fetched")
        img = cs.load_parent(ppath)
    cs.apply_ops(img, ops, seed)
    note = (f"gen={o['gen']} " if "gen" in o else f"parent={o.get('parent')} ") + f"ops={ops} seed={seed}"
    return cs.write(img, base, note)


def existing(out_dir, iid):
    for ext in IMAGE_EXTS:
        p = os.path.join(out_dir, iid + ext)
        if os.path.exists(p):
            return p
    return None


# ── main ────────────────────────────────────────────────────────────────────


def fetch_one(e, args):
    out = args.out
    base = os.path.join(out, e["id"])
    if e["source"] == "alias":
        return "alias", None
    have = existing(out, e["id"])
    if have and not args.force:
        return "skip", have
    if have and args.force:
        os.remove(have)
    # v1 entries already fetched by fetch-test-images.sh: hard-link, do not re-download
    if args.reuse_dir and "v1" in e["sets"] and not args.force:
        src = os.path.join(args.reuse_dir, e["id"] + ".fits")
        if os.path.exists(src):
            try:
                os.link(src, base + ".fits")
            except OSError:
                shutil.copy2(src, base + ".fits")
            return "linked", base + ".fits"
    s = e["source"]
    if s == "hips2fits":
        p = f_hips2fits(e, base)
    elif s == "skyview":
        p = f_skyview(e, base)
    elif s == "legacysurvey":
        p = f_legacysurvey(e, base)
    elif s == "panstarrs":
        p = f_panstarrs(e, base)
    elif s == "sdss":
        p = f_sdss(e, base)
    elif s == "ztf":
        p = f_ztf(e, base)
    elif s == "wise":
        p = f_wise(e, base)
    elif s == "skymapper":
        p = f_skymapper(e, base)
    elif s == "tess":
        p = f_tess(e, base, os.path.join(out, ".cache"))
    elif s == "lco":
        p = f_lco(e, base)
    elif s == "synth":
        p = f_synth(e, base, out)
    else:
        raise ValueError(f"unknown source {s}")
    return "ok", p


def select(entries, args):
    if args.id:
        return [e for e in entries if e["id"] in args.id]
    out = entries
    if args.tier:
        out = [e for e in out if e["tier"] in args.tier]
    if args.source:
        out = [e for e in out if e["source"] in args.source]
    if args.set:
        out = [e for e in out if set(e["sets"]) & set(args.set)]
    if args.dataset:
        out = [e for e in out if any(d.lower() in e["dataset"].lower() for d in args.dataset)]
    return out


def with_parents(sel, entries):
    """Add the parents of selected derived entries so they can be built."""
    byid = {e["id"]: e for e in entries}
    want = {e["id"] for e in sel}
    todo = list(sel)
    while todo:
        e = todo.pop()
        for key in ("parent", "file"):
            p = e["opts"].get(key)
            if p and p in byid and p not in want:
                want.add(p)
                todo.append(byid[p])
    return [e for e in entries if e["id"] in want]


def print_stats(entries):
    from collections import Counter

    def table(title, key):
        c = Counter(key(e) for e in entries)
        print(f"\n{title}")
        for k, v in sorted(c.items(), key=lambda kv: (-kv[1], kv[0])):
            print(f"  {k:<28}{v:>5}")

    def fovb(e):
        f = e["fov"] * min(e["w"], e["h"]) / max(e["w"], e["h"], 1) if e["w"] else e["fov"]
        for lo, hi, n in ((0, .15, "a <0.15"), (.15, .3, "b 0.15-0.3"), (.3, .6, "c 0.3-0.6"),
                          (.6, 1.2, "d 0.6-1.2"), (1.2, 2.5, "e 1.2-2.5"), (2.5, 6, "f 2.5-6"),
                          (6, 20, "g 6-20"), (20, 999, "h >20")):
            if lo <= f < hi:
                return n
        return "?"

    def scaleb(e):
        s = e["fov"] * 3600 / max(e["w"], e["h"], 1)
        for lo, hi, n in ((0, .5, "a <0.5\""), (.5, 1, "b 0.5-1\""), (1, 2, "c 1-2\""),
                          (2, 4, "d 2-4\""), (4, 10, "e 4-10\""), (10, 30, "f 10-30\""),
                          (30, 1e9, "g >30\"")):
            if lo <= s < hi:
                return n
        return "?"

    def decb(e):
        d = e["dec"]
        return f"{int((d + 90) // 30) * 30 - 90:+d}..{int((d + 90) // 30) * 30 - 60:+d}"

    print(f"{len(entries)} entries")
    table("tier", lambda e: e["tier"])
    table("source", lambda e: e["source"])
    table("dataset", lambda e: e["dataset"])
    table("set", lambda e: ",".join(e["sets"]))
    table("truth", lambda e: e["truth"])
    table("FOV (short side, deg; manifest nominal)", fovb)
    table("pixel scale (manifest nominal)", scaleb)
    table("declination band", decb)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", default=os.path.join(HERE, "corpus.tsv"))
    ap.add_argument("--out", default=os.path.join(HERE, "..", "resources", "corpus"))
    ap.add_argument("--reuse-dir", default=os.path.join(HERE, "..", "resources", "testset"),
                    help="hard-link v1 images from here when present ('' to disable)")
    ap.add_argument("--tier", action="append", default=[])
    ap.add_argument("--source", action="append", default=[])
    ap.add_argument("--set", action="append", default=[])
    ap.add_argument("--dataset", action="append", default=[])
    ap.add_argument("--id", action="append", default=[])
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--stats", action="store_true")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--keep-cache", action="store_true",
                    help="keep downloaded TESS FFIs in <out>/.cache after cropping")
    ap.add_argument("--jobs", type=int, default=4)
    args = ap.parse_args()
    args.jobs = max(1, min(args.jobs, MAX_JOBS))
    args.out = os.path.abspath(args.out)
    if args.reuse_dir:
        args.reuse_dir = os.path.abspath(args.reuse_dir)

    entries = load_manifest(args.manifest)
    sel = select(entries, args)
    if args.stats:
        print_stats(sel)
        return 0
    if args.list:
        for e in sel:
            print(f"{e['id']:<24} {e['tier']}  ra={e['ra']:<10.4f} dec={e['dec']:<9.4f} "
                  f"fov={e['fov']:<7g} {e['w']}x{e['h']:<6} {e['source']:<12} {e['dataset']:<18} "
                  f"{','.join(e['sets'])}")
        return 0

    os.makedirs(args.out, exist_ok=True)
    sums_path = os.path.join(args.out, "SHA256SUMS")
    sums = {}
    if os.path.exists(sums_path):
        for line in open(sums_path):
            parts = line.split()
            if len(parts) == 2:
                sums[parts[1]] = parts[0]

    if args.verify:
        bad = 0
        for e in sel:
            p = existing(args.out, e["id"])
            if not p:
                continue
            name = os.path.basename(p)
            if name not in sums:
                print(f"  no-sum  {name}")
                continue
            if sha256(p) != sums[name]:
                print(f"  CHANGED {name}")
                bad += 1
        print(f"verify: {bad} changed")
        return 1 if bad else 0

    todo = with_parents(sel, entries)
    base_entries = [e for e in todo if e["source"] != "synth"]
    # interleave sources so the per-host limit on one archive does not idle the pool
    by_src = {}
    for e in base_entries:
        by_src.setdefault(e["source"], []).append(e)
    base_entries = []
    while any(by_src.values()):
        for src in list(by_src):
            if by_src[src]:
                base_entries.append(by_src[src].pop(0))
    derived = [e for e in todo if e["source"] == "synth"]
    t0 = time.time()
    lock = threading.Lock()
    results = {"ok": 0, "skip": 0, "linked": 0, "alias": 0, "fail": 0}
    failed = []

    def run(e):
        try:
            st, p = fetch_one(e, args)
        except Exception as ex:  # noqa: BLE001 - report and carry on
            with lock:
                results["fail"] += 1
                failed.append((e["id"], str(ex)[:160]))
                print(f"  FAIL    {e['id']:<24} {e['source']:<12} {ex}", flush=True)
            return
        with lock:
            results[st] += 1
            if p and (st in ("ok", "linked") or os.path.basename(p) not in sums):
                sums[os.path.basename(p)] = sha256(p)
            if st in ("ok", "linked"):
                sz = os.path.getsize(p) / 1e6 if p else 0
                print(f"  {st:<7} {e['id']:<24} {e['source']:<12} {sz:8.1f} MB", flush=True)

    # TESS crops sharing an FFI wait on a per-file lock, so the FFI is fetched once.
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as ex:
        list(ex.map(run, base_entries))
    # derived images are CPU-bound pure Python: use processes
    if derived:
        derived.sort(key=lambda e: e["opts"].get("parent", ""))
        # a derived entry may itself be the parent of another: build in dependency order
        done = {e["id"] for e in entries if existing(args.out, e["id"])}
        pending = list(derived)
        while pending:
            ready = [e for e in pending if e["opts"].get("parent") in (None, *done)]
            if not ready:
                for e in pending:
                    failed.append((e["id"], "parent missing"))
                    results["fail"] += 1
                break
            pending = [e for e in pending if e not in ready]
            nproc = max(1, min(os.cpu_count() or 2, 16))
            with concurrent.futures.ProcessPoolExecutor(max_workers=nproc) as px:
                futs = {px.submit(fetch_one, e, args): e for e in ready}
                for fut in concurrent.futures.as_completed(futs):
                    e = futs[fut]
                    try:
                        st, p = fut.result()
                    except Exception as exn:  # noqa: BLE001
                        results["fail"] += 1
                        failed.append((e["id"], str(exn)[:160]))
                        print(f"  FAIL    {e['id']:<24} synth        {exn}", flush=True)
                        continue
                    results[st] += 1
                    done.add(e["id"])
                    if p and st == "ok":
                        sums[os.path.basename(p)] = sha256(p)
                        print(f"  built   {e['id']:<24} synth        "
                              f"{os.path.getsize(p) / 1e6:8.1f} MB", flush=True)

    with open(sums_path, "w") as f:
        for name in sorted(sums):
            f.write(f"{sums[name]}  {name}\n")
    cache = os.path.join(args.out, ".cache")
    if not args.keep_cache and os.path.isdir(cache) and not failed:
        shutil.rmtree(cache, ignore_errors=True)

    total = sum(os.path.getsize(os.path.join(args.out, n)) for n in os.listdir(args.out)
                if os.path.isfile(os.path.join(args.out, n)))
    print("\n" + "=" * 70)
    print(f" {results['ok']} fetched/built, {results['linked']} linked from v1, "
          f"{results['skip']} already present, {results['fail']} failed"
          f"   ({time.time() - t0:.0f} s)")
    print(f" {args.out}: {total / 1e9:.2f} GB")
    if failed:
        print(" failed (re-run to retry transient errors):")
        for i, msg in failed:
            print(f"   {i:<24} {msg}")
    print("=" * 70)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
