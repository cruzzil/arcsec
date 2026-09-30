#!/usr/bin/env python3
"""Run arcsec over the benchmark corpus and score it against the ground-truth WCS.

Truth is always taken from the delivered image header (for tier A that header carries
exactly what we asked hips2fits/SkyView for; for the synthetic tier S it is the parent's
truth carried through the transformation), or from a `<id>.truth` sidecar of FITS cards
when one exists. Errors are measured as angular separation at the image centre and the
four corners, so scale and rotation error are visible separately from pointing error.

Usage:
    scripts/benchmark.py [--arcsec target/release/arcsec] [--db ~/star_database]
                         [--db-name d80 | --auto-db] [--jobs 8] [--threads N]
                         [--images resources/testset] [--manifest scripts/test-images.tsv]
                         [--corpus] [--set v1] [--source ztf] [--dataset DSS2]
                         [--workdir /tmp/arcsec_bench_out] [--radius 5] [--timeout 300]
                         [--offset-hint 0.0] [--method quads|tetra] [--stars N]
                         [--max-corner-err 5] [--tier A] [--id foo] [--csv out.csv]
                         [--astap ~/astap_cli] [--by source,tier,fov]

--corpus switches to the expanded corpus (scripts/corpus.tsv, images in
resources/corpus). --tier, --id, --set, --source and --dataset are repeatable; --id
takes precedence over the rest, which combine with AND. --auto-db omits -D so arcsec
picks the database from the field size. --astap also solves every image with ASTAP (-D
from --db-name) and prints a head-to-head. Use --jobs 1 for meaningful timings.

A reported solve whose worst corner error exceeds --max-corner-err arcsec (or
--max-corner-px pixels, whichever is larger) counts as a false positive. When the truth has SIP or TPV distortion, the linear plate model arcsec
fits cannot reach the corners exactly; the threshold is then raised by the worst corner
error of the best possible linear fit (the "linear floor", reported per image).

Per-entry options in the manifest's `extra` column (key=value;...):
    file=<id>        solve another entry's image (negative controls reuse real fields)
    hint_dra=<deg>   hint offset in RA (true degrees on the sky), added to --offset-hint
    hint_ddec=<deg>  hint offset in Dec
    fov_scale=<x>    pass the FOV multiplied by x (wrong-scale controls)
    radius=<deg>     override -r for this entry
    expect=nosolve   (default for tier D) any reported solution is a false positive
    expect=any       tier D stress case: a correct solve is fine, a wrong one is not

The exit status is 1 if arcsec returned any false positive, 0 otherwise.

No third-party dependencies (no astropy/numpy).
"""

import argparse
import concurrent.futures
import math
import os
import re
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fitslite  # noqa: E402
from fitslite import Wcs, angsep  # noqa: E402,F401  (re-exported for callers)

IMAGE_EXTS = (".fits", ".fit", ".fts", ".fits.fz", ".fits.gz", ".fit.gz", ".xisf")

# Installed-catalogue coverage is decided by FOV (image height, which is what --fov
# carries); see arcsec/src/db_select.rs.
DB_FOV_RANGES = {"d80": (0.15, 6.0), "v50": (0.20, 6.0), "d50": (0.20, 6.0),
                 "d20": (0.30, 6.0), "v05": (0.60, 6.0), "d05": (0.60, 6.0),
                 "g05": (3.0, 20.0), "w08": (20.0, 80.0)}


# ── FITS header parsing ──────────────────────────────────────────────────────


def read_header(path, max_blocks=200):
    """Return {keyword: value} for the first image HDU (the one arcsec reads)."""
    return fitslite.read_header(path, max_blocks=max_blocks)


def read_truth_sidecar(path):
    """A `.truth` sidecar: FITS-style cards, one per line (same syntax as .wcs)."""
    return parse_wcs_file(path)


# ── arcsec invocation ────────────────────────────────────────────────────────


def parse_wcs_file(path):
    """Parse arcsec's .wcs output (80-char FITS cards, one per line)."""
    h = {}
    with open(path, "r", errors="replace") as f:
        for line in f:
            key = line[:8].strip()
            if not key or key == "END" or line[8:10] != "= ":
                continue
            val = line[10:].split("/")[0].strip()
            if val.startswith("'"):
                h[key] = val.strip("'").strip()
            elif val in ("T", "F"):
                h[key] = val == "T"
            else:
                try:
                    h[key] = float(val)
                except ValueError:
                    pass
    return h


def parse_astap_ini(path):
    """Parse astap_cli's .ini output (KEY=value, one per line)."""
    h = {}
    if not os.path.exists(path):
        return h
    with open(path, "r", errors="replace") as f:
        for line in f:
            line = line.strip()
            if "=" not in line:
                continue
            k, v = line.split("=", 1)
            k = k.strip().upper()
            v = v.strip()
            if k == "PLTSOLVD":
                h[k] = v == "T"
            else:
                try:
                    h[k] = float(v)
                except ValueError:
                    h[k] = v
    return h


def score(truth, sol, naxis1, naxis2):
    """Errors of `sol` against `truth` at the centre and the four corners."""
    pts = [((naxis1 + 1) / 2.0, (naxis2 + 1) / 2.0),
           (1.0, 1.0), (naxis1, 1.0), (1.0, naxis2), (naxis1, naxis2)]
    errs = []
    for (x, y) in pts:
        tr = truth.pix2sky(x, y)
        sv = sol.pix2sky(x, y)
        errs.append(angsep(tr[0], tr[1], sv[0], sv[1]))
    return errs


def run_astap(args, path, hint_ra, hint_dec, fov_deg, radius, out_base, truth, naxis1,
              naxis2, limit):
    """Solve with astap_cli and score it the same way. Returns a dict."""
    r = {"astap_status": "", "astap_secs": 0.0, "astap_centre": None,
         "astap_corner": None, "astap_scale_err": None}
    for ext in (".wcs", ".ini"):
        try:
            os.remove(out_base + ext)
        except OSError:
            pass
    cmd = [args.astap, "-f", path, "-d", args.db, "-D", args.db_name,
           "-ra", f"{(hint_ra % 360.0) / 15.0:.9f}",
           "-spd", f"{hint_dec + 90.0:.9f}",
           "-fov", f"{fov_deg:.9f}",
           "-r", str(radius),
           "-o", out_base]
    t0 = time.time()
    try:
        subprocess.run(cmd, capture_output=True, text=True, timeout=args.timeout)
    except subprocess.TimeoutExpired:
        r["astap_status"] = "TIMEOUT"
        r["astap_secs"] = round(time.time() - t0, 2)
        return r
    r["astap_secs"] = round(time.time() - t0, 2)

    ini = parse_astap_ini(out_base + ".ini")
    if not ini.get("PLTSOLVD"):
        r["astap_status"] = "NO_SOLVE"
        return r
    sol = Wcs.from_header(ini)
    if sol is None:
        r["astap_status"] = "BAD_WCS"
        return r
    errs = score(truth, sol, naxis1, naxis2)
    r["astap_centre"] = round(errs[0], 3)
    r["astap_corner"] = round(max(errs[1:]), 3)
    r["astap_scale_err"] = round(abs(sol.pixscale() - truth.pixscale()) / truth.pixscale() * 100.0, 5)
    r["astap_status"] = "WRONG" if r["astap_corner"] > limit else "OK"
    return r


def installed_dbs(db_dir):
    out = set()
    try:
        for n in os.listdir(db_dir):
            p = n.split("_")[0]
            if p in DB_FOV_RANGES:
                out.add(p)
    except OSError:
        pass
    return out


def catalogue_covers(fov_deg, dbs, db_name, auto):
    names = dbs if auto else ({db_name} & dbs)
    return any(DB_FOV_RANGES[n][0] <= fov_deg <= DB_FOV_RANGES[n][1] for n in names)


def run_one(entry, args):
    """Solve one image and score it. Returns a result dict."""
    iid, tier, path = entry["id"], entry["tier"], entry["path"]
    opts = entry.get("opts", {})
    res = {"id": iid, "tier": tier, "status": "", "secs": 0.0,
           "err_centre": None, "err_corner": None, "scale_err_pct": None,
           "rot_err_deg": None, "nstars": None, "nquads": None, "note": "",
           "source": entry.get("source", ""), "dataset": entry.get("dataset", ""),
           "sets": entry.get("sets", ""), "truth_q": entry.get("truth", ""),
           "expect": entry.get("expect", ""), "lin_floor": None, "cat_ok": True}

    side = entry.get("truth_path")
    # XISF has no FITS header to read; its sidecar carries NAXIS1/2 as well as the WCS
    hdr = read_header(path) if not path.endswith(".xisf") else {}
    th = read_truth_sidecar(side) if side else hdr
    if not hdr:
        hdr = th
    truth = Wcs.from_header(th)
    if truth is None:
        res["status"] = "NO_TRUTH"
        return res

    naxis1 = int(hdr.get("NAXIS1", 0) or 0)
    naxis2 = int(hdr.get("NAXIS2", 0) or 0)
    if naxis1 == 0 or naxis2 == 0:
        res["status"] = "NO_TRUTH"
        return res

    ps_deg = truth.pixscale()
    fov_deg = ps_deg * max(naxis1, naxis2)
    # What -fov means to ASTAP (and what N.I.N.A. sends): the image height.
    fov_height_deg = ps_deg * naxis2
    cra, cdec = truth.pix2sky((naxis1 + 1) / 2.0, (naxis2 + 1) / 2.0)
    res["fov_deg"] = round(fov_deg, 4)
    res["pixscale_as"] = round(ps_deg * 3600.0, 3)
    res["cat_ok"] = catalogue_covers(fov_height_deg, args.dbs, args.db_name, args.auto_db)

    # A linear plate cannot follow SIP/TPV distortion to the corners; allow for the
    # best linear fit's own corner error so a correct solve is not called wrong.
    floor = 0.0
    if truth.distorted:
        _, floor = fitslite.best_linear(truth, naxis1, naxis2)
        res["lin_floor"] = round(floor, 3)
    # An absolute 5" means a quarter pixel at TESS's 21"/px: never less than
    # --max-corner-px pixels, so coarse images are not held to a sub-pixel standard.
    limit = max(args.max_corner_err, args.max_corner_px * ps_deg * 3600.0) + floor
    res["limit"] = limit

    # Hint: truth centre, optionally pushed off by N field widths, plus any per-entry
    # offset in true degrees on the sky.
    cosd = max(math.cos(math.radians(cdec)), 1e-6)
    dra = args.offset_hint * fov_deg + float(opts.get("hint_dra", 0.0))
    ddec = args.offset_hint * fov_deg + float(opts.get("hint_ddec", 0.0))
    hint_ra = cra + dra / cosd
    hint_dec = max(-89.9, min(89.9, cdec + ddec))
    fov_hint = fov_height_deg * float(opts.get("fov_scale", 1.0))
    radius = opts.get("radius", args.radius)

    out_base = os.path.join(args.workdir, iid)
    for ext in (".wcs", ".ini", ".log"):
        try:
            os.remove(out_base + ext)
        except OSError:
            pass

    cmd = [args.arcsec, "-f", path, "-d", args.db]
    if not args.auto_db:
        cmd += ["-D", args.db_name]
    cmd += [
        "--ra", f"{(hint_ra % 360.0) / 15.0:.9f}",
        "--spd", f"{hint_dec + 90.0:.9f}",
        "--fov", f"{fov_hint:.9f}",
        "-r", str(radius),
        "-o", out_base,
    ]
    if args.method != "quads":
        cmd += ["--method", args.method]
    if args.stars:
        cmd += ["-s", str(args.stars)]
    if args.threads is not None:
        cmd += ["--threads", str(args.threads)]
    cmd += args.extra_arg

    t0 = time.time()
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=args.timeout)
        rc = proc.returncode
    except subprocess.TimeoutExpired:
        res["status"] = "TIMEOUT"
        res["secs"] = round(time.time() - t0, 2)
        return res
    res["secs"] = round(time.time() - t0, 2)

    def with_astap(r):
        if args.astap:
            r.update(run_astap(args, path, hint_ra, hint_dec, fov_hint, radius,
                               out_base + "_astap", truth, naxis1, naxis2, limit))
        return r

    if rc != 0 or not os.path.exists(out_base + ".wcs"):
        res["status"] = "NO_SOLVE"
        res["note"] = f"exit={rc}"
        return with_astap(res)

    sol_h = parse_wcs_file(out_base + ".wcs")
    # A solution's SIP terms are read as written (see fitslite.Wcs.from_header).
    sol = Wcs.from_header(sol_h, keep_sip_constant=True)
    if sol is None:
        res["status"] = "BAD_WCS"
        return with_astap(res)

    errs = score(truth, sol, naxis1, naxis2)
    res["status"] = "OK"
    res["err_centre"] = round(errs[0], 3)
    res["err_corner"] = round(max(errs[1:]), 3)
    res["scale_err_pct"] = round(abs(sol.pixscale() - ps_deg) / ps_deg * 100.0, 5)
    dr = abs(sol.rotation() - truth.rotation()) % 360.0
    res["rot_err_deg"] = round(min(dr, 360.0 - dr), 4)

    with_astap(res)

    ini = out_base + ".ini"
    if os.path.exists(ini):
        for line in open(ini):
            if line.startswith("NSTARS="):
                res["nstars"] = line.strip().split("=")[1]
            elif line.startswith("NQUADS="):
                res["nquads"] = line.strip().split("=")[1]
            elif line.startswith("RMS="):
                res["note"] = "rms=" + line.strip().split("=")[1]
    return res


# ── manifest ─────────────────────────────────────────────────────────────────


def parse_opts(extra):
    """key=value;key=value -> dict (bare tokens are kept under their own name)."""
    out = {}
    if not extra or "=" not in extra:
        return out
    for part in extra.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def find_image(images_dir, iid):
    for ext in IMAGE_EXTS:
        p = os.path.join(images_dir, iid + ext)
        if os.path.exists(p):
            return p
    return None


def load_manifest(manifest, images_dir):
    """Entries whose image is present. Columns (whitespace-separated):

        id tier ra dec fov_deg width height source extra [dataset sets truth]

    The last three are optional, so scripts/test-images.tsv still loads."""
    rows = []
    with open(manifest) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = re.split(r"[ \t]+", line)
            if len(parts) < 8:
                continue
            rows.append(parts)
    out = []
    for parts in rows:
        iid, tier = parts[0], parts[1]
        extra = parts[8] if len(parts) > 8 else ""
        opts = parse_opts(extra)
        target = opts.get("file", iid)
        p = find_image(images_dir, target)
        if not p:
            continue
        tp = os.path.join(images_dir, target + ".truth")
        e = {"id": iid, "tier": tier, "path": p, "fov_manifest": float(parts[4]),
             "source": parts[7], "extra": extra, "opts": opts,
             "dataset": parts[9] if len(parts) > 9 else "",
             "sets": parts[10] if len(parts) > 10 else "",
             "truth": parts[11] if len(parts) > 11 else "",
             "truth_path": tp if os.path.exists(tp) else None}
        e["expect"] = opts.get("expect", "nosolve" if tier == "D" else "solve")
        out.append(e)
    return out


def fov_band(fov):
    if fov is None:
        return "?"
    for lo, hi, name in ((0, 0.15, "<0.15"), (0.15, 0.3, "0.15-0.3"), (0.3, 0.6, "0.3-0.6"),
                         (0.6, 1.2, "0.6-1.2"), (1.2, 2.5, "1.2-2.5"), (2.5, 6, "2.5-6"),
                         (6, 20, "6-20"), (20, 1000, ">20")):
        if lo <= fov < hi:
            return name
    return "?"


# ── main ─────────────────────────────────────────────────────────────────────


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser()
    ap.add_argument("--arcsec", default=os.path.join(here, "..", "target", "release", "arcsec"))
    ap.add_argument("--manifest", default=None,
                    help="default scripts/test-images.tsv, or scripts/corpus.tsv with --corpus")
    ap.add_argument("--images", default=None,
                    help="default resources/testset, or resources/corpus with --corpus")
    ap.add_argument("--corpus", action="store_true",
                    help="use the expanded corpus (scripts/corpus.tsv, resources/corpus)")
    ap.add_argument("--db", default=os.path.expanduser("~/star_database"))
    ap.add_argument("--db-name", default="d80")
    ap.add_argument("--auto-db", action="store_true",
                    help="omit -D and let arcsec pick the database from the field size "
                         "(D-series 0.15-6 deg, G05 3-20, W08 20-80)")
    ap.add_argument("--workdir", default="/tmp/arcsec_bench_out")
    ap.add_argument("--csv", default=None)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--radius", type=float, default=5.0)
    ap.add_argument("--timeout", type=float, default=300.0)
    ap.add_argument("--offset-hint", type=float, default=0.0,
                    help="push the hint off truth by this many field widths")
    ap.add_argument("--method", default="quads")
    ap.add_argument("--stars", type=int, default=None,
                    help="pass -s to arcsec (max detected stars)")
    ap.add_argument("--threads", type=int, default=None,
                    help="pass --threads to arcsec (0 = one per core)")
    ap.add_argument("--astap", default=None,
                    help="path to astap_cli; when given, solve each image with both and compare")
    ap.add_argument("--max-corner-px", type=float, default=1.0,
                    help="the false-positive corner threshold is at least this many "
                         "pixels (matters only above 5\"/px)")
    ap.add_argument("--max-corner-err", type=float, default=5.0,
                    help="corner error (arcsec) above which a reported solve is "
                         "counted as a FALSE POSITIVE, not a success")
    ap.add_argument("--extra-arg", action="append", default=[],
                    help="pass this argument to arcsec as well (repeatable), e.g. --extra-arg=--sip")
    ap.add_argument("--tier", action="append", default=[])
    ap.add_argument("--id", action="append", default=[])
    ap.add_argument("--set", action="append", default=[],
                    help="only entries tagged with this subset (e.g. v1, core)")
    ap.add_argument("--source", action="append", default=[],
                    help="only entries fetched from this source (e.g. ztf, tess)")
    ap.add_argument("--dataset", action="append", default=[],
                    help="only entries whose dataset contains this text (e.g. DSS2)")
    ap.add_argument("--by", default=None,
                    help="comma list of breakdowns: tier,source,dataset,fov,set "
                         "(default: tier,source,fov when the manifest has those columns)")
    ap.add_argument("--quiet", action="store_true", help="no per-image lines")
    args = ap.parse_args()

    if args.manifest is None:
        args.manifest = os.path.join(here, "corpus.tsv" if args.corpus else "test-images.tsv")
    if args.images is None:
        args.images = os.path.join(here, "..", "resources", "corpus" if args.corpus else "testset")
    args.arcsec = os.path.abspath(args.arcsec)
    args.images = os.path.abspath(args.images)
    args.dbs = installed_dbs(args.db)
    os.makedirs(args.workdir, exist_ok=True)

    entries = load_manifest(args.manifest, args.images)
    if args.id:
        entries = [e for e in entries if e["id"] in args.id]
    else:
        if args.tier:
            entries = [e for e in entries if e["tier"] in args.tier]
        if args.set:
            entries = [e for e in entries
                       if set(e["sets"].split(",")) & set(args.set)]
        if args.source:
            entries = [e for e in entries if e["source"] in args.source]
        if args.dataset:
            entries = [e for e in entries
                       if any(d.lower() in e["dataset"].lower() for d in args.dataset)]
    if not entries:
        print("no images matched", file=sys.stderr)
        return 1

    print(f"arcsec   : {args.arcsec}")
    print(f"database : {args.db} ({'auto' if args.auto_db else args.db_name};"
          f" installed: {','.join(sorted(args.dbs)) or 'none'})")
    print(f"images   : {len(entries)} from {args.images}")
    print(f"hint     : truth centre + {args.offset_hint} field widths, -r {args.radius}")
    print()

    results = []

    def astap_tag(r):
        if not args.astap:
            return ""
        st = r.get("astap_status", "")
        if st == "OK":
            return f"   | ASTAP OK  centre={r['astap_centre']:.3f}\""
        if st == "WRONG":
            return f"   | ASTAP WRONG corner={r['astap_corner']:.1f}\""
        return f"   | ASTAP {st or '-'}"

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = {ex.submit(run_one, e, args): e for e in entries}
        for fut in concurrent.futures.as_completed(futs):
            r = fut.result()
            results.append(r)
            if (r.get("astap_status") == "WRONG" and r["lin_floor"] is not None
                    and r["astap_centre"] <= max(args.max_corner_err, 10.0, 2.0 * r["pixscale_as"])):
                r["astap_status"] = "INEXACT"
            if r["status"] == "OK" and r["err_corner"] > r.get("limit", args.max_corner_err):
                # Distorted truth, right place, corners beyond what a linear plate can
                # reach even allowing for the best linear fit: not a false positive,
                # not a success either.
                if r["lin_floor"] is not None and \
                        r["err_centre"] <= max(args.max_corner_err, 10.0, 2.0 * r["pixscale_as"]):
                    r["status"] = "INEXACT"
                else:
                    r["status"] = "WRONG"
            if args.quiet:
                continue
            if r["status"] == "WRONG":
                print(f"  {r['id']:<22} {r['tier']}  WRONG !! "
                      f"centre={r['err_centre']:>8.3f}\"  corner={r['err_corner']:>9.3f}\"  "
                      f"rot={r['rot_err_deg']:>7.3f}d  <-- FALSE POSITIVE" + astap_tag(r))
            elif r["status"] == "INEXACT":
                print(f"  {r['id']:<22} {r['tier']}  INEXACT  "
                      f"centre={r['err_centre']:>8.3f}\"  corner={r['err_corner']:>9.3f}\"  "
                      f"linear floor={r['lin_floor']:.1f}\"  (distortion)" + astap_tag(r))
            elif r["status"] == "OK":
                print(f"  {r['id']:<22} {r['tier']}  OK       "
                      f"centre={r['err_centre']:>8.3f}\"  corner={r['err_corner']:>9.3f}\"  "
                      f"scale={r['scale_err_pct']:>8.4f}%  {r['secs']:>6.2f}s" + astap_tag(r))
            else:
                print(f"  {r['id']:<22} {r['tier']}  {r['status']:<8} "
                      f"{r['note']:<24} {r['secs']:>6.2f}s" + astap_tag(r))

    results.sort(key=lambda r: (r["tier"], r["id"]))

    def is_fp(r):
        """A wrong answer: a bad solve anywhere, or any solve of a must-fail control."""
        return r["status"] == "WRONG" or (r["expect"] == "nosolve"
                                          and r["status"] in ("OK", "INEXACT"))

    if args.csv:
        cols = ["id", "tier", "source", "dataset", "sets", "truth_q", "status", "fov_deg",
                "pixscale_as", "cat_ok", "lin_floor", "err_centre", "err_corner",
                "scale_err_pct", "rot_err_deg", "nstars", "nquads", "secs", "note",
                "astap_status", "astap_centre", "astap_corner", "astap_scale_err",
                "astap_secs"]
        with open(args.csv, "w") as f:
            f.write(",".join(cols) + "\n")
            for r in results:
                vals = ("" if r.get(c) is None else str(r.get(c, "")) for c in cols)
                f.write(",".join(f'"{v}"' if "," in v else v for v in vals) + "\n")
        print(f"\nCSV: {args.csv}")

    def summarise(tier_filter, label):
        rs = [r for r in results if r["tier"] in tier_filter]
        if not rs:
            return
        ok = [r for r in rs if r["status"] == "OK"]
        wrong = [r for r in rs if r["status"] == "WRONG"]
        inexact = [r for r in rs if r["status"] == "INEXACT"]
        nocat = [r for r in rs if r["status"] not in ("OK", "INEXACT") and not r["cat_ok"]]
        print(f"\n  {label}: {len(ok)}/{len(rs)} correct"
              f"   ({len(wrong)} FALSE POSITIVE{'S' if len(wrong) != 1 else ''}"
              + (f", {len(inexact)} inexact (distortion)" if inexact else "")
              + f", {len(rs) - len(ok) - len(wrong) - len(inexact)} no-solve"
              + (f", of which {len(nocat)} outside the installed catalogues" if nocat else "")
              + ")")
        for r in wrong:
            print(f"    !! {r['id']:<22} corner={r['err_corner']:.1f}\"  "
                  f"rot={r['rot_err_deg']:.2f} deg  nquads={r['nquads']}")
        if ok:
            c = sorted(r["err_centre"] for r in ok)
            k = sorted(r["err_corner"] for r in ok)
            s = sorted(r["scale_err_pct"] for r in ok)
            t = sorted(r["secs"] for r in ok)
            med = lambda v: v[len(v) // 2]
            print(f"    centre err   median {med(c):8.3f}\"   max {c[-1]:9.3f}\"")
            print(f"    corner err   median {med(k):8.3f}\"   max {k[-1]:9.3f}\"")
            print(f"    scale err    median {med(s):8.4f}%   max {s[-1]:9.4f}%")
            print(f"    time         median {med(t):8.2f}s   max {t[-1]:9.2f}s")

    print("\n" + "=" * 78)
    print(" SUMMARY")
    print("=" * 78)
    summarise({"A"}, "Tier A (synthetic, exact truth)")
    summarise({"B"}, "Tier B (survey pixels, header truth)")
    summarise({"C"}, "Tier C (real observing frames, pipeline/reference truth)")
    summarise({"S"}, "Tier S (simulated camera artefacts, derived exact truth)")
    tierd = [r for r in results if r["tier"] == "D"]
    if tierd:
        solved = [r for r in tierd if r["status"] in ("OK", "WRONG", "INEXACT")]
        bad = [r for r in tierd if is_fp(r)]
        print(f"\n  Tier D (negative/stress): {len(solved)}/{len(tierd)} returned a solution,"
              f" {len(bad)} false positive{'s' if len(bad) != 1 else ''}")
        for r in tierd:
            if is_fp(r):
                mark = "FALSE POSITIVE"
            elif r["status"] == "OK":
                mark = "solved correctly (expect=any)"
            else:
                mark = "correctly failed"
            extra = f" centre={r['err_centre']}\"" if r["status"] in ("OK", "WRONG", "INEXACT") else ""
            print(f"    {r['id']:<22} {r['status']:<9} {mark}{extra}")

    # Breakdowns (tiers A/B/C/S; D is reported above).
    has_cols = any(r["source"] for r in results) and args.manifest.endswith("corpus.tsv")
    by = args.by.split(",") if args.by else (["tier", "source", "fov"] if has_cols else [])
    pos = [r for r in results if r["expect"] != "nosolve"]

    def keyf(dim):
        if dim == "fov":
            return lambda r: fov_band(r.get("fov_deg"))
        if dim == "set":
            return lambda r: r["sets"] or "-"
        return lambda r: r.get(dim) or "-"

    for dim in by:
        groups = {}
        for r in pos:
            groups.setdefault(keyf(dim)(r), []).append(r)
        if len(groups) < 2 and dim != "tier":
            continue
        print("\n" + "-" * 78)
        print(f" BY {dim.upper()}   (tiers A/B/C/S and expect=any controls)")
        print(f"   {'':<22}{'n':>5}{'correct':>9}{'rate':>7}{'FP':>5}{'inexact':>8}{'nocat':>7}"
              f"{'ctr med':>9}{'crn med':>9}{'t med':>8}")
        order = sorted(groups)
        if dim == "fov":
            bands = ["<0.15", "0.15-0.3", "0.3-0.6", "0.6-1.2", "1.2-2.5", "2.5-6", "6-20", ">20", "?"]
            order = [b for b in bands if b in groups]
        for g in order:
            rs = groups[g]
            ok = [r for r in rs if r["status"] == "OK"]
            fp = [r for r in rs if r["status"] == "WRONG"]
            ie = [r for r in rs if r["status"] == "INEXACT"]
            nc = [r for r in rs if r["status"] not in ("OK", "INEXACT") and not r["cat_ok"]]
            med = lambda v: sorted(v)[len(v) // 2] if v else float("nan")
            print(f"   {g:<22}{len(rs):>5}{len(ok):>9}{100.0 * len(ok) / len(rs):>6.0f}%{len(fp):>5}{len(ie):>8}"
                  f"{len(nc):>7}{med([r['err_centre'] for r in ok]):>8.2f}\""
                  f"{med([r['err_corner'] for r in ok]):>8.2f}\"{med([r['secs'] for r in ok]):>7.2f}s")

    if args.astap:
        pv = [r for r in results if r["tier"] in ("A", "B", "C", "S")]
        p_ok = [r for r in pv if r["status"] == "OK"]
        p_wr = [r for r in pv if r["status"] == "WRONG"]
        a_ok = [r for r in pv if r.get("astap_status") == "OK"]
        a_wr = [r for r in pv if r.get("astap_status") == "WRONG"]
        both = [r for r in pv if r["status"] == "OK" and r.get("astap_status") == "OK"]
        only_p = [r for r in pv if r["status"] == "OK" and r.get("astap_status") != "OK"]
        only_a = [r for r in pv if r["status"] != "OK" and r.get("astap_status") == "OK"]
        neither = [r for r in pv if r["status"] != "OK" and r.get("astap_status") != "OK"]
        d_astap = [r for r in tierd if r.get("astap_status") == "WRONG"
                   or (r["expect"] == "nosolve" and r.get("astap_status") == "OK")]

        print("\n" + "=" * 78)
        print(" HEAD TO HEAD vs ASTAP  (tiers A/B/C/S, n=%d)" % len(pv))
        print("=" * 78)
        print(f"  {'':<22}{'arcsec':>12}{'ASTAP':>12}")
        print(f"  {'correct':<22}{len(p_ok):>12}{len(a_ok):>12}")
        print(f"  {'false positives':<22}{len(p_wr):>12}{len(a_wr):>12}")
        print(f"  {'no solve':<22}{len(pv)-len(p_ok)-len(p_wr):>12}{len(pv)-len(a_ok)-len(a_wr):>12}")
        if tierd:
            print(f"  {'tier D false pos.':<22}{len([r for r in tierd if is_fp(r)]):>12}{len(d_astap):>12}")

        def stat(rs, kc, kk, ks, kt):
            if not rs:
                return None
            m = lambda k: sorted(x[k] for x in rs if x.get(k) is not None)
            c, kk_, ss, tt = m(kc), m(kk), m(ks), m(kt)
            md = lambda v: v[len(v) // 2] if v else float("nan")
            return md(c), md(kk_), md(ss), md(tt), (sum(tt) if tt else 0.0)

        ps = stat(both, "err_centre", "err_corner", "scale_err_pct", "secs")
        as_ = stat(both, "astap_centre", "astap_corner", "astap_scale_err", "astap_secs")
        if ps and as_ and both:
            print(f"\n  On the {len(both)} images BOTH solved correctly (medians):")
            print(f"  {'':<22}{'arcsec':>12}{'ASTAP':>12}")
            print(f"  {'centre err (\")':<22}{ps[0]:>12.3f}{as_[0]:>12.3f}")
            print(f"  {'corner err (\")':<22}{ps[1]:>12.3f}{as_[1]:>12.3f}")
            print(f"  {'scale err (%)':<22}{ps[2]:>12.4f}{as_[2]:>12.4f}")
            print(f"  {'time (s)':<22}{ps[3]:>12.3f}{as_[3]:>12.3f}")
            if as_[3] > 0:
                print(f"  {'speedup':<22}{as_[3]/max(ps[3],1e-6):>11.2f}x{'':>12}")

        print(f"\n  both correct: {len(both)}   arcsec only: {len(only_p)}   "
              f"ASTAP only: {len(only_a)}   neither: {len(neither)}")
        if only_p and len(only_p) <= 60:
            print("  arcsec solved, ASTAP did not:")
            for r in only_p:
                print(f"    + {r['id']:<22} {r.get('astap_status','-')}")
        if only_a:
            print("  ASTAP solved, arcsec did not:")
            for r in only_a:
                print(f"    - {r['id']:<22} arcsec={r['status']}")
        print("=" * 78)

    # Non-zero when arcsec reported any wrong solution, so a script or CI step can
    # gate on the false-positive count without parsing the summary.
    return 1 if any(is_fp(r) for r in results) else 0


if __name__ == "__main__":
    sys.exit(main())
