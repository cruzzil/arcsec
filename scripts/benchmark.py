#!/usr/bin/env python3
"""Run arcsec over the benchmark corpus and score it against the ground-truth WCS.

Truth is always taken from the delivered FITS header (for tier A that header carries
exactly what we asked hips2fits/SkyView for). Errors are measured as angular separation
at the image centre and the four corners, so scale and rotation error are visible
separately from pointing error.

Usage:
    scripts/benchmark.py [--db ~/star_database] [--db-name d80] [--jobs 8]
                         [--images resources/testset] [--radius 5]
                         [--offset-hint 0.0] [--tier A] [--id foo] [--csv out.csv]

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

# ── FITS header parsing ──────────────────────────────────────────────────────


def read_header(path, max_blocks=200):
    """Return {keyword: value} for the primary HDU. Values are float/str/bool."""
    hdr = {}
    with open(path, "rb") as f:
        for _ in range(max_blocks):
            block = f.read(2880)
            if len(block) < 2880:
                break
            done = False
            for i in range(0, 2880, 80):
                card = block[i:i + 80].decode("latin-1")
                key = card[:8].strip()
                if key == "END":
                    done = True
                    break
                if not key or card[8:10] != "= ":
                    continue
                val = card[10:].split("/")[0].strip()
                if val.startswith("'"):
                    hdr[key] = val.strip("'").strip()
                elif val in ("T", "F"):
                    hdr[key] = val == "T"
                else:
                    try:
                        hdr[key] = float(val.replace("E", "e").replace("D", "e"))
                    except ValueError:
                        hdr[key] = val
            if done:
                break
    return hdr


class Wcs:
    """A TAN WCS: CRPIX/CRVAL plus a 2x2 CD matrix in degrees per pixel."""

    def __init__(self, crpix1, crpix2, crval1, crval2, cd):
        self.crpix1, self.crpix2 = crpix1, crpix2
        self.crval1, self.crval2 = crval1, crval2
        self.cd = cd  # [[cd1_1, cd1_2], [cd2_1, cd2_2]]

    @classmethod
    def from_header(cls, h):
        try:
            crpix1, crpix2 = h["CRPIX1"], h["CRPIX2"]
            crval1, crval2 = h["CRVAL1"], h["CRVAL2"]
        except KeyError:
            return None

        if "CD1_1" in h:
            cd = [[h.get("CD1_1", 0.0), h.get("CD1_2", 0.0)],
                  [h.get("CD2_1", 0.0), h.get("CD2_2", 0.0)]]
        else:
            cdelt1 = h.get("CDELT1")
            cdelt2 = h.get("CDELT2")
            if cdelt1 is None or cdelt2 is None:
                return None
            if "PC1_1" in h:
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
        return cls(crpix1, crpix2, crval1, crval2, cd)

    def pix2sky(self, x, y):
        """1-based FITS pixel coords -> (ra_deg, dec_deg) via the TAN deprojection."""
        u = x - self.crpix1
        v = y - self.crpix2
        xi = math.radians(self.cd[0][0] * u + self.cd[0][1] * v)
        eta = math.radians(self.cd[1][0] * u + self.cd[1][1] * v)

        ra0 = math.radians(self.crval1)
        dec0 = math.radians(self.crval2)
        sd0, cd0 = math.sin(dec0), math.cos(dec0)

        denom = cd0 - eta * sd0
        ra = ra0 + math.atan2(xi, denom)
        dec = math.atan2(sd0 + eta * cd0, math.hypot(denom, xi))
        return math.degrees(ra) % 360.0, math.degrees(dec)

    def pixscale(self):
        """Degrees per pixel along the first image axis."""
        return math.hypot(self.cd[0][0], self.cd[1][0])

    def rotation(self):
        return math.degrees(math.atan2(self.cd[1][0], self.cd[1][1]))


def angsep(ra1, dec1, ra2, dec2):
    """Angular separation in arcseconds."""
    r1, d1 = math.radians(ra1), math.radians(dec1)
    r2, d2 = math.radians(ra2), math.radians(dec2)
    c = math.sin(d1) * math.sin(d2) + math.cos(d1) * math.cos(d2) * math.cos(r1 - r2)
    return math.degrees(math.acos(max(-1.0, min(1.0, c)))) * 3600.0


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


def run_astap(args, path, hint_ra, hint_dec, fov_deg, out_base, truth, naxis1, naxis2):
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
           "-r", str(args.radius),
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
    r["astap_status"] = "WRONG" if r["astap_corner"] > args.max_corner_err else "OK"
    return r


def run_one(entry, args):
    """Solve one image and score it. Returns a result dict."""
    iid, tier, path = entry["id"], entry["tier"], entry["path"]
    res = {"id": iid, "tier": tier, "status": "", "secs": 0.0,
           "err_centre": None, "err_corner": None, "scale_err_pct": None,
           "rot_err_deg": None, "nstars": None, "nquads": None, "note": ""}

    hdr = read_header(path)
    truth = Wcs.from_header(hdr)
    if truth is None:
        res["status"] = "NO_TRUTH"
        return res

    naxis1 = int(hdr.get("NAXIS1", 0))
    naxis2 = int(hdr.get("NAXIS2", 0))
    if naxis1 == 0 or naxis2 == 0:
        res["status"] = "NO_TRUTH"
        return res

    ps_deg = truth.pixscale()
    fov_deg = ps_deg * max(naxis1, naxis2)
    cra, cdec = truth.pix2sky((naxis1 + 1) / 2.0, (naxis2 + 1) / 2.0)
    res["fov_deg"] = round(fov_deg, 4)
    res["pixscale_as"] = round(ps_deg * 3600.0, 3)

    # Hint: truth centre, optionally pushed off by N field widths.
    hint_ra = cra + args.offset_hint * fov_deg / max(math.cos(math.radians(cdec)), 1e-6)
    hint_dec = max(-89.9, min(89.9, cdec + args.offset_hint * fov_deg))

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
        "--fov", f"{fov_deg:.9f}",
        "-r", str(args.radius),
        "-o", out_base,
    ]
    if args.method != "quads":
        cmd += ["--method", args.method]
    if args.stars:
        cmd += ["-s", str(args.stars)]
    if args.threads is not None:
        cmd += ["--threads", str(args.threads)]

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
            r.update(run_astap(args, path, hint_ra, hint_dec, fov_deg,
                               out_base + "_astap", truth, naxis1, naxis2))
        return r

    if rc != 0 or not os.path.exists(out_base + ".wcs"):
        res["status"] = "NO_SOLVE"
        res["note"] = f"exit={rc}"
        return with_astap(res)

    sol_h = parse_wcs_file(out_base + ".wcs")
    sol = Wcs.from_header(sol_h)
    if sol is None:
        res["status"] = "BAD_WCS"
        return with_astap(res)

    # Compare at the centre and the four corners.
    pts = [((naxis1 + 1) / 2.0, (naxis2 + 1) / 2.0),
           (1.0, 1.0), (naxis1, 1.0), (1.0, naxis2), (naxis1, naxis2)]
    errs = []
    for (x, y) in pts:
        tr = truth.pix2sky(x, y)
        sv = sol.pix2sky(x, y)
        errs.append(angsep(tr[0], tr[1], sv[0], sv[1]))

    res["status"] = "OK"
    res["err_centre"] = round(errs[0], 3)
    res["err_corner"] = round(max(errs[1:]), 3)
    res["scale_err_pct"] = round(abs(sol.pixscale() - ps_deg) / ps_deg * 100.0, 5)
    dr = abs(sol.rotation() - truth.rotation()) % 360.0
    res["rot_err_deg"] = round(min(dr, 360.0 - dr), 4)

    if args.astap:
        res.update(run_astap(args, path, hint_ra, hint_dec, fov_deg,
                             out_base + "_astap", truth, naxis1, naxis2))

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


# ── main ─────────────────────────────────────────────────────────────────────


def load_manifest(manifest, images_dir):
    out = []
    with open(manifest) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = re.split(r"[ \t]+", line)
            if len(parts) < 8:
                continue
            iid, tier = parts[0], parts[1]
            p = os.path.join(images_dir, iid + ".fits")
            if os.path.exists(p):
                out.append({"id": iid, "tier": tier, "path": p,
                            "fov_manifest": float(parts[4])})
    return out


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser()
    ap.add_argument("--arcsec", default=os.path.join(here, "..", "target", "release", "arcsec"))
    ap.add_argument("--manifest", default=os.path.join(here, "test-images.tsv"))
    ap.add_argument("--images", default=os.path.join(here, "..", "resources", "testset"))
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
    ap.add_argument("--max-corner-err", type=float, default=5.0,
                    help="corner error (arcsec) above which a reported solve is "
                         "counted as a FALSE POSITIVE, not a success")
    ap.add_argument("--tier", action="append", default=[])
    ap.add_argument("--id", action="append", default=[])
    args = ap.parse_args()

    args.arcsec = os.path.abspath(args.arcsec)
    args.images = os.path.abspath(args.images)
    os.makedirs(args.workdir, exist_ok=True)

    entries = load_manifest(args.manifest, args.images)
    if args.id:
        entries = [e for e in entries if e["id"] in args.id]
    elif args.tier:
        entries = [e for e in entries if e["tier"] in args.tier]
    if not entries:
        print("no images matched", file=sys.stderr)
        return 1

    print(f"arcsec   : {args.arcsec}")
    print(f"database : {args.db} ({args.db_name})")
    print(f"images   : {len(entries)} from {args.images}")
    print(f"hint     : truth centre + {args.offset_hint} field widths, -r {args.radius}")
    print()

    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = {ex.submit(run_one, e, args): e for e in entries}
        for fut in concurrent.futures.as_completed(futs):
            r = fut.result()
            results.append(r)
            def astap_tag(r):
                if not args.astap:
                    return ""
                st = r.get("astap_status", "")
                if st == "OK":
                    return f"   | ASTAP OK  centre={r['astap_centre']:.3f}\""
                if st == "WRONG":
                    return f"   | ASTAP WRONG corner={r['astap_corner']:.1f}\""
                return f"   | ASTAP {st or '-'}"

            if r["status"] == "OK" and r["err_corner"] > args.max_corner_err:
                r["status"] = "WRONG"
                print(f"  {r['id']:<16} {r['tier']}  WRONG !! "
                      f"centre={r['err_centre']:>8.3f}\"  corner={r['err_corner']:>9.3f}\"  "
                      f"rot={r['rot_err_deg']:>7.3f}d  <-- FALSE POSITIVE" + astap_tag(r))
            elif r["status"] == "OK":
                print(f"  {r['id']:<16} {r['tier']}  OK       "
                      f"centre={r['err_centre']:>8.3f}\"  corner={r['err_corner']:>9.3f}\"  "
                      f"scale={r['scale_err_pct']:>8.4f}%  {r['secs']:>6.2f}s" + astap_tag(r))
            else:
                print(f"  {r['id']:<16} {r['tier']}  {r['status']:<8} "
                      f"{r['note']:<24} {r['secs']:>6.2f}s" + astap_tag(r))

    results.sort(key=lambda r: (r["tier"], r["id"]))

    if args.csv:
        cols = ["id", "tier", "status", "fov_deg", "pixscale_as", "err_centre",
                "err_corner", "scale_err_pct", "rot_err_deg", "nstars", "nquads",
                "secs", "note", "astap_status", "astap_centre", "astap_corner",
                "astap_scale_err", "astap_secs"]
        with open(args.csv, "w") as f:
            f.write(",".join(cols) + "\n")
            for r in results:
                f.write(",".join("" if r.get(c) is None else str(r.get(c, "")) for c in cols) + "\n")
        print(f"\nCSV: {args.csv}")

    def summarise(tier_filter, label):
        rs = [r for r in results if r["tier"] in tier_filter]
        if not rs:
            return
        ok = [r for r in rs if r["status"] == "OK"]
        wrong = [r for r in rs if r["status"] == "WRONG"]
        print(f"\n  {label}: {len(ok)}/{len(rs)} correct"
              f"   ({len(wrong)} FALSE POSITIVE{'S' if len(wrong) != 1 else ''}"
              f", {len(rs) - len(ok) - len(wrong)} no-solve)")
        for r in wrong:
            print(f"    !! {r['id']:<16} corner={r['err_corner']:.1f}\"  "
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
    tierd = [r for r in results if r["tier"] == "D"]
    if tierd:
        solved = [r for r in tierd if r["status"] in ("OK", "WRONG")]
        print(f"\n  Tier D (negative/stress): {len(solved)}/{len(tierd)} returned a solution")
        for r in tierd:
            bad = r["status"] in ("OK", "WRONG")
            mark = "FALSE POSITIVE" if bad else "correctly failed"
            extra = f" centre={r['err_centre']}\"" if bad else ""
            print(f"    {r['id']:<18} {r['status']:<9} {mark}{extra}")
    if args.astap:
        pv = [r for r in results if r["tier"] in ("A", "B")]
        p_ok = [r for r in pv if r["status"] == "OK"]
        p_wr = [r for r in pv if r["status"] == "WRONG"]
        a_ok = [r for r in pv if r.get("astap_status") == "OK"]
        a_wr = [r for r in pv if r.get("astap_status") == "WRONG"]
        both = [r for r in pv if r["status"] == "OK" and r.get("astap_status") == "OK"]
        only_p = [r for r in pv if r["status"] == "OK" and r.get("astap_status") != "OK"]
        only_a = [r for r in pv if r["status"] != "OK" and r.get("astap_status") == "OK"]
        neither = [r for r in pv if r["status"] != "OK" and r.get("astap_status") != "OK"]

        print("\n" + "=" * 78)
        print(" HEAD TO HEAD vs ASTAP  (tiers A+B, n=%d)" % len(pv))
        print("=" * 78)
        print(f"  {'':<22}{'arcsec':>12}{'ASTAP':>12}")
        print(f"  {'correct':<22}{len(p_ok):>12}{len(a_ok):>12}")
        print(f"  {'false positives':<22}{len(p_wr):>12}{len(a_wr):>12}")
        print(f"  {'no solve':<22}{len(pv)-len(p_ok)-len(p_wr):>12}{len(pv)-len(a_ok)-len(a_wr):>12}")

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
        if only_p:
            print("  arcsec solved, ASTAP did not:")
            for r in only_p:
                print(f"    + {r['id']:<16} {r.get('astap_status','-')}")
        if only_a:
            print("  ASTAP solved, arcsec did not:")
            for r in only_a:
                print(f"    - {r['id']:<16} arcsec={r['status']}")
        print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main())
