#!/usr/bin/env python3
"""Build scripts/corpus.tsv: the v1 manifest plus the v2 expansion.

This is the provenance of the corpus, kept so the selection can be audited and
re-run - but the committed corpus.tsv is the authoritative artefact. Several steps
query live archives (MocServer coverage, IRSA ZTF/WISE metadata, the LCO archive,
the TESS S3 listing, SkyServer, SkyMapper SIAP), whose answers change over time, so
a re-run with the same seed will not reproduce the file byte for byte.

    scripts/corpus-select.py > scripts/corpus.tsv          # queries live services
    scripts/corpus-select.py --cache /tmp/corpus-q.json    # reuse answers already fetched

Standard library only. Queries are sequential and throttled.
"""

import argparse
import json
import math
import os
import random
import re
import sys
import time
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
UA = "arcsec-corpus-select/1.0 (+https://github.com/cruzzil/arcsec; benchmark corpus, low rate)"
SEED = 20260929
CACHE = {}
CACHE_PATH = None


def get(url, timeout=120, tries=4):
    if url in CACHE:
        return CACHE[url]
    delay = 5
    for attempt in range(tries):
        try:
            time.sleep(0.4)
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                body = r.read().decode("utf-8", "replace")
            CACHE[url] = body
            if CACHE_PATH and len(CACHE) % 10 == 0:
                save_cache()
            return body
        except Exception as ex:  # noqa: BLE001
            code = getattr(ex, "code", None)
            if attempt + 1 == tries or (code and code < 500 and code != 429):
                print(f"# query failed: {url[:120]} ({ex})", file=sys.stderr)
                return ""
            time.sleep(delay)
            delay *= 3
    return ""


def save_cache():
    if CACHE_PATH:
        with open(CACHE_PATH, "w") as f:
            json.dump(CACHE, f)


# ── geometry ────────────────────────────────────────────────────────────────


def gal_b(ra, dec):
    """Galactic latitude (deg) of an ICRS position."""
    ra, dec = math.radians(ra), math.radians(dec)
    ra_gp, dec_gp = math.radians(192.85948), math.radians(27.12825)
    sb = (math.sin(dec) * math.sin(dec_gp)
          + math.cos(dec) * math.cos(dec_gp) * math.cos(ra - ra_gp))
    return math.degrees(math.asin(max(-1, min(1, sb))))


def rand_sky(rng):
    return rng.uniform(0, 360), math.degrees(math.asin(rng.uniform(-1, 1)))


# ── HiPS surveys: id -> (dataset label, native-ish output scale range "/px) ──

HIPS = {
    "CDS/P/DSS2/red": ("DSS2-red", 1.0, 1.7),
    "CDS/P/DSS2/blue": ("DSS2-blue", 1.0, 1.7),
    "CDS/P/DSS2/NIR": ("DSS2-IR", 1.0, 1.7),
    "CDS/P/2MASS/J": ("2MASS-J", 1.0, 2.0),
    "CDS/P/2MASS/H": ("2MASS-H", 1.0, 2.0),
    "CDS/P/2MASS/K": ("2MASS-K", 1.0, 2.0),
    "CDS/P/PanSTARRS/DR1/g": ("PS1-g", 0.5, 1.5),
    "CDS/P/PanSTARRS/DR1/r": ("PS1-r", 0.5, 1.5),
    "CDS/P/PanSTARRS/DR1/i": ("PS1-i", 0.5, 1.5),
    "CDS/P/PanSTARRS/DR1/z": ("PS1-z", 0.5, 1.5),
    "CDS/P/Skymapper/DR4/g": ("SkyMapper-g", 0.8, 1.6),
    "CDS/P/Skymapper/DR4/r": ("SkyMapper-r", 0.8, 1.6),
    "CDS/P/Skymapper/DR4/i": ("SkyMapper-i", 0.8, 1.6),
    "CDS/P/DES-DR2/r": ("DES-r", 0.5, 1.5),
    "CDS/P/DES-DR2/i": ("DES-i", 0.5, 1.5),
    "CDS/P/DESI-Legacy-Surveys/DR10/r": ("LS-DR10-r", 0.5, 1.5),
    "CDS/P/DESI-Legacy-Surveys/DR10/z": ("LS-DR10-z", 0.5, 1.5),
    "CDS/P/unWISE/W1": ("unWISE-W1", 2.75, 2.75),
    "CDS/P/allWISE/W1": ("allWISE-W1", 1.4, 2.75),
    "CDS/P/ZTF/DR7/r": ("ZTF-DR7-r", 1.0, 2.0),
    "CDS/P/ZTF/DR7/g": ("ZTF-DR7-g", 1.0, 2.0),
    "CDS/P/VPHAS/DR4/r": ("VPHAS-r", 0.5, 1.5),
    "CDS/P/IPHAS/DR2/r": ("IPHAS-r", 0.5, 1.5),
    "CDS/P/DECaPS/DR2/r": ("DECaPS-r", 0.5, 1.5),
    "wfau.roe.ac.uk/P/VISTA/VHS/J": ("VHS-J", 1.0, 1.5),
    "wfau.roe.ac.uk/P/VISTA/VHS/K": ("VHS-K", 1.0, 1.5),
    "wfau.roe.ac.uk/P/UKIDSS/LAS/K": ("UKIDSS-LAS-K", 1.0, 1.5),
    "CDS/P/DENIS/I": ("DENIS-I", 1.0, 1.5),
    "CDS/P/GALEXGR6_7/NUV": ("GALEX-NUV", 1.5, 3.0),
    "CDS/P/TESS/2yr": ("TESS-HiPS", 13.0, 25.0),
    "CDS/P/SHASSA/DU": ("SHASSA-cont", 26.0, 50.0),
}
# weights for random draws: all-sky workhorses less often, so the rarer surveys show up
WEIGHT = {"DSS2": 1.0, "2MASS": 0.8, "GALEX": 0.3, "TESS": 0.0, "SHASSA": 0.0}

MAXPIX = 3000


def covering(ra, dec, sr):
    """HiPS ids (of those above) whose coverage encloses a disc."""
    body = get("https://alasky.cds.unistra.fr/MocServer/query?" + urllib.parse.urlencode(
        {"RA": f"{ra:.4f}", "DEC": f"{dec:.4f}", "SR": f"{sr:.3f}", "intersect": "enclosed",
         "ID": "*", "dataproduct_type": "image", "get": "id"}))
    ids = set(body.split())
    return [h for h in HIPS if h in ids]


def hips_row(iid, tier, ra, dec, fov, hips, rng, aspect=1.0, rot=None, sets="v2",
             scale=None):
    label, lo, hi = HIPS[hips]
    s = scale if scale else rng.uniform(lo, hi)
    npx = fov * 3600.0 / s
    if npx > MAXPIX:
        npx = MAXPIX
    w = int(round(npx))
    h = int(round(npx / aspect))
    extra = f"hips={hips}"
    if rot is not None:
        extra += f";rot={rot:.1f}"
    return [iid, tier, f"{ra:.4f}", f"{dec:.4f}", f"{fov:g}", str(w), str(h), "hips2fits",
            extra, label, sets, "exact"]


# ── v1 conversion ───────────────────────────────────────────────────────────


def v1_rows():
    rows = []
    for line in open(os.path.join(HERE, "test-images.tsv")):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        p = re.split(r"[ \t]+", s)
        src, extra = p[7], p[8]
        if src == "hips2fits":
            ds = HIPS.get(extra, (extra.replace("CDS/P/", "").replace("/", "-"),))[0]
            if extra == "CDS/P/SDSS9/color":
                ds = "SDSS9-color"
            if extra == "CDS/P/DSS2/color":
                ds = "DSS2-color"
            truth = "exact"
        elif src == "skyview":
            ds, truth = f"SkyView-{extra}", "resampled"
        elif src == "legacysurvey":
            ds, truth = f"LS-DR10-{extra}", "resampled"
        elif src == "panstarrs":
            ds, truth = f"PS1-{extra}", "header"
        elif src == "sdss":
            ds, truth = f"SDSS-{extra.split('/')[3]}", "header"
        else:
            ds, truth = src, "header"
        rows.append(p[:9] + [ds, "v1", truth])
    return rows


# ── v2 sections ─────────────────────────────────────────────────────────────


def sec_random_hips(rng, n=80):
    rows = []
    wsum_ids = list(HIPS)
    k = 0
    tries = 0
    while k < n and tries < 400:
        tries += 1
        ra, dec = rand_sky(rng)
        fov = math.exp(rng.uniform(math.log(0.2), math.log(2.2)))
        cov = [h for h in covering(ra, dec, fov * 0.75) if h in wsum_ids]
        cov = [h for h in cov if not any(t in h for t in ("TESS", "SHASSA"))]
        if not cov:
            continue
        weights = [next((v for t, v in WEIGHT.items() if t in h), 2.0) for h in cov]
        hips = rng.choices(cov, weights)[0]
        aspect = rng.choice([1.0, 1.0, 1.5, 4 / 3.0])
        rot = rng.uniform(0, 360) if rng.random() < 0.4 else None
        b = gal_b(ra, dec)
        tag = "plane" if abs(b) < 10 else ("highlat" if abs(b) > 50 else "midlat")
        k += 1
        rows.append(hips_row(f"rnd_{k:03d}", "A", ra, dec, fov, hips, rng, aspect, rot,
                             sets=f"v2,random,{tag}"))
    return rows


def sec_wide(rng):
    rows = []
    fovs = [3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 24, 28]
    k = 0
    for fov in fovs:
        for _ in range(40):
            ra, dec = rand_sky(rng)
            if "CDS/P/TESS/2yr" in covering(ra, dec, fov * 0.72):
                break
        else:
            continue
        k += 1
        scale = max(13.0, fov * 3600 / MAXPIX)
        aspect = rng.choice([1.0, 1.5])
        rows.append(hips_row(f"wide_tess_{k:02d}", "A", ra, dec, fov, "CDS/P/TESS/2yr", rng,
                             aspect, rng.uniform(0, 360), sets="v2,wide", scale=scale))
    k = 0
    for fov in [15, 20, 25, 30, 35, 40, 45, 50]:
        for _ in range(60):
            ra, dec = rand_sky(rng)
            # SHASSA covers dec < +15; a MOC "enclosed" query at this radius never
            # succeeds (the survey has small gaps), so test the declination limits only
            if dec + 0.6 * fov < 12 and dec - 0.6 * fov > -88:
                break
        else:
            continue
        k += 1
        scale = max(26.0, fov * 3600 / 2500)
        rows.append(hips_row(f"wide_shassa_{k:02d}", "A", ra, dec, fov, "CDS/P/SHASSA/DU", rng,
                             1.5, None, sets="v2,wide", scale=scale))
    for k, fov in enumerate([4, 5, 6, 8, 10, 12], 1):
        ra, dec = rand_sky(rng)
        rows.append(hips_row(f"wide_dss_{k:02d}", "A", ra, dec, fov, "CDS/P/DSS2/red", rng,
                             1.5, None, sets="v2,wide,coarse", scale=fov * 3600 / MAXPIX))
    return rows


NAMED = [
    # id, ra, dec, fov, hips, note
    ("obj_betelgeuse", 88.793, 7.407, 1.0, "CDS/P/DSS2/red"),
    ("obj_rigel", 78.634, -8.202, 1.0, "CDS/P/DSS2/red"),
    ("obj_canopus", 95.988, -52.696, 1.5, "CDS/P/DSS2/red"),
    ("obj_arcturus", 213.915, 19.182, 1.0, "CDS/P/PanSTARRS/DR1/r"),
    ("obj_antares", 247.352, -26.432, 1.0, "CDS/P/DSS2/red"),
    ("obj_polaris", 37.955, 89.264, 1.0, "CDS/P/DSS2/red"),
    ("obj_rosette", 97.98, 4.94, 1.5, "CDS/P/DSS2/red"),
    ("obj_m16", 274.70, -13.81, 0.8, "CDS/P/DSS2/red"),
    ("obj_m20", 270.62, -23.03, 0.6, "CDS/P/DSS2/red"),
    ("obj_veil_e", 313.0, 31.2, 1.2, "CDS/P/DSS2/red"),
    ("obj_ic1396", 324.7, 57.5, 2.0, "CDS/P/DSS2/red"),
    ("obj_heart", 38.2, 61.45, 1.5, "CDS/P/DSS2/red"),
    ("obj_m33", 23.462, 30.660, 1.0, "CDS/P/DSS2/blue"),
    ("obj_m51", 202.470, 47.195, 0.5, "CDS/P/PanSTARRS/DR1/g"),
    ("obj_m104", 189.998, -11.623, 0.4, "CDS/P/DSS2/red"),
    ("obj_ngc253", 11.888, -25.288, 0.7, "CDS/P/DSS2/blue"),
    ("obj_cena", 201.365, -43.019, 0.6, "CDS/P/DSS2/red"),
    ("obj_lmc_bar", 80.894, -69.756, 1.5, "CDS/P/DSS2/red"),
    ("obj_smc", 13.187, -72.829, 1.5, "CDS/P/DSS2/red"),
    ("obj_m11", 282.77, -6.27, 0.5, "CDS/P/DSS2/red"),
    ("obj_47tuc", 6.024, -72.081, 0.6, "CDS/P/DSS2/red"),
    ("obj_ngc3532", 166.41, -58.75, 1.5, "CDS/P/DSS2/red"),
    ("obj_m7", 268.46, -34.79, 1.5, "CDS/P/DSS2/red"),
    ("obj_m35", 92.27, 24.33, 0.8, "CDS/P/DSS2/blue"),
    ("obj_b68", 261.50, -23.83, 0.3, "CDS/P/DSS2/red"),
    ("obj_sgra", 266.417, -29.008, 0.5, "CDS/P/2MASS/K"),
    ("obj_anticentre", 86.4, 28.9, 1.0, "CDS/P/DSS2/red"),
    ("obj_coalsack", 192.5, -62.8, 1.5, "CDS/P/DSS2/red"),
    ("obj_hyades_core", 66.75, 15.87, 2.0, "CDS/P/DSS2/blue"),
    ("obj_orion_belt", 83.0, -1.0, 2.0, "CDS/P/DSS2/red"),
    ("obj_ophiuchus_rho", 246.4, -24.4, 1.5, "CDS/P/DSS2/red"),
    ("obj_m27", 299.90, 22.72, 0.3, "CDS/P/PanSTARRS/DR1/r"),
    ("obj_m57", 283.40, 33.03, 0.25, "CDS/P/PanSTARRS/DR1/r"),
    ("obj_m1", 83.633, 22.014, 0.3, "CDS/P/PanSTARRS/DR1/r"),
    ("obj_sirius_2mass", 101.287, -16.716, 1.0, "CDS/P/2MASS/J"),
    ("obj_m42_wise", 83.82, -5.39, 1.0, "CDS/P/allWISE/W1"),
]


def sec_named(rng):
    rows = []
    for iid, ra, dec, fov, hips in NAMED:
        aspect = rng.choice([1.0, 1.5])
        rows.append(hips_row(iid, "A", ra, dec, fov, hips, rng, aspect, None, sets="v2,named"))
    return rows


def sec_survey_ladder(rng):
    """One southern field through every survey that covers it."""
    ra, dec = 60.0, -45.0
    rows = []
    for hips in ["CDS/P/DSS2/red", "CDS/P/DSS2/blue", "CDS/P/2MASS/K", "CDS/P/DES-DR2/r",
                 "CDS/P/Skymapper/DR4/r", "wfau.roe.ac.uk/P/VISTA/VHS/J", "CDS/P/allWISE/W1",
                 "CDS/P/unWISE/W1", "CDS/P/GALEXGR6_7/NUV", "CDS/P/DENIS/I",
                 "CDS/P/DESI-Legacy-Surveys/DR10/z"]:
        lab = HIPS[hips][0]
        rows.append(hips_row(f"surv_s_{lab.lower().replace('-', '_')}", "A", ra, dec, 0.8, hips,
                             rng, 1.0, None, sets="v2,survey"))
    return rows


# camera-like parents for tier S: (id, hips, fov, width, height)
PARENTS = [
    ("cam_ps1_a", "CDS/P/PanSTARRS/DR1/r", 1.0, 2400, 1600),
    ("cam_ps1_b", "CDS/P/PanSTARRS/DR1/i", 0.8, 2400, 1600),
    ("cam_dss_a", "CDS/P/DSS2/red", 1.2, 2400, 1600),
    ("cam_dss_b", "CDS/P/DSS2/red", 1.5, 2400, 1600),
    ("cam_dss_c", "CDS/P/DSS2/blue", 1.0, 2000, 1500),
    ("cam_2mass_a", "CDS/P/2MASS/J", 1.0, 2400, 1600),
    ("cam_sm_a", "CDS/P/Skymapper/DR4/r", 1.0, 2400, 1600),
    ("cam_des_a", "CDS/P/DES-DR2/r", 0.7, 2400, 1600),
    ("cam_ztf_a", "CDS/P/ZTF/DR7/r", 1.2, 2400, 1600),
    ("cam_dss_d", "CDS/P/DSS2/red", 2.0, 2400, 1600),
    ("cam_dss_e", "CDS/P/DSS2/red", 0.6, 1600, 1200),
    ("cam_ls_a", "CDS/P/DESI-Legacy-Surveys/DR10/r", 0.6, 2000, 1500),
    ("cam_tess_a", "CDS/P/TESS/2yr", 6.0, 1500, 1000),
    ("cam_tess_b", "CDS/P/TESS/2yr", 9.0, 1500, 1000),
    ("cam_tess_c", "CDS/P/TESS/2yr", 4.0, 1200, 800),
    ("cam_2mass_b", "CDS/P/2MASS/H", 1.5, 2400, 1600),
]

RECIPES = [
    ("amateur", "vignette:0.35,gradient:0.15:30,noise:1.5,hot:400,hotcol:2,u16", "amateur"),
    ("lens", "distort:{K},vignette:0.5,noise:1,u16", "distort"),
    ("pincush", "distort:-{K},noise:1,u16", "distort"),
    ("osc", "bayer:RGGB,noise:1,hot:200,u16", "osc"),
    ("oscgbrg", "bayer:GBRG,vignette:0.3,u16", "osc"),
    ("trail", "trail:6:35,noise:1,u16", "tracking"),
    ("trail12", "trail:12:120,u16", "tracking"),
    ("defocus", "blur:2.5,noise:1,u16", "focus"),
    ("clouds", "clouds:0.6:0.4,noise:1,u16", "clouds"),
    ("sat", "saturate:99.0,u16", "saturation"),
    ("junk", "satellite:3,glow:0.5,noise:2,hot:800,u16", "artefacts"),
    ("flipx", "flipx,u16", "orient"),
    ("flipy", "flipy,f32", "orient"),
    ("rot90", "rot90:1,u16", "orient"),
    ("rot270", "rot90:3,u16", "orient"),
    ("transp", "transpose,u16", "orient"),
    ("u8", "noise:0.5,u8", "format"),
    ("i32", "i32", "format"),
    ("f64", "f64", "format"),
    ("gz", "u16,gzip", "format"),
    ("xisf", "u16,xisf", "format"),
    ("keepwcs", "u16,keepwcs", "format"),
    ("bin2", "bin:2,noise:0.5,u16", "binning"),
    ("worst", "distort:{K},vignette:0.5,gradient:0.3:200,clouds:0.4:0.5,trail:4:60,noise:2,hot:500,bayer:RGGB,u16", "amateur"),
]


def sec_parents_and_synth(rng):
    rows = []
    parents = []
    for iid, hips, fov, w, h in PARENTS:
        for _ in range(60):
            ra, dec = rand_sky(rng)
            if hips in covering(ra, dec, fov * 0.75) and abs(gal_b(ra, dec)) > 5:
                break
        lab = HIPS[hips][0]
        s = fov * 3600.0 / max(w, h)
        rows.append([iid, "A", f"{ra:.4f}", f"{dec:.4f}", f"{fov:g}", str(w), str(h),
                     "hips2fits", f"hips={hips}", lab, "v2,parent", "exact"])
        parents.append((iid, ra, dec, fov, w, h, lab, s))
    # every recipe at least once, then fill to ~70 with random pairings
    pairs = []
    for i, rec in enumerate(RECIPES):
        pairs.append((parents[i % len(parents)], rec))
    while len(pairs) < 72:
        pairs.append((rng.choice(parents), rng.choice(RECIPES)))
    seen = set()
    k = 0
    for (pid, ra, dec, fov, w, h, lab, s), (rname, ops, tag) in pairs:
        if (pid, rname) in seen:
            continue
        seen.add((pid, rname))
        k += 1
        half_diag = 0.5 * math.hypot(w, h)
        disp = rng.uniform(4, 25)  # px of radial displacement at the corner
        kval = disp / half_diag ** 3
        o = ops.replace("{K}", f"{kval:.3e}")
        ww, hh = w, h
        if "rot90:1" in o or "rot90:3" in o or "transpose" in o:
            ww, hh = h, w
        if "bin:2" in o:
            ww, hh = ww // 2, hh // 2
        rows.append([f"s_{pid[4:]}_{rname}", "S", f"{ra:.4f}", f"{dec:.4f}", f"{fov:g}", str(ww),
                     str(hh), "synth", f"parent={pid};ops={o};seed={SEED % 100000 + k}",
                     f"synth-{lab}", f"v2,synth,{tag}", "derived"])
    return rows


def sec_skyview(rng, n=30):
    surveys = [("DSS1R", "all"), ("DSS1B", "all"), ("DSS2R", "all"), ("DSS2B", "all"),
               ("DSS2IR", "all"), ("2MASS-J", "all"), ("2MASS-H", "all"), ("2MASS-K", "all"),
               ("WISE 3.4", "all"), ("WISE 4.6", "all"), ("GALEX Near UV", "all")]
    rows = []
    for k in range(1, n + 1):
        sv, _ = surveys[(k - 1) % len(surveys)]
        ra, dec = rand_sky(rng)
        fov = rng.choice([0.3, 0.5, 0.7, 1.0, 1.2])
        native = {"DSS1R": 1.7, "DSS1B": 1.7, "WISE 3.4": 1.375, "WISE 4.6": 1.375,
                  "GALEX Near UV": 1.5}.get(sv, 1.0)
        npx = min(2000, int(fov * 3600 / native))
        rot = rng.uniform(0, 360) if rng.random() < 0.3 else None
        extra = "survey=" + sv.replace(" ", "%20")
        if rot is not None:
            extra += f";rot={rot:.1f}"
        rows.append([f"sv2_{k:02d}", "B", f"{ra:.4f}", f"{dec:.4f}", f"{fov:g}", str(npx),
                     str(npx), "skyview", extra, "SkyView-" + sv.replace(" ", ""), "v2,random",
                     "resampled"])
    return rows


def sec_legacy(rng, n=30):
    rows = []
    k = 0
    while k < n:
        ra, dec = rand_sky(rng)
        if dec < -60 or dec > 80:
            continue
        if "CDS/P/DESI-Legacy-Surveys/DR10/r" not in covering(ra, dec, 0.4):
            continue
        k += 1
        ps = rng.choice([0.262, 0.262, 0.4, 0.5, 0.75, 1.0])
        size = rng.choice([1200, 1800, 2400, 3000]) if ps < 0.5 else rng.choice([1200, 1800, 2400])
        band = rng.choice(["g", "r", "r", "z", "i"])
        fov = ps * size / 3600.0
        rows.append([f"ls2_{k:02d}", "B", f"{ra:.4f}", f"{dec:.4f}", f"{fov:.4f}", str(size),
                     str(size), "legacysurvey", f"band={band};layer=ls-dr10", f"LS-DR10-{band}",
                     "v2,random", "resampled"])
    return rows


def sec_ps1(rng, n=20):
    rows = []
    k = 0
    while k < n:
        ra, dec = rand_sky(rng)
        if dec < -29:
            continue
        k += 1
        # native 0.25"/px pixels. fitscut cuts one skycell (~0.4 deg), so larger
        # cutouts at random positions mostly run off its edge into NaN padding
        size = rng.choice([2400, 2400, 2800])
        out = rng.choice([1200, 1600, size])
        band = rng.choice(["g", "r", "i", "z", "y"])
        fov = size * 0.25 / 3600
        rows.append([f"ps1v2_{k:02d}", "B", f"{ra:.4f}", f"{dec:.4f}", f"{fov:.4f}", str(out),
                     str(out), "panstarrs", f"band={band};size={size};out={out}", f"PS1-{band}",
                     "v2,random", "header"])
    return rows


def sec_sdss(rng, n=25):
    sql = ("select top 400 run,rerun,camcol,field,ra,dec from Field "
           "where quality=3 and (field % 37) = 5 order by run,camcol,field")
    body = get("https://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch?"
               + urllib.parse.urlencode({"cmd": sql, "format": "csv"}))
    fields = []
    for line in body.splitlines()[2:]:
        p = line.split(",")
        if len(p) >= 6:
            fields.append(p)
    rng.shuffle(fields)
    rows = []
    for k, (run, rerun, camcol, field, ra, dec) in enumerate(fields[:n], 1):
        band = rng.choice(["g", "r", "r", "i", "z"])
        rows.append([f"sdss2_{k:02d}", "B", f"{float(ra):.4f}", f"{float(dec):.4f}", "0.225",
                     "2048", "1489", "sdss", f"frame={rerun}/{run}/{camcol}/{band}/{field}",
                     f"SDSS-{band}", "v2,random", "header"])
    return rows


def sec_ztf(rng, n=50):
    rows = []
    k = 0
    tries = 0
    while k < n and tries < 200:
        tries += 1
        ra, dec = rand_sky(rng)
        if dec < -27:
            continue
        filt = rng.choice(["zg", "zr", "zr", "zi"])
        seeing_hi = rng.choice([2.2, 2.2, 3.0, 4.5])
        where = f"filtercode='{filt}' AND infobits=0 AND seeing<{seeing_hi}"
        body = get("https://irsa.ipac.caltech.edu/ibe/search/ztf/products/sci?" + urllib.parse.urlencode(
            {"POS": f"{ra:.4f},{dec:.4f}", "ct": "csv", "WHERE": where,
             "columns": "filefracday,field,ccdid,qid,filtercode,imgtypecode,seeing,ra,dec"}))
        lines = [ln.split(",") for ln in body.splitlines()[1:] if ln.strip()]
        if not lines:
            continue
        p = rng.choice(lines[:200])
        _, _, _, ffd, field, ccd, qid, fc, itc, seeing, cra, cdec = p[:12]
        fn = f"ztf_{ffd}_{int(field):06d}_{fc}_c{int(ccd):02d}_{itc}_q{qid}_sciimg.fits"
        size = rng.choice([1500, 2000, 2000, 2500, 3000])
        # cutout centred on the quadrant centre so it stays inside the CCD
        k += 1
        b = gal_b(float(cra), float(cdec))
        tag = "plane" if abs(b) < 10 else "offplane"
        rows.append([f"ztf_{k:02d}", "B", f"{float(cra):.4f}", f"{float(cdec):.4f}",
                     f"{size * 1.012 / 3600:.4f}", str(size), str(size), "ztf",
                     f"product={fn};size={size};seeing={float(seeing):.2f}", f"ZTF-{fc[1]}",
                     f"v2,random,{tag}", "header-tpv"])
    return rows


def sec_wise(rng, n=25):
    rows = []
    k = 0
    tries = 0
    while k < n and tries < 150:
        tries += 1
        ra, dec = rand_sky(rng)
        if k < 3:  # make sure the poles are in
            ra, dec = [(0.0, 89.0), (180.0, -89.0), (270.0, 66.5)][k]
        body = get("https://irsa.ipac.caltech.edu/ibe/search/wise/allsky/4band_p1bm_frm?" + urllib.parse.urlencode(
            {"POS": f"{ra:.4f},{dec:.4f}", "ct": "csv", "WHERE": "band=1 AND qual_frame=10",
             "columns": "scan_id,frame_num,crval1,crval2,band"}))
        lines = [ln.split(",") for ln in body.splitlines()[1:] if ln.strip()]
        if not lines:
            continue
        p = rng.choice(lines[:100])
        cols = body.splitlines()[0].split(",")
        d = dict(zip(cols, p))
        band = "2" if rng.random() < 0.25 else "1"
        k += 1
        fr = f"{d['scan_id']}{int(d['frame_num']):03d}"
        rows.append([f"wise_{k:02d}", "B", f"{float(d['crval1']):.4f}", f"{float(d['crval2']):.4f}",
                     "0.776", "1016", "1016", "wise", f"frame={fr};band={band}", f"WISE-L1b-W{band}",
                     "v2,random", "header-sip"])
    return rows


def s3list(prefix, delimiter=None, start_after=None, max_keys=1000):
    q = {"list-type": "2", "prefix": prefix, "max-keys": str(max_keys)}
    if delimiter:
        q["delimiter"] = delimiter
    if start_after:
        q["start-after"] = start_after
    body = get("https://stpubdata.s3.amazonaws.com/?" + urllib.parse.urlencode(q))
    prefixes = re.findall(r"<Prefix>([^<]*)</Prefix>", body)[1:]
    keys = re.findall(r"<Key>([^<]*)</Key>.*?<ETag>&quot;([0-9a-f]+)&quot;</ETag>.*?<Size>(\d+)</Size>",
                      body)
    return prefixes, keys


def sec_tess(rng, n_ffi=18):
    rows = []
    sectors = rng.sample(range(1, 70), n_ffi)
    k = 0
    for sec in sectors:
        pre = f"tess/public/ffi/s{sec:04d}/"
        years, _ = s3list(pre, "/")
        if not years:
            continue
        days, _ = s3list(years[0], "/")
        if not days:
            continue
        day = days[len(days) // 2]
        cam, ccd = rng.randint(1, 4), rng.randint(1, 4)
        cdir = f"{day}{cam}-{ccd}/"
        _, keys = s3list(cdir, None, None, 60)
        ffic = [(kk, et) for kk, et, sz in keys if kk.endswith("_ffic.fits") and "-" not in et]
        if not ffic:
            continue
        key, etag = ffic[len(ffic) // 2]
        rel = key[len("tess/public/ffi/"):]
        # crops: one big and one or two smaller, at varying offsets
        crops = [(2048, 44, 0)] if rng.random() < 0.4 else [(1024, 44 + rng.choice([0, 1024]), rng.choice([0, 1024]))]
        crops.append((512, 44 + rng.randint(0, 1536), rng.randint(0, 1536)))
        if rng.random() < 0.5:
            crops.append((384, 44 + rng.randint(0, 1664), rng.randint(0, 1664)))
        for size, x0, y0 in crops:
            k += 1
            fov = size * 21.0 / 3600
            rows.append([f"tess_{k:02d}", "B", "0", "0", f"{fov:.3f}", str(size), str(size), "tess",
                         f"key={rel};x0={x0};y0={y0};size={size};md5={etag}",
                         f"TESS-FFI-cam{cam}", "v2,random,wide" if size >= 1024 else "v2,random",
                         "header-sip"])
    return rows


def sec_lco(rng, n=45):
    """Public, pipeline-reduced (e91) science frames, spread over years and telescopes."""
    rows = []
    seen_targets = set()
    want = {"0m4": 24, "1m0": 16, "2m0": 5}
    got = {k: 0 for k in want}
    tries = 0
    while sum(got.values()) < sum(want.values()) and tries < 120:
        tries += 1
        cls = rng.choice([c for c in want if got[c] < want[c]] or list(want))
        year = rng.randint(2016, 2024)
        day = rng.randint(1, 360)
        start = time.strftime("%Y-%m-%d", time.gmtime(time.mktime((year, 1, 1, 0, 0, 0, 0, 0, 0)) + day * 86400))
        end = time.strftime("%Y-%m-%d", time.gmtime(time.mktime((year, 1, 1, 0, 0, 0, 0, 0, 0)) + (day + 3) * 86400))
        body = get("https://archive-api.lco.global/frames/?" + urllib.parse.urlencode(
            {"public": "true", "reduction_level": "91", "configuration_type": "EXPOSE",
             "limit": "100", "start": start, "end": end}))
        try:
            d = json.loads(body)
        except ValueError:
            continue
        cands = []
        for r in d.get("results", []):
            tel = r.get("telescope_id", "")
            if not tel.startswith(cls):
                continue
            if r.get("proposal_id") in ("LCOEngineering", "auto_focus", "calibrate") or \
                    r.get("target_name", "").lower() in ("auto_focus", ""):
                continue
            if (r.get("exposure_time") or 0) < 20:
                continue
            if not r.get("area"):
                continue
            key = (r.get("target_name"), r.get("instrument_id"))
            if key in seen_targets:
                continue
            cands.append(r)
        if not cands:
            continue
        r = rng.choice(cands)
        seen_targets.add((r.get("target_name"), r.get("instrument_id")))
        poly = r["area"]["coordinates"][0][:4]
        ras = [p[0] % 360 for p in poly]
        if max(ras) - min(ras) > 180:
            ras = [(x + 180) % 360 - 180 for x in ras]
        ra = (sum(ras) / 4) % 360
        dec = sum(p[1] for p in poly) / 4
        side = max(abs(poly[0][1] - poly[1][1]), abs(poly[1][1] - poly[2][1]),
                   abs(poly[0][1] - poly[2][1]))
        got[cls] += 1
        i = sum(got.values())
        md5 = (r.get("version_set") or [{}])[0].get("md5", "")
        rows.append([f"lco_{i:02d}", "C", f"{ra:.4f}", f"{dec:.4f}", f"{side:.3f}", "0", "0", "lco",
                     f"frame={r['id']};md5={md5};inst={r['instrument_id']};filter={r['primary_optical_element']}"
                     f";exp={r['exposure_time']:.0f}",
                     f"LCO-{cls}-{r['instrument_id'][:2]}", "v2,real", "pipeline"])
    return rows


def sec_skymapper(rng, n=15):
    rows = []
    k = 0
    tries = 0
    while k < n and tries < 80:
        tries += 1
        ra, dec = rand_sky(rng)
        if dec > 0:
            continue
        band = rng.choice(["g", "r", "i"])
        body = get("https://api.skymapper.nci.org.au/public/siap/dr4/query?" + urllib.parse.urlencode(
            {"POS": f"{ra:.4f},{dec:.4f}", "SIZE": "0.17", "BAND": band, "FORMAT": "image/fits",
             "INTERSECT": "covers", "RESPONSEFORMAT": "CSV"}))
        lines = [ln for ln in body.splitlines()[1:] if ln.strip()]
        if not lines:
            continue
        m = re.search(r"image=([0-9-]+)&", rng.choice(lines))
        if not m:
            continue
        k += 1
        rows.append([f"smss_{k:02d}", "B", f"{ra:.4f}", f"{dec:.4f}", "0.17", "1231", "1231",
                     "skymapper", f"image={m.group(1)};size=0.17", f"SkyMapper-native-{band}",
                     "v2,random,narrow", "header-tpv"])
    return rows


def sec_negatives(rng, pool):
    rows = []
    gens = [("neg_noise_a", "noise", 2000, 1500, 1.0), ("neg_noise_b", "noise", 4000, 3000, 2.0),
            ("neg_noise_c", "noise", 1000, 1000, 0.5), ("neg_noise_d", "noise", 3000, 2000, 5.0),
            ("neg_flat_a", "flat", 2400, 1600, 1.2), ("neg_flat_b", "flat", 1600, 1200, 8.0),
            ("neg_fake_100", "fakestars:100", 2400, 1600, 1.0),
            ("neg_fake_300", "fakestars:300", 2400, 1600, 1.5),
            ("neg_fake_600", "fakestars:600", 3000, 2000, 2.0),
            ("neg_fake_1000", "fakestars:1000", 3000, 2000, 1.0),
            ("neg_fake_wide", "fakestars:400", 2000, 1400, 10.0),
            ("neg_fake_narrow", "fakestars:150", 1600, 1200, 0.3)]
    for iid, g, w, h, fov in gens:
        ra, dec = rand_sky(rng)
        rows.append([iid, "D", f"{ra:.4f}", f"{dec:.4f}", f"{fov:g}", str(w), str(h), "synth",
                     f"gen={g};seed={rng.randint(1, 99999)}", "synthetic", "v2,negative", "none"])
    for i, pid in enumerate(["cam_ps1_a", "cam_dss_a", "cam_2mass_a", "cam_sm_a", "cam_dss_d",
                             "cam_tess_a"], 1):
        b = 64 if i % 2 else 128
        rows.append([f"neg_shuffle_{i}", "D", "0", "0", "0", "0", "0", "synth",
                     f"parent={pid};ops=shuffle:{b},u16;seed={rng.randint(1, 99999)}", "synthetic",
                     "v2,negative", "none"])
    # real, solvable fields with a hint 30-60 deg away and a radius that cannot reach
    for i, iid in enumerate(rng.sample(pool, 8), 1):
        dra, ddec = rng.uniform(30, 60) * rng.choice([-1, 1]), rng.uniform(-20, 20)
        rows.append([f"neg_hint_{i}", "D", "0", "0", "0", "0", "0", "alias",
                     f"file={iid};hint_dra={dra:.1f};hint_ddec={ddec:.1f};radius=10", "alias",
                     "v2,negative", "none"])
    # wrong FOV: a correct solve is fine, a wrong one is a false positive
    for i, (iid, sc) in enumerate(zip(rng.sample(pool, 4), [2.0, 0.5, 1.5, 0.67]), 1):
        rows.append([f"stress_fov_{i}", "D", "0", "0", "0", "0", "0", "alias",
                     f"file={iid};fov_scale={sc};expect=any", "alias", "v2,stress", "none"])
    # tiny crops of extended objects with few or no stars
    for iid, ra, dec, hips in [("neg_m17_core", 275.20, -16.17, "CDS/P/DSS2/red"),
                               ("neg_eta_car_core", 161.265, -59.685, "CDS/P/DSS2/red"),
                               ("neg_m31_ps1", 10.6847, 41.2690, "CDS/P/PanSTARRS/DR1/g"),
                               ("neg_orion_2massk", 83.818, -5.390, "CDS/P/2MASS/K")]:
        rows.append([iid, "D", f"{ra:.4f}", f"{dec:.4f}", "0.04", "400", "400", "hips2fits",
                     f"hips={hips}", HIPS[hips][0], "v2,negative", "exact"])
    return rows


HEADER = """# arcsec benchmark corpus manifest (v1 + v2) - see docs/test-images.md
#
# Generated by scripts/corpus-select.py; fetched by scripts/fetch-corpus.py; scored by
# scripts/benchmark.py --corpus. Whitespace-separated columns:
#
#   id tier ra dec fov_deg width height source extra dataset sets truth
#
# tier    A synthetic cutout (exact truth) | B survey pixels (header truth) |
#         C real observing frames (pipeline truth) | S simulated camera artefacts derived
#         from an A/B parent (exact derived truth) | D negative / stress controls
# extra   source-specific key=value;... (v1 rows keep their original bare token)
# sets    comma list: v1 = the original 103 (historical numbers), v2 = the expansion,
#         plus tags (random, named, wide, coarse, plane, parent, synth, negative, ...)
# truth   exact | resampled | header | header-sip | header-tpv | pipeline | derived | none
#
# For sources whose geometry is set by the archive (ztf, wise, tess, lco, sdss) ra/dec/
# fov/width/height are nominal; benchmark.py always takes truth from the delivered file.
"""


def main():
    global CACHE_PATH, CACHE
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", default=None, help="JSON cache of query answers")
    ap.add_argument("--only", default=None, help="comma list of sections to (re)build")
    args = ap.parse_args()
    if args.cache:
        CACHE_PATH = args.cache
        if os.path.exists(args.cache):
            CACHE = json.load(open(args.cache))

    sections = [
        ("random HiPS cutouts (tier A)", sec_random_hips),
        ("wide fields: TESS and SHASSA HiPS, coarse DSS (tier A)", sec_wide),
        ("named objects and hard field types (tier A)", sec_named),
        ("survey ladder: one southern field, every survey (tier A)", sec_survey_ladder),
        ("SkyView resampled surveys (tier B)", sec_skyview),
        ("DESI Legacy Surveys DR10 cutouts (tier B)", sec_legacy),
        ("Pan-STARRS1 stack cutouts, native skycell WCS (tier B)", sec_ps1),
        ("SDSS corrected frames, native (tier B)", sec_sdss),
        ("ZTF science-image cutouts, native pixels, TPV WCS (tier B)", sec_ztf),
        ("WISE L1b single exposures, native pixels, SIN-SIP WCS (tier B)", sec_wise),
        ("SkyMapper DR4 native cutouts, TPV WCS (tier B)", sec_skymapper),
        ("TESS full-frame-image crops, SIP WCS, 21\"/px (tier B)", sec_tess),
        ("Las Cumbres Observatory public frames, BANZAI pipeline WCS (tier C)", sec_lco),
        ("camera-like parents and simulated camera artefacts (tiers A, S)", sec_parents_and_synth),
    ]
    out = [HEADER]
    out.append("# " + "=" * 88)
    out.append("# v1 - the original 103 entries of scripts/test-images.tsv, unchanged")
    out.append("# " + "=" * 88)
    for r in v1_rows():
        out.append("  ".join(r))
    pool = []
    for i, (title, fn) in enumerate(sections):
        rng = random.Random(SEED + i)
        rows = fn(rng)
        out.append("#")
        out.append(f"# --- v2: {title} " + "-" * max(0, 80 - len(title)))
        for r in rows:
            out.append("  ".join(r))
            if r[1] in ("B",) and r[7] in ("ztf", "legacysurvey", "sdss"):
                pool.append(r[0])
        print(f"# {title}: {len(rows)}", file=sys.stderr)
        save_cache()
    rng = random.Random(SEED + 99)
    out.append("#")
    out.append("# --- v2: negative and stress controls (tier D) " + "-" * 40)
    for r in sec_negatives(rng, pool):
        out.append("  ".join(r))
    save_cache()
    print("\n".join(out))


if __name__ == "__main__":
    main()
