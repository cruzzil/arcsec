# Benchmark Image Corpus

A standing set of FITS images with trustworthy ground truth, used to measure arcsec's
solve rate, positional accuracy and false-positive rate across a realistic range of field
sizes, star densities and image quality.

Companion to [plate-solving.md](plate-solving.md). There are two manifests:

* **v1** — the original 103 images ([`scripts/test-images.tsv`](../scripts/test-images.tsv),
  §3.1–3.5, results in §6). Kept unchanged so historical numbers stay comparable.
* **the expanded corpus** — 635 entries from ten archives and two synthetic families
  ([`scripts/corpus.tsv`](../scripts/corpus.tsv), §2.8–2.10 and §3.6, results in §7). It
  contains v1 as the subset `v1`.

§9 runs both against ASTAP and [seiza](https://github.com/theatrus/seiza) as well.

```bash
# v1, as before
scripts/fetch-test-images.sh                  # all tiers → resources/testset/ (~2.2 GB)
scripts/benchmark.py --auto-db                # solve and score everything (§6.8)

# expanded corpus
scripts/fetch-corpus.py                       # → resources/corpus/ (~6.5 GB new, v1 hard-linked)
scripts/fetch-corpus.py --stats               # coverage tables, no download
scripts/benchmark.py --corpus --auto-db       # per-tier, per-source and per-FOV breakdowns
scripts/benchmark.py --corpus --auto-db --set v1   # the historical subset only
```

Images are **not** committed — `resources/` is gitignored — the manifests plus the fetch
scripts are the reproducible artefact.

---

## 1. What "ground truth" means here

Not all truth is equal, and the metric you can trust depends on where the image came from.
The corpus is stratified into four tiers.

| Tier | Ground truth | Truth accuracy | What it exercises |
|---|---|---|---|
| **A — Synthetic cutouts** | We *request* the centre, FOV and projection, so the true WCS is known exactly | exact (< 1 mas) | Correctness of the core solver across FOV and star density |
| **B — Survey pixels** | The archive's own pipeline WCS in the header (TAN, TAN-SIP, SIN-SIP or TPV) | 0.02–0.3″ (TESS ~1″) | Real PSFs, real noise, real artefacts, real distortion; sub-arcsecond accuracy |
| **C — Real observing frames** | The observatory pipeline's astrometric solution (LCO BANZAI, fitted to Gaia) | ~0.2–0.3″ | Individual reduced exposures from 0.4 m–2 m telescopes: seeing, tracking, nebulosity, sparse fields |
| **S — Simulated camera artefacts** | A tier-A/B parent's truth carried exactly through each transformation (§2.10) | as the parent | Vignetting, gradients, hot pixels, Bayer mosaics, trailing, defocus, clouds, lens distortion, flips, sample formats |
| **D — Negative / stress** | Known to be unsolvable, or known-hard | n/a | **False-positive rate** and graceful failure |

Tier D is the one most benchmark suites omit and the one that matters most for us: a solver
that always returns *something* is worse than one that says "no solution". See
[plate-solving.md §11.1](plate-solving.md#111-the-catalogue-path-performs-no-verification--fixed-2026-09-02).

---

## 2. Sources

Every URL below was fetched and the resulting FITS header inspected on 2026-08-31; the
WCS keywords quoted are what actually came back.

### 2.1 CDS HiPS2FITS — tier A (primary)

Renders any HiPS survey into a FITS cutout with a WCS we specify. This is what
`scripts/hips_solve_test.sh` already uses.

```
https://alasky.cds.unistra.fr/hips-image-services/hips2fits
    ?hips=CDS/P/DSS2/blue      # or CDS/P/DSS2/color, CDS/P/2MASS/J,
                               #    CDS/P/PanSTARRS/DR1/r, CDS/P/SDSS9/color
    &width=1000&height=1000
    &fov=1.5                   # degrees, largest dimension
    &projection=TAN&coordsys=icrs
    &ra=83.82&dec=-5.39
    &format=fits
```

Verified header: `CTYPE1 = 'RA---TAN'`, `CRVAL1 = 83.82`, `CRVAL2 = -5.39`,
`CDELT1 = -0.0015000856795217`, `CRPIX = 500.0`. The requested centre appears verbatim in
`CRVAL`, so positional truth is exact.

Free, no key, no rate-limit problems at this volume. Cutouts are capped at 50 Mpixel.

*Caveat 1*: HiPS tiles are resampled, so the PSF is not a real instrumental PSF and faint
stars are smoothed. Tier A tests the geometry, not the detector.

*Caveat 2 — sampling, and it is the important one.* A HiPS rendering must be requested at
close to the **native resolution of the underlying survey** (~1″/px for DSS2). Render the
same field coarser and the stars blend into a few pixels, centroids stop corresponding to
individual catalogue stars, and the solve fails for reasons that have nothing to do with
the solver. Measured on one field (RA 51.40, Dec +49.85, DSS2 blue, 1.5°):

| Requested pixels | Resulting scale | arcsec |
|---|---|---|
| 1000 × 1000 | 5.40″/px | no solve |
| 2000 × 2000 | 2.70″/px | no solve |
| 4000 × 4000 | **1.35″/px** | **solved** |

So the manifest sizes tier-A entries to hold roughly **1.25″/px**, not to a fixed pixel
count. Two consequences:

* Wide fields cannot be tested this way at all. HiPS2FITS caps cutouts at 50 Mpixel, so a
  5° field tops out near 3.6″/px and a 15° field near 10.8″/px — both far coarser than
  DSS2 can support. Entries above ~3° are therefore a **limitation of the corpus, not a
  measurement of arcsec**, and must be reported as such. A genuine wide-field test needs a
  real camera-lens frame (tier C), where the optics sample the stars properly even at
  8″/px.
* This is not an arcsec defect. A real telescope frame at 3″/px has a PSF sampled across
  2–3 pixels; a plate scan resampled to 3″/px has stars smeared below one pixel. The two
  look nothing alike to a centroiding detector.

### 2.2 NASA SkyView — tier A and B

Same idea, many more surveys, and it also returns the survey's own resampled data.

```
https://skyview.gsfc.nasa.gov/current/cgi/pskcall
    ?Survey=DSS2R              # DSS1B, DSS2R, DSS2B, DSS2IR, 2MASS-J/H/K,
                               #   SDSSg/r/i, WISE 3.4, GALEX Near UV, ...
    &Position=180.0,20.0
    &Size=0.5                  # degrees
    &Pixels=1000
    &Return=FITS
```

Verified header: `SIMPLE = T / Written by SkyView`, `CTYPE1 = 'RA---TAN'`,
`CRVAL1 = 180.0`, `CRVAL2 = 20.0`, `CDELT1 = -0.00050`, `RADESYS = 'FK5'`,
`SURVEY = 'DSS2RED'`.

Note `RADESYS = FK5` rather than ICRS — a <0.1″ difference, irrelevant at our tolerances
but worth knowing if we ever chase milliarcseconds.

### 2.3 DESI Legacy Imaging Surveys — tier B

Real ground-based CCD data (DECam/Mosaic/90Prime), astrometrically tied to Gaia.

```
https://www.legacysurvey.org/viewer/cutout.fits
    ?ra=180.0&dec=20.0
    &layer=ls-dr10
    &pixscale=1.0              # arcsec/pixel
    &size=512                  # pixels (max ~3000)
    &bands=r
```

Verified header: `SURVEY = 'LegacySurvey'`, `VERSION = 'DR10'`, `CTYPE1 = 'RA---TAN'`,
`CRVAL1 = 180.`, `CD1_1 = -0.000277777777777778` (= 1″/px), `CRPIX1 = 256.5`.

Coverage is the DESI footprint (roughly `-20° < dec < +80°` away from the galactic plane),
so pick fields accordingly. A clean, sky-subtracted, real-PSF test at any pixel scale you
ask for.

### 2.4 Pan-STARRS1 — tier B

Two-step: ask for the skycell filename, then cut it.

```
# 1. find the stack file
https://ps1images.stsci.edu/cgi-bin/ps1filenames.py?ra=180.0&dec=20.0&size=1024&format=fits&filters=r
# 2. cut it out
https://ps1images.stsci.edu/cgi-bin/fitscut.cgi
    ?red=<filename from step 1>&format=fits&size=1024&ra=180.0&dec=20.0
```

Verified header: `CTYPE1 = 'RA---TAN'`, `CRVAL1 = 179.999999999994`,
`CDELT1 = 6.94444461259988E-05` (0.25″/px), and notably
`CRPIX2 = 29325. / Reference pixel shifted for cutout`.

That large off-image `CRPIX` is a genuinely useful edge case: it exercises whether our WCS
comparison code handles a reference pixel far outside the array. Coverage is `dec > -30°`.

### 2.5 SDSS corrected frames — tier B

Full, uncropped survey frames — 2048 × 1489 at 0.396″/px, ~0.22° × 0.16°.

```
https://data.sdss.org/sas/dr17/eboss/photoObj/frames/301/2505/3/frame-r-002505-3-0038.fits.bz2
                                                    ^rerun ^run ^camcol ^band-run-camcol-field
```

Verified header: `RUN = 2505`, `CAMCOL = 3`, `CTYPE1 = 'RA---TAN'`,
`CRVAL1 = 276.945390655`, `CD1_1 = 1.52163755940E-08`, `CD1_2 = 0.000109988419300` —
note the near-90° rotation encoded in the CD matrix (SDSS scans along great circles). An
excellent rotation-handling test, and bz2-compressed so the fetch script must decompress.

### 2.6 MAST / HST — tier D (distortion stress)

```
https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/<rootname>_drc.fits
```

Verified reachable (HTTP 200), but a single ACS/WFC drizzled product is **295 MB** and the
field is only ~3.4′ — below what any ASTAP database supports (`d80`'s published floor is
0.15°). Include at most one, and expect it to fail with our current catalogue; it is here to
document the limit, not to pass. A purpose-built index would be needed. None is in the
manifest at present.

### 2.7 Real observing frames — tier C

* **Our own NGC 3372 set** — 240 frames, ASI533MC Pro + RedCat 51, ~3.13″/px, ~2.6° field.
  Not in the repository or the manifest; `scripts/bench_all.sh` expects the frames as
  `resources/*.fits`. This is the most valuable tier-C
  data we have: a single night, one target, so it isolates frame-to-frame consistency.
  Pick ~10 frames spanning the session (including frame 0200, which ASTAP failed to solve
  and arcsec solved — see `ARCSEC_VS_ASTAP.md`).
* **Free practice datasets** — several astrophotographers publish raw light frames:
  [Light Vortex Astronomy sample data](https://www.lightvortexastronomy.com/sample-image-data.html)
  (M42, ATIK 383L+ mono, per-filter FITS),
  [AstroBackyard practice data](https://astrobackyard.com/your-astrophoto-skills/),
  [Telescope Live curated datasets](https://telescope.live/datasets) (professional remote
  telescopes, pre-calibrated FITS; free tier requires an account).
  These are the right stress test for gradients, star bloat and imperfect tracking, but
  they arrive **without** a WCS, so the ground truth must be a reference solve.
* **Practical note**: for any tier-C image *without* a pipeline WCS, generate the reference with *two* independent
  solvers (`astap-cli` and `solve-field`) and only admit the image to the corpus if they
  agree to better than 2″. Disagreement means the truth is not trustworthy, not that one
  solver is wrong.

### 2.8 Sources added for the expanded corpus

Every source below was probed on 2026-09-29 with a small fetch whose header was inspected
before anything was designed around it. `scripts/fetch-corpus.py` implements each one.

| Source (fetcher) | Access | Pixels | Truth WCS | Tier | Notes |
|---|---|---|---|---|---|
| **CDS hips2fits**, 22 more HiPS surveys (`hips2fits`) | GET, no key | resampled to our TAN grid | exact (requested) | A | DSS2 r/b/IR, 2MASS J/H/K, PS1 g/r/i/z, SkyMapper DR4 g/r/i, DES DR2 r/i, Legacy DR10 r/z, unWISE/allWISE W1, ZTF DR7 g/r, IPHAS r, VHS J/K and UKIDSS-LAS K (WFAU HiPS), DENIS I, GALEX NUV, **TESS 2-yr (13″/px)** and **SHASSA continuum (26″/px)** for 3°–50° fields |
| **SkyView** (`skyview`) | GET, no key | resampled | requested grid | B (as v1) | DSS1 R/B, DSS2 R/B/IR, 2MASS J/H/K, WISE 3.4/4.6, GALEX NUV; `RADESYS = FK5` |
| **Legacy Surveys DR10** (`legacysurvey`) | GET, no key | coadd, resampled on request | header TAN | B | 0.26–1.0″/px, g/r/i/z |
| **Pan-STARRS1** (`panstarrs`) | two GETs | skycell stack, `output_size` rebinning | header TAN, legacy `PCiiijjj` | B | the cutout must stay inside one ~0.4° skycell; larger random cutouts come back NaN-padded (five dropped) |
| **SDSS DR17 frames** (`sdss`) | GET, bz2 | native, full frames | header TAN | B | random fields from a SkyServer SQL query (`quality = 3`) |
| **ZTF science images** (`ztf`) | IRSA IBE search + cutout | **native**, 1.01″/px | header **TPV** | B | 1500–3000 px cutouts of single exposures; seeing 1.3″–4.5″; galactic plane included |
| **WISE L1b single exposures** (`wise`) | IRSA IBE | **native**, 2.75″/px, 1016² | header **SIN-SIP** | B | poles included; W1 and W2 |
| **SkyMapper DR4 images** (`skymapper`) | SIAP + cutout | **native**, 0.5″/px, 16-bit | header **TPV** | B | cutouts capped at 0.17° by the service — right at d80's 0.15° floor |
| **TESS full-frame images** (`tess`) | AWS open data (`stpubdata`) | **native**, 21″/px | header **TAN-SIP**, strong distortion | B | 384–2048 px crops (2.2°–12°) of calibrated FFIs; the only native wide-field camera data in the corpus |
| **Las Cumbres Observatory** (`lco`) | archive API (anonymous, public frames) | **native** reduced exposures, fpacked | BANZAI pipeline TAN (astrometry.net vs Gaia) | C | 0.4 m + SBIG STL-6303 / QHY600, 1 m Sinistro, 2 m Spectral; frames with `WCSERR ≠ 0` rejected |

Verified headers worth knowing about:

* **ZTF** cutouts arrive gzip-compressed with `CTYPE = 'RA---TPV'` and 30 `PVi_j`
  coefficients. The TPV polynomial moves stars by up to ~0.2″ across a 0.56° cutout.
* **WISE L1b** frames use `CTYPE = 'RA---SIN-SIP'` — a *SIN* projection, which the old
  harness could not read. The SIP terms move the corners by 5.5″.
* **TESS** FFIs have an empty primary HDU; the image and its `TAN-SIP` WCS are in HDU 1,
  whose `A_2_0` of 2×10⁻⁵ displaces the corners of a full CCD by ~20 px (7′). One FFI
  (sector 69) carried no celestial WCS at all and was dropped. The fetcher crops the
  science area to a single-HDU float32 file and shifts `CRPIX`, which leaves SIP exact.
* **LCO** frames are RICE-compressed `BINTABLE`s (`ZIMAGE = T`); arcsec reads them through
  CFITSIO, the harness maps `ZNAXISn` back to `NAXISn`. The pipeline's `WCSERR` flag is
  checked, and three frames with a failed fit were dropped.
* **SkyMapper** cutouts are 16-bit with `CRPIX2 = -962` — another off-image reference
  pixel.

### 2.9 Investigated and not used

| Candidate | Why not |
|---|---|
| **nova.astrometry.net user uploads** — the obvious source of real amateur frames with solved WCS | `robots.txt` disallows all agents (and names AI crawlers specifically), image pages sit behind an "are you human" gate, and the per-image Creative Commons choice is only on those pages. Not scraped. |
| **Palomar Transient Factory** (IRSA) | Works, but headers carry both SIP and `PVi_j` terms under a `TAN-SIP` CTYPE, so the truth is ambiguous; 33 MB and ~100 s per frame. ZTF covers the same instrument class. |
| **Mellinger all-sky mosaic** (HiPS) | "Copyright Axel Mellinger. All rights reserved." Would have been the best 20°–60° source; SHASSA used instead. |
| **Kepler/K2 FFIs** | whole-focal-plane multi-extension files of hundreds of MB; TESS fills the role. |
| **UKIDSS/VISTA native frames** (WSA/VSA) | form-driven multi-extension products; the WFAU HiPS are used instead (tier A). |
| **DASCH / StarGlass plate scans** | large plates and an API that returned 404s during the probe; revisit for very wide fields. |
| **HST/MAST** | as §2.6: fields far below the catalogues' floor. |
| **Practice datasets** (Light Vortex, AstroBackyard, Telescope Live) | no redistribution licence stated, or an account required; no WCS. |
| **ESO archive, photutils data, AAVSO VPhot, Astrometrica samples** | as in v1: complex products, too few images, or accounts. |
| **LCO raw (`e00`) frames** | public and genuinely raw (bias, hot pixels, no flat), but truth would have to be carried over from the reduced frame; deferred. |

No licence-clear, scriptable source of real *amateur* frames with a trustworthy WCS was
found. LCO's 0.4 m telescopes (SBIG STL-6303 and QHY600 cameras) are the closest real
hardware in the corpus; tier S fills the rest synthetically, clearly labelled.

### 2.10 Tier S — simulated camera artefacts

`scripts/corpus_synth.py` (standard library only, seeded, deterministic) applies a recipe
to a *parent*: 16 camera-shaped tier-A cutouts (`cam_*`, 1200–2400 px, 0.6°–9°). Every
operation carries the truth exactly: geometric operations transform the CD matrix and
CRPIX with the pixels, trailing moves CRPIX by the photocentre shift, and radial
distortion `r' = r(1 + K r²)` about CRPIX is written as the exact SIP terms
`A_3_0 = A_1_2 = B_2_1 = B_0_3 = K`. The truth goes to a `<id>.truth` sidecar and the
image header carries no WCS, like a raw camera frame (`keepwcs` keeps it, for one case).

| Family | Operations | Count |
|---|---|---|
| amateur | vignetting, sky gradient, noise, hot pixels and columns, `u16` + `BZERO`; and one "worst case" stacking distortion, gradient, clouds, trailing, noise, hot pixels and a Bayer mosaic | 2 |
| lens distortion | barrel / pincushion, 4–25 px at the corner | 10 |
| colour camera | RGGB / GBRG Bayer mosaic (`BAYERPAT` set) | 5 |
| tracking, focus | trailing 6–12 px, three-box-blur defocus | 6 |
| sky | clouds (smooth random transmission + scattered light), saturation clipping, satellite trails, amp glow | 11 |
| orientation | `flipx`, `flipy`, `rot90`, `rot270`, `transpose` | 16 |
| formats | `u8`, `i32`, `f64`, `u16` gzip, XISF, WCS kept in header, 2×2 binning | 15 |

Tier D gains generated controls from the same module: pure noise, flat gradients,
random Gaussian "stars" (100–1000 of them — star-like, no real asterisms), and 64/128-px
block shuffles of real fields.

---

## 3. The corpus

Field selection is driven by what actually breaks solvers, not by what looks pretty.

### 3.1 Field-size ladder (tier A)

Solve behaviour changes qualitatively with FOV because it sets the spiral step, the
catalogue window and the quad scale. Sample geometrically:

```
    0.25°   0.5°   1.0°   1.5°   2.0°   3.0°   5.0°   10°
      │      │      │      │      │      │      │      │
      └─ small refractor + small sensor        └─ camera lens / all-sky
                                   └─ the regime our NGC 3372 data sits in (2.6°)
```

At the wide end, ASTAP-family databases thin out and distortion dominates; at the narrow
end, `d80` runs out of stars (published floor: 0.15°). Both ends are where we expect to
fail, and knowing *where* the cliff is, is the point.

The manifest carries the whole ladder, plus `stress_narrow` (0.10°) and `stress_wide15`
(15°); with `--auto-db` the wide end is solved from G05 rather than the D-series.

### 3.2 Star-density ladder (tier A/B)

Galactic latitude is the single best proxy for star density. Sample it deliberately:

| Field | RA, Dec (deg) | b (approx) | Character |
|---|---|---|---|
| Cygnus-X | 305.55, +40.73 | +2° | extremely dense, plane |
| Ophiuchus | 261.0, −7.0 | +13° | dense, near plane |
| Carina / NGC 3372 | 161.2, −59.7 | −1° | dense **and** strong nebulosity |
| Lyra (Vega) | 279.23, +38.78 | +19° | dense-ish, one very bright star |
| Perseus | 51.40, +49.85 | −8° | moderate |
| Bootes | 218.0, +35.0 | +64° | sparse |
| Ursa Major | 180.0, +65.0 | +51° | sparse |
| Cepheus (polar) | 340.0, +70.0 | +16° | polar geometry, `cos δ` stress |
| Near the NCP | 0.0, +88.0 | +25° | RA convergence, wrap-around |
| Near the SCP | 0.0, −88.0 | −27° | as above, southern |

The two polar fields matter disproportionately: `solver.rs` contains explicit pole-wrap
handling (`δ > π/2 → π − δ`, flip RA by π) and an RA-offset guard, and no unit test covers
it; the `pole_*` corpus entries are its only coverage.

### 3.3 Field-type ladder

| Type | Example | Why it is hard |
|---|---|---|
| Plain star field | Bootes 218, +35 | the control case — should always solve |
| Open cluster | M45 (56.75, +24.12), M67 | crowding; blended centroids |
| Globular cluster | M13 (250.42, +36.46), ω Cen | detection saturates in the core |
| Emission nebula | M42 (83.82, −5.39), NGC 3372 | background model fights the nebulosity |
| Face-on galaxy | M101 (210.80, +54.35) | extended source rejected as "too large"? |
| Very bright star | Vega, Sirius (101.29, −16.72) | bloom, diffraction spikes, saturation |
| Blank high-latitude | 195.0, +28.0 | few stars — tests the `n < 15` quad modes |

### 3.4 Negative controls (tier D)

These must return "no solution". Any solve reported here is a **false positive** and a
release blocker.

Only the first kind is in the manifest today: five tiny (0.03°–0.05°) HiPS crops of
nebula and galaxy cores (`neg_nebula_core`, `neg_m42_tiny`, `neg_m42_core2`,
`neg_m8_core`, `neg_m31_core`). Items 2–5 are proposed and not yet built.

1. A pure-nebula crop with almost no stars (small HiPS cutout inside M42's core).
2. Random Gaussian noise, no sources — synthesise locally, no download needed.
3. A correct star field, but solved with a **deliberately wrong hint** 40° away and a
   small `-r` — should exhaust the spiral and fail, not invent a solution.
4. A correct star field, but with a **deliberately wrong `--fov`** (2× and 0.5×) —
   currently expected to fail; documents the scale sensitivity of
   [§11.3](plate-solving.md#113-we-cannot-solve-without-a-good-pixel-scale-estimate).
5. A shuffled star field: take a solvable image and randomly permute 8×8 pixel blocks.
   Star-like sources, no real asterisms.

### 3.5 Distortion cases (tier D)

None of these is in the manifest as a tier-D entry; `fov_10p0` and `stress_wide15` are the
closest, and they are tier A and solve (worst corner error 4.5″ on the 15° field).

* An HST drizzled product (SIP in the header) — expected to fail on FOV grounds today.
* A wide-field (≥ 8°) HiPS cutout in TAN — at that width the tangent-plane approximation
  itself introduces several arcseconds of error at the corners, which is exactly the
  regime where a linear-only fit shows its limits.
* If available, a real camera-lens frame from tier C with visible barrel distortion.

### 3.6 The expanded corpus (v2)

v1 is a careful set of ladders around a handful of fields; it cannot say whether a
result generalises. v2 adds 532 entries chosen to cover the space rather than to
illustrate it, with random fields doing most of the work so that nobody — including the
solver's authors — picked them. The selection is recorded in
[`scripts/corpus-select.py`](../scripts/corpus-select.py) (seeded; archive queries are
cached) and its output, [`scripts/corpus.tsv`](../scripts/corpus.tsv), is the artefact.

| Group | Entries | How chosen |
|---|---|---|
| `rnd_*` random HiPS cutouts | 80 | uniform on the sphere; FOV log-uniform 0.2°–2.2°; survey drawn from those whose MOC encloses the field; 40 % rotated; 3:2 and 4:3 aspect ratios |
| `wide_*` | 27 | TESS HiPS 3°–24°, SHASSA continuum 15°–50°, DSS2 4°–12° deliberately under-sampled (`coarse`) |
| `obj_*` named fields | 36 | bright stars (Betelgeuse, Canopus, Polaris …), nebulae (Rosette, M16, Heart, Veil), galaxies (M33, M51, Cen A, LMC, SMC), clusters (47 Tuc, M11, M7), dark clouds (B68, Coalsack), Sgr A* in K |
| `surv_s_*` survey ladder | 11 | one southern field through 11 surveys |
| `sv2_*`, `ls2_*`, `ps1v2_*`, `sdss2_*`, `ztf_*`, `wise_*`, `smss_*`, `tess_*` | 221 | random positions within each archive's footprint (§2.8) |
| `lco_*` | 42 | random public reduced science frames 2016–2024, one per target and camera |
| `cam_*` parents + `s_*` tier S | 16 + 65 | §2.10 |
| negative and stress controls | 34 | §3.7 |

Entries that failed permanently on the first fetch (outside a footprint, a 404, a bad
pipeline WCS) stay in the manifest as commented rows with the reason, so the selection
remains auditable. 21 were dropped that way.

**Coverage** (all 635 entries; FOV and pixel scale measured from the delivered truth,
position from the image centre):

| FOV (long side) | A | B | C | S | D | total |
|---|---|---|---|---|---|---|
| < 0.15° | 1 | 5 |  |  | 9 | 15 |
| 0.15–0.3° | 16 | 92 | 5 |  | 6 | 119 |
| 0.3–0.6° | 33 | 55 | 37 |  | 5 | 130 |
| 0.6–1.2° | 90 | 57 |  | 30 | 9 | 186 |
| 1.2–2.5° | 60 | 12 |  | 16 | 6 | 94 |
| 2.5–6° | 8 | 28 |  | 7 | 1 | 44 |
| 6–20° | 18 | 6 |  | 12 | 3 | 39 |
| > 20° | 8 |  |  |  |  | 8 |

| Pixel scale | A | B | C | S | D | total |
|---|---|---|---|---|---|---|
| < 0.5″ |  | 80 | 17 |  | 14 | 111 |
| 0.5–1″ | 13 | 21 | 23 |  | 2 | 59 |
| 1–2″ | 161 | 81 | 2 | 35 | 16 | 295 |
| 2–4″ | 28 | 31 |  | 10 | 3 | 72 |
| 4–10″ | 5 |  |  | 1 | 1 | 7 |
| 10–30″ | 21 | 42 |  | 19 | 3 | 85 |
| > 30″ | 6 |  |  |  |  | 6 |

| Declination | −90…−60 | −60…−30 | −30…0 | 0…+30 | +30…+60 | +60…+90 |
|---|---|---|---|---|---|---|
| entries | 64 | 109 | 146 | 173 | 118 | 25 |

| Galactic latitude | \|b\| < 5° | 5°–15° | 15°–40° | > 40° |
|---|---|---|---|---|
| entries | 59 | 94 | 238 | 244 |

| Sample type | entries |
|---|---|
| 32-bit float | 376 (+ 42 fpacked float, LCO) |
| 16-bit signed (DSS, SkyMapper) | 128 |
| 16-bit unsigned via `BZERO = 32768` | 60 (+ 1 gzip, 3 XISF) |
| 64-bit float / 32-bit int / 8-bit | 13 / 9 / 3 |

| Orientation | entries |
|---|---|
| mirrored parity (det CD > 0) | 125 — ZTF, TESS, half of LCO, flipped tier S |
| rotated ≥ 45° from north-up | 205 |

| Truth | entries |
|---|---|
| exact (requested grid, tier A) | 243 |
| resampled grid (SkyView, Legacy Surveys) | 77 |
| survey header, linear TAN | 53 |
| survey header with SIP | 67 |
| survey header with TPV | 58 |
| observatory pipeline (LCO) | 42 |
| derived exactly from a parent (tier S) | 65 |
| none (negative controls) | 30 |

Distinct datasets (survey × band, or camera): 66 real ones from ten archives, plus the
synthetic families. Median image size 3.8 Mpixel; the corpus occupies 6.4 GB beyond v1
(8.2 GB including the hard-linked v1 files).

### 3.7 Negative and stress controls in v2

| Entries | What | Expected |
|---|---|---|
| `neg_noise_*` (4), `neg_flat_*` (2) | generated noise / smooth gradient, 0.5°–8° hint | no solution |
| `neg_fake_*` (6) | 100–1000 random Gaussian stars, 0.3°–10° hint | no solution |
| `neg_shuffle_*` (6) | real fields with 64- or 128-px blocks permuted | no solution |
| `neg_hint_*` (8) | real, solvable ZTF / Legacy / SDSS fields, hint 30°–60° away, `-r 10` | no solution |
| `neg_m17_core`, `neg_eta_car_core`, `neg_m31_ps1`, `neg_orion_2massk` | 0.04° crops of extended objects | no solution |
| `stress_fov_*` (4) | real fields told the wrong FOV (×2, ×0.5, ×1.5, ×0.67) | `expect=any`: a correct solve is fine, a wrong one is a false positive |

---

## 4. Metrics

For each image, compare the solved WCS against truth. The comparison must be done at
**several points on the image**, not just the centre — a centre-only check hides scale and
rotation errors, which is precisely the class of error a linear fit gets wrong.

```
    for p in {centre, 4 corners}:
        (α_t, δ_t) = truth_wcs(p)
        (α_s, δ_s) = solved_wcs(p)
        err(p)     = ang_sep((α_t,δ_t), (α_s,δ_s))   in arcsec
```

| Metric | Definition | Target |
|---|---|---|
| **Solve rate** | fraction of tier A–C images solved | ≥ 98% (A/B), ≥ 95% (C) |
| **Centre error** | `ang_sep` at the image centre | < 1″ median, < 3″ max |
| **Corner error** | max `ang_sep` over the four corners | < 2″ median, < 10″ max |
| **Scale error** | `\|CDELT2_solved − CDELT2_truth\| / CDELT2_truth` | < 0.05% |
| **Rotation error** | `\|CROTA2_solved − CROTA2_truth\|` | < 0.05° |
| **False-positive rate** | tier-D images reported as solved | **0** |
| **Wall time** | per image, warm cache | report median + p95 |

Corner error minus centre error is the useful derived quantity: it isolates scale,
rotation and distortion error from pointing error. It is also why SIP fitting
([§12.4](plate-solving.md#124-fit-sip-distortion-and-honour---sip---done)) cannot be
shown to help on this corpus: its images are reprojected and distortion-free, so the best
a `--sip` run can do is leave every solution linear, which it does. `benchmark.py
--extra-arg=--sip` scores SIP solutions through their polynomials.

### 4.1 Reference-frame caveats

* SkyView returns `RADESYS = FK5`; Legacy Survey, PS1 and Gaia are ICRS. The difference is
  < 0.1″ — below our targets, but it should be *stated* rather than silently absorbed, and
  if we ever tighten the centre-error target below 0.5″ it must be corrected for.
* Survey WCS headers are themselves solutions, with their own errors (0.02–0.3″). Tier B
  therefore cannot validate us below ~0.3″; only tier A can.
* Comparing `CRVAL` directly is wrong when the two solutions use different reference
  pixels. Always compare **sky positions of the same pixel coordinates**, as above.

### 4.2 Rules added for the expanded corpus

* **Distorted truth.** A linear plate cannot follow SIP or TPV distortion to the
  corners. For such images `benchmark.py` fits the best linear TAN plate to the truth
  (least squares over a 9 × 9 grid, tangent point at the centre) and reports its worst
  corner error as the **linear floor**; the false-positive threshold becomes
  `threshold + floor`. A solve that is still beyond that but whose centre is right
  (within 10″ or 2 px) is **INEXACT**: the field is identified, the plate is not good
  enough at the edges. INEXACT is neither a success nor a false positive and is counted
  separately. Linear floors in the corpus: ZTF ≈ 0.2″, SkyMapper ≈ 0.1″, WISE 5.5″,
  tier-S lenses 7″–140″, TESS 25″–1000″.
* **Coarse pixels.** 5″ is a quarter of a TESS pixel. The threshold is now
  `max(5″, 1 px)` (`--max-corner-px`); it changes nothing below 5″/px, so v1 results are
  unaffected.
* **Catalogue coverage.** A no-solve whose field (image height, which is what `--fov`
  carries) lies outside every installed database's range is reported as *nocat* rather
  than silently counted as a solver failure. With d80 + g05 + w08 installed that is
  only fields under 0.15°.
* **Per-entry overrides** (`extra` column): `file=` reuses another entry's image,
  `hint_dra` / `hint_ddec` move the hint, `fov_scale` lies about the FOV, `radius`
  overrides `-r`, `expect=any` marks a stress case.

---

## 5. Running the benchmark

| Script | Covers |
|---|---|
| **`scripts/fetch-test-images.sh`** | materialises tiers A, B and D from the manifest into `resources/testset/` |
| **`scripts/benchmark.py`** | solves every manifest entry, scores centre and corner error against the header truth, counts false positives, optional ASTAP head-to-head, CSV output |
| `scripts/hips_solve_test.sh` | blind mode (`-i`), fixed 1.5° FOV, 20 fields, solve rate + accuracy |
| `scripts/hips_extended_test.sh` | blind mode, 20 further fields including expected-failure types |
| `scripts/hips_fov_sweep.sh` | blind mode, 8 fields × 6 FOVs (0.5°–5°) |
| `scripts/bench_all.sh` | tier C, our NGC 3372 set (`resources/*.fits`), arcsec vs `astap_cli` |
| `scripts/compare-solvers.sh` | single-image arcsec vs ASTAP diff |

The three HiPS scripts download their own fields into `/tmp` and test the blind path; they
predate the manifest and are not part of the corpus.

Of the gaps this section originally listed, the corner-error computation and the single
run-everything command are done (`benchmark.py`), and `benchmark.py` exits 1 when arcsec
returns any false positive (including a tier-D solve), so it can gate a CI run on its own.

### 5.1 Manifest format

`scripts/test-images.tsv`, whitespace-separated columns (spaces, despite the extension),
`#` comments:

```
id              tier  ra         dec       fov_deg  width  height  source        extra
fov_1p50        A     51.400     49.850    1.5      4300   4300    hips2fits     CDS/P/DSS2/blue
dens_cygnusx    A     305.550    40.730    1.5      4300   4300    hips2fits     CDS/P/DSS2/red
sv_ngc188       B     11.800     85.240    0.5      1000   1000    skyview       DSS2R
ls_field1       B     180.000    20.000    0.284    1024   1024    legacysurvey  r
ps1_m67         B     132.830    11.810    0.284    4096   4096    panstarrs     r
sdss_2505_38    B     276.945    -0.164    0.225    2048   1489    sdss          301/2505/3/r/38
neg_m42_tiny    D     83.822     -5.391    0.03     400    400     hips2fits     CDS/P/DSS2/red
```

`ra`/`dec`/`fov_deg` are the ground truth for tier A; for tier B they are the *request*
and the truth is read back from the delivered header. `benchmark.py` takes the truth from
the delivered header for every tier.

The manifest holds 103 entries: 64 tier A, 34 tier B, 5 tier D. By source: 69 hips2fits,
11 Legacy Survey, 10 SkyView, 8 SDSS, 5 Pan-STARRS.

### 5.2 The expanded corpus

`scripts/corpus.tsv` extends the v1 columns with three more; `benchmark.py` accepts
either file (the extra columns are optional):

```
id        tier ra        dec      fov   width height source  extra                                    dataset   sets              truth
ztf_01    B    240.6466  13.1505  0.56  2000  2000   ztf     product=ztf_2024..._sciimg.fits;size=2000  ZTF-g     v2,random,midlat  header-tpv
s_ps1_a_lens S 207.8766  24.6735  1     2400  1600   synth   parent=cam_ps1_a;ops=distort:3.9e-09,...   synth-PS1-r v2,synth,distort derived
neg_hint_3 D   0         0        0     0     0      alias   file=sdss2_08;hint_dra=-52.5;radius=10      alias     v2,negative       none
```

For archive-defined geometry (ZTF, WISE, TESS, LCO, SDSS) the position and size columns
are nominal; truth always comes from the delivered file.

| Tool | Does |
|---|---|
| `scripts/fetch-corpus.py` | fetch / build everything; `--list`, `--stats`, `--tier`, `--source`, `--set`, `--dataset`, `--id`, `--force`, `--verify`, `--jobs` (≤ 4), `--reuse-dir` (hard-links v1 images already in `resources/testset`) |
| `scripts/benchmark.py --corpus` | as before, plus `--set`, `--source`, `--dataset`, `--by tier,source,dataset,fov,set` |
| `scripts/corpus_synth.py` | the tier-S / generated tier-D recipes (imported by the fetcher) |
| `scripts/fitslite.py` | stdlib FITS reading/writing and TAN/SIN/SIP/TPV WCS, shared by the three above |
| `scripts/corpus-select.py` | how the v2 fields were chosen (provenance; re-running queries live archives) |

The fetcher is deliberately gentle: at most four requests in flight and two per host,
0.5–1 s between requests to one host, a descriptive User-Agent, exponential backoff
honouring `Retry-After`, `.part` files renamed only after validation, and nothing
re-fetched that is already present. Every file is checked (FITS structure, a WCS the
harness can read, the centre where the manifest says, under 25 % blank pixels, and the
archive's own md5 where one is published — LCO, and the AWS ETag for TESS); a
`SHA256SUMS` of what was kept is written next to the images for `--verify`. Cutout
services stamp dates into headers, so a re-fetch is not byte-identical: the sums are a
local integrity check, not a global one.

A full fetch took 77 minutes on a home connection (465 files; the hips2fits queue
dominates, because it is limited to two concurrent requests), plus a minute of CPU
for tier S.

---

## 6. Results — 103 images, arcsec vs ASTAP

Corpus of 103 images (64 tier A, 34 tier B, 5 tier D), both solvers given the same true
field centre and FOV. ASTAP is `astap_cli` CLI-2026.07.30 with the d80 database; arcsec
runs with `--auto-db`, so it uses d80 up to 6° and G05 beyond. A reported solve counts as
**correct** only if its worst corner error is under 5″; anything else is a false positive.
The earlier rows of the table below used d80 only and a corpus that was one or two images
different in each tier as entries moved between tiers (§6.6).

Last re-measured 2026-09-25 on the 0.1.0 code, at both `-r 3` and the script's default
`-r 5` (identical results).

### 6.1 Where it started and where it stands

| | correct | false positives | tier A | tier B | tier D solved |
|---|---|---|---|---|---|
| Baseline (2026-09-01) | 48 | 4 | 26/62 | 22/34 | 0/7 |
| + 9-NN quad redundancy | 64 | 4 | 38/62 | 26/34 | 0/7 |
| + star-level verification | 86 | 0 | 54/63 | 32/34 | 0/6 |
| + spread check, no star trim | 89 | 0 | 55/63 | 34/34 | 0/6 |
| + `.290`/`.001` catalogues, `--auto-db` | 90 | 0 | 56/64 | 34/34 | 0/5 |
| + catalogue read across tiles (§7.7) | **92** | **0** | **58/64** | **34/34** | 0/5 |
| **ASTAP CLI-2026.07.30** | 47 | 0 | 39/64 | 8/34 | 0/5 |

Accuracy of arcsec's correct solves (2026-09-25):

| Tier | centre (median / max) | corner (median / max) | scale (median) | time (median, `--jobs 1`) | time (median, `--jobs 8`) |
|---|---|---|---|---|---|
| A | 0.740″ / 1.807″ | 0.999″ / 4.473″ | 0.0054% | 0.19 s | 0.52 s |
| B | **0.158″** / 1.123″ | **0.272″** / 1.719″ | 0.0059% | 0.11 s | 0.20 s |

Times depend heavily on `--jobs`: each arcsec process also uses every core, so running
eight images at once inflates every per-image time. Use `--jobs 1` for timing.

On the 45 images both solvers get right, arcsec's corner accuracy is slightly ahead
(0.897″ vs 0.938″) at essentially equal centre accuracy (0.690″ vs 0.660″). ASTAP solves
only two images arcsec does not (`type_m31`, `type_m44`); arcsec solves 45 that ASTAP does
not.

Speed on those 45 is at parity: median 0.15 s for both with `--jobs 1` (1.00×). It was
0.60× ASTAP when 126 quads per star and verification first landed; §6.5 is how it got
back.

### 6.2 What moved the numbers

**Quad redundancy (48 → 64).** For ≥60 stars, quads were built from each star's 3 nearest
neighbours — one quad per star. Which three stars are nearest depends entirely on which
stars are in the list, and the image and catalogue lists never match. Building all C(k,4)
subsets of the k nearest neighbours instead gives overlapping quads that survive a missing
neighbour. Recall against k, with and without verification:

```
  neighbours    5    6    7    8    9   10   12
  no verify    64   72   76   75   76   69    -     false positives  4  5  7 12 16 21
  verify       68   77   81   83   85   83   80     false positives  0  0  0  0  0  0  0
```

Without verification, recall and false positives rise together — redundancy buys recall by
spending precision. With verification, precision is free and 9 neighbours is the peak.

**Star-level verification (64 → 86, 4 FP → 0).** Every accepted position is now checked by
inverting the plate, projecting every catalogue star into pixel space, pairing each with the
nearest detected star, and re-fitting on those pairs (radii 6 → 3 → 2 px). It both rejects
coincidences and replaces a fit built from a handful of quad centroids with one built from
200–375 star positions — which is why it *raised* recall as well as killing false positives.

**Spread check + dropping the star trim (86 → 89).** A count threshold alone is not enough:
M31 at 2° passed with 22 matched stars clustered in the galaxy's core and a rotation error
of 1.56° (154″ at the corners). Requiring the matched stars to span ≥ 0.20 of the image
half-diagonal fixes that. Correct solves sit at 0.48–0.58, so it is not a tight threshold.

Spread alone could not separate M31 (22 stars, spread 0.221, rms 0.65″) from the Dec −88°
field (46 stars, spread 0.207, rms 0.66″) — which is *correct* to 2.3″ and was the
31° false positive before verification existed. The star count separates them, so the
minimum verified count went 12 → 30.

Removing the brightest-half star trim then took tier B to 34/34.

### 6.3 Negative results

Recorded so they are not re-tried:

* **Removing the star trim alone**, before verification: tier A +5, tier B −5, net zero.
  Only worth doing once verification exists.
* **A dense-field pre-threshold** (branch `exp/speed`): pick a higher detection level from a
  sampled pixel histogram so crowded frames measure fewer candidates. Solve count unchanged
  at 89/0, but **slower** — 0.60× → 0.54× (and 0.58× with strided sampling). The extra pass
  costs more than the HFD measurements it saves, and when it yields fewer than `max_stars`
  the normal cascade runs anyway. Reverted.
* **More stars.** `-s` 500 → 1000 → 2000 gives 89 → 85 → 76 correct and 0 → 3 → 2 false
  positives, with median time 0.32 s → 0.84 s → 8.15 s. Denser star lists make the k-NN
  neighbourhoods tighter, so centroid noise grows relative to quad size. 500 is the optimum.
* **Catalogue depth and matching tolerance**, swept earlier (`-s` to 16000, `-t` to 0.05):
  no effect on the failures they were aimed at.
* **The level-1 grid histogram limit** (`stars.rs`, `let upper = 65500.max(2 * background)`).
  The `max` means the `2 * background` term can never lower the limit, so every one of the
  ~169 grid cells histograms the full 16-bit range — 262 kB zeroed per cell. Reading it as
  an inverted `min` is the obvious conclusion, and it is wrong to act on: with `min`, the
  corpus still gives 90/103 and 0 false positives, but **total runtime rose from 56.6 s to
  62.7 s** and accuracy was a wash (15 images marginally worse, 10 marginally better, every
  difference sub-arcsec). Note the median time *falls* from 0.28 s to 0.18 s, which is what
  makes this trap convincing if you look at the median alone. Reverted, with the reason
  recorded at the line.

* **Accepting a scale hypothesis's solution as found** (§6.9). A field found at a
  scale well away from its own has been verified against a catalogue window of the
  wrong size. Measured with `--fov` × 0.5 and × 2 before the second solve existed: the
  corners a median 0.07″ worse than at the true scale, up to 1.7″, and with the ladder at
  × 0.25 two near misses past 5″ (`dens_lyra` 5.04″, `type_m45` 6.5″). Solving again at
  the solved scale (radius 0) removed all three and made the rest identical; it is now
  done whenever the solved scale is 5% off.
* **Trying a database only within √2 of its published range** in the ladder. It
  dropped the true scale of `stress_narrow` (0.1°, which d80, published from 0.15°,
  solves), so a factor 2 is used.
* **Cancelling ladder hypotheses after the winner**: kept, but it saves under 2% on the
  slowest no-scale solves; the earlier hypotheses, which must finish, dominate. Sharing
  detections between hypotheses saved another 5%.

### 6.4 What still fails (8 of 64 tier A, 0 of 34 tier B)

```
  dens_carina  dens_scutum  dens_vela        very dense galactic-plane fields
  type_dbl_clus  type_m8                     dense cluster / emission nebula
  type_m31  type_m44  type_sirius            dominated by one bright or extended object
```

Every remaining failure is a crowded field or a frame dominated by a single extended object,
where detection is saturated by nebulosity or bloom and the detected stars are not catalogue
stars. `fov_5p00`, `dens_bootes` and `pole_ra_wrap` — previously failures — now solve.

### 6.5 Speed — profiled, optimised and parallelised

Two rounds of work, both verified to leave the corpus result at exactly 89/103 correct
with 0 false positives.

**Round 1 — what the flamegraph said.** `cargo flamegraph` on a 1° field put **88% of
runtime in star measurement**:

```
  57.1%  measure_star
  12.7%  sort_by<f64>  (median_f64's, attributed to a neighbouring symbol)
  10.1%  drift::sort<f64, median_f64>
   8.3%  quicksort<f64, median_f64>
```

`measure_star` calls `median_f64` twice per candidate — annulus background, then its MAD —
and a crowded field produces tens of thousands of candidates, so a full `sort_by` for a
*single order statistic* was ~31% of total runtime.

| change | effect |
|---|---|
| `median_f64` by quickselect (`select_nth_unstable_by`) instead of a full sort | sorts 31% → 2% |
| annulus buffer and distance histogram on the stack, not heap-allocated per candidate | with the above, −15.5% |
| `round_sqrt`: exact integer `round(sqrt(n))` replacing `sqrt().round()` in the aperture loop | `round` was 6.3%, gone; −6% |
| quad matching over a compact `f32` array, indexed on `ratios[4]` not `ratios[0]` | −9% |

The last one is worth spelling out. `Quad` is 72 bytes, so scanning the tolerance window
touched far more memory than a five-ratio comparison needs; the ratios as `f32` are 20
bytes, and survivors are re-checked at full `f64` precision — the same trick
`catalog::anet` already used for its code array. And the index key matters: the ratios are
`d2/d1 … d6/d1` with `d` descending, so `ratios[0]` clusters near 1 and gives the widest
search window, while `ratios[4]` is the most spread and the most selective.

**Round 2 — parallelism.** Three independent axes:

* **Spiral positions** are independent, so they run a batch at a time across a thread pool.
  Semantics are unchanged: batches go in spiral order and the lowest index in a batch wins,
  so the position returned is the one the serial loop would have returned. Position 0 is
  tried alone first — it is the hint and usually solves outright, and spawning a pool for it
  made single-position solves *slower*.
* **Detection** scans horizontal bands. The `img_sa` "already detected" map makes the scan
  order-dependent, so each band gets its own marker buffer covering its rows plus
  `BAND_OVERLAP = 90` (the widest a marking can reach: radius `3 × hfd`, `hfd` capped at 30),
  and results are merged with a positional dedup.
* **`build_histogram` and `normalize_for_detection`** are exact reductions — integer
  histogram bins, and min/max/count — so splitting them gives bit-identical results.

```
                          serial 8-image set     median tier A     vs ASTAP
  before any of this           29240 ms              0.30 s          0.60x
  after round 1                20917 ms              0.29 s          0.75x
  + parallel spiral            17975 ms              0.27 s
  + parallel detection         17691 ms              0.23 s          0.94x
  + parallel reductions        17285 ms              0.20 s          1.07x
```

**arcsec is now faster than ASTAP again** (1.07×) while solving 89 images to its 47. The
spiral parallelism is what moved the multi-position cases (`dens_scutum` 7609 → 3795 ms);
detection parallelism is what moved the common single-position case.

Remaining profile of a typical solve: 50% `measure_star`, 10% `detect_pass_banded`, 7% quad
matching, 7% `value_subpixel`, 10% the quickselects. It is detection almost end to end, now
spread across cores, with nothing obviously wasteful left.

**Round 3 — a deep read of `measure_star` and the frame scan.** All four changes verified
**bit-identical**: same 89/103, same 0 false positives, and all 89 fitted positions
unchanged to five decimals.

| change | effect |
|---|---|
| centroid cascade by Chebyshev rings | **−23%** |
| annulus gathered by row-runs instead of testing the enclosing box | −9% |
| bilinear weights and bounds hoisted out of the aperture loops | −7% |
| compile-time `round_sqrt` table (integer sqrt was 3.6% of a profile) | −5% |
| frame scan: absolute threshold, tested before the marker lookup | −1.5% |

The ring decomposition is the interesting one. The centroid loop shrinks its aperture
`14 → 12 → 10 → …` until the star is "boxed", re-scanning the whole `(2rs+1)²` box each
time. But the per-pixel test does not depend on `rs`, so the box of half-width `r` is
exactly the union of Chebyshev rings `0..r`. Accumulate each ring once from a single scan
of the widest box and the cascade becomes prefix sums over at most 15 rings. A star that
shrinks 14 → 6 used to read 841+625+441+289+169 = 2365 pixels; it now reads 841.

The annulus one is the same idea applied to a shape: the background annulus is ~91 pixels
but was gathered by testing all 961 pixels of the enclosing 31×31 box, so 90% of the work
was rejection. Per row, `r1² < i²+j² ≤ r2²` is just `min_abs ≤ |i| ≤ hi` — two contiguous
runs, found with `isqrt`.

```
                          serial 8-image set     median tier A     vs ASTAP
  before any of this           29240 ms              0.30 s          0.60x
  after round 1                20917 ms              0.29 s          0.75x
  + parallel spiral            17975 ms              0.27 s
  + parallel detection         17691 ms              0.23 s
  + parallel reductions        17285 ms              0.20 s
  + round 3                    10337 ms              0.17 s          1.00x
```

**Cumulative: −65%.** arcsec is back at parity with ASTAP per image while solving 89 to its
47 — and it is doing considerably more work per solve (126 quads per star, plus star-level
verification). Measure the ASTAP ratio with `--jobs 1`; under contention it swings between
0.88× and 1.07× and means nothing.

Remaining profile of a typical solve: 38% `measure_star` (including its sampling closure),
11% quad matching, 10% the frame scan, 11% the two quickselects in `median_f64`, 6% the
background histogram. Well distributed, with no single obvious waste left. The quickselects
are the largest remaining single item — two per candidate over ~91 annulus pixels — but the
obvious ways to cut them (subsampling the annulus, a cheaper scale estimator) all change the
background estimate, and everything above was achieved without changing a single result.

The negative and neutral results from this work — including the level-1 grid histogram
limit, whose "obvious" fix is slower — are recorded in §6.3 so they are not re-tried.

Re-measured 2026-09-25 with `--jobs 1`: arcsec and ASTAP both at a 0.15 s median on the 45
shared solves (1.00×), with arcsec solving 90 to ASTAP's 47.

### 6.6 A truth bug found by the comparison

On `ps1_m67`, arcsec and ASTAP both reported corner errors of ~5626″ — and **agreed with
each other to 0.5″**. Two independent solvers converging on the same "wrong" answer is a
signal that the truth is wrong, not the solvers.

Pan-STARRS skycell headers carry the rotation matrix in the **legacy `PCiiijjj` form**
(`PC001001 = -1.0`), not the modern `PC1_1`. The harness read neither, fell back to
`CROTA2 = 0`, and so built a truth WCS with the **RA axis sign flipped**. Combined with the
`CRPIX2 = -24034` that PS1 cutouts carry, that lever arm threw the corners out by more than
a degree.

With `PCiiijjj` support added, `ps1_m67` solves correctly. The lesson for the corpus:
**always sanity check a "failure" that two independent solvers agree on.**

`stress_narrow` was likewise misclassified: at 0.10° it was assumed to be below the d80
floor, but it solves correctly to 0.607″, so it moved from tier D to tier A.

### 6.7 An off-centre hint exposed a rotation error (fixed in 0.1.1)

Every run above gives arcsec the true centre as its hint, and that hid a real bug: the
plate was fitted in the tangent plane of the spiral position that matched, not of the
image centre, and a linear fit there absorbs the projection's curvature as a rotation
that grows with the offset and with declination. The centre stayed right (under 1″) and
the star-level RMS stayed small, so nothing flagged it, but the corners did not:

| `--offset-hint` | tier A correct (0.1.0 → 0.1.1) | tier A false positives | tier B correct |
|---|---|---|---|
| 0 | 56 → 56 | 0 → 0 | 34 → 34 |
| 0.3 fields | 8 → 55 | 47 → 0 | 29 → 31 |
| 0.6 fields | 4 → 51 | 48 → 1 | 23 → 29 |

In 0.1.0 the errors reached 1.7° of rotation at Dec −80 and 1600″ at the corners of the
10° field. The fix (`recentre` in `pipeline/solver.rs`) refits the verified matches in
the image centre's tangent plane. The one remaining 0.6-field false positive is
`fov_5p00` at 5.9″, solved from a spiral position 3° away with 49 verified stars.

Real hints are rarely exact (a mount's position, a blind estimate), so **also run the
benchmark with `--offset-hint 0.3`** when changing matching or the fit.

### 6.8 Reproducing

```bash
scripts/fetch-test-images.sh                       # 103 images, ~2.2 GB
scripts/benchmark.py --auto-db --astap ~/astap_cli --radius 3 --csv results.csv
scripts/benchmark.py --auto-db --astap ~/astap_cli --radius 3 --jobs 1   # for timings
```

`benchmark.py` expects the star databases in `~/star_database` (`--db`) and the binary at
`target/release/arcsec` (`--arcsec`). Without `--auto-db` it passes `-D d80`
(`--db-name`), and `stress_wide15` (15°) then fails: 89/103 rather than 90.

### 6.9 Unknown or wrong pixel scale (2026-10-05)

The scale search and the second solve at the solved scale of
[plate-solving.md §10.3e](plate-solving.md#103e-unknown-or-wrong-pixel-scale-autoscalers-2026-10-05),
measured against 0.5.1 (`main`, e7e1f69). `--no-scale` withholds the scale (no `--fov`,
and any FOCALLEN/XPIXSZ removed from a copy; no corpus image has both, so none needed
it); `--fov-scale X` passes the true field times X; `--fov-search` is the new opt-in
flag, passed with `--extra-arg`. `--auto-db`, true-centre hint, `-r 5`, no index in the
catalogue or database directory (so the hinted path alone), 8 jobs; times are only
comparable within a row.

```bash
scripts/benchmark.py --auto-db --no-scale
scripts/benchmark.py --auto-db --fov-scale 0.5 [--extra-arg=--fov-search]
```

v1, 103 images (95 of them solvable at the true scale):

| mode | 0.5.1 correct / FP | now correct / FP | now, `--fov-search` |
|---|---|---|---|
| true `--fov` | 95 / 0 | 95 / 0 (every `.wcs` scores identically) | – |
| true `--fov`, hint 0.3 fields off | 95 / 0 | 95 / 0 | – |
| no scale | 85 / 0 | **98 / 0** | – |
| no scale, hint 0.3 fields off | 88 / 0 | **97 / 0** | – |
| `--fov` × 0.5 | 79 / **1** (`dens_lyra`, 5.04″) | 80 / 0 | **98 / 0** |
| `--fov` × 2 | 69 / 0 | 69 / 0 | **98 / 0** |
| `--fov` × 4 | 3 / 0 | 3 / 0 | **95 / 0** |
| `--fov` × 0.25 | 0 / 0 | 0 / 0 | **98 / 0** |

98 is more than the 95 that solve at the true scale: `dens_scutum`, `type_dbl_clus` and
`type_m8`, crowded fields of §6.4, solve at a hypothesis smaller than their scale (a
catalogue window narrower than the frame), and then again at their own. Every image that
solves both ways solves to the same answer: against the true-`--fov` run, 92 of 95
identical with the scale withheld and 93 of 95 with `--fov-search` at every factor,
none more than 0.004″ apart at the corners.

The expanded corpus, 635 entries, same protocol:

| | 0.5.1 | now |
|---|---|---|
| true `--fov`: correct / FP | 593 / 1 | 593 / 1 (589 identical; the other four are the `stress_fov_*` controls, which pass the field × 0.5–2 on purpose and are now solved again at their own scale: within 0.04″ of before, about 1 s more each) |
| no scale: correct / FP | 433 / 0 | **585 / 1** |
| no scale: time in failed searches (total) | 692 s | 168 s |
| tier D (39 controls), no scale: total time | 93 s | 114 s |

The one false positive is `wide_shassa_01`, which is the 0.5.1 false positive with the
true field too (a 30″ corner on a wide SHASSA plate); with the scale withheld 0.5.1 did
not solve it at all. No new image is solved wrongly.

**Time.** `--jobs 1`, two interleaved rounds, v1 (s, median / mean / p90 / total):

| | 0.5.1 | now |
|---|---|---|
| true `--fov`, all 103 | 0.086 / 0.244 / 0.56 / 25.1 | 0.087 / 0.242 / 0.56 / 25.0 |
| no scale, all 103 | 0.088 / 0.302 / 0.84 / 31.2 | 0.18 / 0.57 / 1.6 / 58.9 |
| no scale, the 5 negative controls, each | 0.03 | 0.06 |

A solve with the scale known costs what it did. Without one, a solve costs about
twice what it did when 1″/px happened to be close enough, because a solution more than
5% from 1″/px is solved again at its own scale (the tier A frames are 1.24″/px), and
the 13 images 0.5.1 could not solve take 0.9–6 s. The ladder's worst case is bounded:
17 hypotheses of 9 positions, one detection per binning and minimum star size, no
seeded fallback; on the corpus's controls it adds 0.04–1.7 s to a failure (`neg_hint_*`
at `-r 10`: 4.6–8.1 s → 4.8–9.4 s).

**With the blind index present** (287 MB `d80.arcsecix` beside the database), `-r 10`
so that it is consulted: true `--fov` 95 / 0 in both; no scale 88 / 0 (0.5.1: the
index, which searches 0.3–60″/px, rescues three) against **98 / 0**, and the failed
searches' total time 71 s → 0.8 s.

**astap_cli with `-fov 0`** does the same kind of search (`Trying FOV: 9.5`, 6.3, 4.2,
… 0.37°, a full `-r` spiral at each), and prints and writes to its `.ini` the
`Warning scale was inaccurate! Set FOV=…d, scale=…"` that arcsec now writes too. Probed
on `decp20`/`decp40` with `-fov` × 0.94–1.06, it warns beyond about 5%.

## 7. Results — expanded corpus (635 entries)

Measured 2026-10-01 on the code after the catalogue-read fix (§7.7, unreleased), against
the 0.1.2 release binary as "before" (the numbers this section showed until then).
`--auto-db` (d80 0.15°–6°, g05 3°–20°, w08 20°–80° installed), true centre as the hint,
`-r 5`, 12 jobs. Four solver changes since then raise the counts to 540 correct (true
centre) and 493 (offset hint); §7.8 has those results and every status change. The
distortion model (§7.9) takes them to 558 and 537, and the crowded-field work (§7.10)
to 589 and 580.

```bash
scripts/benchmark.py --corpus --auto-db --by tier,source,dataset,fov,set --csv run.csv
scripts/benchmark.py --corpus --auto-db --offset-hint 0.3
scripts/benchmark.py --corpus --auto-db --astap ~/astap_cli
```

A full arcsec run takes 2.5 minutes on 24 cores. These are measurements on a benchmark
built by the same project, from images chosen for coverage rather than typical use; they
are not a statement about how often arcsec solves a user's frames, and the ASTAP column
in §7.5 is a reference point, not a ranking.

### 7.1 Overall

| | n | correct (0.1.2 → now) | false positives | inexact | no solve |
|---|---|---|---|---|---|
| Tier A — synthetic cutouts | 234 | 192 → **206** (88 %) | 1 → 1 | – | 41 → 27 |
| Tier B — survey pixels | 255 | 196 → **206** (81 %) | 2 → 3 | 5 → 4 | 52 → 42 (2 below any catalogue) |
| Tier C — LCO frames | 42 | 37 → **38** (90 %) | 0 → 0 | – | 5 → 4 |
| Tier S — simulated camera artefacts | 65 | 50 → **54** (83 %) | 2 → 0 | 6 → 7 | 7 → 4 |
| **A + B + C + S** | **596** | **475 → 504 (85 %)** | **5 → 4** | **11 → 11** | **105 → 77** |
| Tier D — must-fail controls | 35 | – | **0 → 0** | – | 35 correctly refused |
| Tier D — wrong-FOV stress (`expect=any`) | 4 | 3 → 3 solved correctly | 0 | – | 1 |
| v1 subset | 98 + 5 | 90 → **92** | 0 → 0 | – | 8 → 6 |

Accuracy of the correct solves (median centre / corner, now): A 0.52″ / 0.77″, B 0.16″ /
0.33″, C 0.23″ / 0.63″, S 0.60″ / 1.01″ (0.1.2: A 0.53″ / 0.81″, B 0.15″ / 0.33″, C 0.24″ /
0.62″, S 0.62″ / 1.25″). Of the 475 images both versions solve, 39 moved more than 0.2″
closer to the truth at the worst corner and 9 moved further away. Median time 0.2–0.5 s
per image.

**By field of view** (long side):

| FOV | n | correct (0.1.2 → now) | FP | inexact |
|---|---|---|---|---|
| < 0.15° | 6 | 4 → 4 | 0 → 0 | 0 → 0 |
| 0.15–0.3° | 113 | 103 → 103 (91 %) | 0 → 0 | 0 → 0 |
| 0.3–0.6° | 126 | 116 → 117 (93 %) | 0 → 0 | 0 → 0 |
| 0.6–1.2° | 180 | 156 → **163** (91 %) | 0 → 0 | 6 → 2 |
| 1.2–2.5° | 88 | 62 → **75** (85 %) | 0 → 0 | 2 → 2 |
| 2.5–6° | 43 | 9 → **17** (40 %) | 1 → 0 | 2 → 5 |
| 6–20° | 36 | 23 → 23 (64 %) | 4 → 4 | 1 → 2 |
| > 20° | 8 | 5 → 5 | 0 → 0 | 0 → 0 |

**By source:**

| Source | n | correct (0.1.2 → now) | FP | inexact | centre / corner (median, now) |
|---|---|---|---|---|---|
| hips2fits (22 surveys) | 234 | 192 → **206** (88 %) | 1 → 1 | 0 | 0.52″ / 0.77″ |
| ZTF | 43 | **43** → **43** | 0 | 0 | 0.10″ / 0.35″ |
| SDSS | 33 | **33** → **33** | 0 | 0 | 0.20″ / 0.33″ |
| Pan-STARRS1 | 20 | **20** → **20** | 0 | 0 | 0.20″ / 0.32″ |
| SkyView | 37 | 35 → 36 | 0 | 0 | 0.55″ / 0.68″ |
| Legacy Surveys | 40 | 37 → 37 | 0 | 0 | 0.05″ / 0.10″ |
| LCO (tier C) | 42 | 37 → 38 | 0 | 0 | 0.23″ / 0.63″ |
| WISE L1b | 25 | 16 → **20** | 0 | 4 → 0 | 0.32″ / 8.6″ (5.5″ linear floor) |
| SkyMapper native, 0.17° | 15 | 8 → 8 | 0 | 0 | 0.07″ / 0.23″ |
| TESS FFI crops | 42 | 4 → **9** | 2 → 3 | 1 → 4 | 6.9″ / 49″ |
| tier S | 65 | 50 → **54** | 2 → 0 | 6 → 7 | 0.60″ / 1.01″ |

The gains are where §7.7 says they should be: fields of 0.6°–6°, which straddle a d80 tile
boundary more often than not, and above all the 1.2°–6° bands, where a tile is not much
larger than the field. Among the tier-A HiPS (0.1.2 figures): TESS 2-yr 10/16 (3°–24°),
SHASSA 4/8 (15°–50°), coarse DSS 3/6, named objects 29/36, random galactic-plane fields
17/23; now `cam_ztf_a`, the three 3°–5° TESS-HiPS fields, `obj_heart`, `obj_rosette`,
`obj_orion_belt`, `type_m31`, `type_m44` and six random fields solve as well, and
`obj_47tuc` no longer does (§7.7).

**Tier S by family** (0.1.2): Bayer mosaics 5/5, trailing 3/3, saturation 2/2, satellites
and junk 5/5, binning 3/3, clouds 3/4, formats (u8, i32, f64, gzip, XISF, header WCS)
10/12, orientation 13/16, defocus 2/3, amateur 1/2 (+1 inexact), lens distortion 3/10 (5
inexact, 1 false positive). Every format and orientation miss has a parent that also fails
(`cam_dss_c`, `cam_ztf_a`) or is a coarse TESS-HiPS parent, so none is a format bug. Now
`cam_ztf_a` solves and so do its clouds and rot90 variants, `s_dss_c_flipx`,
`s_tess_c_defocus` and `s_tess_c_flipx`; `s_tess_c_flipy` no longer does.

### 7.2 With the hint 0.3 fields off (`--offset-hint 0.3`)

| | correct (0.1.2 → now) | FP | inexact |
|---|---|---|---|
| A | 183 → **199** | 4 → 1 | 0 → 0 |
| B | 169 → **173** | 3 → 2 | 24 → 30 |
| C | 39 → 39 | 0 → 0 | 0 → 0 |
| S | 38 → **48** | 2 → 1 | 9 → 9 |
| **A + B + C + S** | **429 → 459** | **9 → 4** | **33 → 39** |
| v1 subset | 86 → 86 | 0 → 0 | 0 → 0 |
| tier D | – | 0 → 0 | – |

By field of view the gains are again 0.6°–2.5° (+24) and 2.5°–6° (+4). Tier B's six more
inexact solves are seven TESS FFI crops that used to fail outright and now solve to the
linear plate's distortion floor, less `ztf_34`, which went from inexact to correct.

In 0.1.2 the losses against the true-centre hint concentrated where the image carries
distortion: WISE drops from 16 correct to 1 (19 inexact — fitted from an off-centre spiral
position, the linear plate follows the SIN-SIP distortion differently), TESS, wide
TESS-HiPS fields and tier-S lenses. That is unchanged. Three of 0.1.2's four tier-A false
positives here (`rnd_057` 5.6″, `rnd_074` 8.5″, both galactic-plane fields, and
`wide_tess_02`) are now correct solves.

### 7.3 False positives, all of them (true-centre hint)

| Entry | FOV | Corner error | What it is |
|---|---|---|---|
| `wide_shassa_01` | 15° at 26″/px | 31.0″ (1.2 px) | a correct field just over the one-pixel threshold |
| `tess_03`, `tess_25`, `tess_32` | 12° TESS FFI crops | 2200–2830″, centre 100–180″ | right area, but one linear plate across 1000″ of distortion: genuinely wrong at the edges |

`tess_03` is new with the catalogue fix (it failed with the true-centre hint, and was
already a false positive with the offset one); the two tier-S false positives of 0.1.2
went, `s_tess_b_pincush` (279″, now no solve) and `s_tess_c_flipx` (12.2″, now correct at
3.2″). None in tiers C or D, and none below 6°. The dangerous ones are the three wide
*distorted* cases: once the fit residuals show structure, a linear plate should be
refused — or a distortion model fitted — rather than reported as a solve.

*Update (§7.9):* the distortion model now reports the TESS frames correctly (corners at
the ~1000″ linear floor) and `s_tess_b_pincush` too; only `wide_shassa_01` remains.

### 7.4 Failure categories that point at solver work

1. **Native wide-field camera data: TESS FFIs 9/42** (4/42 in 0.1.2). 21″/px, a PSF of 1–2 pixels, and
   25″–1000″ of SIP distortion. The solver *finds* the field — on `tess_05` (2.2°) the
   first spiral position reports 614 matching quad references — and then accepts no
   solution. The same sky through the TESS HiPS (resampled, undistorted) solves 10/16.
   Candidates: verification radii (6 → 3 → 2 px) smaller than the distortion, and the
   handling of undersampled stars (`Minimum star size: 1.5″`, HFD). This is the closest
   thing in the corpus to a camera-lens frame, and the clearest case for
   [plate-solving.md §12.4](plate-solving.md#124-fit-sip-distortion-and-honour---sip).
   *Update (§7.9):* the distortion model makes every TESS frame that solves correct
   (23/42 centre, 16/42 offset); those still failing are the correspondence problem.
   *Update (§7.10):* it was mostly a detection problem — the bright stars were refused
   as too large for the measuring box — and the catalogue-seeded fallback does the
   rest: 41/42 centre, 37/42 offset.
2. **Lens distortion (tier S): 3/10 correct, 5 inexact, 1 false positive.** arcsec's
   corners come out *worse than the best linear plate* (for example 28″ against a 16.5″
   floor): the star-level refit keeps only matches within 2 px, so the edges, where the
   distortion is largest, drop out and the plate is fitted to the centre. WISE shows the
   same pattern. *Update (§7.9):* fixed — 10/10, and WISE 20/25 with either
   hint.
3. **2.5°–6° is the weakest band (17/43, 9/43 in 0.1.2).** It is where d80 hands over
   to g05 and pixel scales pass 10″/px; most entries are TESS (FFI or HiPS) or coarse DSS.
4. **Real-telescope misses (LCO 37/42).** M43 (bright nebula), two frames with FWHM 6–7″
   (seeing or defocus), a sparse *z*-band field, and a 10′ 2 m frame (below d80's floor).
   ASTAP solves four of the five (§7.5), so they are worth a direct look.
5. **Narrow native frames at the catalogue floor.** SkyMapper at 0.17° solves 8/15;
   ASTAP solves 14/15 of the same images.
6. **Crowded and nebulous fields**, as in v1: named nebulae (Rosette, Heart, M16, Orion's
   belt, Coalsack, B68) fail; random galactic-plane fields solve 17/23. *Update
   (§7.10):* of the crowded DSS fields only `dens_scutum`, `type_dbl_clus` and `type_m8`
   still fail.

### 7.5 ASTAP on the same images (reference only)

Measured against arcsec 0.1.2, before the catalogue-read fix. `astap_cli`
CLI-2026.07.30 with `-D d80` only (so no catalogue above 6°), the same hints and the same
scoring. ASTAP was run with the harness's defaults and no per-image tuning;
its low counts on Legacy Surveys (0/40) and SDSS (5/33) repeat v1's pattern and probably
reflect settings (downsampling, star count) as much as ability. Read this as "where the
two differ", not as a ranking.

| | n | arcsec | ASTAP |
|---|---|---|---|
| correct, fields ≤ 6° | 552 | 447 | 300 |
| false positives (A/B/C/S) | 596 | 5 | 2 |
| tier-D false positives | 35 | 0 | 0 |
| median centre / corner error, the 278 both solve | | 0.27″ / 0.55″ | 0.26″ / 0.70″ |

ASTAP solves 24 images arcsec does not: 6 SkyMapper 0.17° cutouts, 4 of the 5 LCO misses,
8 random tier-A fields (mostly galactic plane; PS1, ZTF, 2MASS, SkyMapper), `obj_rosette`,
`obj_orion_belt`, v1's `type_m31` and `type_m44`, and two tier-S images of `cam_ztf_a`.
Those 24 are the most direct pointers to what arcsec's detection or matching still misses.
Since the catalogue-read fix arcsec solves `obj_rosette`, `obj_orion_belt`, `type_m31`,
`type_m44` and the `cam_ztf_a` images among them too.

### 7.6 Ground-truth checks

* **Cross-solver agreement.** On the 287 images both solvers answer, the two agree with
  each other much better than either agrees with a *resampled* truth: hips2fits 0.11″
  between solvers against 0.37″ to truth, SkyView 0.08″ against 0.48″, LCO 0.03″ against
  0.24″ (ZTF: 0.02″ against 0.09″). So tier A's "exact" truth is exact for the *grid*,
  but the stars in it carry the underlying survey's plate solution (DSS ~0.3–0.5″), and
  the LCO pipeline WCS is good to ~0.25″. Claims below those levels need survey headers
  such as ZTF's.
* **No survey header was found wrong** in the §6.6 sense (both solvers agreeing with each
  other everywhere but not with the truth). The one image flagged, `s_dss_a_pincush`, is
  a synthetic distortion case where both linear solvers fit the same compromise.
* **WISE's `A_0_0` — ignored, by decision.** WISE L1b headers carry SIP constant terms
  (`A_0_0` ≈ 0.72 px, `B_0_0` ≈ −0.08 px), which the SIP convention does not define.
  Applied, every solved WISE frame came out 0.64 px off in the same *pixel* direction
  whatever its orientation on the sky; ignored, ~0.1 px — so the harness ignores them.
  But on the one WISE frame ASTAP solves (`wise_13`), ASTAP agrees with the *applied*
  convention (1.8″ from arcsec). The WISE Explanatory Supplement (§IV.4.d) says the
  per-frame SIP terms absorb a differential-aberration fit, which does not settle it.
  Decision (2026-09-30): keep ignoring them, and treat WISE centre errors as uncertain
  at the 2″ level.
* **Fetch-time checks** dropped three LCO frames whose pipeline flagged its own fit
  (`WCSERR ≠ 0`) and one TESS FFI with no celestial WCS.

### 7.7 The catalogue read across tiles (2026-10-01)

[plate-solving.md §11.11](plate-solving.md#1111-the-1476-catalogue-read-starves-every-tile-but-the-first--fixed-2026-10-01):
the `.1476` reader filled the whole star budget from the first tile a field overlaps, so a
field straddling a tile boundary was matched on one side only; the `.290` reader
under-sampled a tile covering most of the field. Both now return the field's brightest
`max_stars` stars from every tile it overlaps. The numbers in §7.1–7.3 are after the fix.

**Images whose status changed, true-centre hint** (32 newly correct, 3 no longer solved,
net +29):

* Now correct — 26 from no solve: `cam_ztf_a`, `lco_29`, `obj_heart`, `obj_orion_belt`,
  `obj_rosette`, `rnd_026`, `rnd_035`, `rnd_036`, `rnd_057`, `rnd_074`, `rnd_079`,
  `s_dss_c_flipx`, `s_tess_c_defocus`, `s_ztf_a_clouds`, `s_ztf_a_rot90`, `sv2_06`,
  `tess_19`, `tess_22`, `tess_24`, `tess_33`, `tess_37`, `type_m31`, `type_m44`,
  `wide_tess_01`–`03`; 5 from inexact: `tess_13`, `wise_04`, `wise_09`, `wise_10`,
  `wise_24` (the WISE corners fall from 10.5–26″ to 8.3–9.5″, inside their linear
  floor); 1 from false positive: `s_tess_c_flipx` (12.2″ → 3.2″).
* Now inexact rather than no solve: `s_tess_b_lens`, `tess_12`, `tess_15`, `tess_21`,
  `tess_23` — distorted TESS fields that now find the right field and stop at the linear
  plate's floor.
* No longer a false positive: `s_tess_b_pincush` (now no solve).
* Newly a false positive: `tess_03` (12° TESS FFI, 2830″ at the corners). It was a false
  positive in 0.1.2 too, with the offset hint; the cause is the 1000″ of distortion, not
  the catalogue.
* No longer solved: `obj_47tuc`, `s_tess_c_flipy`, `tess_09`. On 47 Tuc (0.6°) the
  field's 500 brightest catalogue stars now come from the cluster, which Gaia resolves and
  the image does not, so the catalogue reaches only mag 13.9 where the image's stars go
  to mag 16; the old read took one tile that held none of the core and happened to match
  the image better. `s_tess_c_flipy` and `tess_09` are coarse (12–21″/px) TESS frames
  that now fail verification at the first position; their siblings (`s_tess_c_flipx`,
  `s_tess_c_rot270`, the other FFI crops) gained.

**With the hint 0.3 fields off** (45 better, 7 worse): 36 new correct solves (among them
`dens_bootes`, `dens_norma`, `obj_ic1396`, `obj_rigel`, `obj_rosette`, `obj_smc`,
`type_m44`, eight tier-S SkyMapper and DSS variants, `wide_tess_10`/`11` and `ztf_34`
from inexact), of which 4 were false positives (`rnd_057`, `rnd_074`, `s_tess_c_defocus`, `wide_tess_02`),
`tess_03` and `tess_25` no longer false positives, and seven TESS FFI crops inexact
rather than unsolved. Lost: `fov_10p0`, `s_dss_c_flipx`, `s_dss_c_xisf`, `type_dbl_clus`,
`type_m31`, `wide_dss_05` (no solve), and `tess_43`, a new 988″ false positive of the
distorted-TESS kind. Most losses were marginal solves before (30–70 matched quads, verified
barely above the 30-star minimum): with the hint off-centre only part of the catalogue
window overlaps the image, and a one-sided read that happened to cover the overlap
matched it better than the whole window's brightest does. `fov_3p00` and `s_ps1_b_i32`
are of the same kind: they solve in this run and failed in the run with the position-hash
cut below.

**Matched-star coverage.** With `--sip`, the v1 solves whose verified stars leave a cell
of the 3×3 grid empty fell from 26 of 90 to 7 of 92. In the remaining seven the gap is in
the image (`ps1_big_c`'s window holds only 297 catalogue stars, all of them read).

**Speed.** The v1 subset with `--jobs 1 --threads 1`, two alternating rounds: total
71.0 s → 62.6 s and 70.0 s → 69.7 s, median 0.19 s → 0.18 s and 0.18 s → 0.19 s. A full
spiral that finds nothing (`neg_shuffle_1`, `--threads 1`) went from 7.5 s to 6.7 s:
stopping every tile at the field's magnitude limit reads fewer records than filling the
budget from one tile.

**Negative result — breaking magnitude ties by position.** The budget is usually cut
part-way through a 0.1 mag group, and within a group a tile's records run south to north,
so the stars kept from the last group favour the south of the field (as they always have,
and as in ASTAP). Choosing them by a hash of position instead makes the cut spatially
even. Over the whole corpus it was a wash — true-centre hint 506 correct / 5 false
positives against the file-order cut's 504 / 4, offset hint 458 / 6 against 459 / 4, the
false-positive differences all in distorted TESS fields — so the simpler file-order cut,
which is also what a single-tile read has always done, was kept.

### 7.8 Solver robustness (2026-10-01)

Four changes after 0.2.0, from a diagnosis of the remaining failures, each measured on
the whole corpus (both hints, `--auto-db`, 12 jobs) against the 0.2.0 release and against
the change before it. Cumulative, in the order they were committed:

| | true-centre hint: correct / FP / inexact | offset hint: correct / FP / inexact | tier D FP | v1 (centre / offset) |
|---|---|---|---|---|
| 0.2.0 | 504 / 4 / 11 | 459 / 4 / 39 | 0 / 0 | 92 / 86 |
| 1. similarity check by singular values | 504 / 4 / 11 | 460 / 4 / 39 | 0 / 0 | 92 / 86 |
| 2. sigma-clip the quad pairs | 525 / 7 / 12 | 475 / 6 / 38 | 0 / 0 | 92 / 89 |
| 3. database limit on image stars | 528 / 7 / 12 | 479 / 6 / 38 | 0 / 0 | 92 / 90 |
| 4. sparse images | **540 / 7 / 12** | **493 / 6 / 38** | **0 / 0** | 92 / 90 |

By tier, true-centre hint, 0.2.0 → now: A 206 → 221, B 206 → 220, C 38 → 42, S 54 → 57;
offset hint: A 199 → 209, B 173 → 187, C 39 → 42, S 48 → 55. No image that 0.2.0 solved
correctly is lost in either run. `stress_fov_3` (tier D, `expect=any`) now solves,
correctly. Of the 507 images both versions solve with the true-centre hint none moved more
than 0.2″ closer to or further from the truth at the worst corner except `tess_24` (49.1″
→ 50.1″, linear floor 42″); with the offset hint 12 moved closer and one further
(`type_m45`, 2.9″ → 3.6″).

**1. Plate similarity by singular values** ([plate-solving.md §7.2](plate-solving.md#72-the-six-plate-constants)).
`solve_plate_constants` compared the two *row* norms of the plate, which a sheared matrix
passes; in 0.1.2 the tier-D control `neg_hint_3` (offset hint) solved to a plate with rows
2.73 and 2.64 long and singular values in the ratio 3.03. It now requires σmax/σmin ≤ 1.08;
the largest on any correct solve is 1.027 (TESS FFIs), every other ≤ 1.0066. No status
change with the true-centre hint; with the offset hint `obj_heart` solves (its first
neighbouring position's quad fit was anisotropic past the old row check but well inside
1.08, and verified). The blind front-end's estimates on v1 (14 images it places, hint two
fields off, `-i` with the 4107–4119 indexes) were unchanged, score for score.

**2. Sigma-clipping the quad pairs.** The triangle path clipped its pattern pairs before
the fit; the quad path did not, so a few wrong quads in the winning vote cell could drag
the fit past the similarity check and the right position was abandoned. Clipping always
and clipping only after a refused fit (or a failed verification) gave identical statuses;
always was kept, as it costs no second fit-and-verify on positions that fail.

* Newly correct, true-centre hint (21): `cam_dss_c`, `s_dss_c_f64`, `s_dss_c_xisf`,
  `cam_tess_c`, `s_tess_c_flipy`, `obj_b68`, `obj_coalsack`, `obj_m16`, `rnd_012`,
  `rnd_049`, `lco_09`, `sv2_07`, `tess_09`, `tess_16`, `tess_27`, `tess_39`,
  `wide_shassa_04`/`05`/`08`, `wide_tess_09`/`10`. `s_tess_b_pincush` inexact (221″,
  synthetic pincushion).
* Newly correct, offset hint (15): `fov_10p0`, `ls_big`, `obj_coalsack`, `obj_veil_e`,
  `s_dss_c_flipx`, `s_dss_c_xisf`, `s_tess_b_xisf`, `s_tess_c_rot270`, `s_tess_c_trail`,
  `s_ztf_a_clouds`, `sv_sdssr`, `tess_22` (from inexact), `wide_dss_05`, `wide_tess_01`,
  and `s_tess_a_flipy` from a false positive.
* **New false positives: 12° TESS FFI crops only** — `tess_10`, `tess_18`, `tess_41`
  (true-centre hint) and `tess_03`, `tess_25`, `tess_32` (offset hint), all from no solve.
  They are the class §7.3 already lists (`tess_03`/`25`/`32` were false positives with
  the true-centre hint before): ~1000″ of SIP distortion, the right field found, a linear
  plate reported 2000–2500″ out at the corners. All have 30–54 verified stars of 500
  (under 11 %) at 1.0–1.3 px rms, but so do correct solves of crowded nebulae
  (`dens_carina`: 6 %, 1.0 px), so no threshold on those separates them; they need the
  structured-residual test or the distortion model of §7.3, not a tighter count.

**3. Database limit.** ASTAP caps the image stars it uses at the database's density times
the field's area. arcsec now does the same (plate-solving.md §10.2): with d80 it binds
below ~0.25° (45 images). Correct: `rnd_080` (973 detections, d80 holds 378 in the
window; capped at 344 it verifies 132), `ls2_02`, `ls2_21` (true-centre hint);
`rnd_080`, `ls2_02`, `ps1v2_07`, `stress_narrow` (offset hint). Nothing lost, no new false
positive. Verifying against every detection rather than the capped list gave identical
results.

**4. Sparse images.** Two parts, which pay mostly together (alone, +5 and +1 with the
true-centre hint; together +12):

* When the image has fewer stars than it may use, the catalogue read is denser than the
  image; the window's brightest `k = n · oversize² · long/short` catalogue stars, the
  image's density, now add their quads to the full-depth ones, provided the read holds at
  least 2.5 k stars. (Long over short: the window is square on the long side, so on a
  landscape frame `H/W` would halve the density.) Every image this solves had a read at
  least 3.7 times `k` (`lco_38`; the SkyMapper frames 4.7–7.3, `rnd_027` 18).
* The verified-star minimum is `min(30, max(10, ⌈0.15 n⌉))`, so unchanged from 194
  detections. Below 30 the plate scale must be within 10 % of the hint's and the
  star-level rms at most 0.5 px. The accepted relaxed solves have 16–28 stars, scales
  within 0.2 % of the hint and 0.16–0.46 px rms. Without the gate `ls2_25` (0.09°, Legacy
  Surveys) is accepted, with 12 stars at 1.36 times the hint's scale and 2.9 px rms, 5°
  from the truth: a false positive in both runs, and the only difference the gate makes.

Correct: `lco_20`, `lco_25`, `lco_38`, `rnd_027`, `rnd_044`, and seven of SkyMapper's
0.17° frames, `smss_01`/`02`/`03`/`05`/`06`/`10`/`15` (with the offset hint also `ls2_21`
and `ps1v2_14`). SkyMapper native now solves 15/15 and LCO 42/42 with the true-centre
hint.

**Speed.** v1 subset, `--jobs 1 --threads 1`, alternating builds, two rounds: total
62.6 s → 60.4 s and 62.6 s → 60.1 s, median 0.18 s → 0.18 s and 0.19 s → 0.18 s.
Full-spiral controls run alone (all cores), two runs each, 0.2.0 → now: `neg_hint_1`
26.6 → 26.7 s, `neg_hint_2` 12.7 → 12.9 s, `neg_hint_5` 6.3 → 6.2 s, `neg_shuffle_1`
0.7 → 0.7 s (6.8 → 6.3 s with `--threads 1`), but `neg_hint_3` 13.8 → 20.3 s: a
90-star image whose catalogue read is 3.4 times the image-density count, so the
density-matched quads are added, and matched, at each of its 6300 spiral positions. That
is the cost of the sparse-image change on a search that fails everywhere; a solve that
succeeds pays it only at the positions it visits. In the 12-job corpus runs the eight
`neg_hint_*` controls take 42–86 s (0.2.0: 45–89 s), well inside the harness's 300 s,
and no image timed out.

**Negative results.**

* *Reading the catalogue at the image's depth instead of adding its quads.* In the
  diagnosis (against 0.1.2), switching the read wholesale to the density-matched depth
  gained 11 images and lost 27, those whose faint detections are real and match the
  deeper catalogue. Hence the additive form.
* *Adding the density-matched quads whenever the image is short of stars.* The same
  statuses as with the 2.5× condition, but the extra quads are matched at every spiral
  position, and on searches that fail everywhere on images of 90–180 stars, where the
  read is only about twice `k`, they made the search 50 % slower: `neg_hint_2` (159
  stars) 6.7 s → 10.0 s with `-r 2 --threads 1`, and 12.7 s → 20.6 s for the whole
  search; matching is 58 % of that profile and scales with the catalogue quad count.
* *A count floor without the scale and rms gate.* A Poisson chance-match floor of 10 in
  the diagnosis let four tier-D controls verify; the relaxed count alone admits `ls2_25`
  (above).
* *Clipping only on failure* (after a refused fit, or also after a failed verification):
  the same statuses as clipping always, with a second fit-and-verify at every failing
  position.

### 7.9 Distortion (2026-10-02)

The solver-robustness changes (§7.8) found more fields, and among them three more
12° TESS frames reported with a linear plate 2000–2800″ out at the corners: six false
positives per run, all of the class §7.3 describes. The fix is a distortion model fitted
inside the solve ([plate-solving.md §10.3b](plate-solving.md#103b-distortion-pipelinedistortionrs)):
once a position verifies, the catalogue is re-read about the image centre, a polynomial
plate (up to cubic) is fitted by re-matching every catalogue star as it improves, and the
linear plate closest to it over the frame is reported; a strong distortion the model
cannot follow over the frame is refused. The WCS written stays linear unless `--sip` is
given. Whole corpus, both hints, `--auto-db`, 12 jobs; correct / false positive /
inexact over tiers A, B, C and S (596 images):

| | true-centre hint | offset hint | `--sip`, true-centre | `--sip`, offset |
|---|---|---|---|---|
| 0.2.0 | 504 / 4 / 11 | 459 / 4 / 39 | 512 / 5 / 2 | 465 / 4 / 33 |
| solver robustness (§7.8) | 540 / 7 / 12 | 493 / 6 / 38 | 549 / 8 / 2 | 499 / 6 / 32 |
| **+ distortion model** | **558 / 1 / 0** | **537 / 1 / 1** | **557 / 2 / 0** | **537 / 1 / 1** |

Tier D: no control solves in any run. v1 subset: 92 / 90 (centre / offset), as §7.8.
No image correct in 0.2.0 or after §7.8 is lost, in either run, with or without `--sip`.

| tier (centre hint) | 0.2.0 | §7.8 | now | | tier (offset hint) | 0.2.0 | §7.8 | now |
|---|---|---|---|---|---|---|---|---|
| A | 206/1/0 | 221/1/0 | 221/1/0 | | A | 199/1/0 | 209/1/0 | 210/1/0 |
| B | 206/3/4 | 220/6/4 | 230/0/0 | | B | 173/2/30 | 187/5/29 | 221/0/1 |
| C | 38/0/0 | 42/0/0 | 42/0/0 | | C | 39/0/0 | 42/0/0 | 42/0/0 |
| S | 54/0/7 | 57/0/8 | 65/0/0 | | S | 48/1/9 | 55/0/9 | 64/0/0 |

By source the changes are all in TESS (centre hint 9/3/4 → 13/6/4 → 23/0/0; offset
1/2/11 → 2/5/10 → 16/0/1), WISE (offset 1/0/19 → 20/0/0; centre 20/25 throughout, worst
corners 7.8–9.5″ → 7.0–7.8″), the synthetic distortion set (3 → 10 of 10 centre, 2 → 10
offset), and with the offset hint `ls2_06` and `wide_shassa_02` (the second chance). By
field size: 2.5–6° 17 → 22 → 27 (offset 9 → 13 → 24), 6–20° 23/4/2 → 25/7/3 → 34/1/0
(offset 20/3/3 → 24/5/3 → 31/1/0), 0.6–1.2° with the offset hint 141/0/22 → 144/0/22 →
166/0/0 (WISE). No band loses anything.

**Status changes against §7.8 (PR #10).** True-centre hint: from false positive to
correct `tess_03`, `tess_10`, `tess_18`, `tess_25`, `tess_32`, `tess_41`; from inexact to
correct `tess_12`, `tess_15`, `tess_21`, `tess_23` and the eight synthetic distortion
frames `s_des_a_worst`, `s_dss_a_pincush`, `s_dss_d_lens`, `s_ps1_b_lens`,
`s_tess_a_lens`, `s_tess_b_lens`, `s_tess_b_pincush`, `s_tess_c_pincush`. Offset hint:
false positive to correct `tess_03`, `tess_25`, `tess_32`, `tess_41`, `tess_43`; inexact
to correct 37 (the same eight synthetic frames and `s_ps1_a_lens`; `tess_07`/`12`/`13`/`15`/`19`/`23`/`24`/`33`/`44`; and 19 WISE
frames, `wise_02`–`25`); no solve to correct
`ls2_06`, `wide_shassa_02`. Nothing moves the other way. Against 0.2.0 the list is §7.8's
plus these: every 0.2.0 false positive and inexact solve except `wide_shassa_01` (both
hints) and `tess_34` (offset, inexact) is now correct.

**The ten false positives.** Of 0.2.0's four with the true-centre hint, `tess_03`,
`tess_25` and `tess_32` are correct (corners 990–1002″ against linear floors of
988–999″); `wide_shassa_01` remains, unchanged and not a distortion case: a 15° SHASSA
HiPS cutout, TAN, every star fitted to 0.6 px, but 16.8″ (0.64 px) from the truth at the
centre and 31″ at a corner against a 26″ limit — the survey's own astrometry or the
cutout, not the solve. Of 0.2.0's offset-hint four, `tess_41` and `tess_43` are correct,
`s_tess_a_flipy` was already fixed by §7.8, and `wide_shassa_01` is as above. The six
§7.8 added (`tess_10`, `tess_18`, `tess_41` centre; `tess_03`, `tess_25`, `tess_32`
offset) are all correct. With `--sip` one more remains in both 0.2.0 and §7.8 and now:
`wide_shassa_07`, a 47° SHASSA field whose model is not significant (F 3.7) and which
`fit_sip`, on its own F-test at 4, fits a cubic to that moves the corners from 29″ to
122″ out.

**`--sip`.** The model's pairs reach the corners, so `fit_sip` (unchanged) now fits
distorted frames over the whole field. Worst corners: 12° TESS 15–41″ (linear floor
~1000″), 5.8° TESS 8–51″, 2.9° TESS 6.5–55″, WISE 1.6–2.8″ (floor 5.5″), the synthetic
set 0.35–10″. Without `--sip` they sit at the linear floor: TESS 12° 990–1004″, WISE
7.0–7.8″.

**What was tried, and what decided the constants.**

* *Fitting from the spiral position's catalogue only.* With the offset hint a third of
  the frame has no catalogue stars, the pairs reach 3–6 of the 9 cells, and the model
  cannot be used: the 12° frames stayed false positives and the WISE frames inexact. The
  full-frame re-read fixed all of them (offset 517 → 536 correct, 6 → 1 false positives).
* *Without the quad seeds* `tess_03` stayed a false positive: its verified stars sit in
  one corner, and the quadratic grown from them never reaches the far third; the quad
  centroids, spread wider, pull the cubic there.
* *A quadratic gate on the cubic.* Testing the cubic only after a quadratic beat the
  linear plate kept every synthetic radial lens linear: radial distortion about the centre
  has no quadratic part.
* *Reporting the model whenever its F-test passes at 4.* Survey images with no
  distortion reach F = 4–24 on the final pairs (catalogue and centroid systematics), and
  the frame-wide plate then moved some corners: `sdss_d` 0.33″ → 3.2″, `ls2_27`
  0.10″ → 1.09″, `sv2_01` 0.86″ → 3.2″, a dozen others by 0.1–0.4″. None became wrong,
  but nothing gained either. At `MIN_REPORT_F` = 30 those keep their plate; the only
  distorted frame below 30 whose status depended on the model is `tess_34` with the
  offset hint (F = 7, inexact as in §7.8). ZTF has real sub-pixel distortion (F up to 47,
  floors 0.2–0.65″): `MIN_DEPARTURE_PX` = 1 leaves those plates alone.
* *The second chance without the full-frame read.* It turned `wide_shassa_02` (offset)
  into a false positive from a model fitted to part of the frame; with the re-read and the
  coverage requirement it is correct.
* *Refusal on the fraction or rms of verified stars* (§7.8's negative result) stays
  rejected; the refusal rule uses the model instead, and no corpus solve comes near it
  (largest wide-pair cubic F 21 among unused models 2 px or more from the plate; the rule
  needs 100 and 3 px).

**Speed.** The machine was shared with other benchmark runs throughout (load 50–170 on
24 cores), so wall-clock times are unreliable; these are CPU times (user + system) per
process. v1 subset, `--threads 1`, each image run with all three builds in turn:
0.2.0 69.9 s, §7.8 66.2 s, now 66.5 s in total (medians 0.21 / 0.22 / 0.22 s) in the
quieter round; 140.9 / 134.7 / 138.3 s in a round under twice the load. The model costs
one catalogue read and a few milliseconds of matching and small Cholesky fits per solve
(a few per cent of a 0.2 s solve in a profile, after replacing the Givens solver with normal
equations for it; with Givens it was about 10 %). Full-spiral no-solve controls (all
threads, CPU seconds, 0.2.0 / §7.8 / now): `neg_hint_1` 700 / 587 / 545, `neg_hint_2`
245 / 234 / 240, `neg_hint_3` 216 / 312 / 307, `neg_hint_5` 99 / 96 / 97,
`neg_shuffle_1` 9.7 / 10.2 / 9.9. A search that never verifies never fits the model,
and no control reached the second chance, so these differ only by noise; `neg_hint_3`'s
increase is §7.8's. No harness timeout (300 s) in any run; the slowest image in the
corpus runs was 186 s (`neg_orion_2massk`, offset hint, at the heaviest load; in the other
runs the slowest took 82–112 s).

**Still failing.** 19 TESS frames (2.2–5.9°) do not solve at all: §7.4's correspondence
problem (few detections have a Gaia counterpart at 21″/px in the TESS band), not
distortion — `tess_05` finds 614 quad references at the right position and keeps too few
consistent ones to fit. `tess_34` (2.2°) with the offset hint is inexact: its distortion
is real but weak on the pairs (F = 7).

### 7.10 Crowded and nebulous fields (2026-10-02)

§9 left 37 of the 596 images in tiers A/B/C/S unsolved with the true-centre hint, 22 of
them fields seiza solves: crowded and nebulous DSS fields, 47 Tuc, the Galactic Centre,
coarse DSS, TESS crops and WISE frames. Three changes and one refinement
([plate-solving.md §10.2, §10.3b, §10.3d](plate-solving.md#103d-the-catalogue-seeded-fallback-quadsseededrs)),
measured on the whole corpus with both hints against `main` (0.3.0 with the faster
failed search of §9.8, whose results are 0.3.0's, image for image), `--auto-db`, 12
jobs:

| | true-centre hint: correct / FP / inexact | offset hint: correct / FP / inexact | tier D FP | v1 (centre / offset) |
|---|---|---|---|---|
| `main` (§9.8) | 558 / 1 / 0 | 537 / 1 / 1 | 0 / 0 | 92 / 90 |
| **now** | **589 / 1 / 1** | **580 / 1 / 1** | **0 / 0** | **95 / 95** |

No image `main` solves correctly is lost in either run; the false positive in both is
`wide_shassa_01`, unchanged (§7.9). The one new inexact solve is `wise_11` (no solve
before; corners 10.8″ against a 5.5″ linear floor, 31 verified stars too few for the
distortion model to cover the frame).

| tier | centre: `main` → now | offset: `main` → now | | source | centre | offset |
|---|---|---|---|---|---|---|
| A | 221/1/0 → 230/1/0 | 210/1/0 → 225/1/0 | | hips2fits | 221 → 230 | 210 → 225 |
| B | 230/0/0 → 252/0/1 | 221/0/1 → 248/0/1 | | TESS FFI crops | 23 → **41** of 42 | 16 → **37** |
| C | 42/0/0 → 42/0/0 | 42/0/0 → 42/0/0 | | WISE | 20 → 23 of 25 | 20 → 23 |
| S | 65/0/0 → 65/0/0 | 64/0/0 → 65/0/0 | | Legacy Surveys | 39 → 40 | 37 → 40 |

By field size (true-centre hint): 2.5–6° 27 → **43 of 43**, 1.2–2.5° 79 → 86, 0.6–1.2°
168 → 173, below 0.15° 5 → 6, 6–20° 34 → 35; with the offset hint 2.5–6° 24 → 40,
1.2–2.5° 77 → 85, 0.6–1.2° 166 → 173, 0.3–0.6° 122 → 125, 6–20° 31 → 35. Accuracy of
what both solve: median worst corner 0.60″ → 0.57″ with the true-centre hint (59 images
more than 0.2″ closer to the truth, 9 further), 0.68″ → 0.54″ with the offset hint (137
closer, 10 further).

**The 22 fields seiza solved and `main` did not.** Now solved with the true-centre hint:
17 — `dens_carina`, `dens_vela`, `obj_47tuc`, `obj_sgra`, `type_sirius`, `rnd_054`,
`wide_dss_01`/`02`/`04`, `tess_05`/`17`/`20`/`31`/`40`, `wise_06`/`08`/`15`; with the
offset hint 16 (`dens_carina` and `obj_47tuc` already solved with it). Still not:
`dens_scutum`, `type_dbl_clus`, `type_m8`, `wise_16`, and `wise_11` (inexact).

**Status changes, true-centre hint** (all from no solve; 31 correct, 1 inexact):
`dens_carina`, `dens_vela`, `obj_47tuc`, `obj_sgra`, `rnd_054`, `type_sirius`,
`wide_dss_01`/`02`/`04`, `ls2_25`, `tess_04`/`05`/`06`/`08`/`11`/`17`/`20`/`26`/`28`/`29`/`30`/`31`/`35`/`36`/`38`/`40`/`42`/`43`,
`wise_06`/`08`/`15`; `wise_11` inexact.

**Status changes, offset hint** (43 correct, 1 inexact; from no solve unless noted):
`cam_dss_c`, `s_dss_c_f64`, `dens_vela`, `obj_m27`, `obj_sgra`, `rnd_012`, `rnd_049`,
`rnd_054`, `type_m31`, `type_sirius`, `stress_wide15`, `wide_dss_01`/`02`/`04`,
`wide_shassa_03`, `wide_tess_13`, `ls2_25`, `ls2_27`, `ls_p14`,
`tess_04`/`05`/`06`/`08`/`09`/`10`/`11`/`16`/`17`/`18`/`20`/`21`/`26`/`29`/`30`/`31`/`36`/`37`/`39`/`42`,
`tess_34` (from inexact), `wise_06`/`08`/`15`; `wise_11` inexact.

**Diagnosis.** For each failing group the detections were compared with the catalogue
through the truth WCS, and the catalogue's brightest stars in the frame were measured
at their true positions with the detection's own rules:

* *Crowded and nebulous DSS fields* (`dens_carina`/`scutum`/`vela`, `type_dbl_clus`,
  `type_m8`, `type_sirius`): 0–2 of the 30 brightest catalogue stars were among the
  500 detections, and only 5–12% of the detections had a catalogue counterpart. Of the
  60 brightest, the detection refused 9–34 as *too large* (saturated discs 20–75 pixels
  across: the 3σ isophote reaches the 14-pixel box) and 2–22 as noise (inside
  saturated nebulosity, where nothing can be measured); the detections were the faint
  stars between them, below the catalogue's depth.
* *Coarse DSS* (`wide_dss_01/02/04`, 4.8–9.6″/px) and `rnd_054`: the brightest
  detections did match (21–23 of the top 30), but 52–53 of the 60 brightest catalogue
  stars were refused as too large, so the image's 500 were mostly fainter than the
  catalogue's 500.
* *TESS FFI crops* (21″/px): 40 of the 60 brightest catalogue stars in `tess_05` were
  refused as too large (a core three pixels wide, wings 3σ above the sky past 16
  pixels), 7 as not a disc. Correspondence 30%: the quads did not match, and a
  catalogue-seeded search (below) without a significance test found wrong plates that
  verified 30–32 stars, chance in a frame of 500 stars on 384 × 384 pixels.
* *WISE*: a passband (W1–W4) far from Gaia's BP ranks stars differently (`wise_15`: 0
  of the 30 highest-SNR detections are catalogue stars), cosmic rays and streaks
  dominate single L1b exposures (`wise_11`: 470 of 500 detections, half-flux diameter
  1.4 px against 3.0 for its stars), saturated cores are masked (`wise_16`), or
  stray light leaves 49 detections (`wise_06`).
* *47 Tuc and the Galactic Centre* (`obj_47tuc`, `obj_sgra`): the catalogue's 500
  brightest come from a core the image does not resolve; correspondence 9%.

**What changed.**

1. *Bright stars measured in a larger box* (`detection::stars::measure_large`). A
   candidate the ordinary measurement refuses as too large or not a disc is re-centred
   on its brightest pixel and measured in a 32-pixel box at 5% of its peak; trails,
   galaxies and knots are refused by a second-moment test. These stars do not count
   towards the `-s` that ends the detection cascade. Alone (with nothing else changed)
   this took the true-centre count to 583 and the offset count to 562, and brought the
   only new false positives of the work (`tess_17`, a wrong field verifying 86 stars at
   3.4 px rms; `wide_shassa_03`, 1.2 px at a corner with the offset hint), which
   changes 2 and 4 removed.
2. *Significance of a verification* (`MIN_SIGNIFICANCE` = 4, and an rms within the last
   match radius). Every correct solve verifies at least 7.9 times the matches expected
   by chance at the frame's density (`tess_40`, 121 against 15.3); the wrong TESS plates
   1.5–2.7 times; two wrong plates (`wise_15`, `obj_m27`, both from the fallback) had
   4.4 px rms after the 2 px pass, where no correct solve exceeds 1.35 px. Without it
   the fallback reported four wrong TESS fields and `wise_15` with the true-centre hint,
   and eight TESS fields, `wise_15` and `obj_m27` with the offset hint.
3. *Catalogue-seeded fallback* (`quads::seeded`), after a spiral that found nothing.
   Eight solves with the true-centre hint come from it (`ls2_25`, `obj_sgra`,
   `tess_05`, `tess_26`, `wide_dss_02`, `wise_08`, `wise_15`, and `wise_11`, inexact),
   twenty with the offset hint; the rest come from the detection change and the spiral.
4. *Linear refit on the model's pairs* when the distortion model's full-frame match has
   1.5 times the verified stars: with the offset hint it made `wide_shassa_03` (a
   1.2-pixel corner) and `type_m45` (5.7″, after change 1) correct and `tess_34` correct
   from inexact, and moved the median worst corner from 0.70″ to 0.59″ on that run;
   with the true-centre hint it changes three plates, no status.

**Speed.** CPU time (user + system) per process, `--threads 1`, alternating the two
builds image by image, two rounds and the lower of each, load average 0.3–1.6, against
`main` after §9.8:

| | `main` | now |
|---|---|---|
| v1 subset (103), total | 30.8 s | 31.6 s, 3 more solved |
| ... the 92 both solve: total / median | 20.7 s / 0.114 s | 21.7 s / 0.116 s (median per-image ratio 1.03) |
| the 8 `neg_hint_*` controls (full spiral, `-r 10`) | 249.3 s | 250.5 s |
| the 27 other tier-D controls | 57.3 s | 58.7 s |

The solved images that slow down are the wide fields full of saturated stars, where the
discs measured as large stars replace the pieces of them the cascade used to count and
it runs a second level (plate-solving.md §10.2): `stress_wide15` 0.27 → 0.52 s,
`fov_10p0` 0.40 → 0.46 s, `fov_3p00`, `fov_5p00` +10%. Every other solve is within a few
milliseconds; the out-of-line measurement and the ring walk took the typical overhead
from 18% (`dec000`, 99 → 118 ms in the first build) to under 4%. A failed search pays
the fallback once: 0.1–0.3 s of CPU (`neg_shuffle_6` 0.21 → 0.46 s, `neg_shuffle_5`
0.24 → 0.47 s), under 1% of a full spiral. With all threads the fallback is still one
thread, so its share of a failure's wall time is larger on a many-core machine.

**Negative results.**

* *Counting the large stars towards the cascade's `-s`.* `wide_dss_05` stopped at the
  first level with 525 stars (157 of them large) instead of going on to the gridded
  level that finds 10 522, and lost its solve.
* *A separate bright-star pass after the cascade*, leaving the cascade's own stars
  exactly as before and replacing only the off-centre pieces of a source: the same speed
  as `main` on solved images (the large measurement inline costs a second cascade level
  on a few wide fields, `stress_wide15` 0.27 → 0.52 s), but 587 / 577 correct rather
  than 589 / 580: it loses `tess_38` and `tess_43`, and with the offset hint `obj_m27`,
  `tess_08`, `tess_21` and `type_m31` (gaining `tess_28`).
* *A 64-pixel box for large stars* (`type_dbl_clus` and `type_m8` have bright stars
  50–75 pixels across): neither solved, and `dens_scutum` and the offset `dens_carina`
  were lost.
* *An elongation filter on every detection* (eigenvalue ratio of the second moments
  above half the peak, ≤ 4), measured on an earlier build of this work: lost the two
  trailed tier-S frames `s_2mass_a_trail` and `s_sm_a_trail12` (all their stars are
  elongated) and made `s_des_a_worst` inexact, with nothing gained; at 2.5 it removed
  half of `wise_11`'s detections (undersampled stars measure as elongated) and lost it.
  The large-star measurement keeps its own test, where trails and galaxy cores are the
  risk.
* *Spatially diverse selection* (the brightest detections capped per cell of a 4 × 4
  grid, at twice the mean), on the same build: no status change in either run; 18
  corners moved, by at most 1.4″. `dens_scutum`, the one crowded field it was aimed at,
  solves in some variants of this work and not in others (it verifies 34 stars at the
  first position when it does): marginal, not a selection effect.
* *Verifying the fallback's candidates against the 2000 deep detections* instead of the
  500 the spiral uses: no new solve (`dens_scutum`, `type_m8`, `type_dbl_clus`, `wise_16`
  stay unsolved); not kept.
* *A larger fallback budget* (10⁸ and 2·10⁸ instead of 4·10⁷): no further solve; the
  remaining failures' brightest catalogue stars are not among the detections at all
  (`type_m8` 7 of the 100 brightest detected, `type_dbl_clus` 9), so their seed quads
  cannot be found. Indexing only the 500 or 1000 brightest detections for the fallback's
  image pairs: a wash (each solved one case the other did not).
* *Trying the fallback's short quads first* (fewer image pairs per quad, so more quads
  within the budget): `tess_05` 1.8·10⁷ → 2.5·10⁶ units, `obj_sgra` 1.8·10⁷ → 7·10⁶, but
  `tess_26` 2.9·10⁶ → 2.9·10⁷ and `wide_dss_02` lost (a quad of near neighbours misses
  more often in a coarse frame). The catalogue's brightness order was kept.
* *A census without the significance floor*, or excluding the quad's own four stars
  from it: the first sent thousands of wrong TESS transforms to the verification (and
  through it, before change 2); the second lost `wise_15`.
* *Fewer stars, or more* (`-s` 40–2000, the cheapest way to change both depths): some
  targets solve at 150 or 1000–2000 stars and others lose; no single setting helps
  overall, which is what pointed to the measurement and the fallback instead.

## 8. Licensing and attribution

The corpus is deliberately not committed: `resources/` is gitignored and the manifests
are the artefact — each user downloads the images from the providers for their own use.
The data are public, but several providers attach conditions, recorded here. If any of
this is ever published or redistributed, credit as the providers ask and check each
one's current terms first.

| Source | Terms as found (2026-09) | Attribution |
|---|---|---|
| CDS hips2fits / MocServer | free service | CDS, Strasbourg; plus the underlying survey |
| DSS (DSS1, DSS2) | STScI/AURA; the plates are copyright Caltech / ROE / AAO — use for research and education with acknowledgement | the standard DSS acknowledgement |
| 2MASS, WISE / allWISE / unWISE, ZTF via IRSA | public | UMass & IPAC/Caltech (2MASS); JPL/UCLA (WISE); ZTF (NSF, Caltech) |
| Pan-STARRS1, GALEX, TESS via MAST / AWS open data | public | PS1 Surveys; GALEX; the TESS mission acknowledgement |
| SDSS DR17 | public | SDSS-IV acknowledgement |
| Legacy Surveys DR10, DES DR2, DECaPS | public | the DESI Legacy Imaging Surveys / DES acknowledgements (NOIRLab, DOE, NSF) |
| SkyMapper DR4 | public release | SkyMapper acknowledgement (ANU) |
| VISTA VHS, UKIDSS LAS (WFAU HiPS) | public ESO / WFAU releases | VHS / UKIDSS acknowledgements |
| DENIS, IPHAS | public | survey acknowledgements |
| SHASSA | "by courtesy of Swarthmore College"; NSF-funded | the SHASSA acknowledgement page |
| TESS 2-yr HiPS | NASA/MIT/TESS and Ethan Kruse (USRA) | as stated |
| SkyView | NASA/GSFC HEASARC | "The SkyView virtual observatory" |
| Las Cumbres Observatory | public archive frames (after the proprietary period), anonymous API | "This work makes use of observations from the Las Cumbres Observatory global telescope network." |
| HST / MAST | STScI | observation programme ID |
| Tier S / generated tier D | derived locally from the above | as the parent |

Two policy notes for whoever maintains this:

* Several archives' `robots.txt` disallow crawlers on the very endpoints their
  documentation tells scripts to use (SkyView `cgi`, PS1 `cgi-bin`, `data.sdss.org`, and
  the LCO archive API). The corpus treats a documented, rate-limited API call as
  permitted and a crawl as not, and fetches nothing that is not listed in the manifest.
  **nova.astrometry.net is the exception**: its robots file names AI crawlers and its
  pages sit behind a human check, so it was not used at all.
* The Mellinger mosaic ("all rights reserved") was excluded on licence grounds.

---

## 9. Results: arcsec vs ASTAP vs seiza

Measured 2026-10-02. [seiza](https://github.com/theatrus/seiza) is a Rust plate solver
(Apache-2.0) with a hinted solver, a blind solver against a prebuilt whole-sky index,
and compatibility modes for N.I.N.A. (ASTAP's command line) and Siril (astrometry.net's
`solve-field`). This section runs it on the same corpus, with the same hints and the
same scoring, next to arcsec and ASTAP. As with §7.5, the corpus was built by this
project to cover arcsec's weak spots, not to be typical of anyone's frames: read the
numbers as where the three differ, not as a ranking.

The arcsec numbers are for `main` after the solver-robustness and distortion work of
§7.8–7.9 (d30e6e7). The comparison was first made against the 0.2.0 release; those
arcsec figures are kept, dated, at the end of §9.3 and §9.5. The ASTAP and seiza runs
were not repeated: neither program changed.

### 9.1 What was run

| | arcsec | ASTAP | seiza |
|---|---|---|---|
| version | `main` d30e6e7 (after 0.2.0); blind: PR #13 `custom-index` 4f0760a, unreleased | `astap_cli` CLI-2026.07.30 | 0.18.20 (`main`, 66bcde6), built with Rust 1.97.1 |
| licence | MIT | MPL-2.0 | Apache-2.0 |
| catalogue | d80 (0.15°–6°), g05 (3°–20°), w08 (20°–80°), chosen by field size (`--auto-db`) | the same directory; no `-D`, so ASTAP chooses by field size itself (§9.2) | `stars-deep-gaia17.bin`, Gaia DR3 to G 17, 154 M stars, 1.54 GB |
| blind index | its own index built from d80 (PR #13): 287 MB (fields 0.3°–30°) or 698 MB (0.15°–30°); astrometry.net 4107–4119 for the historical row (§9.6) | – | `blind-gaia16.idx`, Gaia DR3 to G 16, 1.63 GB |
| how it is driven | `-f -d --ra --spd --fov -r 5 -o` | `-f -d -ra -spd -fov -r 5 -o` | one-shot `seiza worker`, a hinted request (§9.2) |

seiza's data came from its own downloader (`seiza download-data prebuilt --file …`,
bundle `catalog-bundle-v4-2026-10-01`, SHA-256 verified): the two files above plus
`stars-gaia.bin` (G ≤ 15, 367 MB) and `stars-lite-tycho2.bin` (25 MB), 3.6 GB in all,
about 2.7 GB to download (zstd). The bundle states no licence of its own; the stars are
Gaia DR3 (ESA/Gaia/DPAC, CC BY-SA 3.0 IGO) and Tycho-2. Given a directory, seiza uses the
deepest catalogue in it, so the main runs use the G ≤ 17 file; §9.4 repeats the corpus
with the G ≤ 15 file its documentation recommends as the general default.

### 9.2 How each solver was driven, and why

* **The same hint for all three.** `benchmark.py` computes the hint once per image (the
  true centre, or 0.3 fields off with `--offset-hint 0.3`) and gives it to each solver in
  its own units: arcsec and ASTAP get RA in hours, south-pole distance and the FOV of the
  image height; seiza gets RA/Dec in degrees and the pixel scale that FOV implies, with
  its default ±20 % scale tolerance. All three get `-r 5`.
* **seiza's native hinted solver, not its ASTAP mode.** `seiza solve` prints a rounded
  summary (centre to 1e-5°, corners to 1e-4°, rotation to 0.01°), which is not enough to
  score corners on a 15° field, so the harness sends the same request through a one-shot
  `seiza worker` process: the same `seiza::solve::solve` call, with the WCS returned in
  full as JSON. On `fov_1p50` the two give the same solution and take 0.44 s and 0.52 s.
  The worker is asked for 200 stars, which is what `seiza solve` requests. seiza's
  ASTAP-compatible mode is *not* like for like: it clamps the radius to 0.5°–3°, uses
  at least 200 stars, ignores `-d`/`-D`/`-z`, and when the hinted solve fails it runs a
  full blind solve over 0.1–20″/px. §9.4 runs it separately (`--seiza-mode astap`).
* **ASTAP picks its own database.** §7.5 ran ASTAP with `-D d80` only, which left it
  nothing above 6°. Here it gets no `-D`, and `astap_cli` then chooses for itself
  (G05 for a 10° field, D80 for 1.5°, checked with `-progress`). In practice this changes
  little: ASTAP solves 3 of the 36 fields between 6° and 20° and none above 20°.
* **Compressed inputs.** seiza recognises images by file extension and does not open
  `.fits.fz` (the 42 LCO frames, tier C) or `.fits.gz` (`s_dss_b_gz`): "The file
  extension `.fz` was not recognized as an image format". Those 43 are a format gap, not a
  search failure, so the main runs use a mirror of the corpus in which those files are
  decompressed (astropy, same pixels and header) for all three solvers. On the files as
  delivered, seiza solves 503 rather than 546 and ASTAP 303 rather than 304 (it does not
  read the `.gz` file either); arcsec reads both and its result is the same (558).
* **Scoring is identical.** Every solution is scored at the centre and four corners
  (§4), with the 5″-or-one-pixel threshold, the linear floor for SIP/TPV truth and the
  INEXACT class (§4.2), for every solver. A tier-D control that is answered at all counts
  as a false positive, as for arcsec.
* **Parallelism.** The accuracy runs used `--jobs 12` on a machine that another job was
  also loading (load average up to 150), so their times are not used; §9.5 has the
  timing runs.

Commands (the corpus mirror is `resources/corpus` with the 43 compressed files
decompressed in place of the originals):

```bash
scripts/benchmark.py --corpus --images <mirror> --auto-db --astap ~/astap_cli \
    --seiza <seiza> --seiza-data ~/.local/share/seiza-data --jobs 12 --by tier,source,fov
scripts/benchmark.py ... --offset-hint 0.3
scripts/benchmark.py --auto-db --astap ~/astap_cli --seiza <seiza> --seiza-data ...   # v1
```

### 9.3 Accuracy and solve rate

**v1 (103 images), true-centre hint:**

| | arcsec | ASTAP | seiza |
|---|---|---|---|
| correct (of 98) | 92 | 47 | **97** |
| false positives | 0 | 0 | 0 |
| no solve | 6 | 51 | 1 |
| tier D (5) false positives | 0 | 0 | 0 |
| median centre / corner, the 47 all three solve | 0.68″ / 0.86″ | 0.69″ / 0.94″ | 0.79″ / 1.13″ |

seiza solves everything arcsec does except `ls_p14` (a 0.14° Legacy Survey field), plus
six it does not: `dens_carina`, `dens_scutum`, `dens_vela`, `type_dbl_clus`, `type_m8`,
`type_sirius` — the crowded and bright-star fields of §6.4. With the hint 0.3 fields off:
arcsec 86, ASTAP 46, seiza 96 with one false positive (`stress_wide15`, 15°, corner 58″,
5 px). (v1 is where main changed least: arcsec's v1 counts are the same as 0.2.0's.)

**Expanded corpus (635 entries), true-centre hint:**

| | arcsec | ASTAP | seiza |
|---|---|---|---|
| correct (of 596 in tiers A/B/C/S) | **558** | 304 | 546 |
| false positives | **1** | 2 | 22 |
| inexact (distorted truth) | **0** | 9 | 24 |
| no solve | 37 | 281 | **4** |
| tier D (35 must-fail) false positives | 0 | 0 | 3 |
| wrong-FOV stress entries (4, `expect=any`) solved correctly | 4 | 0 | 0 |
| median centre / corner, the 302 all three solve | **0.25″ / 0.49″** | 0.26″ / 0.71″ | 0.38″ / 0.77″ |

| Tier | n | arcsec | ASTAP | seiza |
|---|---|---|---|---|
| A — synthetic cutouts | 234 | 221 (1 FP) | 131 (2 FP) | 223 (11 FP) |
| B — survey pixels | 255 | 230 | 97 (1 inexact) | 224 (8 FP, 19 inexact) |
| C — LCO frames | 42 | 42 | 41 | 42 |
| S — simulated camera artefacts | 65 | 65 | 35 (8 inexact) | 57 (3 FP, 5 inexact) |

| FOV (long side) | n | arcsec | ASTAP | seiza |
|---|---|---|---|---|
| < 0.15° | 6 | 5 | 0 | 3 |
| 0.15–0.3° | 113 | 113 | 39 | 112 |
| 0.3–0.6° | 126 | 125 | 97 | 125 |
| 0.6–1.2° | 180 | 171 | 111 (3 inexact) | 177 |
| 1.2–2.5° | 88 | 79 | 53 (1 FP, 2 inexact) | 86 (2 inexact) |
| 2.5–6° | 43 | 27 | 1 (1 inexact) | 16 (9 FP, 18 inexact) |
| 6–20° | 36 | 34 (1 FP) | 3 (1 FP, 3 inexact) | 20 (12 FP, 4 inexact) |
| > 20° | 8 | 8 | 0 | 7 (1 FP) |

| Source | n | arcsec | ASTAP | seiza |
|---|---|---|---|---|
| hips2fits (22 surveys) | 234 | 221 (1 FP) | 131 (2 FP) | 223 (11 FP) |
| ZTF | 43 | 43 | 42 | 43 |
| SDSS | 33 | 33 | 5 | 33 |
| Pan-STARRS1 | 20 | 20 | 9 | 20 |
| SkyView | 37 | 37 | 27 | 37 |
| Legacy Surveys | 40 | 39 | 0 | 36 |
| LCO (tier C) | 42 | 42 | 41 | 42 |
| WISE L1b | 25 | 20 | 0 (1 inexact) | 25 |
| SkyMapper native, 0.17° | 15 | 15 | 14 | 15 |
| TESS FFI crops | 42 | 23 | 0 | 15 (8 FP, 19 inexact) |
| tier S | 65 | 65 | 35 (8 inexact) | 57 (3 FP, 5 inexact) |

**Accuracy, on the images all three solve** (median centre / worst corner):

| Tier | n | arcsec | ASTAP | seiza |
|---|---|---|---|---|
| A | 129 | 0.33″ / 0.65″ | 0.42″ / 0.79″ | 0.54″ / 0.93″ |
| B | 97 | 0.12″ / 0.34″ | 0.13″ / 0.68″ | 0.31″ / 0.65″ |
| C | 41 | 0.23″ / 0.65″ | 0.22″ / 0.73″ | 0.19″ / 0.72″ |
| S | 35 | 0.29″ / 0.55″ | 0.27″ / 0.56″ | 0.33″ / 0.86″ |

Pairs that arcsec and seiza both solve, by survey (centre / corner): ZTF 0.10″ / 0.35″
against 0.43″ / 0.68″, Legacy Surveys 0.05″ / 0.10″ against 0.14″ / 0.40″, SDSS 0.20″ /
0.33″ against 0.26″ / 0.54″, SkyMapper 0.07″ / 0.26″ against 0.14″ / 0.35″, LCO 0.23″ /
0.65″ against 0.19″ / 0.72″. Where the truth is good to better than 0.3″ (§7.6), arcsec's
plates are closer to it, and on tier A as well; on the LCO frames, where the pipeline
truth is the limit, the two are level. seiza fits its plate to the 200 brightest
detections on an 8-bit stretched copy of the image, arcsec to up to 500 stars measured on
the linear pixels, which is the likely difference; we did not test it directly.

**Who solves what** (correct only):

* arcsec and seiza both: 524. arcsec only: 34. seiza only: 22. ASTAP solves nothing that
  neither of the others does; it solves one image arcsec does not (`rnd_054`) and one
  seiza does not (`wide_tess_08`, which seiza answers 1.2 px out).
* seiza only (22): the crowded and bright-star fields that remain arcsec's standing
  failures (`dens_carina`, `dens_scutum`, `dens_vela`, `obj_47tuc`, `obj_sgra`,
  `type_dbl_clus`, `type_m8`, `type_sirius`); `rnd_054`; three coarse DSS fields
  (`wide_dss_01/02/04`); five TESS FFI crops (`tess_05/17/20/31/40`); five WISE frames
  (`wise_06/08/11/15/16`). *Update (§7.10):* arcsec now solves 17 of them with the
  true-centre hint (16 with the offset hint); `dens_scutum`, `type_dbl_clus`, `type_m8`
  and `wise_16` still fail, and `wise_11` is inexact (10.8″ at the corners against a
  5.5″ linear floor). arcsec's count with the true-centre hint is 589, against seiza's
  546; the seiza-only list is now those five.
* arcsec only (34): ten coarse TESS fields that seiza answers just outside the
  threshold (`cam_tess_c`, `wide_tess_01/02/04/05/07/08/09/12/13`) and three tier-S
  variants of them (`s_tess_c_flipx/flipy/rot270`); thirteen distorted TESS FFI crops
  that seiza answers wrong or inexact (`tess_03/09/10/12/13/18/21/23/25/32/33/41/44`) and
  five tier-S lens and pincushion cases (`s_dss_d_lens`, `s_tess_a_lens`,
  `s_tess_b_lens`, `s_tess_b_pincush`, `s_tess_c_pincush`) — main's distortion model
  (§7.9) at work; and `ls2_07`, `ls2_21`, `ls_p14`, which seiza does not solve.

**False positives, by kind** (true-centre hint):

| Kind | arcsec | ASTAP | seiza |
|---|---|---|---|
| right field, corners 1–3 px out (coarse 12–29″/px TESS and SHASSA HiPS) | 1 (`wide_shassa_01`, 1.2 px) | – | 14 |
| distorted TESS FFI, one linear plate over 500–2800″ of distortion | – | – | 8 |
| plate wrong at the corners (centre within 14″) | – | 2 (`rnd_074`, `wide_tess_06`) | – |
| a must-fail control answered | – | – | 3 |

The three tier-D answers are different in kind. `neg_fake_wide` (400 random stars,
10°) is a wrong field, 5.4° off. `neg_shuffle_6` (`cam_tess_a` with 128-px blocks
permuted) comes out at the right centre with 12″ corners: a solution for a scrambled
image. `neg_orion_2massk` is a 0.04° 2MASS K crop that is a control only because it lies
below every ASTAP database; seiza's G ≤ 17 catalogue reaches it and the answer is
correct (0.07″ / 0.33″).

**Offset hint (0.3 fields), expanded corpus:**

| | arcsec | ASTAP | seiza |
|---|---|---|---|
| correct | **537** | 280 | 531 |
| false positives | **1** | 2 | 25 |
| inexact | 1 | 6 | 32 |
| no solve | 57 | 308 | 8 |
| tier D false positives | 0 | 0 | 2 |
| median centre / corner, the 279 all three solve | 0.27″ / 0.60″ | 0.26″ / 0.69″ | 0.39″ / 0.80″ |

arcsec and seiza both solve 502; arcsec only 35 (again mostly coarse and distorted TESS,
plus `wise_03`, `wise_23`, `sdss2_18`, `sdss2_19`), seiza only 29. arcsec's one false
positive is `wide_shassa_01` again. Nine of seiza's 25 + 2 are worse than near misses:
seven wide fields (`stress_wide15`, `wide_shassa_02/04/05/06`, `wide_tess_12/13`,
15°–42°) whose offset puts the true centre outside `-r 5`, answered with corners 3–11 px
out, and two wrong fields (`wide_shassa_08`, 35° away, and `neg_shuffle_6`, 3.8°). ASTAP
returns nothing for all nine; arcsec solves `wide_tess_12` correctly and returns nothing
for the rest.

**Against 0.2.0** (the first run of this comparison, 2026-10-02): arcsec 0.2.0 scored
504 correct / 4 false positives / 11 inexact with the true-centre hint and 459 / 4 / 39
with the offset hint, against seiza's 546 / 22 / 24 and 531 / 25 / 32; at that point
seiza solved 57 images arcsec did not and arcsec 15 that seiza did not. The gains since
(§7.8–7.9) are mostly the fields where seiza was ahead (SkyMapper, the LCO misses,
random galactic-plane fields, tier S) and the distorted TESS frames.

### 9.4 Variations on seiza

| Run (expanded corpus, true centre) | correct | FP | inexact | no solve | tier D FP |
|---|---|---|---|---|---|
| seiza hinted, G ≤ 17 catalogue (§9.3) | 546 | 22 | 24 | 4 | 3 |
| seiza hinted, G ≤ 15 `stars-gaia.bin` | 484 | 22 | 24 | 66 | 3 |
| seiza hinted, files as delivered (`.fz`/`.gz` unread) | 503 | 22 | 24 | 47 | 3 |
| seiza ASTAP-compatible mode (blind fallback) | 544 | 21 | 27 | 4 | 10 |

* **Catalogue depth matters below 0.3°.** With the G ≤ 15 catalogue seiza solves 54 of
  the 113 fields of 0.15°–0.3° rather than 112, and 1 of 6 below 0.15°; above 0.3° the
  counts differ by at most two per band. Small fields need the 1.5 GB file.
* **ASTAP mode** solves about the same images, but its blind fallback answers seven more
  tier-D controls. Six are `neg_hint_*` — real fields given a hint 30°–60° away, which
  the blind search finds correctly; they are failures only by the rule that a hinted
  solver must not wander that far. One is a wrong answer for pure noise (`neg_noise_c`,
  2.6° off), and `wise_16` becomes a wrong field (1.8° off) where the hinted solver
  found it. A run in this mode took 60–90 s on some images (the blind fallback).

### 9.5 Speed

(For arcsec, superseded by §9.8, which re-times 0.3.0 and the faster failed search
after it; the ASTAP and seiza numbers here stand.)

Timing runs used `--jobs 1`, so one solver ran at a time. Two sets: v1 (103 images) and a
stratified 98-image subset of the expanded corpus (4–16 per FOV band across the sources,
plus 10 tier-D controls, chosen by a fixed rule before any timing). Two thread settings:
each solver's default, and one thread (`--threads 1` for arcsec, `RAYON_NUM_THREADS=1`
for seiza — its only thread pool — and every solver pinned to one core with
`taskset -c 5`; ASTAP has no thread option). Two rounds of each.

ASTAP and seiza were timed together with arcsec 0.2.0, each image going through all
three in turn, the order alternating between rounds. arcsec `main` was timed afterwards
on its own, with the same protocol on the same quiet machine; 0.2.0 and main take the
same time on the images whose result did not change, which cross-checks the two
sessions. The machine was shared with another job that came and went: each run waited
for a load average under 4, and a sampler recorded the CPU used by processes outside the
run's own process tree every 5 s; a run with a mean above 1 core or a peak above 4 was
discarded and repeated (three were). The runs kept had 0.2–0.3 cores of outside load on
average and at most 1.3. Totals agreed between rounds to within 2 % (arcsec main: within
0.5 %), and solve counts were identical. Times are wall-clock per process, start to exit,
warm page cache, averaged over the two rounds.

**Per image, median / mean / p90 / total (s):**

| Set, threads | | arcsec (main) | ASTAP | seiza |
|---|---|---|---|---|
| v1, default | all 103 | 0.15 / 0.31 / 0.63 / 31.5 | 0.55 / 0.87 / 1.31 / 89.8 | 0.18 / 0.24 / 0.36 / 24.7 |
| | solved (n = 92 / 47 / 97) | 0.15 / 0.24 / 0.38 / 22.1 | 0.15 / 0.27 / 0.77 / 12.9 | 0.18 / 0.21 / 0.34 / 20.7 |
| | not solved (n = 11 / 56 / 6) | 0.70 / 0.85 / 2.16 / 9.4 | 0.93 / 1.37 / 2.50 / 76.9 | 0.76 / 0.66 / 1.59 / 3.9 |
| v1, one thread | all 103 | 0.18 / 0.58 / 0.68 / 59.6 | 0.56 / 0.87 / 1.32 / 90.0 | 0.19 / 0.26 / 0.43 / 26.8 |
| | solved | 0.18 / 0.28 / 0.57 / 25.6 | 0.15 / 0.28 / 0.76 / 13.1 | 0.19 / 0.24 / 0.39 / 23.1 |
| | not solved | 3.26 / 3.10 / 7.34 / 34.1 | 0.92 / 1.37 / 2.49 / 77.0 | 0.76 / 0.63 / 1.49 / 3.8 |
| subset, default | all 98 | 0.14 / 0.63 / 0.40 / 62.0 | 0.34 / 1.26 / 1.62 / 123.8 | 0.17 / 0.39 / 0.42 / 38.6 |
| | solved (n = 83 / 33 / 88) | 0.14 / 0.18 / 0.34 / 15.2 | 0.16 / 0.20 / 0.37 / 6.5 | 0.17 / 0.18 / 0.27 / 16.0 |
| | not solved (n = 15 / 65 / 10) | 0.23 / 3.12 / 13.0 / 46.9 | 0.51 / 1.80 / 2.94 / 117.3 | 1.33 / 2.26 / 13.0 / 22.6 |
| subset, one thread | all 98 | 0.16 / 5.66 / 0.74 / 554.7 | 0.34 / 1.24 / 1.62 / 121.9 | 0.17 / 0.39 / 0.47 / 38.2 |
| | solved | 0.16 / 0.20 / 0.40 / 16.8 | 0.16 / 0.20 / 0.38 / 6.4 | 0.16 / 0.18 / 0.26 / 15.9 |
| | not solved | 1.11 / 35.9 / 158 / 537.8 | 0.49 / 1.78 / 2.97 / 115.5 | 1.32 / 2.24 / 12.9 / 22.4 |

"Solved" means the solver returned an answer, right or wrong; the counts are arcsec /
ASTAP / seiza.

**On the same images.** Where arcsec and seiza both solve correctly (91 v1 images, 74
subset images), arcsec's time is 0.85–0.96 of seiza's (geometric mean of the per-image
ratio: 0.92 and 0.85 at default threads, 0.96 and 0.92 on one thread); medians 0.155 s
against 0.174 s on v1 and 0.137 s against 0.167 s on the subset. On the 47 v1 images all
three solve, the medians are 0.155 s, 0.147 s and 0.177 s (arcsec, ASTAP, seiza) at
default threads and 0.182 s, 0.148 s and 0.188 s on one thread. A solved image costs
about the same in all three.

**The totals are decided by the failures, and arcsec's failed search is the expensive
one.** With all cores, arcsec gives up in a median 0.2–0.7 s, but its spiral runs to the
edge of the radius before it does, and the tier-D `neg_hint_*` controls (real fields with
the hint 30°–60° away and `-r 10`) take 6–26 s each. On one thread the same three take
72 s, 158 s and over 300 s (`neg_hint_1`, stopped by the timeout): 532 of the subset's
555 s. seiza gives up sooner (its slowest failure, `neg_fake_300`, took 13 s; its others
1–3 s) and fails least often, so it has the lowest totals: 25 s against arcsec's 32 s
and ASTAP's 90 s on v1, and 39 s against 62 s and 124 s on the subset. ASTAP spends 77
of its 90 s on v1 in the 56 images it does not solve, and its failures on the widest
fields (the SHASSA mosaics, 20°–50°) take 6–17 s each.

**Threads.** seiza and ASTAP take the same time on one pinned core as with the whole
machine (seiza's hinted solve is close to serial; ASTAP's CLI did not use more than one
core here). arcsec uses every core for star detection and the spiral search; on one
thread its v1 total nearly doubles (31.5 → 59.6 s), almost all of it in failed searches
(9.4 → 34.1 s).

**By field of view** (median / total per band, default threads, subset):

| FOV | n | arcsec (main) | ASTAP | seiza |
|---|---|---|---|---|
| < 0.15° | 4 | 0.05 / 0.2 | 0.51 / 21.6 | 0.17 / 1.2 |
| 0.15–0.3° | 13 | 0.09 / 1.3 | 0.55 / 12.3 | 0.15 / 2.4 |
| 0.3–0.6° | 15 | 0.14 / 2.2 | 0.25 / 8.9 | 0.16 / 2.6 |
| 0.6–1.2° | 16 | 0.14 / 3.1 | 0.22 / 4.8 | 0.18 / 2.8 |
| 1.2–2.5° | 14 | 0.18 / 2.6 | 0.24 / 3.3 | 0.17 / 2.5 |
| 2.5–6° | 12 | 0.32 / 3.4 | 0.58 / 8.9 | 0.15 / 2.3 |
| 6–20° | 10 | 0.20 / 3.0 | 0.62 / 8.6 | 0.19 / 1.9 |
| > 20° | 4 | 0.13 / 0.5 | 15.6 / 50.5 | 0.24 / 1.0 |
| tier D | 10 | 0.15 / 45.7 | 0.22 / 4.8 | 1.33 / 21.9 |

**Fixed costs.** seiza's worker reports its own timings: of a typical 0.17 s solve, 28 ms
is reading the image, 15–35 ms star detection and about 115 ms the solve; process start,
catalogue open (memory-mapped) and the JSON exchange add 7–8 ms. `neg_m42_tiny`
(400 × 400 pixels, which every solver refuses at once) gives the floor for a whole
process: arcsec 15–25 ms, seiza 8–18 ms, ASTAP 60 ms. Times are with a warm page cache;
the first solve after a reboot, which has to read the catalogue from disk, was not
measured.

**Against 0.2.0** (2026-10-02, timed in the same session as ASTAP and seiza): v1 totals
32.6 s (default threads) and 62.1 s (one thread), subset 79.2 s and 735.4 s. Main is
faster on the subset because 13 more images solve, which on 0.2.0 were failed searches
(the 0.15°–0.3° band alone went from 13.6 s to 1.3 s, and from 138 s to 1.2 s on one
thread); the `neg_hint_*` controls cost the same in both.

### 9.6 Blind solving

`--blind-index <index>` gives arcsec `-i` and seiza a blind request against its own index. Neither gets
a position. arcsec still needs a hint on its command line, for the ordinary search it
falls back to when the index finds nothing, so the harness gives it the antipode of the
true centre with `-r 0`: only the index can find the field. Both keep the field size
(arcsec `--fov`, ±20 %; seiza a pixel-scale range of ±20 %, where its own default is
0.1–20″/px).

arcsec's blind numbers here use its own pattern index, **from PR #13 (`custom-index`,
4f0760a), not yet merged or released**; see [offline-index.md](offline-index.md). The
Astrometry.net 4100-series results measured with 0.2.0 are kept as a historical row.

| | arcsec, own index (PR #13) | arcsec, astrometry.net (0.2.0, historical) | seiza |
|---|---|---|---|
| index | built from the user's d80: 287 MB (fields 0.3°–30°) or 698 MB (0.15°–30°) | 4107–4119, 340 MB, Tycho-2, quads 22′ and up | `blind-gaia16.idx`, 1.63 GB (Gaia DR3 to G 16) |
| also needs | the d80/g05/w08 databases the hinted solver uses | the same | the 1.54 GB G ≤ 17 catalogue, for verification |
| acceptance | the hinted solver's star-level verification (≥ 30 matched stars, spread over the frame) | the same | ≥ 12 matched stars, RMS < 2 px |

**Expanded corpus (596 in tiers A/B/C/S), no position, scale given:**

| | arcsec, 698 MB index | arcsec, 287 MB index | seiza |
|---|---|---|---|
| correct | 473 | 399 | **494** |
| false positives | **1** | **1** | 58 |
| … of which a different part of the sky (11°–176° away) | 0 | 0 | 42 |
| inexact | 0 | 0 | 18 |
| no solve | 122 | 196 | 26 |
| tier-D controls answered | 4 (all `neg_hint_*`, correct) | 2 (the same kind) | 8 (6 `neg_hint_*` correct, `neg_orion_2massk` 143° away, `neg_shuffle_6`) |
| median centre / corner, the 438 both (698 MB) solve | 0.31″ / 0.66″ | | 0.43″ / 0.90″ |

| FOV (long side) | n | arcsec, 698 MB | arcsec, 287 MB | seiza |
|---|---|---|---|---|
| < 0.15° | 6 | 0 | 0 | 0 (5 FP) |
| 0.15–0.3° | 113 | 54 | 2 | 79 (10 FP) |
| 0.3–0.6° | 126 | 113 | 95 | 114 (11 FP) |
| 0.6–1.2° | 180 | 165 | 161 | 172 (4 FP) |
| 1.2–2.5° | 88 | 76 | 76 | 85 (1 FP, 2 inexact) |
| 2.5–6° | 43 | 24 | 24 | 15 (15 FP, 13 inexact) |
| 6–20° | 36 | 34 (1 FP) | 34 (1 FP) | 21 (12 FP, 3 inexact) |
| > 20° | 8 | 7 | 7 | 8 |

**v1, blind:** arcsec 76 of 98 with the 698 MB index and 64 with the 287 MB one, no false
positives; seiza 86 with 4 false positives (three of them other parts of the sky:
`stress_narrow`, `ls_north`, `ls_p14`); arcsec with the astrometry.net 4100 series
(0.2.0) 12.

* **The blind solves are mostly the hinted ones.** Of the images each solves hinted with
  the true centre (§9.3), arcsec's 698 MB index finds 473 of 558 blind and seiza's 488 of
  546. arcsec's misses are below 0.3° (the 287 MB index stops there, the 698 MB one at
  0.15°, and half the 0.15°–0.3° fields still fail) and in the SDSS, Pan-STARRS and
  Legacy Survey frames (8, 9 and 26 of 33, 20 and 40); seiza's are spread more evenly.
  Both solve 438; arcsec only 35 (distorted and coarse TESS, tier-S lenses, and ten
  fields that seiza places elsewhere on the sky or does not solve), seiza only 56 (SDSS
  and Pan-STARRS frames, SkyMapper, the crowded fields of §9.3).
* **The false positives differ in kind.** arcsec's one is the hinted solver's
  `wide_shassa_01` (1.2 px out). seiza returned 42 images (and one tier-D control) at a
  different part of the sky, between 11° and 176° from the truth: 0.1°–0.8° fields
  (random fields, Legacy Survey, WISE) and 2.9°–5.9° TESS FFI crops; its hinted solver, which only looks near the hint, made none of
  these. arcsec accepts a blind hypothesis only through the hinted solver's star-level
  verification; seiza's blind acceptance (≥ 12 matched stars, RMS < 2 px) is looser.
  The `neg_hint_*` controls are real fields with a hint 30°–60° away; a blind solve ignores
  the hint, so finding them is correct here.
* **The indexes are not the same size or source.** arcsec's are built locally from the
  d80 database the user already has (287 or 698 MB on top of it); seiza's are a 1.63 GB
  index plus a 1.54 GB catalogue from Gaia DR3, downloaded. Below 0.3° the 698 MB index is
  what arcsec needs.

**Timing**, 98-image subset, quiet machine (`--jobs 1`, default threads, outside load
0.2–0.35 cores, one round), median / mean / p90 / total (s):

| | arcsec, 698 MB | arcsec, 287 MB | seiza |
|---|---|---|---|
| correct (of 88) / FP | 66 / 0 | 58 / 0 | 64 / 16 |
| all 98 | 0.31 / 0.55 / 1.25 / 54.1 | 0.27 / 0.52 / 1.25 / 50.8 | 1.40 / 4.26 / 9.57 / 417 |
| answered | 0.37 / 0.67 / 1.38 / 45.9 (n = 68) | 0.54 / 0.71 / 1.46 / 41.4 (n = 58) | 1.38 / 1.65 / 1.60 / 142 (n = 86) |
| not answered | 0.21 / 0.27 / 0.41 / 8.2 (n = 30) | 0.18 / 0.23 / 0.38 / 9.4 (n = 40) | 30.7 / 23.0 / 46.4 / 276 (n = 12) |

arcsec's blind solve costs about what its hinted one does, and gives up in a fraction of
a second (it tries at most six hypotheses and then runs the ordinary search at `-r 0`).
seiza takes 1.4 s for an answer and 30–50 s to give up after its 400 hypotheses; under
load those searches ran past 300 s, and 26 of its corpus runs had to be repeated with a
longer timeout (none then solved; one, `sv2_07`, returned a different part of the sky).

**The historical row.** With 0.2.0 and the astrometry.net 4107–4119 files, arcsec solved
12 of v1 and 21 of the subset blind (seiza 86 and 64), nothing below 0.6° on v1, with no
false positives; its blind solves took a median 6 s on a quiet machine against seiza's
1.5 s, mostly loading and sorting the index files.

Commands: `scripts/benchmark.py --corpus --images <mirror> --auto-db --blind
<index.arcsecix> --arcsec <custom-index build> --seiza <seiza> --seiza-data
~/.local/share/seiza-data` (and without `--seiza` for the second index).

### 9.7 Caveats

* **One machine, one corpus.** A 24-core x86-64 Linux (WSL2) machine, warm page cache.
  The corpus is mostly reprojected survey cutouts and survey frames; tier C (42 LCO
  frames) is the only set of real observing frames, and none are DSLR or one-shot-colour
  camera frames apart from the simulated ones of tier S.
* **The catalogues differ.** arcsec and ASTAP use the same ASTAP `.1476/.290/.001` files;
  seiza uses its own Gaia tiles, to G 17. A field one solver misses for want of catalogue
  stars may be a catalogue result rather than a solver one; the G ≤ 15 run shows how much
  depth moves seiza.
* **Defaults, not tuning.** Every solver ran with its defaults apart from the hint. ASTAP
  in particular might solve more of the Legacy Surveys and SDSS frames with other
  `-z`/`-s` settings (§7.5).
* **The harness was written by arcsec's authors**, its thresholds (5″ or one pixel at a
  corner) were set for arcsec's own work, and the corpus was chosen to find arcsec's
  weaknesses — and arcsec's changes since 0.2.0 were measured on the same corpus, so its
  current numbers have had the benefit of being tuned against it, which seiza's and
  ASTAP's have not. seiza is held to the same corner threshold; many of its "false
  positives" are the right field with a plate one to three pixels out at the corners on
  12–29″/px images, which some users would accept.
* **seiza's own numbers** (its README) were measured differently: real camera frames,
  its own catalogue recommendations for each solver, Windows. They are not contradicted
  or reproduced here.

### 9.8 Faster failed searches (2026-10-02, after 0.3.0)

§9.5 found arcsec as fast as seiza and ASTAP on the images it solves, and far slower on
the ones it does not: a failed search visits every spiral position out to `-r`, and on
one thread the three `neg_hint_*` controls of the subset were 532 of its 555 s. The work
per position is now 4–8× smaller, with the same results (every `.wcs` of the corpus is
byte-identical, true-centre and offset hint, with and without an index installed, and at
`-r 10`); what changed and why is in
[plate-solving.md §10.3c](plate-solving.md#103c-what-a-failed-search-costs-2026-10-02).
With an index installed, a hint far outside `-r` now usually ends the search at once
([offline-index.md §2.7](offline-index.md#27-command-line)).

Same protocol as §9.5: `--jobs 1`, default threads and one thread (`--threads 1
--taskset 5`), two rounds, a wait for a load average under 4 before each run, and a
sampler of the CPU used outside the run; runs with a mean above 1 core or a peak above 4
were discarded and repeated (two of 0.3.0's were; the outside load came from other jobs
on the machine). The runs kept had 0.2–0.3 cores of outside load on average. arcsec
0.3.0 (`main`, 26f36e5) was re-timed in this session, alternating in order with an
intermediate build of this branch; the final build (532b78c) was timed straight after,
on the same protocol. Per image, median / mean / p90 / total (s), averaged
over the two rounds; totals agreed between rounds within 2 % (0.3.0) and 5 % (this
branch):

| Set, threads | | 0.3.0 | this branch |
|---|---|---|---|
| v1, default | all 103 | 0.16 / 0.32 / 0.67 / 33.2 | 0.08 / 0.24 / 0.58 / 24.9 |
| | solved (92) | 0.16 / 0.25 / 0.40 / 23.3 | 0.08 / 0.19 / 0.32 / 17.6 |
| | not solved (11) | 0.81 / 0.90 / 2.25 / 9.9 | 0.48 / 0.66 / 1.93 / 7.3 |
| v1, one thread | all 103 | 0.19 / 0.64 / 0.71 / 65.4 | 0.12 / 0.31 / 0.65 / 31.5 |
| | solved | 0.19 / 0.29 / 0.56 / 26.9 | 0.12 / 0.23 / 0.49 / 21.2 |
| | not solved | 3.67 / 3.50 / 8.20 / 38.5 | 0.81 / 0.93 / 2.24 / 10.3 |
| subset, default | all 98 | 0.15 / 0.68 / 0.40 / 66.3 | 0.09 / 0.21 / 0.34 / 20.2 |
| | solved (83) | 0.15 / 0.19 / 0.36 / 15.9 | 0.09 / 0.14 / 0.29 / 11.4 |
| | not solved (15) | 0.24 / 3.36 / 14.0 / 50.4 | 0.10 / 0.59 / 3.33 / 8.8 |
| subset, one thread | all 98 | 0.17 / 5.93 / 0.82 / 580.8 | 0.11 / 1.17 / 0.36 / 114.6 |
| | solved | 0.16 / 0.21 / 0.41 / 17.6 | 0.11 / 0.16 / 0.34 / 13.3 |
| | not solved | 1.25 / 37.5 / 173 / 563.2 | 0.22 / 6.75 / 40.7 / 101.3 |

The solved and unsolved sets are the same in both columns (0.3.0's `neg_hint_1` on one
thread is a 300 s timeout rather than "no solution"). Next to §9.5's ASTAP and seiza
(timed in an earlier session, so only roughly comparable): on the subset, one thread,
this branch's 115 s total is now level with ASTAP's 122 s and three times seiza's 38 s;
with all cores, 20 s against seiza's 39 s and ASTAP's 124 s. On v1, 25 s and 31 s against
seiza's 25 s and 27 s.

**The `neg_hint_*` controls**, real fields with the hint 30°–60° away and `-r 10`, where
the whole radius has to be searched (s; no index installed, as above):

| | 0.3.0, default | this branch, default | 0.3.0, one thread | this branch, one thread |
|---|---|---|---|---|
| `neg_hint_1` (`ls2_22`, 0.25°) | 27.9 | 3.33 | > 300 (timeout) | 45.3 |
| `neg_hint_4` (`sdss2_20`, 0.23°) | 14.0 | 3.56 | 173 | 40.7 |
| `neg_hint_7` (`sdss2_10`, 0.23°) | 6.7 | 1.00 | 81.2 | 13.5 |

Solved images are faster too (median 0.16 → 0.08 s on v1 at default threads, 0.19 →
0.12 s on one thread), since even position 0 matches its quads the new way.

**With the 698 MB index installed** (`d80_015.arcsecix` beside the database), all eight
`neg_hint_*` controls, one round:

| | 0.3.0, default | this branch, default | 0.3.0, one thread | this branch, one thread |
|---|---|---|---|---|
| `neg_hint_1` | 28.3 | 0.42 | 376 | 1.38 |
| `neg_hint_2` | 13.1 | 3.57 | 163 | 40.7 |
| `neg_hint_3` | 19.6 | 8.04 | 211 | 74.9 |
| `neg_hint_4` | 14.1 | 3.56 | 177 | 42.4 |
| `neg_hint_5` | 7.4 | 0.57 | 82.4 | 1.69 |
| `neg_hint_6` | 6.7 | 1.42 | 78.7 | 18.9 |
| `neg_hint_7` | 6.8 | 0.24 | 81.9 | 0.60 |
| `neg_hint_8` | 7.2 | 0.51 | 86.8 | 1.53 |
| total | 103 | 18.3 | 1256 | 182 |

All sixteen answers are "no solution", as they must be. Four (`_1`, `_5`, `_7`, `_8`) are
fields the index verifies elsewhere on the sky, so the search stops at once; the other
four are fields of 0.16°–0.23° that the index does not find (§9.6: half of the
0.15°–0.3° band fails blind), and they still pay for the whole spiral. (0.3.0's one-thread
run had 0.3 cores of outside load on average, with a peak of 4.25 from a compile; its
two repeats were busier and are not used.)

**N.I.N.A.-style, `-r 180`, hint at the antipode** (the eight images of
[offline-index.md §7.3](offline-index.md#73-ninastyle--r-180-no-useful-hint), default
threads, one at a time):

No index (the spiral over the whole sky): wall time and total CPU (user + system) of
each run, the two builds back to back on each image. These ran with 24 threads each on
a machine another job was also using (load average 9–75 at the starts), so the CPU
column is the fairer comparison:

| image (field) | answer, both builds | 0.3.0 wall / CPU | this branch wall / CPU |
|---|---|---|---|
| `lco_12` (0.48°) | correct | > 900 s (timeout) / > 13 100 s | 146 s / 2 630 s |
| `ztf_06` (0.56°) | correct | 793 s / 12 460 s | 86 s / 2 000 s |
| `type_m101` (1.0°) | correct | 436 s / 4 220 s | 29 s / 673 s |
| `decp20` (1.0°) | correct | 251 s / 4 075 s | 28 s / 660 s |
| `dens_blank` (1.5°) | correct | 113 s / 1 853 s | 13 s / 304 s |
| `fov_3p00` (3.0°) | correct | 32 s / 449 s | 3.7 s / 79 s |
| `rnd_022` (0.35°) | **wrong**, the same position 149° off in both (as in offline-index §7.3) | 129 s / 2 182 s | 15 s / 348 s |
| `wide_dss_01` (4.0°) | no solution | 20 s / 255 s | 2.0 s / 44 s |

CPU per image is 5–7× less. With the 698 MB index installed, the automatic path finds
the field before the spiral, as in offline-index §7.3; quiet machine, default threads,
one image at a time, wall time:

| | 0.3.0 | this branch |
|---|---|---|
| the seven that solve (all correct, `rnd_022` included) | 0.94–3.36 s | 0.36–2.70 s |
| `wide_dss_01` (no solution) | 22.2 s | 3.9 s |
| total, eight images | 33.3 s | 10.5 s |
