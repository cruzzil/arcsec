# Benchmark Image Corpus

A standing set of FITS images with trustworthy ground truth, used to measure arcsec's
solve rate, positional accuracy and false-positive rate across a realistic range of field
sizes, star densities and image quality.

Companion to [plate-solving.md](plate-solving.md). Fetch everything with:

```bash
scripts/fetch-test-images.sh                  # all tiers → resources/testset/
scripts/fetch-test-images.sh --tier A         # just the synthetic-truth tier
scripts/fetch-test-images.sh --list           # show the manifest without downloading
```

The manifest lives at [`scripts/test-images.tsv`](../scripts/test-images.tsv). Images are
**not** committed — `resources/` is gitignored — the manifest plus the fetch script are the
reproducible artefact.

---

## 1. What "ground truth" means here

Not all truth is equal, and the metric you can trust depends on where the image came from.
The corpus is stratified into four tiers.

| Tier | Ground truth | Truth accuracy | What it exercises |
|---|---|---|---|
| **A — Synthetic cutouts** | We *request* the centre, FOV and projection, so the true WCS is known exactly | exact (< 1 mas) | Correctness of the core solver across FOV and star density |
| **B — Survey cutouts** | The archive's own pipeline WCS in the header | 0.02–0.3″ | Real PSFs, real noise, real artefacts; sub-arcsecond accuracy |
| **C — Real observing frames** | A reference solution from ASTAP and/or astrometry.net, cross-checked | ~0.5–1″ | Gradients, trailing, Bayer mosaics, hot pixels, clouds |
| **D — Negative / stress** | Known to be unsolvable, or known-hard | n/a | **False-positive rate** and graceful failure |

Tier D is the one most benchmark suites omit and the one that matters most for us: a solver
that always returns *something* is worse than one that says "no solution". See
[plate-solving.md §11.1](plate-solving.md#111-the-catalogue-path-performs-no-verification).

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
field is only ~3.4′ — below what the ASTAP `d80` database supports (documented for FOV ≥
0.6°). Include at most one, and expect it to fail with our current catalogue; it is here to
document the limit, not to pass. A denser database (`d50`/`v50`) or a purpose-built index
would be needed.

### 2.7 Real observing frames — tier C

* **Our own NGC 3372 set** — 240 frames, ASI533MC Pro + RedCat 51, ~3.13″/px, already in
  `resources/` and already used by `scripts/bench_all.sh`. This is the most valuable tier-C
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
* **Practical note**: for any tier-C image, generate the reference with *two* independent
  solvers (`astap-cli` and `solve-field`) and only admit the image to the corpus if they
  agree to better than 2″. Disagreement means the truth is not trustworthy, not that one
  solver is wrong.

### 2.8 Rejected sources, and why

* **ESO Science Archive / Phase 3** — excellent data, but products are large multi-extension
  mosaics and the query interface is not a simple parameterised GET. Not worth the
  complexity for our purposes.
* **astropy / photutils example data** — small, convenient, but only a handful of images
  and several lack a usable WCS or enough stars.
* **AAVSO VPhot, Astrometrica samples** — access requires accounts.

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
                     └─ the regime our NGC 3372 data sits in (4.6°)
```

At the wide end, ASTAP-family databases thin out and distortion dominates; at the narrow
end, `d80` runs out of stars (documented floor: 0.6°). Both ends are where we expect to
fail, and knowing *where* the cliff is, is the point.

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
handling (`δ > π/2 → π − δ`, flip RA by π) and an RA-offset guard, and none of it is
covered by a test today.

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

* An HST drizzled product (SIP in the header) — expected to fail on FOV grounds today.
* A wide-field (≥ 8°) HiPS cutout in TAN — at that width the tangent-plane approximation
  itself introduces several arcseconds of error at the corners, which is exactly the
  regime where a linear-only fit shows its limits.
* If available, a real camera-lens frame from tier C with visible barrel distortion.

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
rotation and distortion error from pointing error, and is the number that will move when
[§12.4](plate-solving.md#124-fit-sip-distortion-and-honour---sip) (SIP fitting) lands.

### 4.1 Reference-frame caveats

* SkyView returns `RADESYS = FK5`; Legacy Survey, PS1 and Gaia are ICRS. The difference is
  < 0.1″ — below our targets, but it should be *stated* rather than silently absorbed, and
  if we ever tighten the centre-error target below 0.5″ it must be corrected for.
* Survey WCS headers are themselves solutions, with their own errors (0.02–0.3″). Tier B
  therefore cannot validate us below ~0.3″; only tier A can.
* Comparing `CRVAL` directly is wrong when the two solutions use different reference
  pixels. Always compare **sky positions of the same pixel coordinates**, as above.

---

## 5. Running the benchmark

The existing scripts already cover part of this and should be extended rather than
replaced:

| Script | Covers |
|---|---|
| `scripts/hips_solve_test.sh` | tier A, fixed 1.5° FOV, ~12 fields, solve rate + accuracy |
| `scripts/hips_extended_test.sh` | tier A, 20 further fields including expected-failure types |
| `scripts/hips_fov_sweep.sh` | tier A, the FOV ladder of §3.1 |
| `scripts/bench_all.sh` | tier C, our NGC 3372 set, arcsec vs `astap-cli` |
| `scripts/compare-solvers.sh` | single-image arcsec vs ASTAP diff |
| **`scripts/fetch-test-images.sh`** | **new** — materialises tiers A/B from the manifest |

What is missing and should be added:

1. A **corner-error** computation (all current scripts compare centres only).
2. A **tier-D harness** that asserts failure, and fails the run if any negative control
   solves.
3. A single `scripts/benchmark.sh` that runs every tier and emits one CSV plus a summary,
   so a regression is one command and one diff.

### 5.1 Manifest format

`scripts/test-images.tsv`, tab-separated, `#` comments:

```
id      tier  ra        dec      fov_deg  width  height  source        extra
m42     A     83.82     -5.39    1.5      1000   1000    hips2fits     CDS/P/DSS2/blue
cygx    A     305.55    40.73    1.0      1000   1000    hips2fits     CDS/P/DSS2/red
ngc188  B     11.80     85.24    0.5      1024   1024    skyview       DSS2R
ls_sparse B   180.00    20.00    0.28     1024   1024    legacysurvey  r
ps1_m67 B     132.83    11.81    0.07     1024   1024    panstarrs     r
sdss_1  B     276.95    -0.16    0.22     -      -       sdss          301/2505/3/r/38
```

`ra`/`dec`/`fov_deg` are the ground truth for tier A; for tier B they are the *request*
and the truth is read back from the delivered header.

---

## 6. Results — 103 images, arcsec vs ASTAP

Corpus of 103 images (63 tier A, 34 tier B, 6 tier D), d80 database, both solvers given the
same true field centre and FOV, `-r 3`. ASTAP is `astap_cli` CLI-2026.07.30. A reported
solve counts as **correct** only if its worst corner error is under 5″; anything else is a
false positive.

### 6.1 Where it started and where it ended (2026-09-02)

| | correct | false positives | tier A | tier B | tier D solved |
|---|---|---|---|---|---|
| Baseline (2026-09-01) | 48 | 4 | 26/62 | 22/34 | 0/7 |
| + 9-NN quad redundancy | 64 | 4 | 38/62 | 26/34 | 0/7 |
| + star-level verification | 86 | 0 | 54/63 | 32/34 | 0/6 |
| + spread check, no star trim | 89 | 0 | 55/63 | 34/34 | 0/6 |
| + `.290`/`.001` catalogues, `--auto-db` | **90** | **0** | **56/64** | **34/34** | 0/5 |
| **ASTAP CLI-2026.07.30** | 47 | 0 | 39/63 | 8/34 | 0/6 |

Accuracy of arcsec's correct solves:

| Tier | centre (median / max) | corner (median / max) | scale (median) | time (median) |
|---|---|---|---|---|
| A | 0.733″ / 1.691″ | 0.958″ / 4.137″ | 0.0053% | 0.33 s |
| B | **0.161″** / 1.123″ | **0.272″** / 1.719″ | 0.0059% | 0.15 s |

On the 45 images both solvers get right, arcsec's corner accuracy is now slightly ahead
(0.897″ vs 0.938″) at essentially equal centre accuracy (0.691″ vs 0.660″). ASTAP solves
only two images arcsec does not (`type_m31`, `type_m44`); arcsec solves 44 that ASTAP does
not.

The cost is speed: **0.60× ASTAP** (median 0.33 s vs 0.18 s on tier A). 126 quads per star
is not free. See §6.5.

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

### 6.4 What still fails (8 of 63 tier A, 0 of 34 tier B)

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

**Negative and neutral results, so they are not re-tried:**### 6.6 A truth bug found by the comparison

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

### 6.7 Reproducing

```bash
scripts/fetch-test-images.sh                       # 103 images, ~2.2 GB
scripts/benchmark.py --astap ~/astap_cli --jobs 10 --radius 3 --csv results.csv
```

## 7. Licensing and attribution

All tier A/B sources are public-domain or freely redistributable *data*, but the corpus is
deliberately not committed: `resources/` is gitignored and the manifest is the artefact.
If any of this is ever published, credit as the providers ask:

| Source | Attribution |
|---|---|
| HiPS2FITS | CDS, Strasbourg; plus the underlying survey (e.g. DSS2 — STScI/AURA) |
| SkyView | "The SkyView virtual observatory", NASA/GSFC HEASARC |
| Legacy Surveys | DESI Legacy Imaging Surveys — NOIRLab/DOE/NSF, DR10 |
| Pan-STARRS1 | PS1 Surveys — University of Hawaii IfA / STScI |
| SDSS | SDSS-IV/DR17 acknowledgement text |
| HST / MAST | STScI; observation programme ID |
| Tier C amateur data | as specified by the individual publisher |
