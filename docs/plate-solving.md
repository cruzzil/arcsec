# Plate Solving: Methods, Mathematics, and Where Arcsec Stands

A survey of astrometric plate-solving algorithms, the mathematics behind them, a precise
description of what **arcsec** currently implements, and an assessment of our gaps.

Written 2026-08-31; §10–§13 revised 2026-09-25 against the 0.1.0 code. Sources are
linked inline and collected in [§14 References](#14-references).

**Contents**

1. [The problem](#1-the-problem)
2. [The canonical pipeline](#2-the-canonical-pipeline)
3. [Stage 1 — Source detection](#3-stage-1--source-detection)
4. [Stage 2 — Asterism construction (the feature families)](#4-stage-2--asterism-construction-the-feature-families)
5. [Stage 3 — Matching and search structures](#5-stage-3--matching-and-search-structures)
6. [Stage 4 — Consensus and outlier rejection](#6-stage-4--consensus-and-outlier-rejection)
7. [Stage 5 — Fitting the WCS](#7-stage-5--fitting-the-wcs)
8. [Stage 6 — Verification](#8-stage-6--verification)
9. [Reference catalogues and indexes](#9-reference-catalogues-and-indexes)
10. [What arcsec does today](#10-what-arcsec-does-today)
11. [Shortcomings](#11-shortcomings)
12. [Improvement roadmap](#12-improvement-roadmap)
13. [Benchmarking](#13-benchmarking)
14. [References](#14-references)

---

## 1. The problem

Given an image of the night sky and a catalogue of stars with known celestial
coordinates, find the mapping

```
    (x, y)  pixel coordinates   ⟶   (α, δ)  right ascension / declination
```

that best explains the image. The output is a **World Coordinate System** (WCS) — in
practice a FITS header carrying `CTYPE1/2`, `CRPIX1/2`, `CRVAL1/2` and either a `CD`
matrix or `CDELT` + `CROTA2`, optionally extended with a distortion polynomial.

Two regimes, and almost every design decision follows from which one you are in:

| Regime | Known a priori | Search space | Typical use |
|---|---|---|---|
| **Hinted** ("nearby") | approximate α, δ **and** pixel scale | a few square degrees | mount-attached solving, autoguiding, plate-solve-and-sync |
| **Blind** ("lost in space") | nothing, or only the pixel scale | the whole sky × all scales | arbitrary uploaded images, star trackers, archival plates |

The blind problem is a generalisation of the star-tracker *lost-in-space* problem, in which
"nothing — not even the image scale — is known" ([Lang et al. 2010][lang2010]).

The essential difficulty is that neither list of stars is a subset of the other: the image
contains hot pixels, cosmic rays, satellites, galaxies and asteroids that are not in the
catalogue; the catalogue contains stars too faint, too saturated or too close to a
neighbour to be detected. Any usable algorithm must therefore match **partial, noisy,
unordered point sets under an unknown similarity transform**, and must do so without ever
being able to enumerate correspondences directly (there are `n!/(n-k)!` of them).

The universal trick is **geometric invariance**: build small local features out of `k`
stars whose descriptor is invariant under translation, rotation, scaling and (usually)
reflection. Matching descriptors reduces an exponential correspondence problem to a
nearest-neighbour lookup.

---

## 2. The canonical pipeline

Essentially every solver — `astrometry.net`, ASTAP, Watney, tetra3, PinPoint,
PlateSolve, Siril, StellarSolver, arcsec — is a specialisation of this pipeline:

```
   ┌────────────────────────────────────────────────────────────────────────┐
   │  FITS / TIFF / raw image                                               │
   └───────────────┬────────────────────────────────────────────────────────┘
                   │
      ┌────────────▼─────────────┐
      │ 1. SOURCE DETECTION      │   background model, threshold, connected
      │    → star list (x,y,flux)│   components, sub-pixel centroid, HFD/FWHM
      └────────────┬─────────────┘   filter: too-small, too-large, saturated
                   │
      ┌────────────▼─────────────┐
      │ 2. ASTERISM BUILDING     │   triangles / quads / n-star patterns
      │    → invariant descriptor│   descriptor ⟂ {translate, rotate, scale, flip}
      └────────────┬─────────────┘
                   │
      ┌────────────▼─────────────┐   ┌──────────────────────────────────────┐
      │ 3. MATCHING              │◀──┤ REFERENCE SIDE                       │
      │    descriptor lookup     │   │  (a) catalogue read + same asterism  │
      │                          │   │      construction, per trial position│
      └────────────┬─────────────┘   │  (b) pre-built offline index of all  │
                   │                 │      catalogue descriptors           │
      ┌────────────▼─────────────┐   └──────────────────────────────────────┘
      │ 4. CONSENSUS / OUTLIERS  │   vote accumulator, median-scale filter,
      │    → inlier match set    │   RANSAC, bijective filter, sigma clipping
      └────────────┬─────────────┘
                   │
      ┌────────────▼─────────────┐
      │ 5. WCS FIT               │   least squares for the 6 plate constants,
      │    → CRVAL/CRPIX/CD      │   then optionally SIP/TPV distortion terms
      └────────────┬─────────────┘
                   │
      ┌────────────▼─────────────┐
      │ 6. VERIFICATION          │   project catalogue → image, count agreeing
      │    accept / reject       │   stars; Bayesian odds or a match threshold
      └────────────┬─────────────┘
                   │ accept                     │ reject
                   ▼                            ▼
             write WCS                    next hypothesis / next sky position
```

The families of solvers differ mainly in **stage 2** (what the descriptor is) and
**stage 3** (how the reference side is organised). Stages 1, 4, 5 are broadly shared.

---

## 3. Stage 1 — Source detection

Garbage in, garbage out: a solver's robustness ceiling is set here. Three lineages matter.

### 3.1 Background estimation

Star detection needs a *local* background and a noise estimate, because real frames have
gradients (light pollution, vignetting, amp glow, moonlight).

**SExtractor / SEP** ([Bertin & Arnouts 1996][bertin1996]) is the reference implementation:
tile the frame into meshes (typically 32×32 or 64×64 px), compute a σ-clipped estimate of
the mode in each mesh, median-filter the mesh grid to suppress bright-object contamination,
then bicubic-spline interpolate back to full resolution. The mode is estimated from the
clipped mean and median as

```
    mode ≈ 2.5 × median − 1.5 × mean
```

**ASTAP / arcsec** use a cheaper global scheme: build a 16-bit histogram of the pixel
values and take its peak (the mode) as the background, falling back to the histogram mean
when that is more than 1.5× the mode. The noise is the RMS deviation from that background
over a sparse grid of sample pixels, iteratively clipped at 3σ until it changes by less
than 5% (at most 7 passes). Two star levels are then read off the histogram's upper tail,
and the detection threshold is

```
    star_level = background + k · noise         (ASTAP uses a multi-pass k)
```

Only the last-resort detection pass estimates a *local* background: it splits the frame
into a ~12×12 grid and runs a histogram σ-clipped mean (upper clip 2σ) in each cell.

The trade-off is explicit: no spatial background model means strong gradients bias the
threshold, but the estimate costs one pass over the pixels instead of a mesh + spline.

### 3.2 Detection and centroiding

* **Thresholding + connected components** (SExtractor, ASTAP): pixels above
  `background + k·σ` are grouped into connected blobs. SExtractor additionally
  **deblends** by re-detecting each object at 30 exponentially-spaced thresholds between
  the peak and the detection threshold, splitting whenever a branch carries ≥ 0.1% (a
  tunable `DEBLEND_MINCONT`) of the total flux.
* **Matched filtering** (DAOFIND, Stetson 1987): convolve with a Gaussian kernel of the
  expected PSF width and find local maxima of the convolved image. Optimal for isolated
  point sources at low S/N.
* **Sub-pixel centroid**: first-moment (centre of gravity) over the background-subtracted
  pixels, or a Gaussian fit. `astrometry.net` fits a Gaussian to a 3×3 grid around each
  peak pixel ([Lang et al. 2010][lang2010]).

Centroid precision drives the final astrometric residual. For a Gaussian PSF of width σ
and total S/N, the centroid error is approximately

```
    σ_centroid ≈ σ_PSF / (S/N)
```

so a 2 px PSF at S/N = 50 gives ~0.04 px — well below the ~0.2–1 px systematics from
distortion and catalogue errors in a typical amateur setup.

### 3.3 HFD vs FWHM

ASTAP (and therefore arcsec) measures **half-flux diameter** rather than FWHM. HFD is the
diameter of the circle centred on the source that contains half the total flux. For a
perfect Gaussian ([Wikipedia: HFD][hfd]):

```
    HFD  = 2.50663 · σ        (= √(2π)·σ)
    FWHM = 2.35482 · σ
    ⟹ HFD ≈ 1.064 · FWHM       for an in-focus Gaussian
```

HFD is preferred because it is an *integral* rather than a *point* measurement: it degrades
gracefully for defocused, trailed, or badly-seeing-blurred stars, where FWHM becomes
double-peaked and meaningless. Its use as a quality filter (`-m/--hfd-min`) is how ASTAP
and arcsec reject hot pixels (HFD too small) and galaxies/nebulae knots (HFD too large).

---

## 4. Stage 2 — Asterism construction (the feature families)

### 4.1 Triangles — 2 invariants

The oldest family ([Groth 1986][groth1986], [Valdes et al. 1995][valdes1995]). Take three
stars, sort the side lengths `s₁ ≥ s₂ ≥ s₃`, and use

```
    descriptor = ( s₂/s₁ , s₃/s₁ )              (arcsec's tetra.rs)
```

or, equivalently, astroalign's ([Beroiz et al. 2020][astroalign]) sorted-ascending form
with `L₂ ≥ L₁ ≥ L₀`:

```
    descriptor = ( L₂/L₁ , L₁/L₀ )
```

Both are invariant to translation, rotation, uniform scale and reflection. The admissible
region is bounded by the triangle inequality (`x ≤ 1 + 1/y` in astroalign's parametrisation).

```
    C                 sides:  a = |AB|, b = |BC|, c = |CA|
    /\                sort:   s₁ ≥ s₂ ≥ s₃
   /  \               code:   (s₂/s₁, s₃/s₁)  ∈ (0,1]²
  /    \
 A──────B             3 stars ⟹ C(n,3) = n(n−1)(n−2)/6 triangles
```

**Strength**: only 3 stars needed, so it works in sparse fields and survives a high
missing-star rate — with `n` catalogue stars and detection completeness `p`, the chance
that a given feature survives is `p³` for triangles vs `p⁴` for quads.

**Weakness**: only **2** invariants. The false-match rate is catastrophic compared with
quads (see §4.4). Practical triangle matchers therefore *always* pair the descriptor match
with a strong consensus stage (RANSAC in astroalign, a bijectivity constraint plus σ-clipping
in arcsec's `tetra.rs`, a voting histogram in FOCAS).

### 4.2 Six-distance quads, 5 ratios — the ASTAP family

Used by ASTAP, [Watney][watney] and arcsec. Take four stars; the six pairwise distances
form an irregular tetrahedron when drawn as a graph. Sort them descending,
`d₁ ≥ d₂ ≥ … ≥ d₆`, and normalise by the longest:

```
    descriptor = ( d₂/d₁ , d₃/d₁ , d₄/d₁ , d₅/d₁ , d₆/d₁ )   ∈ (0,1]⁵
    scale      = d₁                                          (kept separately)
```

```
        A ──────────── B        six distances:
        │╲            ╱│          AB, AC, AD, BC, BD, CD
        │ ╲          ╱ │        sort descending → d₁…d₆
        │  ╲        ╱  │        code = d₂/d₁ … d₆/d₁   (5 numbers)
        │   ╲      ╱   │        d₁ = the "scale" of the quad
        │    ╲    ╱    │
        │     ╲  ╱     │        ⟂ translation, rotation, scale, reflection
        D ──────────── C        (sorting kills the labelling ambiguity)
```

"The five ratios are sufficient for the search. Rotation, scaling, and flipping of the
image have no influence on the five ratios" ([ASTAP algorithm page][astap-alg]).

**Strength**: 5 invariants from a plain sort — no coordinate frame, no symmetry-breaking
case analysis, trivially reflection-invariant. Very low false-match rate.

**Weakness**: the descriptor discards *all* orientation information, so a matched pair of
quads tells you the two quads are congruent but not which star corresponds to which. ASTAP
therefore fits the plate using **quad centroids** rather than star correspondences (see
§7.3). arcsec does the same for its first fit, then recovers individual star pairs by
projecting the catalogue through that fit and re-fitting (§11.4). The descriptor is also
degenerate for symmetric configurations (e.g. a square, where several ratios coincide).

### 4.3 Code-space quads — the astrometry.net family

[Lang et al. 2010][lang2010] use a different, orientation-*preserving* descriptor. Of the
four stars, pick the most widely separated pair `A`, `B`. Define a local frame in which
`A → (0,0)` and `B → (1,1)`; then express `C` and `D` in that frame:

```
    code = (x_C, y_C, x_D, y_D)  ∈ ℝ⁴
```

```
      (0,1)          (1,1)=B      • A and B define the frame
        ┌───────────────┐         • C, D must lie inside the circle
        │      ∘ D      │           with AB as diameter
        │   ∘ C         │         • code = the 4 coordinates of C and D
        │               │
        └───────────────┘         symmetry breaking (Lang et al. §2.2):
      A=(0,0)        (1,0)           x_C ≤ x_D   and   x_C + x_D ≤ 1
```

Concretely, with `ab = B − A`, `s = |ab|²`, and the 45°-rotated basis actually used in the
index files (this is the form arcsec's `blind.rs` implements):

```
    cosθ = (ab_y + ab_x) / s
    sinθ = (ab_y − ab_x) / s

    code_x(P) = −(P−A)_x · sinθ + (P−A)_y · cosθ
    code_y(P) =  (P−A)_x · cosθ + (P−A)_y · sinθ
```

The two symmetry-breaking conditions (`A↔B` swap inverts all codes; `C↔D` swap exchanges
the pairs) reduce the 4! labellings to one canonical representative.

**Strength**: the code is a point in a bounded 4-D space with a *metric* — "nearby code"
is meaningful — so the reference side can be a **kd-tree** and matching is a ball query
rather than a scan. Crucially, a code match also yields the star correspondence `A↔A'`,
`B↔B'`, `C↔C'`, `D↔D'` directly, which is what makes a precise WCS available from a
*single* matched quad.

**Weakness**: not reflection-invariant. Both image parities must be tried (arcsec's
`blind.rs` builds every image quad twice, once per parity). The construction is also
fiddly — the canonical-form rules must match the index builder exactly or nothing matches.

A `DIMQUADS=3` variant exists (just `C` in the `AB` frame, a 2-D code) and is recommended
for wide-angle images; `DIMQUADS=5` exists but is "probably not useful"
([build-index docs][build-index]).

### 4.4 How many invariants do you need? — the false-positive arithmetic

This is the single most useful piece of arithmetic for reasoning about a solver's design.

Let the descriptor have `N` dimensions, each matched to a tolerance `±t` on a quantity
roughly uniform on a unit interval. The probability that a *random* reference feature
matches a given image feature is approximately

```
    P_random ≈ (2t)^N
```

and with `M` image features and `R` reference features in play, the expected number of
random matches is `M · R · (2t)^N`. Working numbers, using arcsec's defaults
(`t = 0.007`, ~500 image quads, ~5000 catalogue quads per trial position, `M·R = 2.5×10⁶`):

| Family | N | t | `(2t)^N` | Expected false matches |
|---|---|---|---|---|
| Triangle (`s₂/s₁, s₃/s₁`) | 2 | 0.007 | 2.0×10⁻⁴ | **~490** |
| Triangle, tightened | 2 | 0.0021 | 1.8×10⁻⁵ | ~44 |
| astrometry.net code | 4 | 0.007 | 3.8×10⁻⁸ | ~0.1 |
| ASTAP 5-ratio quad | 5 | 0.007 | 5.4×10⁻¹¹ | **~1.3×10⁻⁴** |

This is exactly why arcsec's `tetra.rs` carries the constant

```rust
pub const TETRA_TOL_FACTOR: f64 = 0.3;   // 0.007^0.4 ≈ 0.3
```

— an attempt to equalise the false-positive rate by shrinking the triangle tolerance, at
the cost of rejecting true matches whose ratios are perturbed by centroid noise. The
fundamental asymmetry does not go away: a 2-D descriptor cannot be made as selective as a
5-D one without also becoming intolerant of noise. **Triangles must be backed by consensus;
quads can nearly stand alone.**

(The estimate is an order-of-magnitude one: the ratios are not uniformly distributed —
`d₂/d₁` clusters near 1 — so the true rate is somewhat higher, but the ratio *between*
families is right.)

The working numbers above assume one quad per star, as ASTAP builds. arcsec now builds all
C(9,4) = 126 quads from each star's neighbourhood (§10.3), so `M` and `R` are each up to
two orders of magnitude larger and the expected number of random 5-ratio matches per
position rises towards ~0.1–1 (overlapping quads are not independent, so this overstates
it). That is still far below the several agreeing matches a solve needs, and the vote
filter and star-level verification (§6.2, §11.1) are what keep it from mattering.

### 4.5 Star-tracker patterns — Pyramid, TETRA, tetra3

The spacecraft attitude-determination literature solves the same lost-in-space problem
under much harder constraints (milliwatts, milliseconds, no filesystem).

* **Pyramid** ([Mortari et al. 2004][mortari2004]) matches inter-star *angles* using the
  **k-vector** — a search-free range query that returns, in O(1), every catalogue pair whose
  angular separation lies within the measurement uncertainty of an observed pair. It then
  builds up a four-star "pyramid" whose mutual consistency makes a false identification
  essentially impossible, and is explicitly designed to survive many spurious detections.
* **TETRA** ([Brown, Stubis & Cahoy 2017][tetra2017]) replaces search entirely with a
  **directly-addressed hash table**: a 4-star pattern's five normalised distances are
  quantised into bins, the bin tuple is hashed, and the correct catalogue pattern is
  retrieved in **a single database access**. This is the same 5-ratio descriptor as ASTAP's,
  used with a completely different lookup strategy.
* **tetra3** (ESA / Gustav Pettersson, and the `cedar-solve` fork) is the modern
  open-source implementation. Its database is generated for a declared FOV range
  (`min_fov`/`max_fov`) with `pattern_stars_per_fov` (default 10) and
  `verification_stars_per_fov` (default 30); matching uses `pattern_max_error` (default
  0.005) and accepts a solution when the probability of the observed number of verification
  matches arising by chance falls below `match_threshold` (default 1e-3).

The measured payoff is dramatic: on a Raspberry Pi 5, [AstroKeith's comparison][astrokeith]
reports ~12 ms for cedar-solve and ~200 ms for tetra3 versus <1 s for astrometry.net —
but tetra3 "required meticulous database customisation" and cedar-detect is "not as robust"
at star detection, notably failing on defocused stars. That is the trade the whole field
makes: **directly-addressed hashing is unbeatably fast and unforgivingly narrow.**

### 4.6 Family comparison

| | Triangles | 5-ratio quads | Code-space quads | tetra3 hash |
|---|---|---|---|---|
| Stars per feature | 3 | 4 | 4 | 4 |
| Invariant dims | 2 | 5 | 4 | 5 (binned) |
| Reflection-invariant | yes | yes | **no** (2 parities) | yes |
| Gives star correspondence | yes (with ordering) | **no** (centroids only) | **yes** | yes |
| Features from `n` stars | C(n,3) | C(n,4), or n (3-NN, ASTAP), or ≤126·n (9-NN, arcsec) | C(n,4) | limited set |
| Reference structure | kd-tree in 2-D | sorted array / hash bins | kd-tree in 4-D | direct hash |
| Robust to missing stars | best (`p³`) | `p⁴` | `p⁴` | `p⁴`, narrow FOV band |
| Used by | FOCAS, astroalign, arcsec `--method tetra` | ASTAP, Watney, arcsec (default) | astrometry.net, arcsec blind mode | tetra3, cedar-solve |

---

## 5. Stage 3 — Matching and search structures

Two orthogonal decisions: **how the reference descriptors are organised**, and **how the
sky is searched**.

### 5.1 Reference-side data structures

```
 (a) BRUTE FORCE                    for each image feature:
     O(M · R)                         for each reference feature:
                                        compare all N ratios
     ASTAP ≤ 120 quads; arcsec's find_matches()

 (b) SORTED + BINARY SEARCH        sort reference by one ratio;
     O(M · (log R + hits))           binary-search the ±t window on that ratio;
                                     check the remaining N−1 ratios on the survivors
     arcsec's find_matches_sorted(), keyed on ratio[4] = d₆/d₁
     (INDEX_RATIO): the most widely spread ratio gives the
     narrowest window

 (c) HASH BINS                     bin each ratio to width ~2t; hash the bin tuple;
     O(M) expected                   probe the 2^N neighbouring bins for edge cases
     ASTAP ≥ 120 quads: hash_bins = round(1/tolerance) + 2

 (d) kd-TREE                       ball query in N-D code space, radius = code tolerance
     O(M · log R)
     astrometry.net (codes are metric); astroalign (2-D invariant space)

 (e) DIRECTLY-ADDRESSED HASH       quantise → hash → one lookup, no probing
     O(M)                            requires a narrow, pre-declared FOV band
     TETRA / tetra3
```

Note that (c), (d) and (e) all require the reference descriptors to *exist* before the
query. That is only possible if the reference features are built **offline**, which forces
the second decision.

### 5.2 Sky-search strategy: spiral vs pre-indexed

```
 SPIRAL SEARCH (ASTAP, Watney "nearby", arcsec default)
 ─────────────────────────────────────────────────────
     hint (α₀,δ₀) ──► trial position 0 ──► read catalogue stars in a FOV-sized box
                          │                 build reference quads NOW (online)
                          │                 match; if enough matches → fit → done
                          ▼ else
                      trial position 1 (one FOV step away)  ... spiral outward

     +  index = the raw star catalogue (~50 MB), no build step
     +  1–2 positions when the hint is good: very fast
     −  cost grows as O((radius/FOV)²): a 30° radius at 1° FOV ≈ 2800 positions
     −  reference quads are rebuilt from scratch at every position
     −  needs a pixel-scale estimate to choose the step size

 PRE-INDEXED (astrometry.net, tetra3, Watney blind)
 ──────────────────────────────────────────────────
     OFFLINE:  for every catalogue region, at every scale band,
               enumerate quads → compute codes → store in kd-tree / hash

     ONLINE:   for each image quad: one lookup → candidate (code, sky position)
               → build a WCS hypothesis → verify → accept or next quad

     +  cost is independent of search radius — blind solving is *cheap*
     +  no pixel-scale hint needed if you sweep the scale bands
     −  index build takes hours and 10–40 GB
     −  a hinted solve is no faster than a blind one
```

`astrometry.net` implements the pre-indexed side thoroughly. The sky is tiled with
**HEALPix** ([Górski et al. 2005][healpix]), whose defining property is that all
`N_pix = 12·N_side²` pixels have exactly equal area `A = π/(3 N_side²)` — critical because
quad density must be uniform for the code statistics to hold. Grid cells are chosen "about
a third of the size of the query images", 10 stars are selected per cell, quads are built
with diameters in a band of `1` to `√2` cell side-lengths, each star may appear in up to 8
quads, and roughly 16 passes are made over the grid to guarantee redundancy
([Lang et al. 2010][lang2010]).

The result is the published index series: one file per scale band, geometrically spaced
by √2. The 4200-series (built from 2MASS) covers skymark diameters from 2.0′ to 2000′
across 20 scale steps ([astrometry.net README][anet-readme]):

| Index | Skymark diameter | | Index | Skymark diameter |
|---|---|---|---|---|
| 4200 | 2.0′–2.8′ | | 4210 | 60′–85′ |
| 4201 | 2.8′–4.0′ | | 4211 | 85′–120′ |
| 4202 | 4.0′–5.6′ | | 4212 | 120′–170′ |
| 4203 | 5.6′–8.0′ | | 4213 | 170′–240′ |
| 4204 | 8′–11′ | | 4214 | 240′–340′ |
| 4205 | 11′–16′ | | 4215 | 340′–480′ |
| 4206 | 16′–22′ | | 4216 | 480′–680′ |
| 4207 | 22′–30′ | | 4217 | 680′–1000′ |
| 4208 | 30′–42′ | | 4218 | 1000′–1400′ |
| 4209 | 42′–60′ | | 4219 | 1400′–2000′ |

The 4100-series is built from Tycho-2 and is the one to use for wide-angle images; the
5200-series ("LITE") is the modern Gaia-based replacement. Guidance is to install indexes
whose skymarks are "10% to 100%" of your field size.

Note the deep structural difference: astrometry.net's index encodes *which quads exist*,
committing to a specific quad selection at build time. ASTAP's approach re-derives quads at
query time from whatever stars it read, which is more flexible but means the image and
catalogue quad-selection rules must agree — a fragility we return to in §11.2.

---

## 6. Stage 4 — Consensus and outlier rejection

Even a 5-D descriptor produces some false matches, and triangles produce a flood. Every
solver therefore has a consensus stage that asks: *do these candidate matches agree on a
single global transform?*

### 6.1 Median-scale filter

The cheapest test. Each match implies a plate scale `ρ = d₁_image / d₁_catalogue`. True
matches share one `ρ`; false ones scatter. Take the median and keep matches within a
relative tolerance:

```
    ρ̃ = median{ρᵢ},      keep i  if  |ρᵢ − ρ̃| ≤ t · ρ̃
```

This is ASTAP's original filter and arcsec's `filter_by_scale`. It is 1-D and therefore
weak: false matches that happen to have the right scale survive.

### 6.2 Vote accumulators (a Hough transform in transform space)

Stronger: each match implies *both* a scale and a rotation. Bin `(scale, angle)` and take
the peak cell. arcsec's `vote_filter` implements this with 5% scale bins and 10° angle
bins, and — because the pixel→sky transform is usually reflected (`CDELT1 < 0`) — maintains
two accumulators simultaneously:

```
    direct    (det > 0):   key_angle = (φ_img − φ_cat) mod π
    reflected (det < 0):   key_angle = (φ_img + φ_cat) mod π
```

For a pure rotation the *difference* of the two quads' principal-axis angles is constant;
under a reflection the difference scatters but the *sum* is constant. Building both grids
and taking the larger peak handles either parity without knowing it in advance. This is a
genuinely nice trick and, as far as I can find, not something ASTAP does.

### 6.3 RANSAC

The general-purpose answer. Repeatedly sample a minimal set (3 correspondences determine a
similarity/affine transform), fit, count inliers within a pixel threshold, keep the best
model. astroalign uses a modified RANSAC with the Sampson distance and accepts a
transform matching "80% of the triangle matches or 10, whichever is lower", with individual
correspondences accepted at Sampson distance < 3 ([Beroiz et al. 2020][astroalign]).

RANSAC's advantage over a vote grid is that it works directly in *correspondence* space and
so tolerates an arbitrary inlier fraction; its cost is `O(k)` model fits where
`k ≈ log(1−p)/log(1−w³)` for inlier fraction `w`.

### 6.4 Bijectivity and sigma clipping

Two cheap add-ons arcsec uses: the bijective filter on the triangle path, sigma clipping
on both:

* **Bijective filter** — a true correspondence set is one-to-one. Reject any image feature
  matching several catalogue features and vice versa; keep only mutual best matches.
* **Iterative sigma clipping** — fit, compute residuals, drop everything beyond `k·RMS`,
  refit, repeat until stable. arcsec's `sigma_clip_pairs` uses a generous absolute cut
  (10 px in catalogue arcsec) on the first pass to kill gross outliers, then `3σ`
  thereafter, for at most 10 iterations. The first fit is unchecked (`fit_affine`): the
  outliers it exists to remove can skew that fit past the similarity check. Until
  2026-10 only the triangle path clipped; on the quad path a few wrong quads in the
  winning vote cell could make the plate fit fail the similarity check and abandon the
  right position (`obj_coalsack`: the plain fit of the first position's quads was
  refused; clipped, it verified 91 stars).

---

## 7. Stage 5 — Fitting the WCS

### 7.1 The gnomonic (TAN) projection

Nearly all solvers work on the **tangent plane**, because a gnomonic projection maps great
circles to straight lines and, crucially, makes an ideal telescope's mapping *exactly
linear* in the projected coordinates. The FITS convention is defined by
[Calabretta & Greisen 2002][wcs2].

Forward, from sky `(α, δ)` to standard coordinates `(ξ, η)` about tangent point `(α₀, δ₀)`
(this is arcsec's `equatorial_standard`, with `Δα = α − α₀`):

```
    D  = sin δ₀ · sin δ + cos δ₀ · cos δ · cos Δα          (the "projection factor")

    ξ  = −cos δ · sin Δα / D
    η  = −(sin δ₀ · cos δ · cos Δα − cos δ₀ · sin δ) / D
```

Dividing by `cdelt` (in radians/pixel) converts standard coordinates to pixels; arcsec
passes `cdelt = 1` and works in **arcseconds** throughout the matching stage.

Inverse, from `(ξ, η)` back to the sky (arcsec's `standard_equatorial`):

```
    α = α₀ + atan2( −ξ , cos δ₀ − η · sin δ₀ )

    δ = asin( (sin δ₀ + η · cos δ₀) / √(1 + ξ² + η²) )
```

The leading minus signs encode the FITS handedness convention: RA increases to the *left*
in a normally-oriented image, hence `CDELT1 < 0`.

The projection diverges at 90° from the tangent point, which is why a solver must
re-project to a tangent point near the field rather than using a single global one.

### 7.2 The six plate constants

With the catalogue in tangent-plane coordinates, the ideal mapping from pixels is affine:

```
    ξ = a·x + b·y + c
    η = d·x + e·y + f
```

Six unknowns, so three correspondences suffice, and more are solved by least squares.
This is ASTAP's `Xref := a*Xtest + b*Ytest + c` / `Yref := d*Xtest + e*Ytest + f`, and
arcsec's `PlateConstants`.

The `(a,b,d,e)` sub-matrix carries scale, rotation and parity:

```
    scale_x = √(a² + b²)        scale_y = √(d² + e²)         [arcsec/pixel]
    rotation = atan2(d, e)                                    [CROTA2]
    parity   = sign(a·e − b·d)                                [flipped if < 0]
```

arcsec sanity-checks the solution by requiring it to be a similarity transform: the
two singular values of the `(a,b,d,e)` matrix must agree.

```
    q = ½·√((a+e)² + (d−b)²)      r = ½·√((a−e)² + (d+b)²)
    σmax / σmin = (q + r) / |q − r|  ≤  1.08      else  ArcsecError::BadSolution
```

(`math::lsq::plate_anisotropy`, `MAX_PLATE_ANISOTROPY`.) Until 2026-10 the check
compared the two *row* norms, `0.9 ≤ (a² + b²)/(d² + e²) ≤ 1.1`, which a sheared matrix
passes: the tier-D false positive `neg_hint_3` (offset hint, 0.1.2) had
`a, b, d, e = 0.78, −2.62, −0.90, −2.48`, rows 2.73 and 2.64 long but columns 1.19
and 3.61, a singular-value ratio of 3.03. On the corpus the largest ratio of any correct
solve is 1.027 (TESS FFIs, whose distortion the linear plate absorbs as a little
anisotropy); every other correct solve is at most 1.0066.

The FITS `CD` matrix is the same thing in degrees:

```
    CD1_1 = −a/3600    CD1_2 = −b/3600
    CD2_1 = +d/3600    CD2_2 = +e/3600

    CDELT1 = −√(CD1_1² + CD1_2²)     CDELT2 = +√(CD2_1² + CD2_2²)
    CROTA2 = atan2(CD2_1, CD2_2)     [degrees]
```

### 7.3 Solving the least-squares system

The normal equations `AᵀA x = Aᵀb` are simple but square the condition number. Both ASTAP
and arcsec instead use **Givens rotations** — an incremental QR factorisation
(Montenbruck & Pfleger, *Astronomy on the Personal Computer*). Each new observation is
rotated into an upper-triangular accumulator:

```
    for each new row (x_i, y_i, 1 | b_i):
        for each existing column j:
            c = R_jj/h,  s = R_ji/h  with h = √(R_jj² + R_ji²)
            apply the 2×2 rotation to eliminate the new row's j-th entry
    back-substitute
```

This is numerically stable, needs only `O(k²)` storage for `k = 3` unknowns regardless of
the number of observations, and streams — you never form `A`.

The crucial detail for accuracy is **what the correspondences are**:

| Solver | Correspondence used for the fit | Count |
|---|---|---|
| ASTAP | quad **centroids** | `n_matched` |
| arcsec (quads) | quad centroids for the first fit, then individual **stars** matched by projecting the catalogue (`verify_and_refit`) | ≥ 30, typically 200–375 |
| astrometry.net | individual **stars** (A,B,C,D of each matched quad, then all verified stars) | up to hundreds |
| astroalign, SCAMP | individual stars after cross-match | hundreds–thousands |

Fitting on centroids is a real limitation. A quad centroid averages four centroid errors —
so it is *individually* more precise by √4 — but you get one constraint per quad instead of
four, and (worse) the centroid is insensitive to any distortion that is antisymmetric about
the quad, so distortion signal is partially cancelled rather than measured. arcsec uses the
centroid fit only as a starting point and replaces it with a star-level fit (§11.4).

### 7.4 Distortion

A real optical train is not affine. The two conventions in use:

* **SIP** — Simple Imaging Polynomial ([Shupe et al. 2005][sip]). Corrections are applied
  in *pixel* space before the linear `CD` term:

  ```
      u = x − CRPIX1,   v = y − CRPIX2

      f(u,v) = Σ_{p+q ≤ A_ORDER}  A_pq · u^p v^q
      g(u,v) = Σ_{p+q ≤ B_ORDER}  B_pq · u^p v^q

      [ξ]   [CD1_1 CD1_2] [u + f(u,v)]
      [η] = [CD2_1 CD2_2] [v + g(u,v)]
  ```

  with `AP_pq`, `BP_pq` giving the approximate inverse. `CTYPE` becomes `RA---TAN-SIP`.
  Order 2–3 is normal; higher orders overfit and extrapolate disastrously outside the
  fitted region.
* **TPV / PV** — polynomial coefficients applied in *intermediate world* space, used by
  SCAMP and the SExtractor lineage ([Bertin, SCAMP][scamp]). Convertible to SIP.

astrometry.net's `tweak2` fits SIP after the linear solve using an **annealing** strategy:
match image stars to catalogue stars within a radius that is small near the region of
confidence and grows with distance from it, so distant stars are picked up but down-weighted;
then progressively increase both the SIP order and the annealing scale in nested loops.
This is the standard way to bootstrap from "roughly right" to "sub-pixel everywhere".

For amateur refractors at 1–3″/px, distortion is usually below the noise; for fast
astrographs, camera lenses, and any field wider than ~3°, it is the dominant residual.

---

## 8. Stage 6 — Verification

A pattern match is a *hypothesis*. Verification decides whether to believe it, and is what
separates "usually works" from "never lies".

### 8.1 astrometry.net — Bayesian decision theory

The gold standard. Given a candidate WCS, project the index stars into the image and ask
whether the observed source positions are better explained by

* **F (foreground)**: a mixture of a uniform "anywhere in the image" term and a Gaussian
  blob around each projected index star, with width set by the combined positional
  variances; or
* **B (background)**: uniform probability everywhere.

The Bayes factor is the ratio of marginal likelihoods, accumulated over query stars:

```
    K = p(D | F) / p(D | B)
```

With the utility table `u(TP)=+1, u(FP)=−1999, u(FN)=−1, u(TN)=+1` and a prior odds of
`10⁻⁶`, Bayesian decision theory yields the acceptance threshold

```
    K > 10⁹          (log-odds ≈ 20.7)
```

Query stars are added one at a time until `K` crosses the threshold or the stars run out;
a typical solve reports something like `log-odds 35.95 (4.1e15), 31 match, 0 conflict,
70 distractors, 123 index`. The payoff: **99.9% success on survey data with no false
positives**, and remaining failures traceable to catalogue incompleteness rather than the
algorithm ([Lang et al. 2010][lang2010]). Notably, the documented false positives are not
random — they are images with genuine linear features (satellite trails, the ISS) matching
linear *artefacts* in the USNO-B plate scans.

### 8.2 tetra3 — probability of chance matches

Project the catalogue, count stars matching within `match_radius` (default 1% of the FOV),
and compute the probability that this many matches would arise by chance given the star
densities. Accept if below `match_threshold` (default 1e-3).

### 8.3 Match-count thresholds

The pragmatic version, and what ASTAP and arcsec use: accept if the number of consistent
matches exceeds a fixed or star-count-scaled threshold. arcsec uses

```
    catalogue path:  min_quads = 3 + n_stars_image / 140 agreeing quads to attempt a fit,
                     then ≥ MIN_VERIFIED_STARS = 30 individually matched stars spanning
                     ≥ MIN_VERIFY_SPREAD = 0.20 of the image half-diagonal to accept
    blind path:      MIN_VERIFY_SCORE = 18 verified stars, EARLY_STOP_SCORE = 20
```

Both paths are now genuinely verified — they project catalogue or index stars into the
image and count agreement. The catalogue path gained this on 2026-09-02; before that it
accepted the first spiral position that yielded enough quad matches. See §11.1.

---

## 9. Reference catalogues and indexes

| Catalogue | Stars | Mag limit | Epoch | Notes |
|---|---|---|---|---|
| **Gaia DR3** | 1.8 billion | ~21 | **2016.0** | The modern reference frame; source of ASTAP's databases and the 5200-series indexes |
| Tycho-2 | 2.5 M | ~12 | 2000.0 | astrometry.net 4100-series; right density for wide fields |
| 2MASS | 470 M | ~17 (IR) | ~1999 | astrometry.net 4200-series; poor for blue images |
| USNO-B | 1 B | ~21 | ~1950–1990 | The original astrometry.net catalogue; plate artefacts cause the known FP mode |
| UCAC4 | 113 M | ~16 | 2000.0 | Good proper motions |
| HYG | 120 k | ~9 | — | Tiny; useful for unit tests |

### 9.1 Density-limited vs magnitude-limited

ASTAP's databases are notable for being **sorted and truncated on star density, not
magnitude** — D50 means "up to 5000 stars per square degree", D20 = 2000, D05 = 500, D80 =
8000. This is the right choice for a solver: a magnitude cut gives you thousands of stars
in the galactic plane and a dozen at the pole, whereas a density cut gives a roughly
uniform number of quads per field everywhere on the sky, which is what the matching
statistics need. (The equal-area property of HEALPix serves the same goal in
astrometry.net's index.)

### 9.1b The three ASTAP database formats

arcsec reads all three. They share nothing but a naming convention, and which one a
database uses is dictated by the field sizes it targets.

| Format | Databases | Layout | Field range |
|---|---|---|---|
| `.1476` | D80, D50, D20, D05, V50 | 1476 tiles, 36 **equal-declination** rings, 5-byte packed records | 0.15°–6° |
| `.290` | G05, V05 | 290 tiles, 18 **equal-area** rings, *identical* 5-byte records | 3°–20° |
| `.001` | W08 | one all-sky file of f32 triples | 20°–80° |

`.1476` and `.290` differ only in the sky grid — same 110-byte ASCII header, same
packed records — so one record reader serves both. The 290 grid's boundaries are not
documented anywhere I could find; they were derived from the shipped G05 files and
satisfy

```
    sin(dec_k) = -1 + 2 · (cumulative ring weight up to k) / 289
```

with ring RA counts `1, 4, 8, 12, 16, 20, 24, 28, 32, 32, 28, 24, 20, 16, 12, 8, 4, 1`
and the two polar caps counting half a cell each. That reproduces the observed per-ring
declination ranges to better than 0.001°, and `dec_boundaries_are_equal_area` asserts it.

One consequence matters for correctness: the 1476 reader finds its areas by sampling
the field's **four corners**, which is only valid when the field fits inside one
declination ring — hence its hard cap at 5.14°. A 20° field spans many rings and many
RA cells, so the 290 path enumerates every overlapping area instead. Both then read
their areas together, one magnitude step at a time, and keep the field's brightest
`max_stars` stars (§11.11). Filling the budget area by area would take every star from
the southern edge of a 60° field and none from the north, because areas are visited in
declination order.

`.001` is a different thing again: a `u32` star count followed by
`{f32 magnitude × 10, f32 RA radians, f32 Dec radians}` triples, brightest first, no
tiling at all. Its layout is undocumented too; it was decoded from the file and
confirmed by the first two records resolving to Sirius (mag −1.5, RA 101.28°,
Dec −16.72°) and Canopus (mag −0.6, RA 95.99°, Dec −52.70°).

Because the ranges do not overlap much, `-D` is now optional: with no explicit
database arcsec picks the densest installed one whose published range contains the
field, falling back to the nearest range if none does. `-d` is optional too — it
defaults to the directory `arcsec catalog install` writes to. See
[catalogues.md](catalogues.md).

Which layout a database uses is not predictable from its name: V05 is `.290` despite
covering the same field range as the `.1476` D-series, so the installer probes for
whichever extension is present rather than assuming one.

### 9.2 Sky tiling: dec-rings vs HEALPix

ASTAP's `.1476` format divides the sky into 1476 tiles as 36 declination rings, each split
into a ring-dependent number of RA cells — the scheme arcsec reads in
`catalog/areas.rs`. It is simple, seek-friendly, and gives cells of *approximately* equal
area (the ring boundaries are chosen to compensate for `cos δ`).

HEALPix gives *exactly* equal areas and a hierarchical index (`N_pix = 12 N_side²`), at
the cost of a more complex pixel↔coordinate transform. For arcsec's purposes the dec-ring
scheme is adequate; the reason to care about HEALPix would be building our own index
(§12.5).

### 9.3 Epoch and proper motion — a real, silent error source

Gaia DR3 positions are at **epoch 2016.0**. An image taken in 2026 is 10 years later.
Proper motions are not applied by ASTAP's `.1476` format (it stores only RA, Dec and
magnitude in 5 bytes/star) and are not applied by arcsec.

Order of magnitude: the median stellar proper motion is a few mas/yr — negligible. But the
distribution has a long tail: Barnard's Star is 10.4″/yr (104″ over a decade), and there are
thousands of stars above 200 mas/yr (2″/decade). At 1–3″/pixel this shifts a handful of
stars by ~1 pixel, which is well inside typical quad tolerances and simply adds a little
noise. It matters when (a) chasing sub-arcsecond residuals, (b) solving in a sparse field
where one of your few stars is a high-PM object, or (c) using an older catalogue
(USNO-B at epoch ~1970 gives 56 years of drift).

---

## 10. What arcsec does today

### 10.1 Overview

arcsec implements the 5-ratio quad family (§4.2) introduced by ASTAP, with an online
spiral search (§5.2), plus additions of our own: star-level verification and re-fitting of
every candidate solution, an astrometry.net-index **blind** front-end (§4.3), and an
alternative **triangle** matcher (§4.1).

```
                            ┌─────────────────────────────┐
                            │  arcsec -f image.fits …     │
                            └──────────────┬──────────────┘
                                           │
                      ┌────────────────────▼─────────────────────┐
                      │ read FITS / XISF / ASDF                  │
                      │ normalise: replace NaN/Inf; rescale      │
                      │   float data spanning < 4096 counts      │
                      │ RA/Dec hint from --ra/--spd, else header │
                      │ pixel scale from --fov (image height),   │
                      │   else header optics, else 1″/px         │
                      │ FOV = scale × max(width, height)         │
                      └────────────────────┬─────────────────────┘
                                           │
                      ┌────────────────────▼─────────────────────┐
                      │ -z absent or 0: auto-bin if arcsec/px < 1│
                      │   factor = round(1 / arcsec_per_px), ≤16 │
                      │ -D absent: pick the densest installed    │
                      │   database whose FOV range fits          │
                      └────────────────────┬─────────────────────┘
                                           │
                     ┌─────────────────────┴──────────────────────┐
                     │ -i/--index given?                          │
                     └──────┬───────────────────────────┬─────────┘
                         yes│                        no │
              ┌─────────────▼──────────────┐            │
              │ BLIND (§10.4)              │            │
              │ rank index files by scale  │            │
              │ run best 2 in parallel     │            │
              │ → (α,δ) estimate           │            │
              │ narrow search radius to    │            │
              │   max(2·FOV, 5°)           │            │
              │ (on failure: keep the hint)│            │
              └─────────────┬──────────────┘            │
                            └────────────┬──────────────┘
                                         │
                      ┌──────────────────▼──────────────────┐
                      │ CATALOG SPIRAL SOLVE (§10.3)        │
                      └──────────────────┬──────────────────┘
                                         │
                      ┌──────────────────▼──────────────────┐
                      │ write <base>.wcs, <base>.ini,       │
                      │ optionally --update the FITS header │
                      └─────────────────────────────────────┘
```

### 10.2 Star detection (`detection/`)

```
   image ──► get_background()                    16-bit histogram over the frame
             ├─ background = histogram mode (or the mean if > 1.5 × mode)
             ├─ noise = RMS about it on a sparse pixel grid, 3σ-clipped,
             │          iterated until it changes < 5% (≤ 7 passes)
             └─ star_level, star_level2 = histogram tail thresholds
                    │
                    ▼
             find_stars_with_background()      detection cascade, stops once
             │                                 `-s` stars have been found
             ├─ level 4: threshold star_level   (if > 30 σ)
             ├─ level 3: threshold star_level2  (if > 30 σ)
             ├─ level 2: threshold 30 σ
             ├─ level 1: ~12×12 grid of cells, local σ-clipped background,
             │           threshold 7 σ_local
             │    each level scans horizontal bands in parallel
             │    (BAND_OVERLAP = 90 rows), merged with a positional dedup
             ├─ measure_star(): HFD in a 14-px annulus, sub-pixel bilinear
             │     centroid, SNR, flux; reject if not "boxed", single hot
             │     pixel, too large, or any result non-finite
             └─ sort by SNR, keep the top `-s` (default 500)
```

The detected list is not trimmed further. A brightest-half trim (`max(max_stars/2, 50)`)
used to guard the 3-nearest-neighbour quads against faint stars missing from the
catalogue; with 9-NN redundancy and star-level verification it only halved the quad
count, and removing it took tier B to 34/34 (see [test-images.md §6.2](test-images.md#62-what-moved-the-numbers)).
The blind front-end still applies it.

### 10.3 The catalogue spiral solve (`pipeline/solver.rs`)

```
  ┌──────────────────────────────────────────────────────────────────────────┐
  │ A. DETECT              up to `-s` stars (default 500), no further trim   │
  ├──────────────────────────────────────────────────────────────────────────┤
  │ B. BUILD IMAGE QUADS   build_quads(): each star plus its nearest         │
  │                        neighbours, all 4-subsets of that group           │
  │      n < 15  → group of 7 (star + 6 NN), C(7,4) =  35 quads per star     │
  │      n < 30  → group of 6 (star + 5 NN), C(6,4) =  15 quads per star     │
  │      else    → group of 9 (star + 8 NN), C(9,4) = 126 quads per star     │
  │                (IMAGE_NEIGHBOURS = CATALOG_NEIGHBOURS = 9)               │
  │      dedup: reject a quad whose centroid is within 1 px of an existing   │
  │             one (hash grid, 5-px cells, bucket capacity 10)              │
  ├──────────────────────────────────────────────────────────────────────────┤
  │ C. SPIRAL              max_distance = search_radius/FOV + 2              │
  │    positions from SpiralSearch, evaluated in batches of `--threads`      │
  │    (position 0 alone first); the lowest-index position that verifies     │
  │    wins, so the result equals the serial search                          │
  │        (0,0),(1,0),(1,1),(0,1),(−1,1),(−1,0),(−1,−1),(0,−1),(1,−1),…     │
  │                                                                          │
  │        δ_db = δ_hint + FOV·sy         (pole wrap → flip RA by π)         │
  │        α_db = α_hint + FOV·sx / cos(δ_db ∓ FOV/2)                        │
  │        skip if ang_sep(hint, trial) > radius + FOV/2                     │
  │                                                                          │
  │        read_catalog_stars(α_db, δ_db, FOV·oversize, N_required)          │
  │            oversize = 2.0 (n<35), 2·√(35/n), 1.0 (n>140)                 │
  │            N_required = max_stars · oversize²   (max_stars = `-s`)       │
  │            → mmap'd .1476/.290 area files or the .001 file               │
  │                                                                          │
  │        project catalogue → tangent plane (arcsec), sort by x             │
  │        build_quads_presorted() with the *image* star count               │
  │        sort_catalog_quads(): sort by ratios[INDEX_RATIO = 4]             │
  │                                                                          │
  │        find_matches_sorted(): binary-search ±t on ratios[4] over a       │
  │                               compact f32 copy, then check all five      │
  │        vote_filter(): 2-D (scale × angle) accumulator, both parities     │
  │        if fewer than min_quads survive, use filter_by_scale() on the     │
  │           raw matches instead when it keeps more                         │
  │                                                                          │
  │        if matches < min_quads = 3 + n/140:  next position                │
  │                                                                          │
  │        extract_star_pairs(): (image quad centroid, catalogue centroid)   │
  │        sigma_clip_pairs(): drop pairs off the consensus (first pass 10 px │
  │           or 3 × 1.48 MAD, then 3σ); if fewer than min_quads remain:     │
  │           next position                                                  │
  │        solve_plate_constants(): Givens-rotation LSQ, 6 constants;        │
  │           refused unless a similarity (singular-value ratio ≤ 1.08)      │
  │        verify_and_refit(): project every catalogue star through the      │
  │           plate, pair it with the nearest unused detected star within    │
  │           6 → 3 → 2 px, re-fit on those pairs at each radius             │
  │        reject unless ≥ 30 stars matched and their spread ≥ 0.20 of the   │
  │           image half-diagonal:  next position                            │
  ├──────────────────────────────────────────────────────────────────────────┤
  │ D. OUTPUT              derive_wcs(): tangent-plane inverse at the image  │
  │                        centre; un-scale CRPIX and CD/CDELT for binning   │
  └──────────────────────────────────────────────────────────────────────────┘
```

The `--method tetra` variant replaces build/match/filter with `build_triangles`,
`find_triangle_matches` (tolerance × 0.3), `bijective_filter`,
and `filter_triangles_by_scale`; the clipping, the plate fit and `verify_and_refit`
are shared.

The reported `RMS` is the per-star residual of the final `verify_and_refit` pass, and the
count written as `NQUADS` in the `.ini` (and "`n` of `m` quads selected" on stdout) is the
number of verified stars, not quads.

### 10.4 The blind solve (`pipeline/blind.rs`)

```
   detect stars (trimmed to max(max_stars/2, 50) when detection overflows)
                    ──► build image entries from the N_ENTRY_STARS = 30 brightest,
                        in astrometry.net code space (DIMQUADS 3 or 4, canonical
                        form: triangles CX ≤ 0.5; quads CX + DX ≤ 1 and CX ≤ DX)
                          │
                          ▼
   for each parity (normal, then flipped; skip the second if the first ≥ 20):
     scale filter — keep entries whose |AB| in pixels lies within the index's
                    angular scale band mapped through the image scale (×0.8, ×1.2)
                          │
                          ▼
     find_code_matches_into() over the compact f32 code array
     (9 MB, L3-resident) rather than the 70 MB entry array
                          │
                          ▼
     for each (image, index) code match:
        solve_plate_constants on the 3–4 star correspondences
        derive_wcs → an (α, δ) hypothesis
        vote into 0.1° sky bins                  ← VOTE_STEP
                          │
                          ▼
     for every vote cell, in descending vote order:
        take the cell's first hypothesis and run verify_score():
        project the index stars in a ±1.1·FOV declination band into the
        image and count those landing within MATCH_PX = 5 px of a detected
        star (excluding the quad's own stars, which would always match and
        would inflate false positives just as much as true ones)
        stop early once a score reaches EARLY_STOP_SCORE = 20
                          │
                          ▼
   accept if best score ≥ MIN_VERIFY_SCORE = 18
   → feed (α, δ) to the catalogue solver with radius = max(2·FOV, 5°)
```

The `arcsec` binary ranks candidate index files (`collect_index_files`) by how close their
scale band's midpoint is to `FOV/2`, discards any with no overlap with
`[0.2·FOV, 1.5·FOV]`, and runs the best `BLIND_MAX_INDEXES = 2` on separate threads
(fewer if `--threads` is lower), taking the highest-scoring result. If every index fails,
it prints a warning and runs the catalogue solve from the original hint and radius.

### 10.5 Constants and defaults

| Constant | Value | Where | Meaning |
|---|---|---|---|
| `quad_tolerance` | 0.007 | `-t` | per-ratio match tolerance |
| `max_stars` | 500 | `-s` | detection cap; also sets catalogue depth |
| `hfd_min` | 1.5″ | `-m` | minimum star size; converted to binned pixels, floor 0.8 px |
| `search_radius` | 180° | `-r` | spiral radius |
| binning | `round(1/arcsec_per_px)`, ≤ 16 | `-z` absent or 0 | auto downsample when `arcsec/px < 1`; any factor is capped so the binned image keeps ≥ 2 px a side |
| `IMAGE_NEIGHBOURS` / `CATALOG_NEIGHBOURS` | 9 | `quads/build.rs` | quad group size for ≥ 30 stars |
| `INDEX_RATIO` | 4 | `quads/match.rs` | ratio the catalogue quads are sorted and searched on |
| `min_quads` | `3 + n/140` | `solver.rs` | agreeing quads needed to attempt a fit |
| `oversize` | 2.0 → 1.0 | `solver.rs` | catalogue window vs FOV |
| `VERIFY_RADII` | 6, 3, 2 px | `solver.rs` | star-level verification match radii |
| `MIN_VERIFIED_STARS` | 30 | `solver.rs` | stars that must agree to accept a position |
| `MIN_VERIFY_SPREAD` | 0.20 | `solver.rs` | spread of those stars, fraction of the half-diagonal |
| `MAX_PLATE_ANISOTROPY` | 1.08 | `math/lsq.rs` | largest singular-value ratio σmax/σmin of a plate fit's linear part; every fit through `solve_plate_constants` (quad, star-level, blind) must be this close to a similarity |
| `TETRA_TOL_FACTOR` | 0.3 | `quads/tetra.rs` | triangle tolerance scaling |
| `SCALE_STEP` / `ANGLE_STEP` | 0.05 / 10° | `quads/vote.rs` | vote bin sizes |
| `BAND_OVERLAP` | 90 rows | `detection/stars.rs` | overlap between parallel detection bands |
| `MIN_VERIFY_SCORE` | 18 | `blind.rs` | blind acceptance |
| `EARLY_STOP_SCORE` | 20 | `blind.rs` | blind early exit |
| `MATCH_PX` | 5.0 | `blind.rs` | verification match radius |
| `VOTE_STEP` | 0.1° | `blind.rs` | sky vote bin |
| `N_ENTRY_STARS` | 30 | `blind.rs` | brightest detected stars used to build blind quads |
| `BLIND_MAX_INDEXES` | 2 | `arcsec` binary | parallel index files |

### 10.6 Measured performance

On the 103-image benchmark corpus (`scripts/benchmark.py --auto-db --astap`, 2026-09-25):
90 correct and 0 false positives, against ASTAP's 47 and 0. On the 45 images both solve
correctly, median centre error is 0.690″ (ASTAP 0.660″), median corner error 0.897″ (ASTAP
0.938″), and median time per image is the same, 0.15 s each, when run one image at a time
(`--jobs 1`). Details in [test-images.md §6](test-images.md#6-results--103-images-arcsec-vs-astap).

Earlier, from `ARCSEC_VS_ASTAP.md`, on a 240-frame NGC 3372 set (ASI533MC Pro, RedCat 51,
3.13″/px): positional agreement with ASTAP better than 0.3″ on every frame checked, pixel
scale agreement better than 0.005%, and 3–4× faster than ASTAP (0.07–0.10 s vs 0.25–0.35 s
per frame). That predates the 9-NN quads and star-level verification, which made each
solve slower (0.60× ASTAP) until the optimisations in test-images.md §6.5 brought it back
to parity.
For context, [AstroKeith][astrokeith] measures astrometry.net at <1 s on a Pi 5 with an
accurate scale hint, tetra3 at ~200 ms and cedar-solve at ~12 ms. arcsec is firmly in the
"fast hinted solver" class.

---

## 11. Shortcomings

Ordered roughly by how much they limit us.

### 11.1 The catalogue path performs no verification — FIXED 2026-09-02

The problem as first written, before the fix below:

`solve_image` returns the **first** spiral position that yields `min_quads = 3 + n/140`
consistent quad matches. With `n = 250` that is 4 quads. There is no projection of the
catalogue back into the image, no count of agreeing stars, no odds ratio — nothing that
could distinguish a real solution from four coincidences.

The 5-ratio descriptor plus the vote filter makes this *fairly* safe (§4.4 puts the
expected random-match count at ~10⁻⁴ per position), but "fairly safe" compounds over a
spiral: a 30° radius at 1° FOV visits ~2800 positions, so the per-solve false-positive
probability scales with the search area. And `filter_by_scale` is used as a *fallback*
when the vote peak is below `min_quads` — a path with markedly weaker guarantees that we
take precisely when evidence is weakest.

By contrast the blind path *is* verified (`MIN_VERIFY_SCORE = 18`), and the git history
records that this threshold was introduced specifically to eliminate observed false
positives. The same reasoning applies to the catalogue path and has not been applied.

**Fixed**: `verify_and_refit` in `solver.rs` now projects every catalogue star through the
candidate plate, pairs each with the nearest detected star and re-fits on those pairs
(radii 6 → 3 → 2 px), accepting a position only when at least 30 stars agree *and* they span
at least 0.20 of the image half-diagonal. False positives across the 103-image corpus went
from 4 to **0**, and recall rose at the same time (64 → 86 correct) because the re-fit
rescues marginal positions. See
[test-images.md §6.2](test-images.md#62-what-moved-the-numbers).

Historical note, which is why the fix was worth doing. Measured on the 103-image corpus,
**ASTAP returned zero wrong answers and arcsec returned four** — and on two of those ASTAP
solved the same image correctly. The 2026-09-01 run demonstrates the mechanism: 3 of 29
reported tier-A solves were **wrong**, the worst by 0.57° at the image corners
(`pole_scp`, `CROTA2` off by 31.23°). Near the pole the spiral's RA step is `FOV / cos δ`
— 43° per step at Dec −88° — and because the 5-ratio descriptor is rotation invariant, a
ring of circumpolar stars matches itself under that rotation. Six quads agreed, the vote
filter found a perfectly coherent peak, the plate fit converged, and a half-degree error was
returned as a success. That field (`pole_scp`) now solves correctly to 2.3″; see
[test-images.md §6.1](test-images.md#61-where-it-started-and-where-it-stands).


### 11.1b Any FITS containing NaN pixels crashes the solver — FIXED 2026-09-01

Confirmed on a Pan-STARRS cutout (1058 NaN pixels out of 1 048 576):

```
thread 'main' panicked at library/core/src/slice/sort/shared/smallsort.rs:854:
user-provided comparison function does not correctly implement a total order
  ...
  10: arcsec_core::detection::stars::measure_star
  11: arcsec_core::detection::stars::detect_pass
```

The culprit was `median_f64` in `detection/stars.rs` (it has since been rewritten as a
quickselect):

```rust
v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
```

`partial_cmp` returns `None` for NaN, and mapping that to `Equal` breaks transitivity.
Rust's current `driftsort` detects the violation and **panics** (exit 101) rather than
returning garbage. Replacing the NaNs with zeros in the same file turns the crash into a
clean "No solution found", confirming the diagnosis.

NaN is routine in real data — drizzled and reprojected images, chip gaps, masked columns,
cutouts that extend past the survey edge.

**Fix applied**, in three layers:

1. Every float comparator in production code now uses `f64::total_cmp`, which *is* a total
   order for NaN, so no sort can panic regardless of what reaches it. This replaced the
   `partial_cmp(...).unwrap_or(Equal)` pattern in `stars.rs`, `solver.rs`, `build.rs`,
   `match.rs`, `tetra.rs`, `blind.rs`, `anet.rs` and `main.rs`.
2. `measure_star` now returns `None` when any of the resulting `x`, `y`, `snr` or `hfd` is
   non-finite — a NaN patch is not a star, and should not become one.
3. `build_histogram` skips non-finite pixels rather than letting `as usize` saturate them
   to bin 0, which used to bias the background estimate downward.

Non-finite pixels are also replaced up front (see §11.1c). Regression tests
`measure_star_rejects_nan_neighbourhood`, `find_stars_survives_nan_pixels` and
`median_f64_ignores_nothing_but_never_panics_on_nan` cover this.

### 11.1c Detection assumed 16-bit-integer pixel values — FIXED 2026-09-01

`get_background` builds a 65536-bin histogram with

```rust
let v = img.get(x, y) as usize;      // background.rs:22
```

The `f32 → usize` cast truncates toward zero and saturates negatives to 0. Modern survey
data is float-valued in physical units — SDSS frames are in nanomaggies spanning roughly
`[-0.15, 4.5]`, DESI Legacy Survey cutouts `[-0.008, 27]` — so every pixel lands in bins
0–4 and the estimator collapses:

```
Legacy Survey cutout, as delivered:
  0 stars found of the requested 500. Background value is 0.
  Detection level used 1 above background. Noise level is 0.
```

Multiplying the same file's pixels by 3000 and adding an offset of 1000 restores it:

```
same file, rescaled:
  982 stars found of the requested 500. Background value is 1000.
  Detection level used 422 above background. Noise level is 4.
```

and an SDSS frame that returned "insufficient stars" as delivered **solves correctly** once
rescaled. So this is purely a dynamic-range assumption, not a limitation of the algorithm.

The failure was silent and misreported: the user saw `exit 2, insufficient stars`, which
suggests a sparse field rather than a units problem.

**Fix applied**: `ImageBuffer::normalize_for_detection` (`types.rs`), called by the CLI
immediately after the image is read and before binning. It replaces non-finite pixels with the
finite minimum, then:

* leaves the data **exactly** as-is when it already spans ≥ 4096 counts, so camera output
  and 16-bit survey images keep their existing behaviour bit-for-bit;
* otherwise maps the 99.9th percentile to ~20000 counts (the percentile rather than the
  maximum, so a single saturated star cannot swallow the dynamic range), leaving background
  and noise several hundred counts wide.

Measured effect on the benchmark corpus: the three DESI Legacy Survey cutouts and one SDSS
frame went from `exit 2, insufficient stars` to solving, at **0.058″–0.225″ centre error** —
the most accurate results in the whole corpus. Five unit tests in `types.rs` cover the
ADU-passthrough, small-range, NaN, all-NaN and flat-frame cases.

Worth noting: **ASTAP fails on these same files** (`Only 0 stars found in image. Abort`), so
this is now a capability arcsec has and its reference implementation does not.

### 11.2 The image and catalogue quad sets are built by incompatible rules — MITIGATED 2026-09-02

For `n ≥ 60` stars, `build_quads` uses only each star's **3 nearest neighbours** — one quad
per star. Which three stars are nearest depends entirely on *which stars are in the list*.
The image list is the detected stars; the catalogue list is the `n_max · oversize²`
brightest catalogue stars in a window. These sets differ — the image has undetected faint
stars, the catalogue has stars below the detection limit — so a star's 3-NN neighbourhood
in one list is frequently not its neighbourhood in the other, and the two quads simply do
not correspond.

The codebase was visibly fighting this. The `oversize` heuristic, the "trim to the brightest
half" rule (since removed), and the comment

> When detection finds many more stars than requested the extra stars are faint enough to
> be absent from the catalogue, which corrupts 3-NN quads.

are all compensations for a structural mismatch. astrometry.net avoids it entirely by
enumerating quads in a **scale band** with deliberate redundancy (each star in up to 8
quads, ~16 passes), so a matching quad survives even when several stars are missing.

Practical consequence: solve reliability is unusually sensitive to the interaction of
`-s`, the detection threshold and the catalogue density — exactly the parameters a user is
least equipped to tune.

**Mitigated** by building all C(9,4) subsets of each star's 9 nearest neighbours instead of
a single 3-NN quad, so a correspondence survives when a neighbour is missing from one of the
two lists. Tier A went 26 → 55 of 63 and tier B 22 → 34 of 34. The underlying mismatch is
still there — the fix buys redundancy rather than removing the dependence on star selection —
and the 8 remaining tier-A failures are all crowded or extended-object fields.

Historical note. Measured against ASTAP before the fix, this was the single biggest
deficit: ASTAP solved 39 of 62 tier-A images to arcsec's 26. Comparing the two solvers'
logs shows arcsec consistently building about **half** ASTAP's quads, because of the
trim in `solve_image`:

| image | arcsec stars used | arcsec quads | ASTAP stars | ASTAP quads |
|---|---|---|---|---|
| dens_bootes | 250 | 203 | 498 | 384 |
| type_m44 | 250 | 192 | 501 | 402 |

Raising `-s` to 1000 so the trim leaves 500 stars matches ASTAP's counts and recovers 2 of
9 tested failures — but 7 still fail at matched star *and* quad counts, so the trim is a
real cost on top of, not instead of, the correspondence problem described above. The trim
was removed once verification existed. See
[test-images.md §6.2](test-images.md#62-what-moved-the-numbers).

### 11.3 We cannot solve without a good pixel-scale estimate

The spiral **step size is the FOV**. If the FOV estimate is wrong by 2×, the steps are wrong
by 2× and the catalogue window is wrong by 2×, so the correct position is stepped over.
Worse, when neither `--fov` nor `FOCALLEN`/`XPIXSZ` is available, the CLI silently
assumes **1 arcsec/pixel**:

```rust
let ps = image_io::read_pixel_scale(file).unwrap_or(1.0);
```

For a 3000-px frame that asserts a 0.83° field. A DSLR-and-lens frame is 20°+; the solve
cannot succeed and the failure mode gives the user no hint that the scale was the problem.
There is no scale search, no scale refinement, and no warning. astrometry.net sweeps scale
bands by design; tetra3 takes an `fov_estimate` with an explicit `fov_max_error`.

Still open as of 0.1.0. A related trap, now fixed: `--fov` used to be read as the
**larger** image dimension (pixel scale `--fov / max(width, height)`), but ASTAP defines
`-fov` as the image *height* and N.I.N.A. sends exactly that (`FoVH`), so on a landscape
frame the scale came out low by the aspect ratio. The pixel scale is now
`--fov / height`; the internal field size used for database selection and the search
window is still the larger dimension.

### 11.4 The plate fit uses quad centroids, not stars — FIXED 2026-09-02

`extract_star_pairs` returns quad **centre** positions. Consequences:

* With `n_matched` quads we get `n_matched` constraints instead of `4·n_matched`.
* The reported `RMS` is a *quad-centroid* residual, which understates the real per-star
  residual and is not comparable with any other solver's reported RMS.
* Distortion that is antisymmetric within a quad cancels in the centroid, so the residual
  cannot detect it — which in turn means we cannot *fit* distortion (§11.5) without first
  changing this.
* Accuracy is capped well above what the centroid precision would allow.

The fix is standard and not large: after the linear fit, project catalogue stars into the
image, match each to the nearest detected star within a few pixels, and refit on those
pairs. That single step would also supply the verification statistic missing in §11.1 —
the same computation serves both purposes.

**Fixed** by the same change as §11.1: `verify_and_refit` re-fits on individual star pairs
(typically 200–375), and the reported `RMS` is now that per-star residual. The first fit
at each spiral position is still built from quad centroids by `extract_star_pairs`; it is
only a starting point.

### 11.5 No distortion model — FIXED for `--sip`

Every solution used to be a pure 6-parameter affine map. That is fine at 3″/px on a 4.6°
field with a well-corrected refractor; it is not fine for fast astrographs, camera lenses,
or anything wider than a few degrees, where field curvature and barrel distortion produce
radial residuals of several pixels at the corners.

**Fixed** for `--sip` (and `--extract2`, which implies it): `wcs::sip::fit_sip` fits
third-order SIP polynomials to the verified star pairs, which `solve_image` now returns in
`WcsSolution::matched_stars`. See §12.4 for the method and the measurements. The default
solve is still linear, as ASTAP's is.

### 11.6 Several CLI flags are accepted and ignored — FIXED

This used to read "accepted and ignored": a script that passed `--analyse 10` got a full
(slow) solve and no CSV, silently. The flags were then refused with an error until they
were implemented; now every ASTAP option does what `astap_cli` does (§12.6):

| Flag | Status |
|---|---|
| `--sip` | fits SIP distortion, when significant (§12.4) |
| `--check` | the Bayer check-pattern filter, with `-check y` as in ASTAP |
| `--analyse` | median HFD and star count, no solve |
| `--extract` / `--extract2` | the star list as CSV; `--extract2` solves first and adds RA/Dec |
| `--speed` | `auto` or `slow` (a catalogue window of twice the field at every position) |
| `--wcs` | accepted and ignored; the help text says the `.wcs` file is `{always written}` |
| `-f` help text | ~~claims "fits, tiff, png, pbm, jpg"~~ — **fixed**: reads FITS, XISF and ASDF, and the help text says so |

### 11.7 Documentation drift

`ARCSEC_VS_ASTAP.md` states that "the `-z` CLI flag … is parsed but not applied — the image
is always passed to the solver at full resolution", and `FUTURE_IMPROVEMENTS.md` repeats
it. That is no longer true: the CLI implements binning, including an auto mode, and
`solver.rs` correctly un-scales `CRPIX`/`CD`/`CDELT` afterwards. The accuracy table in that
document attributes a ~0.2″ offset to ASTAP binning and arcsec not binning, an explanation
that no longer holds.

Relatedly, the auto-binning comment in the CLI source used to describe a rule the code did
not implement ("Bin when height > 2500 px OR pixel scale < 1 arcsec/px"; only the
pixel-scale test exists). That comment has since been replaced by an accurate one on
`choose_binning`.

### 11.8 Catalogue and epoch limitations

* We depend on ASTAP's star databases (`.1476`, `.290`, `.001`). The format is
  documented by GPL source and the data is Gaia, so there is no legal barrier, but there
  is a practical dependency on a third-party download. `arcsec catalog install` now
  fetches them from ASTAP's distribution and puts them where the solver looks, which
  removes the friction but not the dependency.
* No proper-motion correction (§9.3). The `.1476` record has no room for it.
* No colour/magnitude information is used, so we cannot weight the fit by expected
  detectability or reject a match on implausible photometry.

### 11.9 Search-cost scaling and thread usage

The spiral is `O((r/FOV)²)` positions and each position rebuilds catalogue quads from
scratch. The blind front-end exists precisely to avoid this, but it needs astrometry.net
index files — so a user with only the ASTAP database and no position hint has no fast path.

The thread-usage half is **fixed**: spiral positions are evaluated a batch at a time across
`--threads` workers (lowest spiral index wins, so the result matches the serial search),
and detection, the background histogram and the pixel-range scan are parallel too. The
blind stage still runs at most `BLIND_MAX_INDEXES = 2` index files concurrently.

### 11.10 Smaller items

* `find_matches` (brute force) is retained and exported alongside `find_matches_sorted`;
  `FUTURE_IMPROVEMENTS.md` describes an unimplemented hash path with ASTAP's tuned
  constants (`nrquads ≥ 120`, `hash_bins = round(1/tolerance) + 2`).
* `.ini` output omits `CRPIX1/2` and the `CD` matrix that ASTAP writes; the `.wcs` file
  contains only WCS cards where ASTAP dumps the full original header.
* `CROTA1` is written equal to `CROTA2`, which is conventional but not strictly correct for
  a skewed CD matrix.
* `run_blind_pass` verifies **every** vote cell (sorted by vote count, with an early stop
  once a score reaches 20), not a top-K subset as `blind.rs`'s header comment once
  claimed (the comment has since been corrected). On an image that will not solve, every
  cell is verified, which
  is why blind failures are much slower than blind successes — `bench_all.sh` already
  reduces concurrency to 4 for `quads+blind` because of this.
* Blind verification tests only `hyps[0]`, the first hypothesis deposited in each vote
  cell, rather than the cell's best or a consensus of its members. A cell can therefore
  hold the right answer and be scored on the wrong member.
* No unit test covers an end-to-end solve — there is no `tests/` directory, and the in-file
  unit tests cover components (spiral order, LSQ, areas, coordinate round-trips) but never
  the pipeline. `scripts/benchmark.py` over the 103-image corpus is the integration
  coverage, and it needs data that is not in the repository (fetched by
  `scripts/fetch-test-images.sh`, plus the star databases).

---

### 11.11 The 1476 catalogue read starves every tile but the first — FIXED 2026-10-01

`read_catalog_stars_1476` read the (up to four) database tiles a field overlaps one after
another, each up to the whole star budget, and stopped once the budget was full. The first
tile usually filled it, so a field that straddles a tile boundary got catalogue stars on
one side only. The `.290` reader had the opposite fault: it gave every tile a fixed share,
`(max_stars / n_tiles).max(16) * 2`, so a tile covering most of a G05 field was
under-sampled and the tiles clipping its edges filled the budget with fainter stars.

It rarely stopped a solve outright, since half a field of stars is plenty to match, but
the fit then rested on half the frame and extrapolated to the other. Found while testing
`--sip`: in 26 of the 90 v1 solves the verified stars left at least one cell of a 3×3
grid over the frame empty (`decp80`, `type_ngc7000`: the right half or more). `fit_sip`
refuses to fit a cubic to such a set for that reason.

**Fix.** Both layouts now go through one reader, `read_brightest`, which returns the
brightest `max_stars` stars inside the field window whatever the tile layout. Every tile
the field overlaps is memory-mapped, and since each tile is stored brightest first in
groups of one magnitude step (0.1 mag), the tiles are read *in step*: the brightest unread
group of every tile, then the next, until the window holds `max_stars` stars; the union
is then sorted by magnitude and cut. Each tile is read only as deep as the field's own
magnitude limit, so it costs no more than the single-tile read did — a full-spiral
no-solve run got about 10% faster, because stopping at the field's limit reads fewer
records than filling the budget from one tile. A field inside one tile gets exactly the
stars it got before.

ASTAP instead shares the budget in proportion to each tile's share of the field
(`frac1..frac4` from `find_areas`). That gives the field's brightest stars only where the
sky is uniformly dense; reading by magnitude gives them everywhere, including across the
steep density gradients of the galactic plane. `find_areas_1476` still computes the
fractions; nothing uses them now.

After the fix 7 of the 92 v1 solves leave a grid cell empty, against 26 of 90 before. The
remaining seven are not catalogue-read gaps: in `ps1_big_c`, for example, the whole
window holds 297 catalogue stars, fewer than the budget, so every one is read and the
empty cells are empty in the image. On the expanded corpus the fix is worth +29 correct
with the true centre as hint and +30 with the hint 0.3 fields off, with fewer false
positives in both (5 → 4 and 9 → 4); [test-images.md §7](test-images.md#7-results--expanded-corpus-635-entries)
has the breakdown, and §7.7 the images that changed.

## 12. Improvement roadmap

Ordered by (value ÷ effort) as originally written. Items 1–3 are the ones I would do first,
and they interlock: one piece of machinery — *project the catalogue and match individual
stars* — fixes the verification gap, the accuracy cap and the distortion blocker at once.

Status as of 0.1.0: §12.1 and §12.7 are done, §12.5 option 1 and §12.6 are partly done,
§12.10 is partly done; the rest are open.

### 12.1 Add a star-level refit and verification pass ★ highest value — DONE 2026-09-02

After `solve_plate_constants` succeeds at a spiral position:

```
    1. project all catalogue stars in the window through the candidate WCS
    2. for each, find the nearest detected star; keep pairs within ~2–3 px
    3. refit the 6 plate constants on those pairs (Givens LSQ, already written)
    4. iterate 1–3 twice with a shrinking radius
    5. n_verified = number of pairs at convergence
    6. accept only if n_verified ≥ threshold, else continue the spiral
```

Delivers, in one change: a real verification statistic (§11.1), a 4–10× increase in fit
constraints and a genuine per-star RMS (§11.4), and the residual field needed to fit
distortion (§11.5). Reuses `equatorial_standard`, `solve_plate_constants` and the star list
already in hand; the only new code is a spatial lookup over detected stars (a uniform grid
is plenty at ~500 stars).

The acceptance threshold should be set the way `MIN_VERIFY_SCORE` was: measured against
the HiPS test fields, tuned so that no known-wrong field passes.

**Done** as `verify_and_refit` in `solver.rs`, with three passes at 6, 3 and 2 px. The
thresholds were tuned on the benchmark corpus rather than the HiPS blind fields: at least
`MIN_VERIFIED_STARS = 30` matched stars, spanning at least `MIN_VERIFY_SPREAD = 0.20` of
the image half-diagonal. See §11.1 and §11.4 for the effect.

### 12.2 Report and use a proper odds ratio

With §12.1 in place, replace the bare count with the astrometry.net-style Bayes factor
(§8.1). Even a simplified version — comparing a Gaussian-blob foreground against a uniform
background, accumulated over detected stars — makes the accept/reject decision principled
and scale-free, and gives a `LOGODDS` value to write into the `.ini` for downstream tools.

### 12.3 Handle an unknown or wrong pixel scale

Three levels, in increasing effort:

1. **Warn.** When neither `--fov` nor `FOCALLEN`/`XPIXSZ` is available, say so on stderr
   instead of silently assuming 1″/px. Cheap, and removes a whole class of confusing
   failures.
2. **Refine.** After the first successful fit, compare the fitted `CDELT` with the assumed
   one; if they differ by more than a few percent, re-run with the fitted scale.
3. **Search.** Sweep a geometric ladder of scale hypotheses (√2 steps, as astrometry.net's
   index scales do) around the assumed value, ordering by likelihood. Naturally
   parallelisable across the ladder.

### 12.4 Fit SIP distortion and honour `--sip` — DONE

Once star-level correspondences exist (§12.1), fitting SIP is a linear least-squares
problem in the polynomial coefficients — the same `lsq_fit` with more columns.

**Done** (`arcsec-core/src/wcs/sip.rs`), following ASTAP's `add_sip` so the keywords mean
what `astap_cli -sip`'s do: a full cubic (all ten terms including the constant and linear
ones, `A_ORDER = 3`) from pixel offsets to the offsets the linear WCS puts the catalogue
stars at, with the inverse (`AP`, `BP`) fitted separately from the same pairs; `CTYPE`
becomes `RA---TAN-SIP`; the keywords go in the `.wcs` file and, with `--update`, the
header, in ASTAP's order. Where ASTAP fits quad centroids, arcsec fits the individual
verified stars (typically 100–300). One round of 3σ clipping precedes the final fit.

Unlike ASTAP, the fit is kept only if it is warranted:

* the matched stars must reach every cell of a 3×3 grid over the frame, since a cubic
  extrapolates wildly beyond its stars (before the §11.11 fix a third of the v1 solves
  fell short of this; now 7 of 92 do);
* the 14 extra terms must pass an F-test against a linear fit (F ≥ 4);
* no corner may move by more than 5% of the half-diagonal.

The F-test is what the benchmark demanded. The corpus is survey data, reprojected and so
distortion-free; there, an unconditional cubic fits only centroid noise, which is largest
in the corners. Fitted regardless (with only the coverage test), `--sip` made the worst
corner worse on 57 of 89 solves (median 0.85″ → 1.01″, max 2.96″ → 4.58″, one image past
the 5″ false-positive line). `astap_cli -sip` does the same to its own solutions: median
0.97″ → 1.40″ over its 39 tier-A solves, three past 5″. With the F-test every corpus image
stays linear (F from 0.4 to 3.7) and the `--sip` benchmark is identical to the default one.
On frames with real distortion — corpus images warped by a known radial distortion, the
truth then being TAN plus that SIP — it is found with F from 18 to 225:

| Frame | Distortion at corners | arcsec | arcsec `--sip` | astap_cli | astap_cli `-sip` |
|---|---|---|---|---|---|
| `ra065` | 3 px | 2.73″ | 1.02″ | 2.63″ | 1.21″ |
| `ra065` | 8 px | 6.13″ | 1.23″ | 6.92″ | 8.95″ |
| `type_m101` | 20 px | 15.53″ | 2.21″ | 15.76″ | 1.52″ |

(worst-corner error against the truth; grid RMS falls from 1.2–6.6″ to 0.8″.) Still open:
astrometry.net's `tweak2` shape — re-match with the improved WCS, then raise the order —
would find more stars in a strongly distorted field's corners, where the linear WCS
misses them by more than the 2 px verification radius.

### 12.5 Make quad selection robust to differing star sets

The deeper fix for §11.2. Options, cheapest first:

1. **Build image quads from several neighbour counts** (3-NN *and* 4-NN/5-NN combinations)
   so that a quad survives when one neighbour is missing from the catalogue. Costs more
   quads to match, but `find_matches_sorted` is `O(M log R)` so this is affordable.
   **Done** in the form of all C(9,4) subsets of each star's 9-star neighbourhood (§11.2).
2. **Match the star sets by count, not by magnitude cut.** Choose the catalogue depth so
   that the *number* of catalogue stars in the field equals the number of detected stars,
   rather than over-reading by `oversize²` and hoping.
3. **Build a proper offline index** in the astrometry.net style — scale-banded quads with
   redundancy — which removes the online quad-construction mismatch entirely and makes
   hinted and blind solving the same code path. Fully designed and costed in
   **[offline-index.md](offline-index.md)**: format, builder algorithm, user-facing CLI,
   sizing and staged plan.

   Two things changed since this was first proposed. The size is far less alarming than
   `FUTURE_IMPROVEMENTS.md`'s 10–40 GB — sizing grid cells to the quad diameter rather
   than a third of the field puts a typical rig at **87–350 MB**, built from the ASTAP
   database the user already has, with no new download. But the *motivation* is weaker,
   because item 1 above (implemented as 9-NN redundancy) plus star-level verification
   already absorbed most of the damage, and all 8 remaining tier-A failures
   are detection problems rather than indexing problems. See
   [offline-index.md §1](offline-index.md#1-read-this-first-the-case-is-weaker-than-it-was)
   and its §9, a one-day experiment that settles it before committing two weeks.

### 12.6 Finish or remove the ignored CLI flags — DONE

`--wcs` is accepted and the `.wcs` file always written; the rest now behave as in
`astap_cli` (read from its source, `astap_command_line.lpr`, and checked against the
binary):

* `--analyse` / `--extract` (`detection::analyse`): ASTAP's `analyse_image` — up to four
  detection passes at falling thresholds on the full-resolution image, each starting
  afresh, stopping at the first that finds `-s` stars. stdout is `HFD_MEDIAN=` (one
  decimal) and `STARS=`; `--extract` writes `<image>.csv` (never moved by `-o`). No
  `.ini`, no solve; on Windows `--analyse`'s exit code is
  `round(HFD × 100) × 10⁶ + stars`, as ASTAP's. The per-star measurement is the solver's,
  with one switch: ASTAP's disc test (35% of `(2r − 2)²`, where the solver keeps its
  stricter `(2r)²`). On `type_m101` at SNR 20, 624 of ASTAP's 769 rows are identical to the
  last digit and the counts differ by 2–3% (733 vs 713 at SNR 30).
* `--extract2`: solve (with SIP, as ASTAP forces), then the same CSV with RA and Dec
  through the solution, whether or not the solve succeeded (the header's WCS, if any, when
  it did not). RA/Dec agree with ASTAP's to 0.12″ median.
* `--speed slow`: a catalogue window of twice the field at every position
  (`SearchSpeed::Slow`), capped at one database tile, and `Speed: slow` on stdout.
* `--check y`: `ImageBuffer::check_pattern_filter`, which scales the four Bayer phases to
  the brightest's mean over the central quarter, skipped for a colour image. As in ASTAP it
  needs the `y`; unlike ASTAP a bare `--check` also turns it on (ASTAP silently ignores a
  bare `-check`, which its own help text shows as the way to use it).

Two deliberate differences: `--extract` alone prints the real median HFD where ASTAP
prints its "none" value, 21.5; and `--sip`/`--extract2` add SIP only when significant
(§12.4).

### 12.7 Parallelise the spiral — DONE

Spiral positions are independent up to the "first match wins" rule. Evaluate a batch of
positions across a thread pool, and take the *best-verified* (not the first) result — which
is strictly better than first-wins once §12.1 gives a score to compare. With 8 threads this
turns a wide-radius hinted solve from seconds into hundreds of milliseconds.

**Done**, but deliberately keeping first-wins semantics: positions run in batches of
`--threads`, and the lowest spiral index that verifies wins, so the answer is identical to
the serial search. `dens_scutum` went from 7.6 s to 3.8 s. Choosing the best-verified
position across a batch is still open.

### 12.8 Proper motion

If we ever build our own index (§12.5.3), store `pmra`/`pmdec` (2 bytes each at ~0.5
mas/yr resolution covers ±16″/yr) and propagate to the image epoch from `DATE-OBS`:

```
    α(t) = α₀ + μ_α* · (t − 2016.0) / cos δ₀
    δ(t) = δ₀ + μ_δ  · (t − 2016.0)
```

Low priority for a 3″/px system, but it is the difference between "matches ASTAP" and
"matches Gaia".

### 12.9 Hash-based quad matching

`FUTURE_IMPROVEMENTS.md` records ASTAP's tuned constants (switch to hashing at
`nrquads ≥ 120`; `hash_bins = round(1/quad_tolerance) + 2`). Worth doing *after*
measuring: `find_matches_sorted` already reduced this to `O(M(log R + hits))`, so the
remaining win may be small compared with §12.7.

### 12.10 Testing

* Add an end-to-end test that builds a synthetic star field with a known WCS (the machinery
  already exists in `solver.rs`'s `make_test_scene`), writes it to a temporary FITS, solves
  it and asserts sub-arcsecond agreement. This is the missing regression net. Still open.
* A corpus of FITS frames with reference solutions and a pass/fail harness — **done** as
  `scripts/test-images.tsv`, `scripts/fetch-test-images.sh` and `scripts/benchmark.py`;
  see [test-images.md](test-images.md). The images are fetched rather than checked in.
* Add a false-positive test: solve a field with a deliberately wrong hint far from the
  truth and assert the solver *fails* rather than inventing a solution. Tier D of the
  corpus covers star-poor fields; the wrong-hint case is still open
  (`benchmark.py --offset-hint` moves the hint, but does not assert failure).

---

## 13. Benchmarking

The accuracy and robustness claims above only mean something against a fixed corpus. See
**[test-images.md](test-images.md)** for the benchmark set: where the FITS images with
trustworthy ground truth come from, what the corpus spans (field size, star density, image
quality), the metrics and pass criteria, and the measured results — currently 92 of 103
correct with 0 false positives, against ASTAP's 47.

---

## 14. References

**Primary algorithm papers**

- [lang2010]: Lang, Hogg, Mierle, Blanton & Roweis (2010), *Astrometry.net: Blind astrometric calibration of arbitrary astronomical images*, AJ 139, 1782. [arXiv:0910.2233](https://arxiv.org/abs/0910.2233) · [ar5iv HTML](https://ar5iv.labs.arxiv.org/html/0910.2233) · [IOP](https://iopscience.iop.org/article/10.1088/0004-6256/139/5/1782)
- [groth1986]: Groth, E. J. (1986), *A pattern-matching algorithm for two-dimensional coordinate lists*, AJ 91, 1244. [ADS](https://ui.adsabs.harvard.edu/abs/1986AJ.....91.1244G/abstract)
- [valdes1995]: Valdes, Campusano, Velasquez & Stetson (1995), *FOCAS Automatic Catalog Matching Algorithms*, PASP 107, 1119. [IOP](https://iopscience.iop.org/article/10.1086/133667) · [ADS](https://ui.adsabs.harvard.edu/abs/1995PASP..107.1119V/abstract)
- [astroalign]: Beroiz, Cabral & Sanchez (2020), *Astroalign: A Python module for astronomical image registration*. [arXiv:1909.02946](https://ar5iv.labs.arxiv.org/html/1909.02946)
- [mortari2004]: Mortari, Samaan, Bruccoleri & Junkins (2004), *The Pyramid Star Identification Technique*, NAVIGATION 51(3). [Wiley](https://onlinelibrary.wiley.com/doi/abs/10.1002/j.2161-4296.2004.tb00349.x)
- [tetra2017]: Brown, Stubis & Cahoy (2017), *TETRA: Star Identification with Hash Tables*, AIAA/USU SmallSat. [Semantic Scholar](https://www.semanticscholar.org/paper/TETRA:-Star-Identification-with-Hash-Tables-Brown-Stubis/1889aadccdcfe0e8b19a2e1c0861083131142e41)
- [bertin1996]: Bertin & Arnouts (1996), *SExtractor: Software for source extraction*, A&AS 117, 393.

**Standards and conventions**

- [wcs2]: Calabretta & Greisen (2002), *Representations of celestial coordinates in FITS* (WCS Paper II), A&A 395, 1077. [PDF](https://www.aanda.org/articles/aa/pdf/2002/45/aah3860.pdf) · [arXiv](https://arxiv.org/pdf/astro-ph/0207413)
- [sip]: Shupe, Moshir, Li, Makovoz, Narron & Hook (2005), *The SIP Convention for Representing Distortion in FITS Image Headers*, ADASS XIV. [FITS registry PDF](https://fits.gsfc.nasa.gov/registry/sip/SIP_distortion_v1_0.pdf)
- [healpix]: Górski et al. (2005), *HEALPix: A Framework for High-Resolution Discretization … on the Sphere*, ApJ 622, 759. [IOP](https://iopscience.iop.org/article/10.1086/427976)
- [hfd]: [Half flux diameter — Wikipedia](https://en.wikipedia.org/wiki/Half_flux_diameter) · [AAVSO: Using HFD instead of FWHM](https://www.aavso.org/using-half-flux-diameter-hfd-instead-fwhm)

**Implementations**

- [astap-alg]: Han Kleijn, [*ASTAP star pattern recognition algorithm and astrometric (plate) solving*](https://www.hnsky.org/astap_astrometric_solving.htm) — the algorithm family arcsec implements
- [ASTAP main page](https://www.hnsky.org/astap.htm) — star database variants (D80/D50/D20/D05/V50…), density-limited design
- [anet-readme]: [Astrometry.net README](https://astrometrynet.readthedocs.io/en/latest/readme.html) — index series and skymark table
- [build-index]: [Building index files for Astrometry.net](https://astrometrynet.readthedocs.io/en/latest/build-index.html) — presets, `dimquads`, HEALPix `Nside`
- [astrometry.net source](https://github.com/dstndstn/astrometry.net) — `solver/tweak2.c` (SIP annealing), `util/fit-wcs.c`
- [watney]: [Watney Astrometry Engine](https://github.com/Jusas/WatneyAstrometry) — .NET solver, ASTAP-family algorithm with its own Gaia quad database
- [tetra3](https://tetra3.readthedocs.io/) · [cedar-solve](https://github.com/smroid/cedar-solve) — star-tracker-lineage solvers
- [StellarSolver](https://github.com/rlancaste/stellarsolver) — SEP + astrometry.net as an embeddable library (KStars)
- [twirl](https://github.com/lgrcia/twirl) — pure-Python asterism solver following Lang et al.
- [scamp]: Bertin, [*SCAMP user's guide*](https://faun.rc.fas.harvard.edu/ameisner/scamp.pdf) — cross-matching and PV/TPV distortion
- [Siril platesolving docs](https://siril.readthedocs.io/en/latest/astrometry/platesolving.html)

**Comparisons and background**

- [astrokeith]: AstroKeith, [*Plate-solvers compared*](https://astrokeith.com/blogs/latest-blogs/plate-solvers-compared.html) — astrometry.net vs tetra3 vs cedar-solve on Raspberry Pi
- [Astrometric solving — Wikipedia](https://en.wikipedia.org/wiki/Astrometric_solving)
- [Oleg Ignat, *How astronomic plate-solving works*](https://olegignat.com/how-plate-solving-works/)

<!-- Link reference definitions for the inline citations above. -->

[lang2010]: https://arxiv.org/abs/0910.2233
[groth1986]: https://ui.adsabs.harvard.edu/abs/1986AJ.....91.1244G/abstract
[valdes1995]: https://iopscience.iop.org/article/10.1086/133667
[astroalign]: https://ar5iv.labs.arxiv.org/html/1909.02946
[mortari2004]: https://onlinelibrary.wiley.com/doi/abs/10.1002/j.2161-4296.2004.tb00349.x
[tetra2017]: https://www.semanticscholar.org/paper/TETRA:-Star-Identification-with-Hash-Tables-Brown-Stubis/1889aadccdcfe0e8b19a2e1c0861083131142e41
[bertin1996]: https://ui.adsabs.harvard.edu/abs/1996A%26AS..117..393B/abstract
[wcs2]: https://www.aanda.org/articles/aa/pdf/2002/45/aah3860.pdf
[sip]: https://fits.gsfc.nasa.gov/registry/sip/SIP_distortion_v1_0.pdf
[healpix]: https://iopscience.iop.org/article/10.1086/427976
[hfd]: https://en.wikipedia.org/wiki/Half_flux_diameter
[astap-alg]: https://www.hnsky.org/astap_astrometric_solving.htm
[anet-readme]: https://astrometrynet.readthedocs.io/en/latest/readme.html
[build-index]: https://astrometrynet.readthedocs.io/en/latest/build-index.html
[watney]: https://github.com/Jusas/WatneyAstrometry
[scamp]: https://faun.rc.fas.harvard.edu/ameisner/scamp.pdf
[astrokeith]: https://astrokeith.com/blogs/latest-blogs/plate-solvers-compared.html
