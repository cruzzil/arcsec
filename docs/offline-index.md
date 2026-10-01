# arcsec's Blind Index

A pre-computed pattern index, built from the ASTAP star database the user already has,
that finds an image's field anywhere on the sky without a position hint and — with no
extra cost — without a pixel scale either.

Companion to [plate-solving.md](plate-solving.md) §5.2 (spiral vs pre-indexed), §10.4
(the Astrometry.net blind path) and §12.3/§12.5. First written 2026-09-03 as a plan;
rewritten 2026-10-02 when the index was built. The plan's original sections (sizing by
cell, an astrometry.net-style code index) are summarised in §10, with why the built
design differs.

**Status (2026-10-02): built, behind `--index` and an opt-in automatic fallback.**
`arcsec catalog index build` writes `<db>.arcsecix` into the catalogue directory;
`arcsec -i <file|dir>` solves blind with it; and when an index is installed, a search
wider than five fields round the hint consults it after the first five fields of the
spiral. On the 596-image corpus (tiers A/B/C/S) a blind solve finds 436 images
(87 % of what the hinted solver finds with the true centre as its hint), in a
median 0.4–3 s, with **no false position**: every reported solution passes the
hinted solver's own star-level verification, and the only wrong answers are the four
distorted wide fields the hinted solver also gets wrong. Astrometry.net's 4107–4119 on
the same images: 70 of the 254 fields ≥ 0.6° against the index's 188,
at 25 s against 3 s.

---

## 1. Why build it — the case, revisited

The plan argued (and §9's experiment confirmed) that an index would not rescue the
hinted solver's remaining failures: those are detection problems — saturated, blended or
nebulous fields — and a pre-built index mismatches them just as the online quads do.
That stands. What made the index worth building is the other three arguments, and the
evidence for each came out stronger than the plan expected:

1. **No Astrometry.net dependency for blind solving.** `-i` used to need index files on
   top of the ASTAP database (`anet-4100` is 355 MB and covers only fields ≥ ~0.7°). The
   arcsec index is built from D80 (or whichever database is installed) in 2–13 minutes,
   with nothing downloaded.
2. **Wide-radius and hint-free speed.** The spiral is `O((r/FOV)²)`: N.I.N.A.'s blind
   mode (`-r 180`, no position) took 76–418 s per image. The index answers in a second or
   two whatever the radius (§7.3).
3. **The pixel scale becomes optional.** Every index match implies a scale, so a sweep
   over 0.3–60″/px costs about twice a scale-hinted solve (§7.4). plate-solving.md
   §11.3 listed "cannot solve without a good pixel-scale estimate" as a shortcoming.

---

## 2. Design

### 2.1 Disc-anchored patterns

The idea, which comes from seiza's blind index (Apache-2.0; ideas only, the code here is
an independent implementation), is to index exactly the patterns an image can rebuild
without knowing where it is:

* For each **tier** — a disc radius `r` and a magnitude cap — a star **anchors** patterns
  only if it is the brightest star (no fainter than the cap) within `r` of itself.
* Its **group** is itself plus the next-brightest stars of its disc: 6 stars in the
  wider tiers, 5 in the two deepest.
* Every **4-subset** of the group is a pattern: 15 per anchor (6-star groups) or 5.

An image that contains a disc sees the same locally-brightest stars, so it can rebuild
the group from its own brightest detections. The 4-subsets give redundancy: a group
member that is undetected or ranks differently in the image costs some patterns, not
all of them.

No two anchors share a pattern (a pattern contains its anchor, and one anchor cannot lie
in another's disc — one of them would be the brighter), so nothing needs de-duplicating.

### 2.2 Descriptor, key and correspondence

* **Descriptor**: the quad's six pairwise distances, sorted, divided by the largest —
  the five ratios ASTAP's own quads use. Invariant to translation, rotation, scale and
  reflection, so one entry serves both image parities.
* **Key**: each ratio quantised into 128 bins (1/128 ≈ 0.0078, close to the hinted
  solver's 0.007 tolerance), packed eight bits per dimension into a `u64`. A lookup is a
  binary search for an exact key. A measured ratio within 0.0025 of a bin edge also
  probes the neighbouring bin, which costs one to four lookups per quad instead of 3⁵.
* **Correspondence**: the descriptor does not say which star is which, so both sides
  order a quad's vertices canonically — ascending total distance to the other three.
  Where two totals are within 1.5 %, the image side also tries the swapped order. A
  wrong correspondence fails the affine shape check, never produces a match.

The plan (§10) chose astrometry.net's code space instead, which carries the
correspondence in the code. The quantised 5-ratio hash won on three counts: an exact-key
binary search over a flat sorted array needs no kd-tree or window scan; one entry serves
both parities; and the ratios were already arcsec's descriptor.

### 2.3 Tiers

| Disc radius | Mag cap | Group | Fields served (short side) | Anchors | Patterns | Size |
|---|---|---|---|---|---|---|
| 12° | 4.6 | 6 | 30°–144° | — | — | (not built by default) |
| 6° | 6.1 | 6 | 15°–72° | — | — | (not built by default) |
| 3° | 7.6 | 6 | 7.5°–36° | 1 436 | 21 404 | 0.5 MB |
| 1.5° | 9.2 | 6 | 3.75°–18° | 5 824 | 87 162 | 2.1 MB |
| 0.75° | 10.7 | 6 | 1.9°–9° | 23 387 | 349 669 | 8.4 MB |
| 0.4° | 11.8 | 6 | 1.0°–4.8° | 80 980 | 1 185 086 | 28 MB |
| 0.2° | 12.7 | 6 | 0.5°–2.4° | 285 240 | 3 639 264 | 87 MB |
| 0.1° | 14.2 | 5 | 0.25°–1.2° | 986 604 | 4 428 404 | 106 MB |
| 0.06° | 16.0 | 5 | 0.15°–0.72° | 2 824 154 | 12 900 826 | 310 MB |

Each cap is where a disc of that radius holds 15–20 stars on average (D80, measured:
mag ≤ 6.1 → 0.1 stars/deg², ≤ 9.2 → 7, ≤ 12.7 → 150, ≤ 14.2 → 600, ≤ 16 → 1900). Radii
step by about 2×, so a field normally sees two tiers. "Size" is the tier's keys and quads
(24 bytes a pattern); the star table comes on top (12 bytes a star).

`arcsec catalog index build --min-fov F --max-fov G` builds the tiers that cover fields
F–G: every tier overlapping the span, but of the tiers reaching below F only the widest,
and of those reaching above G only the narrowest. The defaults, 0.3°–30°, give the six
tiers 3°–0.1°: **287 MB**, built in 2 minutes. `--min-fov 0.15` adds the 0.06° tier for
D80's narrowest fields: **698 MB**, 13 minutes (both on 24 threads, under a load average
of 50–80 from other jobs). For comparison D80 itself is 1.3 GB, G05 102 MB.

### 2.4 File format: `ARCSECIX` version 1

One file, memory-mapped, little-endian, every section 8-byte aligned
(`arcsec-core/src/index/format.rs` has the byte layout):

```
header (256 bytes)  magic "ARCSECIX", version, byte-order marker, header length,
                    descriptor bins, tier/star/pattern counts, build time, source
                    database, a table of 5 sections {offset, length, CRC-32},
                    and a CRC-32 of the header itself
tiers               40 bytes each: radius, mag cap, group size, pattern range, anchors
stars               12 bytes each: RA, Dec (f32 radians), mag ×100 (i16), widest tier
star directory      first star of each of 720 quarter-degree declination bands
keys                u64 per pattern, sorted within each tier
quads               4 × u32 star indices per pattern, canonical vertex order
```

* **Instant open, lazy validation.** Opening checks magic, version, byte order, the
  header CRC, and that every section lies inside the file with the size its counts imply
  — a few hundred bytes of work, so the 698 MB index opens in 9 ms. The section CRCs are
  checked only by `arcsec catalog verify` (4 s for 698 MB). Every star reference a lookup
  follows is bounds-checked, so a corrupt body can fail a solve but never read out of
  bounds.
* **Versioned.** Any change to the tier table semantics, the bin count, the descriptor or
  the canonical order changes which patterns exist or how they hash, and must bump
  `VERSION`; an old file is then refused with "rebuild it" rather than silently matching
  nothing.
* **Written atomically**: to `<file>.arcsecix.part`, then renamed.
* **Sorted by star band and RA**, the star table doubles as the verification catalogue:
  `stars_near` finds a field's stars with two binary searches per band.

### 2.5 Builder

`arcsec-core/src/index/build.rs`. For each tier, the sky is processed in 5° declination
strips (one strip for the shallow tiers), reading the strip plus an `r` margin from the
database with `for_each_star_in_dec_band` — each area file only down to the tier's cap,
since files are sorted brightest first. Stars are sorted by (magnitude, RA, Dec), a total
order, so "brighter" means the same thing in every strip. A grid of `r`-sized cells
answers the disc queries; anchors are found in parallel; the result is identical for any
thread count (tested). Stars reached from two strips or two tiers are merged by exact
position at the end, keeping the widest tier.

Build cost on D80, 24 threads, machine load 50–80:

| Index | Tiers | Patterns | Stars | Size | Time | Peak RSS |
|---|---|---|---|---|---|---|
| fields 1°–30° | 3°…0.4° | 1.64 M | 0.43 M | 45 MB | 19 s | 179 MB |
| fields 0.3°–30° (default) | 3°…0.1° | 9.71 M | 4.52 M | 287 MB | 122 s | 847 MB |
| fields 0.15°–30° | 3°…0.06° | 22.6 M | 12.9 M | 698 MB | 782 s | 1.55 GB |

The deepest tier dominates everything: it reads 52 M stars, and its strip reading is
serial. Reading strips in parallel would cut the 13 minutes several-fold if it matters.

### 2.6 Solver

`arcsec-core/src/pipeline/index_solve.rs`:

1. **Detect** every star (not the `-s` brightest by SNR) and **re-rank by aperture
   flux**. This was the decisive fix of the first version: SNR is not a brightness order
   at the bright end — a saturated star's flat top and wide aperture give it a lower SNR
   than a fainter, sharper star — and the index is built from each region's *brightest*
   stars. With SNR order, no index group was ever rebuilt on DSS fields; with flux order,
   the true field ranked first.
2. **Image patterns**: every 4-subset of the 20 brightest stars, at most two per cell of
   a 6×6 grid; plus, round each of the 150 brightest stars, the five brightest within
   seven window radii (short side ÷ 16 … ÷ 2) — the image's version of the disc groups.
   About 7 000 patterns.
3. **Lookup** in the tiers whose field range suits the image. A candidate's implied
   pixel scale (longest edge on the sky ÷ in the image) must lie in the allowed range; a
   4-point affine fit must be within 6 % of a similarity transform, and its scale within
   6 % of the edge ratio. Each survivor is a field hypothesis (centre, scale, affine).
4. **Vote** in (RA, Dec, ln scale) buckets — 5 % of the field on the sky, RA width ÷
   cos δ; 5 % in scale — with the neighbour smoothing, strongest-bucket representative
   and non-maximum suppression of `pipeline/sky_votes.rs` (shared with the Astrometry.net
   path, §9.1). Up to 3 000 regions.
5. **Rank** each region's medoid hypothesis by projecting the index's stars through it
   onto the detections, refitting once on the matches. The score is a **significance**,
   `(matches − E)/√(E + 1)` with `E` the chance matches for that many projected stars,
   plus 2 per agreeing vote beyond the first. Raw match counts let a hypothesis at a far
   too coarse scale, which projects hundreds of stars into the frame, outrank the truth;
   that was the one fix the scale-free mode needed.
6. **Accept** through the hinted solver: the best hypotheses (score ≥ 8, at most six) go
   to `solve_image` with the hypothesis as hint, its scale as field size and no search
   radius. Its star-level verification (≥ 30 matched stars, spread over the frame) is the
   only acceptance test, so a blind solve is held to exactly the hinted solver's
   standard. The index never decides on its own that a field is found — seiza's
   blind-only acceptance (≥ 12 matches, RMS < 2 px) produced wrong-field solves there.

### 2.7 Command line

```bash
arcsec catalog index build                    # deepest installed database, fields 0.3°–30°
arcsec catalog index build --min-fov 0.15     # down to D80's floor (698 MB)
arcsec catalog index build --db ~/star_database -D d80 -o idx.arcsecix
arcsec catalog index info                     # tiers, sizes, source
arcsec catalog list                           # lists built indexes too
arcsec catalog verify                         # checks every section CRC too

arcsec -f image.fits -i ~/.local/share/arcsec/catalogs        # blind: any position
arcsec -f image.fits -i idx.arcsecix --fov 1.2                # blind, scale known
```

* **`-i` names an arcsec index** (a file, or a directory holding one; recognised by
  magic, so Astrometry.net files in the same directory are ignored): the index is tried
  first, blind, whatever the hint and radius; if nothing verifies, the ordinary search
  runs from the hint, as with Astrometry.net files.
* **Automatic, no `-i`**: when an index is installed in the catalogue directory (or
  beside the star database) and `-r` reaches past five fields round the hint, the
  spiral first searches those five fields (minimum 1°); only if that fails is the index
  consulted, restricted to `-r` (plus a field) round the hint unless `-r` covers the
  sky; if that fails too, the full spiral runs as before. Inside five fields the result
  is therefore exactly the spiral's; see §7.5 for the measured effect.
* **N.I.N.A.** in blind mode passes `-r 180` and no `-ra`/`-spd`: with an index installed
  that is the automatic path, and the five-field first stage is skipped when there is no
  hint at all (no `-ra`/`-spd` and no RA/Dec in the header).
* **Pixel scale**: from `--fov` or FOCALLEN/XPIXSZ, ±20 %; with neither, 0.3–60″/px.

---

## 3–6. (Superseded plan sections)

The plan's format (`PLOVIDX`, astrometry.net codes sorted by one dimension), its
cell-grid builder and its staged integration through a shared trait with `AnetIndex`
were not built as written; §2 is what was built and §10 records why. The
`--like`/`--dec-range`/resumable-shard/cost-prompt features of the plan's §5 are
deferred (§8).

---

## 7. Results

All on the expanded corpus (docs/test-images.md §7), `scripts/benchmark.py --corpus
--auto-db`, star databases D80, G05 and W08 in `~/star_database`, 24 cores shared with
other benchmark jobs (load averages 50–150 throughout — **times are upper bounds**, with
8 concurrent jobs). Blind runs use `--blind`, which puts the hint at the antipode of the
truth with `-r 0`, so only the index can find the field.

### 7.1 Blind solve rate by field

Correct solves (median time of the correct ones). "Hinted" is the ordinary solver with
the *true* centre as hint and `-r 5`, the best case; the blind columns get no position.
The 0.15° index has all seven tiers (698 MB), the 0.3° index the default six (287 MB).

| FOV (long side) | n | hinted, true centre | blind, 0.15° index | blind, 0.3° index | blind, no scale (0.15° index) |
|---|---|---|---|---|---|
| < 0.15° | 6 | 4 | 0 | 0 | 0 |
| 0.15–0.3° | 113 | 103 | 51 (0.5 s) | 2 | 43 |
| 0.3–0.6° | 125 | 116 | 109 (0.6 s) | 94 | 102 |
| 0.6–1.2° | 177 | 161 | 160 (1.2 s) | 156 | 157 |
| 1.2–2.5° | 88 | 75 | 73 (2.8 s) | 73 | 74 |
| 2.5–6° | 43 | 17 | 16 (0.5 s) | 16 | 13 |
| 6–20° | 36 | 23 | 23 (1.3 s) | 23 | 21 |
| > 20° | 8 | 5 | 4 (2.5 s) | 4 | 0 |
| **all (A+B+C+S)** | **596** | **504** | **436** (1.0 s) | **368** | **410** |
| tier A / B / C / S | | 206 / 206 / 38 / 54 | 192 / 154 / 36 / 54 | 177 / 107 / 31 / 53 | 183 / 140 / 34 / 53 |

From 0.3° up the blind solve finds 385 of the hinted solver's 397. Every accepted
hypothesis was the top-ranked one in 452 of 455 index solves (rank ≤ 2 in all), and none
needed a second hinted solve. The index stage itself takes a median 0.9 s; the slowest
solves (30–60 s) are images where *detection* takes 20–40 s, which the hinted
verification then repeats.

### 7.2 Against Astrometry.net 4107–4119

The same 254 fields of ≥ 0.6° (tiers A/B/S), scale given, no position. The Astrometry.net
runs were made under load averages of 60–200, so their times are inflated, but even
halved they are an order of magnitude slower: each process loads and sorts its index
files (~35 s each, two in parallel) before matching.

| FOV | n | anet 4107–4119 before §9.1 | anet after §9.1 | arcsec index |
|---|---|---|---|---|
| 0.6–1.2° | 79 | 15 (71 s) | 16 | 72 (3.3 s) |
| 1.2–2.5° | 88 | 24 (18 s) | 25 | 73 (2.8 s) |
| 2.5–6° | 43 | 12 (52 s) | 10* | 16 (0.5 s) |
| 6–20° | 36 | 17 (15 s) | 16 | 23 (1.3 s) |
| > 20° | 8 | 2 | 3 | 4 |
| **all** | **254** | **70** (26 s) | **70*** | **188** (2.7 s) |

\* Two of those were 300 s timeouts under a load average near 200; re-run with a longer
timeout they solve, giving 72 (§9.1). On the original 103-image set the 4100 series finds
12 either way: most of it is narrower than 0.7°, below that series' range.

### 7.3 N.I.N.A.-style: `-r 180`, no useful hint

N.I.N.A. hands ASTAP (and so arcsec) `-r 180` and no position when it blind-solves. Eight
corpus images, hint at the antipode, `-r 180`, two at a time:

| | correct | time |
|---|---|---|
| spiral only (0.2.0) | 1 of 8, **plus one false position** (`rnd_022`, 155° off), 5 timeouts at 900 s | 269 s for the one |
| index installed, automatic (no `-i`) | 7 of 8 | 1.0–4.1 s each |

(The spiral run overlapped a load average of 100–200; the Siril work measured 76–418 s
per image for the same mode on a quiet machine.) The eighth, `wide_dss_01`, fails both
ways. With the antipodal hint the automatic path first spends ~1 s on the five-field
spiral; with no hint at all (N.I.N.A.'s real case when the header has no RA/Dec) that
stage is skipped. The spiral's false position is worth noting on its own: a whole-sky
spiral verifies tens of thousands of positions, and one passed.

### 7.4 No pixel scale at all

`--blind --no-fov`: no `--fov`, and the corpus headers carry no FOCALLEN/XPIXSZ, so the
index searches 0.3–60″/px. 410 correct against 436 with the scale given (table in §7.1),
median 0.9 s: it misses 35 that the scale-given run finds (mostly narrow survey frames,
where the wider search admits more competing hypotheses) and finds 9 that it misses. Fields above
20° fail because their scale is above 60″/px. Before the chance-corrected score (§2.6
step 5) this mode found almost nothing: a hypothesis at 20–50″/px projected hundreds of
stars into the frame and outscored the truth on raw counts.

### 7.5 Hinted solves with an index installed

The automatic path changes nothing a five-field spiral would solve. Corpus, true-centre
hint, `-r 5`, index installed: **504 → 504, identical statuses on all 635 entries**. With
the hint 0.3 fields off: **459 → 462** (`ls_big`, `ps1v2_07`, `sv_sdssr` gained, none
lost, same false positives).

The cost is on images that solve nowhere: the five-field spiral, the index, then the full
spiral, which repeats the first stage. Eight such images, `--jobs 1`, quiet machine (load
4), two rounds: 1.1–6.2 s → 1.9–10.5 s, about 1.7× (solved images unchanged, 0.13–0.15 s).
Before the index checked its leading hypotheses against the star database (§2.6 step 6)
a wrong hypothesis cost a full hinted solve, and some no-solves grew by 15–30 s.

### 7.6 False positives

* Tiers A/B/C/S: the same four as the hinted solver — `wide_shassa_01` (31″ at the
  corners, just over the 1-pixel threshold) and `tess_03`, `tess_25`, `tess_32` (12° TESS
  FFIs with 1000″ of distortion, where any linear plate is wrong at the edges). They are
  the hinted solver's acceptance, reached from a different start.
* Tier D: `neg_hint_1`, `_5`, `_7`, `_8` are reported, correctly placed (0.04–0.24″).
  Those controls are real images given a *wrong hint and a small radius*; they must fail
  only because the field lies outside `-r`. A blind solve ignores the hint by design, so
  in a blind run they are not false positives. The automatic path keeps to `-r`, and in
  the hinted runs with an index installed they are refused as before. Every
  must-fail image with no real field behind it (noise, flats, fake star fields, nebula
  cores) is refused in every mode.

### 7.7 Where blind solving still fails

* **Deep survey frames, 0.15–0.3°** (SDSS 8/33, Pan-STARRS 9/20, SkyMapper 5/15, Legacy
  Survey 26/40, against 33, 20, 8 and 37 hinted): the
  0.06° tier's anchors are mag 12–15 stars, which saturate or are masked in these
  surveys, so the groups cannot be rebuilt. A tier that skips stars brighter than a
  saturation limit, or 6-star groups in the deep tiers (×3 the size), would address it.
* **Below 0.15°**: no tier.
* **The hinted solver's own failures** (dense, nebulous, bright-object fields; distorted
  TESS frames): the index finds some of these fields, but acceptance is the hinted
  solver's, so they fail as before.
* The 0.3° default index loses most of the 0.15–0.3° band (2 of 113 against 51); that is
  the 410 MB the 0.06° tier costs.

---

## 8. What was built, what was deferred

Built: the format and reader (lazy validation, CRCs, versioning), the builder (strips,
parallel, deterministic), `catalog index build|info`, `catalog list|verify` awareness,
the solver, `-i` support, the automatic fallback, the scale-free mode, the Astrometry.net
path's vote fix (§9.1), benchmark flags `--blind` and `--no-fov`, and unit tests on
synthetic databases (format round trip and corruption, builder invariants and thread
independence, an end-to-end blind solve of a rendered field in both parities).

Deferred:

* **Saturated anchors in deep survey frames** (§7.7) — the main remaining loss.
* **Fields below 0.15°**: no tier; a 0.035° tier would need a catalogue deeper than D80's
  density limit allows to be useful.
* **Database choice in the scale-free mode**: the hinted verification uses the database
  chosen for the assumed 1″/px; a hypothesis many times wider than D80's range is then
  verified against the wrong catalogue.
* **Builder conveniences** from the plan: `--like <image>`, `--dec-range`, resumable
  shards, a size/time prompt before building, parallel strip reading.
* **Replacing the spiral's online quads** with the index (plan §6.3) — not attempted;
  the hinted path is untouched.
* **Distribution** of pre-built files (plan phase 6): unnecessary while a build takes
  minutes.

---

## 9. Phase 0 and the go decision

The plan's §9 experiment was adapted: seiza's results already showed disc-anchored
indexing working end to end, so instead of a throwaway in-memory index the experiment
was run on the first version of the real one, on whole-sky indexes, with a debug mode
(`ARCSEC_INDEX_DEBUG=<true WCS>`) that counts, for the true field, how many index
patterns lie in the frame, how many of their stars are detected, how many hash to the
same key, and how many the image side generated.

* First result: **0 of 5 fields** found. The debug counts located the failure exactly:
  index stars detected but at SNR ranks 22–3 000 (DSS fields; bright stars saturated),
  so no group was ever rebuilt. Re-ranking detections by aperture flux: **4 of 5**, all
  at hypothesis rank 0, in 0.2–1 s; the fifth (M101, rank 6) ranked first once votes
  joined the score.
* The plan's interpretation rule — "several of the failing images solve → build it" —
  does not apply as written: the question was never the hinted failures but blind
  solving, and on that the answer was decisive. **Go.**
* The plan's tier-A failure set: the blind index does not solve the dense or nebulous
  fields the hinted solver fails either, as §1 predicted.

### 9.1 Step 0: the Astrometry.net path's vote ranking

`blind.rs` verified only `hyps[0]` of each 0.1° RA/Dec vote cell, binned raw RA (which
splits one field's votes near the poles) and had no scale axis. It now uses
`sky_votes.rs`: (RA, Dec, ln scale) buckets with RA ÷ cos δ, 3×3×3 smoothing, the
medoid of the strongest bucket as representative, and non-maximum suppression.

Measured on the 254 corpus fields of ≥ 0.6° with index files 4107–4119 (`--blind -r 0
--fov`): 70 → 72 correct, the same one false positive (`wide_shassa_01`, a hinted-solver
acceptance). Gained `rnd_052`, `rnd_065`, `wide_shassa_06`; lost `wide_tess_10`.
(`tess_13` and `tess_24` timed out at 300 s in the after run under a load average near
200 and solve with a longer timeout.) Timing could not be compared: the two runs saw
load averages of 60–90 and 160–200. On the 103-image v1 set: 12 → 12. So the bug was
real but cost little here; the 4100 series' weakness is coverage, not ranking. The
calibrated `MIN_VERIFY_SCORE`/`EARLY_STOP_SCORE` were not changed, and no new false
positive appeared; the HiPS scripts (`hips_extended_test.sh`) were not re-run, since
they download from a network service — the corpus run stands in for them.

---

## 10. Alternatives considered

**Astrometry.net code space with a sorted-dimension window scan** (the plan's design).
Not built: see §2.2. Its one advantage, correspondence from the code, is recovered by
the canonical vertex order.

**Cell-grid quads sized to the quad diameter** (the plan's builder and size estimates).
The disc-anchored rule replaces it: anchoring on local brightness maxima gives the
image a way to rebuild the same groups without knowing the grid, where a cell grid
needs passes on offset grids and many quads per star for the same robustness.

**Writing astrometry.net's own format.** Rejected for the same reason as before: their
kd-trees are much more work to write than a sorted flat array.

**HEALPix**: not needed; the builder's strip-and-grid scheme and the star directory's
declination bands are enough.

**Gaia directly** (proper motions): orthogonal; the header records the source database.
