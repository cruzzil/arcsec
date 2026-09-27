# Building an Offline Quad Index

A design and staged plan for pre-computing arcsec's own quad index, including how a
user would build one without a large download.

Companion to [plate-solving.md](plate-solving.md) §5.2 (spiral vs pre-indexed) and
§12.5 (quad-selection robustness). Written 2026-09-03.

**Status (0.1.0, 2026-09-25): not started.** None of the phases below, including the
Phase 0 experiment, has been done; there is no `arcsec index` subcommand and no
`PLOVIDX` reader. The numbers in §1 still hold on the current code (90/103, 0 false
positives, the same 8 tier-A failures).

---

## 1. Read this first: the case is weaker than it was

`FUTURE_IMPROVEMENTS.md` proposed an offline index to fix the quad-correspondence
problem, and [plate-solving.md §12.5](plate-solving.md#125-make-quad-selection-robust-to-differing-star-sets)
listed it as the deep fix. Since then two things landed that took most of that
motivation away:

* **9-nearest-neighbour quad redundancy** — tier A went 26/62 → 38/62 on its own.
* **Star-level verification** — took the corpus to 90/103 with **zero** false
  positives, against ASTAP's 47.

The original argument was "online quad construction is fragile because the image and
catalogue star sets differ". That is still true, but redundancy has absorbed most of
the damage. Of the 8 tier-A images still failing, every one is a dense galactic-plane
field or a frame dominated by one bright or extended object — cases where the
*detected* stars are not catalogue stars at all (saturation, blends, nebulosity). A
pre-built index does not fix that; it changes which fixed star selection we mismatch
against.

So the honest position: **an offline index is unlikely to move the 8 remaining
failures.** Do not build it for that reason. There are three arguments left that do
stand up:

1. **Drop the astrometry.net dependency for blind solving.** Blind mode currently
   needs their index files — a separate download on top of the ASTAP database the user
   already has (`arcsec catalog install anet-4100` is ~160 MB; the Gaia-based
   `anet-5200` LITE series is ~8.8 GB). Building from the ASTAP database means **no new
   data download at all**.
2. **Wide-radius and blind speed.** The spiral costs `O((r/FOV)²)` positions and
   rebuilds catalogue quads at every one. Profiling put 36% of a failing dense solve
   in `find_matches_indexed`, called once per position. An index turns that into one
   lookup per image quad, independent of search radius.
3. **Very sparse fields.** Where only a handful of stars exist, a pre-enumerated quad
   set with deliberate redundancy beats whatever the online 9-NN happens to pick.

§9 proposes a cheap experiment that tests all three before committing to the
subsystem.

---

## 2. Why this is affordable now

`FUTURE_IMPROVEMENTS.md` estimated 10–40 GB, which is what killed the idea. That
figure assumed astrometry.net's parameters: grid cells one third of the field size and
each star in up to 8 quads across ~16 passes. Their 4200-series is ~35 GB for exactly
that reason.

We do not need their robustness envelope, because we have a hint and a verifier. Sizing
a cell to the **quad diameter** rather than a third of the field, with 12 stars per
cell and 6 quads per star at 32 bytes per quad:

| Quad diameter | Cells | Stars | Quads | Band size |
|---|---|---|---|---|
| 0.25° | 660,048 | 7.9 M | 47.5 M | 1584 MB |
| 0.35° | 330,024 | 4.0 M | 23.8 M | 792 MB |
| 0.50° | 165,012 | 2.0 M | 11.9 M | 396 MB |
| 0.71° | 82,506 | 990 k | 5.9 M | 198 MB |
| 1.00° | 41,253 | 495 k | 3.0 M | 99 MB |
| 1.41° | 20,626 | 248 k | 1.5 M | 50 MB |
| 2.00° | 10,313 | 124 k | 743 k | 25 MB |
| 2.83° | 5,157 | 62 k | 371 k | 12 MB |
| 4.00° | 2,578 | 31 k | 186 k | 6 MB |
| 5.66°–22.6° | ≤1,289 | ≤15 k | ≤93 k | ≤3 MB each |

**Full ladder, 0.25° to 24°, whole sky: 3.2 GB.** Per rig it is far less, because a
field of view `F` only needs quads roughly 0.45 F to 1.0 F across:

| Rig | Bands | Size |
|---|---|---|
| NGC 3372 set (3.13″/px, 2.6° field) | 3 | **87 MB** |
| 1.5° field | 3 | **347 MB** |
| 0.5° field | 3 | 2.8 GB |

That is the headline for usability: **a typical user builds 100–400 MB from data they
already have.** Only sub-0.5° fields get expensive, and those are the fields where a
positional hint is easiest to come by anyway.

The knobs, in order of effect on size:

* **Bands built** — quadratic in `1/d`. Building 3 bands instead of 6 is the single
  biggest saving.
* **Declination range** — most people image a limited slice of sky. `--dec-range
  -30:+60` is 0.65 of the sphere; `-10:+70` is 0.55.
* **Cell offset passes** — one pass leaves quads straddling cell boundaries
  unindexed. A second pass on a half-cell-offset grid removes the blind spot and
  doubles the size. Start with one pass and measure whether it matters.
* **Stars per cell and quads per star** — linear, and the least safe to cut.

---

## 3. Format

One file per scale band, memory-mapped, binary-searchable, no new dependencies —
mirroring how the existing readers work. Little-endian throughout.

```
  ┌─ header, 256 bytes ─────────────────────────────────────────────┐
  │ magic          [u8; 8]  "PLOVIDX\0"                             │
  │ version        u32      = 1                                     │
  │ dim_quads      u32      4 (quads) or 3 (triangles, wide fields)  │
  │ n_quads        u64                                              │
  │ n_stars        u64                                              │
  │ band_lo        f64      min quad diameter, radians               │
  │ band_hi        f64      max quad diameter, radians               │
  │ cell_side      f64      grid cell side used, radians             │
  │ dec_lo, dec_hi f64      sky coverage, radians                    │
  │ epoch_jyear    f32      catalogue epoch (ASTAP headers say 2025) │
  │ mag_limit      f32      faintest star included                   │
  │ source         [u8; 32] e.g. "ASTAP d80"                        │
  │ built_unix     i64                                              │
  │ codes_crc32    u32      integrity, cheap to check on load        │
  │ reserved       ...      zero-filled to 256                       │
  ├─ section 1: codes ──────────────────────────────────────────────┤
  │ n_quads × [f32; 4]   (CX, CY, DX, DY), sorted ascending by DY    │
  ├─ section 2: quad stars ─────────────────────────────────────────┤
  │ n_quads × [u32; 4]   indices into section 3, in A,B,C,D order    │
  ├─ section 3: stars ──────────────────────────────────────────────┤
  │ n_stars × [f32; 2]   (RA, Dec) radians                           │
  └─────────────────────────────────────────────────────────────────┘
```

Three choices worth justifying:

**Astrometry.net code space, not our 5-ratio descriptor.** `blind.rs` already
implements the canonical A–B frame construction (for quads `CX + DX ≤ 1` and
`CX ≤ DX`; for triangles `CX ≤ 0.5`) and its matcher, so the online side largely
exists. More importantly a code match yields the **star correspondence** (A↔A′, B↔B′,
…), so a single matched quad gives a full WCS — that is what makes blind solving cheap.
Our 5-ratio descriptor is reflection-invariant and gives no correspondence, which is why
the hinted path's first fit at each position is built from quad *centroids*, and star
pairs only appear afterwards in `verify_and_refit`. The 5 ratios can always be recomputed from the four star positions if
something needs them.

**Sorted by the last code dimension, not the first.** Exactly the finding from the
matcher optimisation: the tolerance window on a tightly-clustered key catches a large
slice of the file. Measure the spread of each code dimension over a real band and sort
on the widest, as `INDEX_RATIO` does for the ratio matcher.

**Star indices, not inline positions.** A star appears in ~6 quads, so indices cost
16 B where inline `(f32, f32)` pairs would cost 32 B, and the star table doubles as
the verification catalogue that `verify_and_refit` needs.

`f32` for positions gives ~0.05″ resolution at these magnitudes — an order of
magnitude finer than our 0.69″ median centre error, and the codes are matched to a
tolerance three orders of magnitude coarser.

---

## 4. Builder algorithm

```
for each band (lo, hi):
    cell = hi                                  # a maximal quad fits inside one cell
    grid = equal_area_cells(cell, dec_lo, dec_hi)

    parallel for each cell in grid:            # embarrassingly parallel
        stars = read_catalog_stars(db, centre_of(cell), cell * 1.2, spc_limit)
        keep the brightest STARS_PER_CELL
        for each star s:
            for each 4-subset of s's k nearest neighbours:
                d = max pairwise separation
                if d not in [lo, hi]: continue
                code = canonical_code(A, B, C, D)      # reuse blind.rs
                emit (code, [A,B,C,D]) , at most QUADS_PER_STAR per s
        write a per-cell shard to a temp file    # resumability

    merge shards, dedup stars, sort by the chosen code dimension, write the band file
```

Notes that matter:

* **Reuse `find_many_quads`** for the neighbour enumeration — it already does k-NN
  plus all `C(k,4)` subsets with hash-grid dedup, and it is well tested.
* **Read through `read_catalog_stars`**, so the builder gets `.1476`, `.290` and
  `.001` support for free and inherits the FOV-aware database selection.
* **Equal-area cells** can come from the same generating rule `areas_290.rs` already
  uses, evaluated at arbitrary resolution: `sin(dec_k) = -1 + 2·k/N_rings` with RA
  counts set to keep cells near-square. No HEALPix dependency.
* **Shard then merge.** The 0.5° band is 380 MB of quads — fine in RAM, but 0.25° is
  1.6 GB and sharding makes the build resumable, which matters when it takes an hour.
* **One parity in the index.** The solver already tries both image parities, so
  indexing both would double the file for nothing.
* **Proper motion.** The ASTAP databases have no proper motions, so the index inherits
  their epoch. Record it in the header (§3) so a future Gaia-sourced builder can
  propagate and the solver can tell the difference.

---

## 5. What the user does

The builder is a subcommand of the existing binary, not a second tool to install:

```bash
# The common case: size the bands from an image you already have.
arcsec index build --like ~/lights/M42_0001.fits

# Or state the field directly.
arcsec index build --fov 1.5

# Only the sky you can actually see, which is most of the saving.
arcsec index build --fov 1.5 --dec-range -30:+60

# Explicit control.
arcsec index build --bands 0.7,1.0,1.5 --db ~/star_database -D d80 --out ~/arcsec-index

arcsec index list ~/arcsec-index      # bands, coverage, size, source, epoch
arcsec index verify ~/arcsec-index    # CRC + solve a synthetic field per band
```

Behaviour that makes it "easy" rather than merely possible:

* **No new data download.** It builds from the ASTAP database the user already has for
  normal solving. This is the whole point.
* **Defaults that need no thought.** `--out` defaults to the directory
  `arcsec catalog` already manages (`catalog_cmd::default_dir`, e.g.
  `~/.local/share/arcsec/catalogs` on Linux, overridable with `ARCSEC_CATALOG_DIR`),
  threads to `max_threads()`, bands to 0.45–1.0 × the field.
* **Says what it will cost before doing it.** `Will build 3 bands (0.71°, 1.00°,
  1.41°), 347 MB, ~12 min on 24 threads. Continue? [y/N]` — and `--yes` for scripts.
* **Resumable.** Interrupt it and re-run; completed cell shards are reused.
* **Progress that means something** — cells done, quads emitted, ETA.
* **The solver finds it without being told.** `-i` already takes a directory and ranks
  astrometry.net files by scale; extend that to recognise `PLOVIDX` by magic and rank
  both kinds together. Then check the catalogue directory when `-i` is absent, so a
  built index is simply used. (Today `-i` is always required for blind mode, even though
  `arcsec catalog install anet-4100` puts the astrometry.net files in that directory.)

For users who would rather not build at all, §8 phase 4 covers publishing pre-built
bands.

---

## 6. Integration with the solver

Three levels, each independently shippable:

1. **Blind, instead of astrometry.net.** `blind.rs`'s `hyp_from_entry` and
   `verify_score`, and `AnetIndex::find_code_matches_into`, work on `AnetIndex`. Introduce a small trait —
   codes, quad stars, star list, scale range — implement it for both `AnetIndex` and
   the new format, and blind mode reads either. Lowest risk: it touches no path that
   the current 90/103 depends on.
2. **Hinted fast path for wide radii.** When the search radius is large enough that
   the spiral would visit many positions, look the image quads up in the index
   directly, vote positions, and hand the best to `verify_and_refit`. Keep the spiral
   as the fallback. Gate on measured wall time, not on a guess.
3. **Replace online catalogue quads entirely.** Only if 1 and 2 show the index
   matching or beating the spiral on the corpus. This is the change that could regress
   90/103, so it goes last and only on evidence.

---

## 7. Sizing and validation gates

Every phase is measured against the existing 103-image corpus with
`scripts/benchmark.py`, and must not regress **90 correct / 0 false positives**.

| Gate | Requirement |
|---|---|
| Reader | Round-trips a synthetic index; rejects truncated and bad-magic files |
| Builder | For a known sky cell, the index quads contain the quads the online path builds for the same stars |
| Blind parity | On the HiPS blind test fields, our index matches or beats the astrometry.net 4100 series on solve rate |
| No regression | Corpus stays at ≥ 90 correct, 0 false positives |
| Speed | Wide-radius solve (`-r 30`) faster than the spiral, measured, not assumed |
| Build cost | 3 bands for a 1.5° field in under 30 min on 24 threads, under 400 MB |

---

## 8. Phases

**Phase 0 — test the premise (≈1 day).** Described in §9. Do this first.

**Phase 1 — format and reader (≈2 days).** Header, three sections, magic and CRC
checks, memory-mapped loader, code-window binary search. Unit tests on a synthetic
index built in memory, as `format_001.rs`'s tests do. No builder yet.

**Phase 2 — builder, single band (≈3 days).** `arcsec index build --bands X`,
serial, one band, writing shards then merging. Validate against the builder gate
above. Correctness before speed.

**Phase 3 — builder, usable (≈3 days).** Parallel over cells, resumable, `--like`,
`--fov`, `--dec-range`, the cost prompt, `index list`, `index verify`, progress.

**Phase 4 — blind integration (≈2 days).** The trait in §6.1 so blind mode reads
either format; run the blind test scripts against both and compare.

**Phase 5 — hinted fast path (≈3 days, optional).** §6.2, gated on measurement.

**Phase 6 — distribution (optional).** Publish pre-built bands as tarballs with
checksums, and a `arcsec index fetch --fov 1.5` that downloads instead of building.
Only worth it once the format has stopped changing.

Roughly two weeks of focused work to the end of phase 4, which is the point where the
astrometry.net dependency goes away.

---

## 9. Phase 0: the experiment that decides whether to build any of this

The subsystem above is a fortnight of work and a new on-disk format to maintain
forever. Before committing, spend a day testing whether an index would actually help,
using a throwaway script rather than production code:

1. Pick the 8 tier-A images that still fail, plus 4 that succeed as controls.
2. For each, build a **tiny** index in memory covering only that field's sky and one
   band at the field size — a few thousand quads, no file format, no CLI.
3. Run the existing `blind.rs` matcher against it and record: does it find the field,
   how many code matches, what does `verify_score` say?
4. Separately, time a `-r 30` hinted solve against the same tiny index versus the
   spiral, to test argument 2 from §1.

Interpretation, decided in advance so the result cannot be rationalised:

* **Several of the 8 solve** → the correspondence argument is alive after all; build
  the full thing.
* **None solve but the controls do** → confirms §1: the remaining failures are
  detection problems, not indexing problems. Then build the index only if the
  astrometry.net-dependency and wide-radius-speed arguments are worth two weeks on
  their own — and fix detection first.
* **Even the controls fail** → the code-space builder disagrees with `blind.rs`'s
  canonical form somewhere. Fix that before drawing any conclusion; it is the most
  likely outcome of a first attempt and the cheapest possible place to discover it.

---

## 10. Alternatives considered

**Write astrometry.net's own format instead of a new one.** Then `anet.rs` reads it
unchanged and their tooling works on our files. Rejected for the builder: their format
carries kd-trees and range tables that are considerably more work to *write* correctly
than a sorted flat array is to write and read. Worth revisiting if we ever want
interoperability rather than independence.

**Index triangles rather than quads for wide fields.** `DIMQUADS=3` is what
astrometry.net recommends above ~10°, and `blind.rs` already handles it. The header has
`dim_quads` for this; treat it as a phase-5 refinement.

**Re-encode Gaia DR3 directly**, per `FUTURE_IMPROVEMENTS.md`. That buys proper
motions and independence from ASTAP's databases, but it is a much larger project
(terabytes of source data) and it is orthogonal: the index format above does not care
where the stars came from, and `source`/`epoch_jyear` in the header record it. Do the
index first, from ASTAP data, and swap the source later if proper motions ever matter.

**HEALPix for the cell grid.** Genuinely better — exactly equal areas and a
hierarchical index — but it is a new dependency and the equal-area rule from
`areas_290.rs` is good enough at the resolutions involved. Reconsider if cell shape
turns out to bias quad density.
