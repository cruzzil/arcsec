# Installing Star Catalogues

arcsec cannot solve anything without a star catalogue, and cannot do photometric
colour calibration without a *photometric* one. Both are one command away:

```bash
arcsec catalog recommend --like ~/lights/M42_0001.fits --photometry
arcsec catalog install d50 v05
```

Everything lands in one place that the solver already reads, so after installing you
can just run `arcsec -f image.fits` with no `-d` and no `-D`.

---

## 1. Which catalogue do I need?

Ask arcsec, and give it either a field size or an image to read one from:

```bash
arcsec catalog recommend --fov 1.5                 # field size in degrees
arcsec catalog recommend --like image.fits         # read FOCALLEN/XPIXSZ from the header
arcsec catalog recommend --fov 1.5 --photometry    # include colour calibration
```

It prints the **smallest** download that covers your field, because there is no point
pulling 1.2 GB for a rig that a 102 MB catalogue serves perfectly well.

The full set:

| Name | Purpose | Download | Fields | Notes |
|---|---|---|---|---|
| `d05` | solving | 102 MB | 0.6°–6° | Gaia DR3, 500 stars/deg². Smallest that works. |
| `d20` | solving | 400 MB | 0.3°–6° | 2000 stars/deg². |
| `d50` | solving | 901 MB | 0.2°–6° | 5000 stars/deg². The usual choice. |
| `d80` | solving | 1.2 GB | 0.15°–6° | 8000 stars/deg². Needed below ~0.2°. |
| `g05` | solving | 102 MB | 3°–20° | Wide fields. The D-series stops at 6°. |
| `w08` | solving | 330 kB | 20°–80° | Very wide fields, to magnitude 8. |
| `v05` | photometry | 117 MB | 0.6°–6° | Johnson-V + Gaia BP-RP colour, 500 stars/deg². |
| `v50` | photometry | 1.0 GB | 0.2°–6° | Johnson-V + BP-RP, 5000 stars/deg². |
| `anet-4100` | blind | 355 MB | 0.7°+ | Astrometry.net Tycho-2 indexes 4107–4119. Solve with no hint. |
| `anet-5200` | blind | 8.8 GB | 0.1°–2° | Astrometry.net Gaia LITE indexes 5200–5202. Large. |

Rules of thumb:

* **Solving**: one D-series database is enough for a telescope. Add `g05` if you also
  shoot with a camera lens, and `w08` if you shoot all-sky. They are small.
* **Colour calibration**: `v05` unless you need the depth of `v50`.
* **Blind solving** (`-i`): only if you cannot supply an approximate position. With a
  mount that reports where it is pointing, you do not need these at all. Build
  arcsec's own index from the star database you already have (§5) rather than
  downloading these: it is faster and finds more fields. The Astrometry.net indexes
  only find the rough position; the solve is then finished against a star database,
  so you need one of those as well.

---

## 2. Installing

```bash
arcsec catalog install d50               # asks for confirmation, shows the size first
arcsec catalog install d50 v05 g05 --yes # several at once, no prompt
arcsec catalog list                      # everything, with what is installed
arcsec catalog verify                    # check for missing or truncated files
arcsec catalog remove d80                # free the space again
arcsec catalog path                      # where they live
```

`install` and `remove` both list what they will do and ask first; `--yes` (or `-y`)
skips the prompt. `--keep-archive` keeps the downloaded `.zip` or `.deb` after
unpacking.

Installs are **resumable**: interrupt one and re-run the same command, and the
download continues from where it stopped rather than starting again, provided the
server supports range requests. Files are written
to a temporary name and only moved into place once complete, so a half-finished
download can never look like an installed catalogue.

---

## 3. Where they go

| Platform | Default location |
|---|---|
| Linux | `$XDG_DATA_HOME/arcsec/catalogs`, else `~/.local/share/arcsec/catalogs` |
| macOS | `~/Library/Application Support/arcsec/catalogs` |
| Windows | `%LOCALAPPDATA%\arcsec\catalogs` |

If none of those variables is set, `~/.arcsec/catalogs` is used. `arcsec catalog path`
prints the directory in effect.

Override for a single command with `--dir`, or permanently with the
`ARCSEC_CATALOG_DIR` environment variable — useful when the catalogues live on a
different disk, or are shared with an ASTAP install:

```bash
export ARCSEC_CATALOG_DIR=/mnt/data/star_databases
```

The solver reads this directory automatically, as long as a star database is
installed there; otherwise it falls back to the current working directory, as ASTAP
does. `-d` is only needed to point somewhere else for one run, and `-D` only to force a
particular database — otherwise arcsec picks the densest installed one whose field
range covers the image.

Astrometry.net indexes install into the same directory, as `index-*.fits` files. The
solver does not use them unless asked, so blind solving means passing that directory
(or one index file) to `-i`:

```bash
arcsec -f image.fits -i "$(arcsec catalog path)" --fov 3
```

---

## 4. Formats, and why there are three

The ASTAP databases come in three on-disk layouts, which arcsec reads transparently.
You do not need to care, but `catalog verify` checks each database's file count
against its layout, so:

| Layout | Grid | Used by |
|---|---|---|
| `.1476` | 1476 tiles, 36 equal-declination rings | D80, D50, D20, D05, V50 |
| `.290` | 290 tiles, 18 equal-area rings | G05, V05 |
| `.001` | one all-sky file | W08 |

Which database uses which is not predictable from its name — `v05` is `.290` even
though it covers the same field range as the `.1476` D-series. Both `catalog` and the
solver probe for whichever is present rather than assuming.

See [plate-solving.md §9.1b](plate-solving.md#91b-the-three-astap-database-formats)
for the format details.

---

## 5. Building a blind index

arcsec can build its own blind-solving index from an installed star database — no
download:

```bash
arcsec catalog index build                  # deepest installed database, fields 0.3°–30°
arcsec catalog index build --min-fov 0.15   # down to D80's narrowest fields
arcsec catalog index info                   # what was built
```

From D80, fields 0.3°–30° take about 2 minutes and 290 MB; down to 0.15° about 13
minutes and 700 MB (`--max-fov` and `--min-fov` choose the range; the build lists the
tiers it will make first). The file, `<db>.arcsecix`, goes into the catalogue directory,
and `catalog list` and `catalog verify` include it.

With it installed:

* `arcsec -f image.fits -i "$(arcsec catalog path)"` solves with no position at all,
  typically in a second or two, and with no pixel scale either if `--fov` and the
  header's FOCALLEN/XPIXSZ are missing.
* Without `-i`, a search of `-r` 10° or more that reaches past five fields round the hint
  (N.I.N.A.'s blind mode sends `-r 180`) tries the index once the first five fields have
  failed, instead of spiralling over the whole sky.

Design and measurements: [offline-index.md](offline-index.md). Generating the star
catalogue itself from Gaia is not possible; `arcsec catalog install` is the supported
route, and it needs no external tools: downloads and unpacking (both `.zip` and `.deb`)
are built in, so nothing has to be on `PATH`.

---

## 6. If you already have ASTAP databases

Point arcsec at them and skip the download entirely:

```bash
export ARCSEC_CATALOG_DIR=~/star_database    # or wherever ASTAP keeps them
arcsec catalog list                          # confirms what it found
```

arcsec reads ASTAP's files directly and the solver never modifies them, so the two can
share a directory. The catalogue manager does act on that directory, though: in a
shared directory `arcsec catalog remove d80` deletes ASTAP's D80 files too (it lists
what it will delete and asks first), and `install` replaces files of the same name.
