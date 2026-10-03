# Installing Star Catalogues

arcsec cannot solve anything without a star catalogue, and cannot do photometric
colour calibration without a *photometric* one. Both are one command away:

```bash
arcsec catalog recommend --like ~/lights/M42_0001.fits --photometry
arcsec catalog install d50 v05
```

Everything lands in one place that the solver already reads, so after installing you
can just run `arcsec -f image.fits` with no `-d` and no `-D`. Installing a solving
database also builds arcsec's **blind index** from it (§5), so a solve needs no position
hint either; the confirmation shows its size, build time and memory before anything
starts.

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
* **Blind solving** (no position hint, or N.I.N.A.'s blind mode): nothing to download.
  `install` builds arcsec's own index from the solving database (§5). The
  Astrometry.net sets (`anet-4100`, `anet-5200`) still work with `-i`, but arcsec's
  index is faster and finds more fields; they only find the rough position, and the
  solve is then finished against a star database, so you need one of those as well.

---

## 2. Installing

```bash
arcsec catalog install d50               # asks for confirmation, shows the size first
arcsec catalog install d50 v05 g05 --yes # several at once, no prompt
arcsec catalog install d80 --no-index    # download only, no blind index
arcsec catalog list                      # everything, with what is installed
arcsec catalog verify                    # check for missing or truncated files
arcsec catalog remove d80                # free the space again
arcsec catalog path                      # where they live
```

`install` and `remove` both list what they will do and ask first; `--yes` (or `-y`)
skips the prompt. `--keep-archive` keeps the downloaded `.zip` or `.deb` after
unpacking. With no terminal to ask (a script, a scheduler) and no `--yes`, both cancel
rather than guess.

Installing a solving database (`d05`, `d20`, `d50`, `d80`, `g05`, `w08`) then builds the
blind index, and the confirmation says what that costs alongside the download:

```text
$ arcsec catalog install d50
Installing into /home/me/.local/share/arcsec/catalogs

  d50          901.3 MB  Gaia DR3 to 5000 stars/deg². The usual choice for solving.

Total download: 901.3 MB
Then build a blind index from D50, fields 0.3°–30° (no blind index is installed yet):
  ~287.4 MB on disk, ~25 s, ~0.6 GB memory
  (--no-index to skip it)

Continue? [y/N]
```

Free disk space is checked before anything is downloaded: the archives, their unpacked
files and the index must fit, with 10 % and 100 MB to spare, or the install stops with
the numbers. §5 has the index's sizes, the options and when you are asked separately.

Installs are **resumable**: interrupt one and re-run the same command, and the
download continues from where it stopped rather than starting again, provided the
server supports range requests. Files are written
to a temporary name and only moved into place once complete, so a half-finished
download can never look like an installed catalogue. The same goes for the index: an
interrupted build leaves nothing behind, and re-running the install, which finds the
database already there, offers the index again.

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

The blind index arcsec builds sits there too, as `<db>.arcsecix` (`d50.arcsecix`), and
the solver uses it without being asked when a search is wide (§5).

Astrometry.net indexes install into the same directory, as `index-*.fits` files. The
solver does not use them unless asked, so blind solving with them means passing that
directory (or one index file) to `-i`:

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

## 5. The blind index

arcsec finds a field with no position hint (and no pixel scale, if need be) using its
own pattern index, built from a star database you already have rather than downloaded.
`catalog install` builds it after installing a solving database; for databases that
were installed before that (arcsec 0.4 and earlier) or that came from ASTAP, build it
once by hand:

```bash
arcsec catalog index build                  # from the installed databases, as install would
arcsec catalog index build --min-fov 0.15   # down to D80's narrowest fields
arcsec catalog index info                   # what was built, from what
```

With it installed:

* `arcsec -f image.fits -i "$(arcsec catalog path)"` solves with no position at all,
  typically in a second or two, and with no pixel scale either if `--fov` and the
  header's FOCALLEN/XPIXSZ are missing.
* Without `-i`, a search of `-r` 10° or more that reaches past five fields round the hint
  (N.I.N.A.'s blind mode sends `-r 180`) tries the index once the first five fields have
  failed, instead of spiralling over the whole sky. Without an index, such a search
  prints a one-line hint on stderr with the command and its cost; it never builds one
  during a solve.

### What gets built

One index per catalogue directory, built from the **deepest** solving database there
and covering the fields of **every** one:

| Database | Fields covered by default | Index | Build time* | Peak memory |
|---|---|---|---|---|
| `d80`, `d50`, `d20` | 0.3°–30° | 287 MB | 22–30 s | 0.5–0.6 GB |
| `d80` with `--index-min-fov 0.15` | 0.15°–30° | 698 MB | 1.5–2 min | 1.3 GB |
| `d05` | 0.6°–30° | 145 MB | 7–10 s | 0.3 GB |
| `g05` | 3°–30° | 12.5 MB | 1 s | 0.1 GB |
| `w08` | 10°–80° | 0.8 MB | instant | — |

\* Measured on a 24-thread desktop; four threads take about 1.7× as long, one thread
about 3.4×. The figures in the prompt come from a cost model fitted to these builds,
scaled to your thread count.

The ranges follow what each database can verify: the index only proposes a position,
and the ordinary solver then has to confirm it against the database, so an index for
fields D05 cannot solve would be wasted. D80 and D50 stop at 0.3° rather than their
0.15°/0.2° floor because the extra tier for 0.15–0.3° fields more than doubles the index
and finds about half of those fields blind ([offline-index.md §7.1](offline-index.md));
`--index-min-fov 0.15` adds it, and `catalog recommend` suggests it for fields narrower
than 0.3°. The tiers down to 0.3° fields are practically the same whichever database
they are built from, which is why D05, D20 and D80 give the same index there.

Install `w08` beside `d50` and the index widens to 80°; install `d80` where a `g05` index
was and it is rebuilt from D80, replacing the old one. A rebuild keeps tiers the
existing index already had, so a `--min-fov 0.15` index stays one.

### When you are asked separately

An ordinary build rides on the install's confirmation. One that is large gets a notice
and its own question, so the download can be accepted and the index declined:

* the index would be over **1 GB** on disk, or
* the build is expected to take over **5 minutes**, or
* it needs more than **half the memory available** (relevant on a Raspberry Pi or a
  small mini-PC: the default index needs about 0.6 GB while it builds).

```text
Note: this blind index is a large job:
  - it needs ~1.3 GB of memory, more than half of the 1.9 GB available

Continue? [y/N] y
Build the blind index as well? [y/N]
```

`--yes` answers both. `catalog index build` asks the same question only when it is run
at a terminal; from a script it prints the notice and builds, as it always has. A build
that would not fit on the disk is refused outright.

### Options

| Command | Option | Effect |
|---|---|---|
| `install` | `--no-index` | Download only. |
| `install` | `--index-min-fov <deg>` | Smallest field (short side) the index serves, e.g. `0.15`. |
| `install` | `--index-max-fov <deg>` | Largest field the index serves. |
| `index build` | `--min-fov`, `--max-fov` | The same, for a manual build; default from the installed databases. |
| `index build` | `-D <db>`, `--db <dir>`, `-o <file>` | Build from another database or directory, or to another file. |
| `index build` | `--threads <n>`, `--yes` | Limit threads; do not ask before a large build. |
| `remove` | `--keep-index` | Keep the index built from a removed database. |

### Keeping it current

* `catalog list` shows the index, its source and the fields it serves; where there is
  none, it shows the command to build one and what it would cost.
* `catalog verify` checks every section's checksum, and whether the database it was
  built from has **changed** since (a new ASTAP release, a re-download, another copy):
  the index records a fingerprint of the source database's files (their names, sizes
  and first 4 kB, not modification times, so copying a database elsewhere does not make
  its index stale). A stale index is reported as a problem; `catalog install <db>` or
  `catalog index build` rebuilds it. Indexes built by arcsec 0.4 carry no fingerprint and
  are reported as unchecked, not stale.
* `catalog remove d80` lists `d80.arcsecix` with the database and removes both, unless
  `--keep-index`; the index works on its own, but nothing would keep it current.

Builds are atomic: nothing is written until the whole index is built in memory, and then
it goes to a `.part` file renamed into place. Ctrl-C during the build leaves nothing;
during the write, the `.part` is removed (on Windows, by the next build). `catalog list`
mentions any `.part` left over.

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
arcsec catalog index build                   # once: the blind index, ~287 MB from D80
```

`catalog list` says when there is no blind index and what building one costs. The
index is the one file arcsec adds to the directory (`d80.arcsecix`); ASTAP ignores it.

arcsec reads ASTAP's files directly and the solver never modifies them, so the two can
share a directory. The catalogue manager does act on that directory, though: in a
shared directory `arcsec catalog remove d80` deletes ASTAP's D80 files too (it lists
what it will delete and asks first), and `install` replaces files of the same name.
