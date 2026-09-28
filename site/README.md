# arcsec website

The arcsec website: a home page and the user documentation, built with
[Astro](https://astro.build/) and [Starlight](https://starlight.astro.build/) and
published to GitHub Pages at <https://cruzzil.github.io/arcsec/> by
[`.github/workflows/site.yml`](../.github/workflows/site.yml).

## Working on it

Needs Node.js 22.12 or newer (CI uses the current LTS). From this directory:

```sh
npm ci                  # install exactly what package-lock.json pins
npm run dev             # live-reloading server at http://localhost:4321/arcsec/
npm run build           # build the static site into dist/
npm run preview         # serve dist/ at http://localhost:4321/arcsec/
npm run check           # type-check the .astro and .ts files
npm run check-links     # after a build: every internal link and #anchor resolves
```

CI runs `check`, `build` and `check-links` on every pull request that touches `site/`,
and deploys on every push to `main` that does.

## Layout

| Path | What it is |
|---|---|
| `src/content/docs/` | The pages, in MDX. `index.mdx` is the home page. |
| `src/components/` | The home page's hero, download table and proof points; the catalogue picker. |
| `src/lib/catalogues.ts` | The catalogue list and the `catalog recommend` rule, for the picker. |
| `src/lib/release.ts` | Fetches the latest release for the download links. |
| `src/assets/hero/` | The home page's sky image and its RA/Dec grid (generated, see below). |
| `src/styles/custom.css` | Colour changes to Starlight's default theme. |
| `astro.config.mjs` | Site URL and base path, sidebar, Starlight settings. |

Links between pages in `src/content/docs/` are relative (`../nina/`), so they work
whatever the base path. Components build links with `url()` from `src/lib/url.ts`.

## Keeping it in step with the code

The site states facts about the CLI; when those change, the site must too.

- **Catalogues.** `src/lib/catalogues.ts` copies the registry in
  `arcsec/src/catalog_cmd/registry.rs` and the rule in `cmd_recommend`. After changing
  either, update it and compare a few answers with `arcsec catalog recommend --fov <deg>`.
- **Options, exit codes and output files** are described in
  `src/content/docs/reference/cli.mdx`; compare with `arcsec --help`.
- **Benchmark figures** on the home page and in the FAQ come from
  `docs/test-images.md` §6.

## Download links

Release archive names include the version (`arcsec-v0.1.2-x86_64-unknown-linux-gnu.tar.gz`),
so the links cannot be written by hand. At build time `src/lib/release.ts` asks the
GitHub API for the latest release (`repos/cruzzil/arcsec/releases/latest`) and picks out
the archive for each platform, its `.sha256` file and `SHA256SUMS`.

- The workflow rebuilds the site after the Release workflow succeeds, so new versions
  appear without a site change.
- Set `GITHUB_TOKEN` to avoid the API's anonymous rate limit; CI passes its token.
- If the request fails, the build still succeeds: the pages link to the releases page
  instead. `ARCSEC_SITE_OFFLINE=1` skips the request, for working offline.

In the browser, the home page's Download button picks the archive for the visitor's
operating system (from `navigator.userAgentData` or the user agent string). Without
JavaScript, or on a platform with no build, it links to the full table.

## The home page image

`src/assets/hero/sky.png` and `sky-grid.json` are made by
`scripts/make-hero-image.py` (pure Python, no packages needed) from an SDSS frame in
the benchmark corpus and the `.wcs` file arcsec wrote when solving it. The script's
docstring has the exact commands. The page credits the image to the Sloan Digital Sky Survey.

## Moving to a custom domain

In `astro.config.mjs`, set `site` to the new origin and `base` to `'/'`; add
`public/CNAME` containing the domain; and set the domain under the repository's
**Settings > Pages**. Nothing else needs to change.
