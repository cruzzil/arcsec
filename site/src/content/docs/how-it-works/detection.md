---
title: Finding the stars
description: How arcsec estimates the sky background and noise, decides what counts as a star, and measures each star's centre and size.
---

Everything later depends on this step. A star that is missed cannot be matched, and a
noise spike taken for a star is one more thing that matches nothing; the solve can
survive plenty of both, but its accuracy is set by how well each star's centre is
measured.

## Background and noise

The sky in an image is not black: it has a background level (light pollution, moonlight,
the camera's offset) and noise about it. arcsec estimates both from a histogram of the
pixel values, as ASTAP does.

- **Background.** In a star field almost every pixel is sky, so the most common pixel
  value (the peak of the histogram, its *mode*) is the background level $B$. If the
  histogram's mean is more than 1.5 times its mode, the mean is used instead.
- **Noise.** The noise $\sigma$ is the root-mean-square difference from the background
  over a sparse grid of sample pixels:

  $$
  \sigma = \sqrt{\frac{1}{N}\sum_i (p_i - B)^2}
  $$

  Stars would inflate that, so the estimate is **sigma-clipped**: samples more than
  $3\sigma$ from the background are dropped and $\sigma$ is computed again, until it
  changes by less than 5 % (at most seven times).

One histogram for the whole frame is cheap but cannot follow a gradient across it, so the
last detection pass (below) estimates the background locally instead, in a grid of about
12 × 12 cells.

## What counts as a star

A pixel is a candidate if it stands well above the background:

$$
p \;>\; B + k\,\sigma
$$

arcsec searches with a cascade of thresholds, brightest first, and stops as soon as it
has enough stars (`-s`, 500 by default): two levels read from the bright tail of the
histogram, then $30\sigma$, and finally $7\sigma$ above the local background of each grid
cell. A bright, clean field stops at the first level; a faint or nebulous one goes on to
the last. The frame is split into horizontal bands that are scanned in parallel.

Each candidate is then measured, and rejected if it is not star-like: a single hot pixel,
a spot that is too small (`-m`, 1.5″ by default), or one that is too large or not round.
Bright stars that overflow the normal measuring box (saturated discs on photographic
plates, or the wide wings of stars on very coarse cameras) are measured again in a bigger
box, because they are the stars a catalogue is most sure of. Elongated sources (trails,
galaxies, knots of nebula) are refused.

## Centre and size

A star's centre is its centre of light, the brightness-weighted mean position of its
pixels above the background:

$$
\bar x = \frac{\sum_i (p_i - B)\,x_i}{\sum_i (p_i - B)}, \qquad
\bar y = \frac{\sum_i (p_i - B)\,y_i}{\sum_i (p_i - B)}
$$

The sums run over a small aperture, which is re-centred on the result and shrunk until
the centre settles. This finds the centre to a small fraction of a pixel. For a star
whose image has a width $\sigma_\text{PSF}$ and a signal-to-noise ratio $S/N$, the
centre is uncertain by roughly

$$
\sigma_\text{centre} \approx \frac{\sigma_\text{PSF}}{S/N}
$$

so a star 2 pixels wide at $S/N = 50$ is placed to about 0.04 pixel. The signal-to-noise
ratio itself compares the star's total light $F$ with the sky noise over its aperture of
radius $r$:

$$
S/N = \frac{F}{\sqrt{F + \pi r^2 \sigma^2}}
$$

The star's size is its **half-flux diameter** (HFD), measured as ASTAP measures it: twice
the light-weighted mean distance of its pixels from the centre,

$$
\text{HFD} = 2\,\frac{\sum_i (p_i - B)\,r_i}{\sum_i (p_i - B)}
$$

For an ideal Gaussian star of width $\sigma$ that comes to

$$
\text{HFD} = \sqrt{2\pi}\,\sigma \approx 2.507\,\sigma, \qquad
\text{FWHM} = 2\sqrt{2\ln 2}\,\sigma \approx 2.355\,\sigma
$$

so HFD is about 1.06 × FWHM. It is used rather than FWHM because it is an average over
the whole star rather than a reading at one height, so it stays meaningful for
out-of-focus, trailed or blurred stars, whose profile may have no clear half-maximum.

## Keeping the best stars

The stars are sorted by signal-to-noise ratio and the best `-s` are kept. The catalogue
cannot help with more stars than it holds: ASTAP's databases are limited by star
*density* (D50 holds up to 5000 stars per square degree, D80 8000), so arcsec also keeps
no more than

$$
\text{density} \times \text{field area}
$$

image stars, the most the catalogue can have in the field. Fainter image stars would only
make patterns the catalogue cannot contain.

## Binning

Very finely sampled images (under 1″ per pixel) have more pixels than the solve needs, and
their stars are many pixels wide. arcsec then bins the image first, averaging blocks of
$n \times n$ pixels with

$$
n = \operatorname{round}\!\left(\frac{1''}{\text{scale}}\right), \quad n \le 16
$$

which makes detection faster and the stars more compact. The solution is converted back
to the original pixels at the end. `-z` sets the factor by hand.

Next: [Star patterns](../patterns/).
