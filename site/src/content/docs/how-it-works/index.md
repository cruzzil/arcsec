---
title: How plate solving works
description: The process and the mathematics behind an arcsec solve, from finding the stars in an image to the World Coordinate System it writes.
---

A plate solver looks at the stars in an image and works out exactly which part of the
sky it shows. This section explains how arcsec does that, one step at a time, with the
mathematics behind each step. Nothing here is needed to *use* arcsec; it is for anyone
who wants to know what happens between pressing "solve" and getting an answer.

## The problem

An image is a grid of pixels, each with a position $(x, y)$. The sky is a sphere, where
each star has a right ascension $\alpha$ and a declination $\delta$. Solving means
finding the function that takes one to the other:

$$
(x,\ y)\ \text{pixels} \quad\longrightarrow\quad (\alpha,\ \delta)\ \text{on the sky}
$$

For a telescope this function is simple: up to a small distortion, the camera takes a
flat projection of the sky, and the mapping is fixed by six numbers (a centre, a scale,
a rotation, and whether the image is mirrored). arcsec writes those numbers as a
**World Coordinate System** (WCS), the FITS standard for it: the sky position of a
reference pixel (`CRVAL1`, `CRVAL2`), that pixel (`CRPIX1`, `CRPIX2`), and a 2 × 2
matrix (`CD1_1` … `CD2_2`) that turns pixel offsets into angles on the sky.

What makes it hard is that the two lists of stars do not agree. The image has hot
pixels, satellites, galaxies and noise that are not in any star catalogue; the
catalogue has stars too faint, too bright or too crowded to be found in the image. And
nothing says which image star is which catalogue star. With a few hundred stars on each
side there are more ways to pair them than atoms in the universe, so the pairing cannot
simply be tried.

## The idea that makes it possible

The trick every plate solver uses is to describe small groups of stars by their
**shape**, in numbers that do not change when the image is moved, rotated, enlarged or
mirrored. Four stars that make a particular shape on the sky make the same shape in the
image, whatever the camera did. So instead of pairing stars, the solver pairs shapes:
it looks up each image shape among the catalogue's shapes, and a handful of agreeing
lookups is enough to fix the whole mapping.

arcsec uses the four-star shapes ("quads") described by Han Kleijn for
[ASTAP](https://www.hnsky.org/astap.htm), and adds steps of its own around them: it
checks every candidate solution star by star, re-fits it on the individual stars, models
lens distortion, and can solve with no position hint at all.

## The steps

```
image ──► 1 find the stars ──► 2 describe them as quads
                                         │
              catalogue near the hint ──►│
                                         ▼
          3 search the sky: read the catalogue at a trial position,
            build its quads, look up each one among the image's
                                         │
                                         ▼
          4 fit the plate: the six numbers that best map one onto the other
                                         │
                                         ▼
          5 check it: do the catalogue's stars land on the image's?
                     │                              │
                  yes▼                            no▼
               write the WCS               next trial position
```

1. **[Finding the stars](detection/).** Estimate the sky background and its noise, find
   the spots that rise clearly above it, and measure each one's centre to a fraction of
   a pixel.
2. **[Star patterns](patterns/).** Group each star with its neighbours into quads, and
   describe each quad by five ratios of its distances. Why five numbers is enough, and
   three-star triangles are not.
3. **[Searching the sky](search/).** Project the catalogue onto a flat plane around a
   trial position, build the same quads there, match them against the image's, and
   keep only the matches that agree on one scale and rotation. If the trial position
   is wrong, spiral outwards.
4. **[Fitting the plate](fitting/).** Solve for the six plate constants by least
   squares, turn them into the WCS, and model lens distortion where the field needs it.
5. **[Checking the answer](verification/).** Project every catalogue star through the
   fit and count how many land on a detected star, compared with how many would by
   chance. A solution is reported only if the count is convincing.

[Blind solving](blind/) is a different first step for when there is no position hint:
it looks shapes up in a pre-built index of the whole sky, and hands the position it
finds to the steps above.

## How much is known in advance

Two things make a search fast: knowing roughly where the telescope points, and knowing
the pixel scale (how many arcseconds of sky each pixel covers, which follows from the
focal length and the camera's pixel size). N.I.N.A. and most capture programs pass both,
and then the answer is usually found at the first trial position, in a tenth of a
second. arcsec can manage without either: [without a scale](search/#when-the-scale-is-unknown)
it tries a range of them, and [without a position](blind/) it uses a blind index.

## Further reading

The pages here simplify. The full engineering notes, with every constant, the reasons
for each design choice, measurements, and the ideas that were tried and dropped, are in
[`docs/plate-solving.md`](https://github.com/cruzzil/arcsec/blob/main/docs/plate-solving.md)
in the repository. The other main sources are:

- Han Kleijn's [description of ASTAP's algorithm](https://www.hnsky.org/astap_astrometric_solving.htm),
  where the quad method comes from.
- D. Lang et al., *Astrometry.net: Blind astrometric calibration of arbitrary
  astronomical images*, Astronomical Journal 139 (2010), 1782, for the index-based blind
  solve and its verification.
- M. R. Calabretta and E. W. Greisen, *Representations of celestial coordinates in
  FITS*, Astronomy & Astrophysics 395 (2002), 1077, for the projection and the WCS.
