---
title: Checking the answer
description: How arcsec decides whether a fitted plate is real - counting the catalogue stars that land on detected stars, against how many would by chance - and what it does when the ordinary search finds nothing.
---

A fit always produces *some* answer, even from coincidences. Before arcsec reports a
solution it tests it, by asking whether the sky it claims to show actually matches the
image. A wrong answer is worse than none: a mount that syncs to a false solution slews
to the wrong place with confidence. On the test images arcsec reports one wrong answer in
596 (see [arcsec compared with ASTAP](../../compared/)).

## Counting stars

The test reuses the star pairing from [the fit](../fitting/#from-quads-to-stars): every
catalogue star in the field is projected into the image through the fitted plate and
paired with a detected star within 2 pixels, one-to-one. A correct solution typically
pairs 200 to 375 stars. A wrong one, which matched on a few coincidental quads, pairs
almost none beyond them.

A plate is accepted only if all of these hold:

- **Enough stars.** At least 30 stars paired (fewer for images with very few stars, see
  below).
- **Spread over the frame.** The paired stars' root-mean-square distance from their own
  centre is at least 0.2 of the image's half-diagonal. A fit pinned down by stars in one
  corner says little about the rest.
- **A tight fit.** The stars' rms residual is within the 2-pixel match radius.
- **More than chance.** See the next section.

## How many would match by chance

In a crowded image, a catalogue star projected to a random spot will often find *some*
detection within 2 pixels. If the image has $\rho$ detections per square pixel, scattered
at random, the chance that at least one lies within radius $r$ of a given point is

$$
P = 1 - e^{-\rho\,\pi r^2}
$$

(the Poisson probability of not finding zero). A plate that puts $N$ catalogue stars in
the frame should therefore pair about

$$
\mathbb{E}[\text{chance pairs}] = N \left(1 - e^{-\rho\,\pi r^2}\right)
$$

stars by luck alone. arcsec accepts a plate only if it pairs at least **four times** that
many. In a sparse field this is never the limit; in a dense one, it is what stops the
30-star rule from passing a wrong plate. On a 2.2° TESS frame with 500 detections on
384 × 384 pixels, a random point has a detection within 2 pixels 4 % of the time, so a
wrong plate putting 450 catalogue stars in the frame finds 18 matches by chance and can
creep to 30 as the fit follows them. The wrong plates seen there had 1.5 to 2.7 times
the chance count; every correct solve on the test images has at least 7.9.

## Images with few stars

An image with 50 detected stars cannot pair 30 of them reliably, even when correct, so the
minimum relaxes to

$$
\min\!\left(30,\ \max\!\left(10,\ \lceil 0.15\,n\rceil\right)\right)
$$

for $n$ detected stars. A handful of pairs is easy to fit by accident, though, so below 30
two more conditions apply: the solved pixel scale must be within 10 % of the one given
(from `-fov` or the header), and the fit's rms must be under half a pixel. A few
coincidences can be fitted by *some* plate at *some* scale; they are not fitted at the
scale the optics give, to half a pixel. When the scale is not known and arcsec is trying
a ladder of them, the relaxation is not used and 30 stars are required.

## When nothing is found

The quad search relies on the image and the catalogue agreeing, roughly, on which stars
are the brightest in each neighbourhood. In fields crowded with nebulosity, or where the
camera's colour band ranks stars differently from the catalogue's, they may not: in some
Milky Way fields, none of the 30 brightest catalogue stars is among the image's 500
detections, because they are saturated, blended or buried in nebula.

So if the whole spiral finds nothing, arcsec makes one more search about the hint, based
on the rank-robust fallback of [seiza](https://github.com/theatrus/seiza). It uses
brightness on the catalogue side only, where it can be trusted:

1. Quads are made from the catalogue's 150 brightest stars in a window 1.5 fields across.
2. The image's 2000 brightest detections, *without* regard to their order, go into a
   lookup table of positions and a table of every pair sorted by length.
3. For each catalogue quad, every image pair as long as the quad's longest side (at the
   known scale, ± 5 %) is tried as that side, in both orders and both mirrorings. Each
   choice fixes a similarity transform, and the quad's other two stars must then land
   within 2.5 pixels of detections.
4. A transform that passes puts all 150 bright catalogue stars on the image, and must find
   detections for at least 10 of them, and at least twice as many as chance would give.
5. Survivors are fitted and checked exactly as above.

The search has a fixed work budget of about a third of a second, so an image that cannot
be solved still fails in bounded time.

## The other direction: what the checks cost

Every rule here was added after a wrong answer on the test images, and each was measured
against the correct solves so that it does not lose them. The full record, including the
images that motivated each rule, is in the repository's
[engineering notes](https://github.com/cruzzil/arcsec/blob/main/docs/plate-solving.md).

Next: [Blind solving](../blind/).
