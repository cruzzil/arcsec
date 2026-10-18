---
title: Blind solving
description: How arcsec finds a field anywhere on the sky with no position hint, using a pre-built index of star patterns, and why every blind answer is still checked by the ordinary solver.
---

The spiral search is fast when the hint is close and slow when it is not: its cost grows
with the square of the search radius in fields, so searching the whole sky for a
0.5° field means tens of thousands of positions. A **blind** solve turns the search
around. Instead of building the catalogue's patterns at each place it looks, it builds
them once, for the whole sky, in advance, so that each image pattern can be looked up
directly.

## The index

`arcsec catalog install` builds the index after downloading a star database, from that
database, typically in a minute or less (a few hundred MB on disk). Building it
means deciding, ahead of time and without seeing any image, which patterns an image is
likely to contain.

The answer arcsec uses (an idea from [seiza](https://github.com/theatrus/seiza)'s blind
index) is to anchor patterns on stars that are **locally the brightest**. For a disc of
radius $r$ on the sky, a star is an anchor if no brighter star lies within $r$ of it.
Its group is itself and the next brightest stars within its disc (six in all), and every
four of those six make a pattern:

$$
\binom{6}{4} = 15 \text{ patterns per anchor}
$$

An image that contains that disc sees the same locally brightest stars, so it can rebuild
the same group from its own brightest detections, even without knowing where it is. The
15 subsets mean a group member that is missing or ranked differently in the image loses
some patterns, not all.

One disc size cannot suit every field, so the index has several **tiers**: disc radii from
3° down to 0.1°, each with a magnitude limit chosen so that a disc holds 15 to 20 stars on
average. Each tier serves fields from about 2.5 to 12 times its radius, and a field
usually sees two.

## The key

Each pattern is described by the same five ratios as the [quads](../patterns/) of the
ordinary search. For a quick lookup, each ratio is rounded into one of 128 bins and the
five bin numbers packed into one 64-bit key,

$$
\text{key} = \sum_{k=1}^{5} \left\lfloor 128\,\frac{d_{k+1}}{d_1} \right\rfloor \cdot 256^{\,5-k}
$$

and the keys are stored sorted, so a lookup is a binary search for an exact value. A bin is
$1/128 \approx 0.0078$ wide, close to the ordinary search's tolerance. A ratio that falls
within 0.0025 of a bin edge also tries the neighbouring bin, so a measurement that lands
just across a boundary is not lost. Because the ratios do not change under mirroring, one
entry serves both normal and mirrored images.

The ratios do not say which star of a pattern is which, so both sides order a pattern's
stars the same way: by the total of each star's distances to the other three. That gives
a star-to-star correspondence, and from it a transform.

## Looking up an image

1. **Detect** every star and rank them by total light, not signal-to-noise: the index is
   built from each region's brightest stars, and a saturated star has a lower
   signal-to-noise ratio than a fainter, sharper one.
2. **Build patterns** from the image's brightest stars, and from the five brightest within
   discs of several sizes around each of the 150 brightest: the image's version of the
   anchor groups. About 7000 patterns.
3. **Look up** each pattern in the tiers that suit the image's size. A hit fixes a
   transform from the four star pairs, which must be within 6 % of a similarity (no
   stretch or shear), at a pixel scale in the allowed range (within 20 % of the given
   scale, or 0.3″ to 60″ per pixel if none is given). Each hit is a hypothesis: a centre
   on the sky, a scale and a rotation.
4. **Vote.** Hypotheses fall into buckets of right ascension, declination and log of
   scale: 5 % of the field on the sky (the right-ascension width divided by
   $\cos\delta$) and 5 % in scale. Real hits pile up in one bucket; chance hits scatter
   across the sky.
5. **Rank** the strongest buckets by projecting the index's stars through each hypothesis
   onto the image and counting how many land on detections. The score is a
   significance rather than a raw count,

   $$
   \text{score} = \frac{m - E}{\sqrt{E + 1}}
   $$

   for $m$ matches where $E$ would be expected by chance, so that a hypothesis at far too
   coarse a scale, which projects hundreds of stars into the frame and catches many by
   luck, cannot outrank the truth.

## The answer is checked like any other

The index never decides on its own that a field is found. The best few hypotheses are
handed to the ordinary solver as hints, each with its scale as the field size and no
search radius, and the solve has to pass [the same checks](../verification/) as a hinted
one: at least 30 stars matched, spread over the frame, well above chance. A blind solve
is therefore held to exactly the hinted solver's standard. On the test images it finds 85
% of the fields the hinted solver finds, in a median half a second, with no wrong
position.

## When the index is used

- With `-i`, the index is tried first, whatever the hint.
- Without `-i`, an installed index is used automatically for wide searches: when the
  search radius `-r` is at least 10° and reaches more than five fields from the hint. The
  spiral searches those first five fields itself, so a good hint is unaffected. This is
  what makes N.I.N.A.'s blind mode (`-r 180` and no position) take seconds rather than
  minutes.
- If the index finds nothing, the ordinary spiral search runs as before.

## Astrometry.net indexes

`-i` also accepts the index files of [Astrometry.net](https://astrometry.net/) (for
example `arcsec catalog install anet-4100`). These describe a quad differently: the two
stars furthest apart, $A$ and $B$, define a frame in which $A$ is at $(0, 0)$ and $B$ at
$(1, 1)$, and the code is the position of the other two stars in that frame,

$$
\text{code} = (x_C,\ y_C,\ x_D,\ y_D)
$$

With $\mathbf{ab} = B - A$ and $s = |\mathbf{ab}|^2$, a point $P$ maps to

$$
\cos\theta = \frac{ab_y + ab_x}{s}, \quad \sin\theta = \frac{ab_y - ab_x}{s}, \qquad
\begin{aligned}
x &= -(P - A)_x \sin\theta + (P - A)_y \cos\theta \\
y &= \phantom{-}(P - A)_x \cos\theta + (P - A)_y \sin\theta
\end{aligned}
$$

Swapping $A$ with $B$, or $C$ with $D$, gives a different code for the same quad, so one
labelling is chosen as canonical: $x_C \le x_D$ and $x_C + x_D \le 1$. Unlike the ratios, this
code carries the star correspondence with it, but it changes under mirroring, so arcsec
builds every image quad twice, once per handedness. Matches vote into 0.1° sky bins and
the best are verified as above. arcsec's own index is the default because it is built from
the star database already installed, covers narrower fields, and finds more of them.
