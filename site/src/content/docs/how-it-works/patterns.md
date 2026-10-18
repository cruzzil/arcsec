---
title: Star patterns
description: How arcsec describes groups of four stars by five distance ratios that survive rotation, scaling and mirroring, and why four stars work far better than three.
---

The solver cannot pair individual stars, because one star looks like any other. It pairs
*patterns* of stars instead, described in a way the camera cannot change.

## What the description must survive

Between the sky and the image, a star pattern is moved (the telescope points somewhere),
rotated (the camera is turned), scaled (the focal length and pixel size), and possibly
mirrored (a diagonal or a flipped readout). Together these make a **similarity
transform**:

$$
\begin{pmatrix} x' \\ y' \end{pmatrix} =
s \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}
\begin{pmatrix} \pm x \\ y \end{pmatrix} +
\begin{pmatrix} t_x \\ t_y \end{pmatrix}
$$

with scale $s$, rotation $\theta$, shift $(t_x, t_y)$ and a sign for the mirror. A
description of a pattern is useful only if it gives the same numbers before and after
any such transform: it must be **invariant**.

Distances between stars are unchanged by moving, rotating and mirroring, and all of them
are multiplied by the same $s$ when the image is scaled. So any **ratio** of two distances
is invariant to all four.

## The quad and its five ratios

Take four stars. Between them there are six distances. Sort them from longest to
shortest, $d_1 \ge d_2 \ge \dots \ge d_6$, and divide each by the longest:

$$
\text{code} = \left(\frac{d_2}{d_1},\ \frac{d_3}{d_1},\ \frac{d_4}{d_1},\ \frac{d_5}{d_1},\ \frac{d_6}{d_1}\right) \in (0, 1]^5
$$

<figure class="hiw-figure">
<svg viewBox="0 0 300 200" role="img" aria-labelledby="quad-title" width="300" height="200">
<title id="quad-title">Four stars A, B, C and D joined by their six distances; the longest, A to C, is d1.</title>
<g stroke="currentColor" stroke-opacity="0.45" stroke-width="1.5" fill="none">
<line x1="40" y1="40" x2="250" y2="55"/>
<line x1="40" y1="40" x2="70" y2="160"/>
<line x1="250" y1="55" x2="215" y2="170"/>
<line x1="70" y1="160" x2="215" y2="170"/>
<line x1="250" y1="55" x2="70" y2="160"/>
</g>
<line x1="40" y1="40" x2="215" y2="170" stroke="var(--sl-color-accent)" stroke-width="3"/>
<g fill="currentColor">
<circle cx="40" cy="40" r="6"/><circle cx="250" cy="55" r="5"/>
<circle cx="215" cy="170" r="5.5"/><circle cx="70" cy="160" r="4"/>
</g>
<g font-size="14" fill="currentColor" font-family="sans-serif">
<text x="22" y="34">A</text><text x="258" y="50">B</text>
<text x="222" y="188">C</text><text x="50" y="180">D</text>
<text x="140" y="98" fill="var(--sl-color-accent)" font-weight="700">d₁</text>
</g>
</svg>
<figcaption>Four stars, six distances. Dividing the other five by the longest, d₁, gives
five numbers that do not depend on where the stars are, how the image is turned, how
large it is, or whether it is mirrored.</figcaption>
</figure>

Sorting is what makes it work: the code does not depend on which star is called A, so
there is no labelling to get right. Two more numbers are kept with each quad, though not
in its code: the longest distance $d_1$ itself, and the direction of that longest side.
Two matching quads then also say what scale and rotation turn one into the other, which
the next step uses.

Quads are matched with a tolerance $t$ on each ratio (`-t`, 0.007 by default): an image
quad and a catalogue quad match if all five ratios agree within $\pm t$.

## Why four stars and not three

Three stars give three distances and so only **two** ratios. That seems nearly as good,
but it is not, and the reason is the most useful piece of arithmetic in plate solving.

Suppose each ratio is spread roughly evenly over its range, and matched to $\pm t$. The
chance that an unrelated pattern matches on one ratio is about $2t$, and on all $N$
ratios

$$
P_\text{chance} \approx (2t)^N
$$

With $M$ image patterns and $R$ catalogue patterns compared, about $M \cdot R \cdot (2t)^N$
false matches turn up by chance. With $t = 0.007$, 500 image patterns and 5000 catalogue
patterns ($M R = 2.5 \times 10^6$):

| Pattern | Ratios $N$ | $(2t)^N$ | False matches expected |
|---|---|---|---|
| Triangle | 2 | $2 \times 10^{-4}$ | about 490 |
| Quad | 5 | $5 \times 10^{-11}$ | about 0.0001 |

Each extra ratio divides the false matches by about 70. With triangles, the true
matches are buried in hundreds of false ones; with quads, almost every match is real.
(The ratios are not really spread evenly, so these are orders of magnitude, but the gap
between the two is real.)

Triangles have one advantage: a pattern of three stars survives more often when stars
are missing. If a fraction $p$ of the catalogue's stars are detected, a triangle survives
with probability $p^3$ and a quad with $p^4$. arcsec has an experimental triangle mode,
`--method tetra`, which tightens the tolerance to $0.3\,t$ to bring the false matches
down and relies on a stricter consensus step; on the test images it solves far fewer
fields than quads.

## Which quads are built

A quad made of four stars from opposite corners of the image is fragile: a little lens
distortion changes it, and any one of the four missing from the catalogue loses it. So
quads are made from **neighbourhoods**. For each star, arcsec takes it and its eight
nearest neighbours and builds every quad from those nine:

$$
\binom{9}{4} = \frac{9!}{4!\,5!} = 126 \text{ quads per star}
$$

(With fewer than 30 stars the group is smaller.) Many of these are near-copies of each
other, and quads whose centres fall within one pixel of an existing quad's are dropped.

The redundancy matters because the image and the catalogue never agree on which stars
are present. ASTAP builds one quad per star, from its three nearest neighbours; if one of
those neighbours is missing from the catalogue, that quad is lost. With 126 quads from
each neighbourhood, many survive whichever star is missing. It is one of the main reasons
arcsec solves fields ASTAP does not, at the cost of more patterns to match.

The catalogue side builds its quads by exactly the same rule, at each place in the sky
the search tries.

Next: [Searching the sky](../search/).
