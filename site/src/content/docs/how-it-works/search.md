---
title: Searching the sky
description: How arcsec projects the catalogue onto the image plane, matches its quads against the image's, finds the matches that agree, and spirals outward from the hint.
---

With the image's quads built, the search compares them with quads built from the star
catalogue around a trial position, starting at the position hint (from `-ra`/`-spd`, or
the image header). If no consistent match turns up there, it moves on.

## Flattening the sky

The catalogue gives stars as $(\alpha, \delta)$ on a sphere; the image is flat. Before
anything can be compared, the catalogue around the trial position $(\alpha_0, \delta_0)$
is projected onto a plane touching the sphere at that point, the **gnomonic** or
**tangent-plane** projection, which is what an ideal camera does. With
$\Delta\alpha = \alpha - \alpha_0$,

$$
D = \sin\delta_0 \sin\delta + \cos\delta_0 \cos\delta \cos\Delta\alpha
$$

$$
\xi = \frac{\cos\delta \,\sin\Delta\alpha}{D}
$$

$$
\eta = \frac{\cos\delta_0 \sin\delta - \sin\delta_0 \cos\delta \cos\Delta\alpha}{D}
$$

$(\xi, \eta)$ are the **standard coordinates**: angles on the plane, measured from the
tangent point, which arcsec keeps in arcseconds. (In the code $\xi$ carries a minus sign,
the FITS convention that right ascension increases to the left.) The projection is
exact for an undistorted lens, maps every great circle to a straight line, and becomes
unusable only towards 90° from the tangent point, far outside any image.

Going back, from the plane to the sky:

$$
\alpha = \alpha_0 + \operatorname{atan2}\!\left(\xi,\ \cos\delta_0 - \eta\sin\delta_0\right)
$$

$$
\delta = \arcsin\frac{\sin\delta_0 + \eta\cos\delta_0}{\sqrt{1 + \xi^2 + \eta^2}}
$$

On this plane, the image is just a moved, turned, scaled and perhaps mirrored copy of the
projected catalogue, exactly the transform the quads are invariant to.

## Reading the catalogue

At each trial position arcsec reads the brightest catalogue stars in a square around it,
somewhat larger than the image: twice the field for sparse images (under 35 stars),
shrinking to the field itself for rich ones (over 140), so that a hint some way off still
puts the image inside the window. It reads as many stars as the image will have in that
area, projects them, and builds their quads by the same rule as the image's.

The star databases are tiled so that a read touches only the tiles that overlap the
window, and the files are memory-mapped rather than loaded, so a read costs little.

## Matching quads

Each catalogue quad is now looked up among the image quads. Comparing every pair would
be slow (tens of thousands on each side), so the image quads are put once, at the start
of the solve, into a grid on two of their ratios, with cells one tolerance $t$ wide. A
catalogue quad can only match image quads in the 3 × 3 cells around its own, and only
those are compared on all five ratios. This turns a scan of every pair into a few
comparisons per quad.

## Keeping the matches that agree

Even five ratios let a few chance matches through, and at a wrong position they are all
there is. The test that separates them is that true matches all describe **the same
transform**. Each match of an image quad with a catalogue quad implies

$$
\text{scale } \rho = \frac{d_1^\text{image}}{d_1^\text{catalogue}}, \qquad
\text{rotation } = \varphi^\text{image} - \varphi^\text{catalogue}
$$

where $\varphi$ is the direction of a quad's longest side. True matches share one scale
and one rotation; chance matches scatter. So each match casts a **vote** into a grid of
scale (bins 5 % wide) by rotation (bins 10° wide), and the fullest cell wins: a Hough
transform, in the space of transforms.

A mirrored image reverses angles, so for a mirrored image it is the **sum**
$\varphi^\text{image} + \varphi^\text{catalogue}$ that is constant, not the difference.
arcsec keeps a second grid for the sum and takes whichever peak is higher, so it need not
know in advance whether the image is mirrored. Directions are taken modulo 180°, because
a side has no preferred end.

The position goes on to a fit only if at least

$$
3 + \frac{n}{140}
$$

quads agree, for $n$ image stars. Otherwise the search moves to the next position.

## The spiral

If the hint is right, the first position is the answer. If not, the search walks
outward in a square spiral, one field width per step:

<figure class="hiw-figure">
<svg viewBox="0 0 300 200" role="img" aria-labelledby="spiral-title" width="300" height="200">
<title id="spiral-title">A square spiral of trial positions: 0 at the hint, then 1 to the east, 2 north-east, 3 north, and round, one field per step.</title>
<g stroke="currentColor" stroke-opacity="0.25" fill="none">
<rect x="60" y="10" width="180" height="180"/>
<line x1="120" y1="10" x2="120" y2="190"/><line x1="180" y1="10" x2="180" y2="190"/>
<line x1="60" y1="70" x2="240" y2="70"/><line x1="60" y1="130" x2="240" y2="130"/>
</g>
<polyline points="150,100 210,100 210,40 150,40 90,40 90,100 90,160 150,160 210,160 255,160" fill="none" stroke="var(--sl-color-accent)" stroke-width="2"/>
<g font-size="15" font-family="sans-serif" fill="currentColor" text-anchor="middle" dominant-baseline="central">
<circle cx="150" cy="100" r="13" fill="var(--sl-color-accent)"/>
<text x="150" y="100" fill="var(--sl-color-black)" font-weight="700">0</text>
<text x="210" y="88">1</text><text x="210" y="28">2</text><text x="150" y="28">3</text>
<text x="90" y="28">4</text><text x="78" y="100">5</text><text x="90" y="172">6</text>
<text x="150" y="172">7</text><text x="210" y="172">8</text>
</g>
</svg>
<figcaption>Trial positions, one field apart. The search goes round in rings until a
position verifies or it has covered the search radius.</figcaption>
</figure>

Position $(s_x, s_y)$ in the spiral is the sky position

$$
\delta = \delta_\text{hint} + s_y \cdot \text{FOV}, \qquad
\alpha = \alpha_\text{hint} + \frac{s_x \cdot \text{FOV}}{\cos\left(\delta \mp \tfrac{1}{2}\text{FOV}\right)}
$$

The cosine is there because lines of right ascension crowd together towards the poles:
one degree of sky east–west is $1/\cos\delta$ degrees of right ascension. It is taken at
the edge of the field nearer the pole, where the crowding is worst, so that neighbouring
positions overlap rather than leave gaps. Positions
further from the hint than the search radius (`-r`) are skipped.

Positions are tried on several threads at once, but the answer is always the
lowest-numbered position that passes, so the result is the same whatever the number of
threads.

A search that succeeds usually stops at position 0, in a tenth of a second. A search
that fails has to visit every position out to the radius, roughly

$$
\pi \left(\frac{r}{\text{FOV}}\right)^2
$$

of them: about 6000 for a 0.25° field searched to 10°. This is why a good hint and a
modest `-r` matter for narrow fields, and why giving up on an image that cannot be solved
is the slow case.

## When the scale is unknown

The quads are independent of scale, but the search is not: the step size and the size of
the catalogue window both come from the field of view. Measured on the test images, a
scale that is wrong by a factor of up to about 1.2 still solves nearly as often, but a
factor of 2 loses a quarter of them.

With no scale given (no `-fov`, and nothing in the header), arcsec tries a ladder of
scales, each $\sqrt 2$ times the last, from 0.25″ to 64″ per pixel. Any true scale is
within $2^{1/4} \approx 19\,\%$ of one of the rungs, inside that tolerance. Each rung
searches only the hint's position and the eight around it, nearest 1″/px first. If the
scale it solves at is more than 5 % from the one it searched with, the image is solved
once more at the solved scale, so the final check uses a catalogue window of the right
size, and arcsec prints the same "scale was inaccurate" warning ASTAP does.
`--fov-search` uses the same ladder, from a quarter to four times a given scale, when a
search at the given scale fails.

Next: [Fitting the plate](../fitting/).
