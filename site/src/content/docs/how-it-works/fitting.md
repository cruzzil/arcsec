---
title: Fitting the plate
description: The six plate constants, how arcsec solves for them by least squares, how they become the WCS that N.I.N.A. and other programs read, and how lens distortion is handled.
---

The agreeing quad matches say *roughly* how the image sits on the sky. This step makes it
exact: it finds the mapping from pixels to the sky that best fits all the matched stars at
once.

## The six plate constants

On the tangent plane (see [Searching the sky](../search/#flattening-the-sky)), an ideal
camera's mapping from pixel $(x, y)$ to standard coordinates $(\xi, \eta)$ is **affine**:

$$
\begin{aligned}
\xi  &= a\,x + b\,y + c \\
\eta &= d\,x + e\,y + f
\end{aligned}
$$

These six numbers are the **plate constants**. $c$ and $f$ are where pixel $(0, 0)$
lands; the matrix

$$
\begin{pmatrix} a & b \\ d & e \end{pmatrix}
$$

carries the scale, the rotation and the mirroring:

- the pixel scale, in arcseconds per pixel, is $\sqrt{a^2 + b^2} = \sqrt{d^2 + e^2}$;
- the rotation is $\operatorname{atan2}(d,\ e)$;
- the image is mirrored if the determinant $ae - bd$ is negative.

Each matched pair of an image position and a catalogue position gives one equation for
$\xi$ and one for $\eta$, so three pairs fix all six constants. There are usually many
more, so the constants are chosen to make the total squared miss as small as possible:

$$
\min_{a,b,c} \sum_i \left(a x_i + b y_i + c - \xi_i\right)^2, \qquad
\min_{d,e,f} \sum_i \left(d x_i + e y_i + f - \eta_i\right)^2
$$

These are two independent least-squares problems, one per axis, built from the same
pixel positions.

## Solving the least squares

Written as matrices, each problem is $A\mathbf{p} \approx \mathbf{b}$, with one row
$(x_i,\ y_i,\ 1)$ of $A$ per star. The textbook answer, the normal equations
$A^\mathsf{T}A\,\mathbf{p} = A^\mathsf{T}\mathbf{b}$, squares the problem's sensitivity to
rounding error. arcsec, like ASTAP, uses **Givens rotations** instead (the method in
Montenbruck and Pfleger's *Astronomy on the Personal Computer*): each new row is rotated
into a small upper-triangular matrix $R$ by a sequence of plane rotations

$$
\begin{pmatrix} c & s \\ -s & c \end{pmatrix}, \qquad
c = \frac{R_{jj}}{h},\quad s = \frac{a_j}{h},\quad h = \sqrt{R_{jj}^2 + a_j^2}
$$

each of which zeroes one entry of the new row. When all rows are in, $R\mathbf{p} =
\mathbf{z}$ is solved by back-substitution. It is numerically stable, needs only a 3 × 3
matrix however many stars there are, and never forms $A$ at all.

## From quads to stars

The first fit uses the quads: the centre of each matched image quad paired with the
centre of its catalogue quad. A quad's five ratios say that two quads have the same shape,
but not which corner is which, so centres are all that can be paired at this stage. A
few pairs that do not fit the rest (shown by residuals much larger than the others) are
removed by **sigma clipping**: fit, drop pairs more than three standard deviations out,
fit again, until nothing changes.

That fit is then used to pair individual stars, which is far more precise. Every catalogue
star is projected into the image through the fit, by inverting it,

$$
\begin{pmatrix} x \\ y \end{pmatrix} =
\begin{pmatrix} a & b \\ d & e \end{pmatrix}^{-1}
\begin{pmatrix} \xi - c \\ \eta - f \end{pmatrix}
$$

and paired with the nearest detected star within 6 pixels, each detection used at most
once. The plate is fitted again on those pairs, and the whole thing repeated within 3
and then 2 pixels. The final fit rests on hundreds of individual stars rather than tens
of quad centres; the same pairs are what [checking the answer](../verification/) counts.

### A sanity check on the shape

A camera pixel is square, so the matrix should be a pure **similarity**: rotation, one
scale, and perhaps a mirror, with no stretch or shear. Its two singular values (the
largest and smallest factors by which it stretches any direction) should then be equal.
For the matrix above they are $q + r$ and $|q - r|$, with

$$
q = \tfrac12\sqrt{(a+e)^2 + (d-b)^2}, \qquad
r = \tfrac12\sqrt{(a-e)^2 + (d+b)^2}
$$

and any fit with

$$
\frac{\sigma_\text{max}}{\sigma_\text{min}} = \frac{q + r}{|q - r|} > 1.08
$$

is refused. A set of chance matches can be fitted by *some* affine map, but seldom by
one that keeps squares square. On the test images no correct solution has a ratio above
1.03.

## Moving the tangent point

The fit is made in the plane that touches the sky at the trial position, which may be
some way from the image's centre. The tangent plane is only exact near its tangent point,
so after a position verifies, arcsec moves the tangent point to the centre of the image,
projects the catalogue again, and refits. Without this, with a hint a third of a field
off, most test images came out several arcseconds to many arcminutes wrong at the
corners, though the fit looked good in its own plane.

## Writing the WCS

FITS describes the same mapping as a reference pixel `CRPIX`, its sky position `CRVAL`,
and a matrix `CD` in degrees per pixel. With the tangent point at the image centre,
`CRPIX` is the centre pixel, `CRVAL` the sky position there, and

$$
\begin{pmatrix} \text{CD1\_1} & \text{CD1\_2} \\ \text{CD2\_1} & \text{CD2\_2} \end{pmatrix}
= \frac{1}{3600}
\begin{pmatrix} -a & -b \\ d & e \end{pmatrix}
$$

(the minus signs because FITS has right ascension increasing to the left). The older
keywords that many programs still read are derived from it as ASTAP derives them: with
$f = -1$ for the usual handedness ($\det \text{CD} < 0$) and $+1$ for a mirrored image,

$$
\text{CDELT1} = f\sqrt{\text{CD1\_1}^2 + \text{CD1\_2}^2}, \qquad
\text{CDELT2} = \sqrt{\text{CD2\_1}^2 + \text{CD2\_2}^2}
$$

$$
\text{CROTA2} = \operatorname{atan2}\left(f \cdot \text{CD2\_1},\ \text{CD2\_2}\right)
$$

If the image was binned for the solve, `CRPIX` and the scale terms are converted back to
the original pixels. These go into the `.wcs` and `.ini` files that N.I.N.A. and others
read (see [the command-line reference](../../reference/cli/)).

## Lens distortion

A real lens is not a perfect projection: towards the edges of the field the stars drift
from where an affine map puts them, a little for a telescope and a lot for a camera lens.
If that drift exceeds the 2-pixel matching radius, the star-level fit keeps only the part
of the field where the linear map happens to fit, and is accurate there alone.

So once a position verifies, arcsec fits a **polynomial** plate as well, in the style of
astrometry.net's `tweak`. Each standard coordinate becomes a polynomial in the scaled
pixel offsets $u, v \in [-1, 1]$ from the centre,

$$
\xi = \sum_{p + q \le 3} A_{pq}\,u^p v^q, \qquad
\eta = \sum_{p + q \le 3} B_{pq}\,u^p v^q
$$

ten terms per axis at third order. The catalogue is read again for the whole frame,
stars are re-paired within a wider radius under the current model, and the model is
refitted, up to ten times, until it stops moving. A higher order is used only if it
earns its extra terms, by an **F-test** against the best lower order:

$$
F = \frac{\left(\text{RSS}_\text{lower} - \text{RSS}_\text{higher}\right) / \Delta k}
         {\text{RSS}_\text{higher} / (n - k_\text{higher})} \;\ge\; 4
$$

(RSS is the sum of squared residuals, $k$ the number of terms, $n$ the number of stars),
and only if its stars cover the frame, so it is not extrapolating into empty corners.

The WCS that arcsec writes stays linear, as ASTAP's is, because that is what N.I.N.A. and
most other programs read. When the model shows real distortion, the linear plate reported
is the one closest to the model **over the whole frame** (a least-squares fit to it on a
9 × 9 grid), not the one that is exact at the centre and twice as far out at the corners.
With `--sip`, the distortion itself is written as SIP polynomial terms too, for programs
that read them.

Next: [Checking the answer](../verification/).
