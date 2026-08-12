# Kerr metric interpolation convergence test

This directory contains a standalone convergence experiment for the
[`hermite3d`](../hermite3d/) value-only quintic Hermite library. It samples the
Kerr four-metric on uniform cell-centered Cartesian grids, interpolates all ten
independent components at the same randomly chosen points, and compares the
interpolated values and Cartesian gradients with analytic expressions.

The test fixes

```text
M = 1
a = 1/2
s = r_+/4 = (2 + sqrt(3))/8
```

where `s` is the LES coordinate radius of the horizon throat. The active cube
is `[-5s,5s]^3`, so its side length is five coordinate horizon diameters.

## Build and run

From this directory, build the convergence executable with:

```sh
make
```

Run it by supplying one or more active grid resolutions and a positive number
of random points:

```sh
./build/kerr_convergence \
    --resolutions 32,48,64,96 \
    --points 10000
```

Equivalently:

```sh
make run ARGS='--resolutions 32,48,64,96 --points 10000'
```

The optional seed accepts decimal or hexadecimal notation:

```sh
./build/kerr_convergence \
    --resolutions 64,96,128 \
    --points 20000 \
    --seed 0x123456789abcdef0
```

The default seed is `0x4b4552524845524d`. Resolution values must be distinct,
strictly increasing, and even. A single resolution is allowed and reports
errors without convergence orders.

## Cell-centered grid

For an active resolution `N`, the grid spacing is

```text
h = 10s/N.
```

The program stores `N+6` points in every direction. Three analytic ghost
layers on each side allow every point in the active cube to use the complete
six-point interpolation stencil. Stored coordinate `i` is

```text
x(i) = -5s + (i - 3 + 1/2) h,
```

and likewise for `y` and `z`. Because `N` is even, coordinates around the
origin are `...,-3h/2,-h/2,h/2,3h/2,...`; the singular puncture is never a
grid sample.

The ten functions use a structure-of-arrays allocation and the contiguous
logical strides

```text
stride0 = 1
stride1 = N+6
stride2 = (N+6)^2.
```

At each query, one call to `hermite3d_interpolate_value_gradient` interpolates
all ten functions, so the coordinate weights are computed once and reused.

## Random points and reported errors

SplitMix64 produces one uniform, unfiltered cloud in the open active cube. The
generator is reset to the same seed at each resolution, so every level uses
identical physical points. Its binary64 midpoint mapping cannot produce the
origin and does not reject, replace, or specialize any query.

Component order is:

```text
tt tx ty tz xx xy xz yy yz zz
```

For every component, the program reports values and derivatives with respect
to `x`, `y`, and `z`. Each table contains the sampled RMS error (labeled
sampled `L2`), sampled maximum error (sampled `Linf`), and adjacent-level
orders computed with the actual spacing ratio. The `tz` component and all of
its derivatives are exactly zero and therefore have undefined (`--`) orders.

These are finite-cloud sampled norms, not continuum norms on the whole cube.
The cube contains the excluded puncture in its continuum domain; spatial
metric components diverge roughly as `r^-4` and their gradients as `r^-5`.
Consequently a close random point can dominate the reported errors and cause
nonmonotone or negative apparent orders. The program intentionally retains
that behavior for later analysis. In a smooth neighborhood separated from the
puncture, the expected local orders are five for values and four for first
derivatives, but the executable enforces no convergence threshold.

## Analytic evaluator and verification

`generate_kerr_exact.py` uses exact SymPy constants and the stable,
axis-regular Cartesian formulas documented in [`kerr_notes`](../kerr_notes/).
It differentiates all ten components symbolically and applies deterministic
common-subexpression elimination. It emits separate value-only and
value-gradient kernels: grid initialization does not pay to compute unused
derivatives.

The generated C source and high-precision fixture header are committed, so a
normal build requires only a C99 compiler, `make`, the math library, and the
neighboring `hermite3d` source. Regenerate or verify them with:

```sh
make regen
make check-generated
```

Regeneration requires Python and SymPy. Run all symbolic, evaluator, library,
and end-to-end smoke checks with:

```sh
make check
```

The evaluator tests cover high-precision references on the inner and outer
coordinate sheets, both rotation axes, the throat, and points immediately on
both sides of it. They also compare the analytic gradient to five-point finite
differences of the separate value-only kernel and test component ordering,
parity, exact-zero channels, invalid inputs, nonfinite results, and unchanged
outputs on failure.

Remove reproducible build products without touching the generated sources:

```sh
make clean
```
