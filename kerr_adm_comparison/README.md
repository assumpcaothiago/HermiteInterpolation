# Kerr ADM Hermite–Lagrange comparison

This standalone test compares two ways of obtaining ten Cartesian ADM fields
and their first derivatives at one point:

- `hermite3d` returns values and analytic gradients of its value-only quintic
  interpolating polynomial from a `6 x 6 x 6` stencil.
- The copied production Lagrange routine returns values from a `7 x 7 x 7`
  stencil. Its gradients are then formed by sampling the Lagrange interpolant
  at four displaced points per Cartesian direction and applying a fourth-order
  centered finite difference.

The fields are ordered as

```text
alpha, betax, betay, betaz,
gammaxx, gammaxy, gammaxz, gammayy, gammayz, gammazz
```

All formulas use fixed parameters `M=1`, `a=1/2`, and LES throat radius

```text
s = (2 + sqrt(3))/8.
```

## Build and run

Supply one or more strictly increasing even active resolutions and a positive
number of random points:

```sh
make -C kerr_adm_comparison
make -C kerr_adm_comparison run \
    ARGS='--resolutions 32,48,64,96 --points 10000'
```

The optional seed accepts decimal or hexadecimal syntax:

```sh
make -C kerr_adm_comparison run \
    ARGS='--resolutions 64,96 --points 20000 --seed 0x12345678'
```

Grid population uses OpenMP. The interpolation loops are intentionally serial
because their reported wall times estimate latency for obtaining all ten ADM
values and gradients at one point. Control grid initialization with the
standard OpenMP environment variables, for example:

```sh
OMP_NUM_THREADS=8 OMP_DYNAMIC=FALSE OMP_PROC_BIND=spread OMP_PLACES=cores \
make -C kerr_adm_comparison run \
    ARGS='--resolutions 64,96 --points 10000'
```

## Geometry and padding

The active cell-centered cube is `[-5s,5s]^3`, with spacing `H=10s/N`.
Random centers are uniform on its exterior part `r>=s` and are reused at every
resolution. Six analytic layers are stored on every face:

```text
stored dimension = N + 12
memory = 10 (N + 12)^3 sizeof(double).
```

Two layers accommodate the `+/-2H` derivative samples and three more
accommodate the radius of the Lagrange stencil. A sixth layer covers the
binary64 case where a point immediately below the upper active face has its
displaced coordinate round onto the outward midpoint and the production
nearest-node rule selects the outer node. All stored samples are
evaluated analytically, including the inner sheet `0<r<s` and the padding
outside the active cube. Even `N` prevents a cell center from coinciding with
the puncture `r=0`.

The storage remains substantial at high resolution: `N=800` requires about
39.9 GiB for the ten field arrays alone. OpenMP accelerates their population
but does not reduce this memory requirement.

The lapse is the nonnegative LES lapse on both sheets. It is continuous but
not differentiable at `r=s`. Centers are never filtered merely because an
interpolation or finite-difference stencil enters the inner sheet; the program
reports how many do. Such stencils can exhibit reduced or erratic convergence,
which is a property of the sampled lapse rather than a hidden test failure.

## Errors and timing

For every field and for values, `d/dx`, `d/dy`, and `d/dz`, the report contains
sampled absolute RMS and maximum errors, scale-normalized versions, and
pairwise orders. Relative RMS is

```text
sqrt(sum(error^2) / sum(exact^2)),
```

and relative maximum error is `max(abs(error))/max(abs(exact))`. They are
undefined for identically zero reference channels. An `ANALYTIC-FD` baseline
applies the same difference stencil directly to exact values, isolating its
truncation error from Lagrange interpolation error.

The finite-difference spacing is always the local grid spacing `H`:

```text
(F(x-2H) - 8F(x-H) + 8F(x+H) - F(x+2H)) / (12H).
```

The latency table excludes allocation, grid construction, analytic-reference
evaluation, norm accumulation, and printing. It reports one Hermite
value-and-gradient call, one center-only Lagrange value call, and the complete
13-call Lagrange value-plus-FD-gradient path. Expected smooth-region orders are
five for Hermite values, seven for Lagrange values, and four for both gradient
paths.

## Generated and copied sources

`generate_kerr_adm_exact.py` uses exact SymPy arithmetic and emits the committed
C evaluator and high-precision fixtures. Ordinary builds require no Python or
SymPy. Regenerate or verify them with:

```sh
make -C kerr_adm_comparison regen
make -C kerr_adm_comparison check-generated
```

The Lagrange source, helper header, and SIMD header under
`vendor/bhah_lagrange/` are an unchanged snapshot of the read-only reference
tree. Their commit and hashes are recorded in `VENDOR.md`. Check the committed
snapshot or compare it with a currently available reference tree using:

```sh
make -C kerr_adm_comparison check-vendor
make -C kerr_adm_comparison check-upstream
```

The second command accepts an alternative location through
`BHAH_REFERENCE=/path/to/bhah_lib`. The build never modifies that directory.

Run all symbolic, library, evaluator, copied-Lagrange, and end-to-end checks
with:

```sh
make -C kerr_adm_comparison check
```
