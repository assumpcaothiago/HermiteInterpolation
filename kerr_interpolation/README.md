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

## OpenMP parallel execution

The convergence executable uses OpenMP for both Cartesian grid initialization
and random-point interpolation. Set the team size with the standard OpenMP
environment variable:

```sh
OMP_NUM_THREADS=8 \
OMP_DYNAMIC=FALSE \
OMP_PROC_BIND=spread \
OMP_PLACES=cores \
./build/kerr_convergence --resolutions 64,96 --points 10000
```

`OMP_NUM_THREADS` is the usual control. `OMP_DYNAMIC=FALSE` requests the full
team, while `OMP_PROC_BIND=spread` and `OMP_PLACES=cores` can improve placement
on multi-core or multi-socket machines. If these variables are unset, the
OpenMP runtime chooses its normal defaults. The program reports the actual
grid-initialization and query team sizes for every resolution.

Grid rows are statically divided among threads while retaining contiguous
`xx0` traversal. This also provides parallel first-touch placement for the
large field allocation. Random queries use static scheduling and private norm
accumulators that are merged after the parallel region. Separate resolutions,
z-axis profile export, CSV output, and final reporting remain serial.

OpenMP reduces computation time but not storage. The ten fields require

```text
10 * (N+6)^3 * sizeof(double)
```

bytes. In particular, `N=800` requires about 39.0 GiB. A machine without
enough physical memory may swap heavily, and a sufficiently large run may be
limited by memory bandwidth or NUMA placement rather than arithmetic.

## Plot a variable on the z axis

The convergence executable can additionally sample one metric component or
first derivative on the rotation axis `x=y=0`. This diagnostic does not alter
the random cloud or its norms. Select a component and quantity with
`COMPONENT:QUANTITY`, where the components are

```text
tt tx ty tz xx xy xz yy yz zz
```

and the quantities are `value`, `dx`, `dy`, and `dz`. For example, export
`d g_xy/dz` at 2000 fixed physical positions for every requested resolution:

```sh
./build/kerr_convergence \
    --resolutions 32,48,64,96 \
    --points 10000 \
    --z-profile xy:dz \
    --z-samples 2000 \
    --z-output build/g_xy_dz.csv
```

The CSV columns are:

```text
component,quantity,resolution,z,analytic,interpolated,error
```

The profile sample count must be even. Half of its midpoint samples lie on
each exterior branch, `-5s<z<-s` and `s<z<5s`. Thus neither the lower LES
sheet `|z|<s` nor the horizon itself is used as a profile query.

Render the CSV with the accompanying Matplotlib script:

```sh
python3 plot_z_profile.py \
    build/g_xy_dz.csv \
    --output build/g_xy_dz.png
```

The upper panel overlays the analytical curve and one Hermite curve per grid
resolution. The lower panel shows absolute error, normally on a logarithmic
scale. The two exterior curves are not connected across the excluded interval
`|z|<s`, and vertical markers identify `z=0` and the throats `z=+/-s`.

For convenience, both export and plotting can be performed with one target:

```sh
make plot \
    PROFILE=xy:dz \
    Z_SAMPLES=2000 \
    ARGS='--resolutions 32,48,64,96 --points 10000'
```

Override `PROFILE_CSV` or `PROFILE_PNG` to choose different output paths. The
plot target requires Matplotlib; the C build and ordinary checks do not.

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

The rectangular grid still contains values on both LES exterior sheets. This
is necessary because the six-point Cartesian stencil of a query close to the
throat can cross the sphere `r=s`. Query selection, rather than grid storage,
is restricted to the chosen exterior sheet.

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

SplitMix64 first generates candidates uniformly in the open active cube and
rejects those with `r<s`. The accepted cloud is therefore uniform on the cube
portion `r>=s` and contains no queries on the lower LES sheet. The complete
accepted cloud is generated once, before any OpenMP region, and its immutable
coordinates are reused at every resolution. Thus thread scheduling cannot
change the sampled points.

Component order is:

```text
tt tx ty tz xx xy xz yy yz zz
```

For every component, the program reports values and derivatives with respect
to `x`, `y`, and `z`. Each table contains the sampled RMS error (labeled
sampled `L2`), sampled maximum error (sampled `Linf`), and adjacent-level
orders computed with the actual spacing ratio. The `tz` component and all of
its derivatives are exactly zero and therefore have undefined (`--`) orders.

These are finite-cloud sampled norms rather than exact continuum norms. On the
sampled exterior domain the metric components and Cartesian derivatives are
finite, including their smooth spatial extension to the throat. The expected
orders are five for values and four for first derivatives, but individual
two-level estimates can still oscillate with interpolation phase, and the
executable intentionally enforces no convergence threshold.

## Analytic evaluator and verification

`generate_kerr_exact.py` uses exact SymPy constants and the stable,
axis-regular Cartesian formulas documented in [`kerr_notes`](../kerr_notes/).
It differentiates all ten components symbolically and applies deterministic
common-subexpression elimination. It emits separate value-only and
value-gradient kernels: grid initialization does not pay to compute unused
derivatives.

The generated C source and high-precision fixture header are committed, so a
normal convergence build requires a C99 compiler with OpenMP support, `make`,
the math library, and the neighboring `hermite3d` source. The exact-evaluator
unit test and `hermite3d` themselves do not acquire an OpenMP dependency.
Regenerate or verify the symbolic output with:

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
