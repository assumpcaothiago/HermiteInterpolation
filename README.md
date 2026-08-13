# Quintic Hermite Interpolation

This repository contains a standalone C99 library and self-contained
mathematical notes for tensor-product quintic Hermite interpolation. The main
implementation is the **value-only** construction: callers provide samples on
a uniform, logically Cartesian grid, while the derivative reconstruction is
incorporated directly into the interpolation weights.

In three dimensions, one query uses a separable `6 x 6 x 6` stencil. The same
coordinate-dependent weights are reused when several scalar functions are
interpolated at that point.

## Repository contents

- [`hermite3d/`](hermite3d/) contains the standalone interpolation library,
  its public C API, tests, and detailed usage documentation.
- [`notes/`](notes/) contains the LaTeX article, exact SymPy derivations, and
  the generated formula appendix.
- [`kerr_notes/`](kerr_notes/) contains checked, axis-regular formulas for the
  subextremal Kerr metric in modified quasi-isotropic Cartesian coordinates.
- [`kerr_interpolation/`](kerr_interpolation/) contains a generated analytic
  Kerr evaluator and a C convergence experiment for metric values and
  Cartesian gradients.
- [`kerr_adm_comparison/`](kerr_adm_comparison/) compares Hermite ADM values
  and direct gradients with the copied seven-point Lagrange interpolator and
  fourth-order finite differences of its interpolated values.
- [`references/`](references/) contains the source material used when deriving
  and checking the mathematics.
- [`mathematica_notebooks/`](mathematica_notebooks/) contains earlier symbolic
  exploration notes.

The library has no dependency on the reference documents or symbolic tools.
They are retained to document and independently verify its mathematical
foundation.

## Library capabilities

`hermite3d` supports:

- Uniform three-dimensional grids with independent positive spacings
  `dxx0`, `dxx1`, and `dxx2`.
- Shifted grid origins and caller-defined element strides.
- Any positive number of scalar functions evaluated at the same point.
- An optimized value-and-gradient interface.
- A full `3 x 3 x 3` derivative jet for testing and analysis, containing
  every derivative
  `d^(q0+q1+q2)/(dxx0^q0 dxx1^q1 dxx2^q2)` with each `qa` in `{0,1,2}`.
- C99 implementation and a C-compatible API that can also be called from C++.
- No internal allocation and no mutable global state.

All returned derivatives come from the same piecewise polynomial as the
interpolated value. Neighboring cell polynomials agree through second normal
derivatives at shared grid planes.

See the [library README](hermite3d/README.md) for the grid descriptor, memory
layout, valid interpolation domain, result indexing, and example calls.

## Build and test the C library

From the repository root, compile the static library and run the complete test
suite with:

```sh
make -C hermite3d check
```

This builds:

- `hermite3d/build/libhermite3d.a`
- The dependency-free C mathematical and API test program.
- A C++ header and linkage smoke test.

The tests cover tensor-polynomial reproduction, all 27 jet entries,
value-gradient consistency, multiple-function calls, nontrivial strides,
reflection symmetry, endpoint finite-difference data, inter-cell continuity,
domain boundaries, and error handling.

To build only the static library or remove build products:

```sh
make -C hermite3d
make -C hermite3d clean
```

The library build requires a C99 compiler, an archiver, `make`, and the
standard mathematics library. The compatibility smoke test additionally
requires a C++11 compiler.

## Build and verify the mathematical notes

The notes distinguish the exact endpoint-jet interpolant from the uniform-grid
value-only construction used by the C library. They derive the one-dimensional
basis, six direct grid weights, physical derivative scaling, error terms,
continuity properties, and multidimensional tensor products.

Run the exact symbolic checks and compile the PDF with:

```sh
make -C notes check
```

This requires Python with SymPy, `latexmk`, and a suitable LaTeX installation.
The command verifies the deterministic generated formula fragment before
building the article.

Useful individual targets include:

```sh
make -C notes formulas  # regenerate only the formula appendix
make -C notes pdf       # build the article
make -C notes clean     # remove auxiliary LaTeX files
```

The exact symbolic verification can also be run without compiling LaTeX:

```sh
python3 notes/derive_hermite.py --check
```

The independent Kerr notes and their symbolic checks can be built with:

```sh
make -C kerr_notes check
```

## Compare Hermite and Lagrange ADM interpolation

The ADM comparison samples the lapse, contravariant shift, and six independent
spatial-metric components on the same cell-centered Kerr grids. Hermite returns
values and gradients in one call; Lagrange supplies values at the center and at
twelve displaced points used by a fourth-order centered difference with step
equal to the grid spacing. The comparison uses the signed LES lapse, which is
smooth through the throat and reconstructs the same four-metric as the
nonnegative convention.

```sh
make -C kerr_adm_comparison check
make -C kerr_adm_comparison run \
    ARGS='--resolutions 32,48,64,96 --points 10000'
```

See the [ADM comparison README](kerr_adm_comparison/README.md) for field order,
analytic padding, error definitions, timing boundaries, and copied-source
provenance.

To compare sensitivity to deterministic multiplicative grid noise while
retaining clean analytic references:

```sh
OMP_NUM_THREADS=8 make -C kerr_adm_comparison run \
    ARGS='--resolutions 200,400 --points 10000 --noise-epsilon 1e-8'
```

## Run the Kerr convergence experiment

The Kerr experiment initializes all ten independent four-metric components on
cell-centered grids and interpolates them together at one deterministic random
point cloud restricted to the LES exterior sheet `r>=s`. Build it and select
the resolutions and point count at run time:

```sh
make -C kerr_interpolation
make -C kerr_interpolation run \
    ARGS='--resolutions 32,48,64,96 --points 10000'
```

Grid initialization and random-point interpolation use OpenMP. For example:

```sh
OMP_NUM_THREADS=8 \
OMP_DYNAMIC=FALSE \
OMP_PROC_BIND=spread \
OMP_PLACES=cores \
make -C kerr_interpolation run \
    ARGS='--resolutions 64,96 --points 10000'
```

The thread count affects execution time but not the grid allocation. At
`N=800`, the ten stored fields still require about 39.0 GiB.

An optional z-axis diagnostic overlays any metric component or first
derivative with its analytical expression:

```sh
make -C kerr_interpolation plot \
    PROFILE=xx:value \
    ARGS='--resolutions 32,48,64,96 --points 10000'
```

Run its generated-expression, exact-evaluator, library, symbolic, and
end-to-end checks with:

```sh
make -C kerr_interpolation check
```

See the [experiment README](kerr_interpolation/README.md) for its exact metric
API, cell-centered ghost grid, deterministic random cloud, and interpretation
of the sampled exterior-domain error norms.

## Numerical scope

The value-only interpolant assumes a uniform grid along each logical
coordinate and requires the complete six-point stencil in every direction.
It does not clamp queries, extrapolate, select ghost-zone policies, or provide
one-sided boundary closures. Those responsibilities remain with the calling
application.

For sufficiently smooth data, the construction is generically fifth-order for
interpolated values, fourth-order for first derivatives, and third-order for
second derivatives taken twice in the same coordinate; mixed first-first
derivatives remain fourth-order. The weights are not positivity preserving, so
the method is intended for smooth interpolation rather than monotonicity- or
shape-constrained reconstruction.
