# Kerr geodesic interpolation test

This package compares the accumulated dynamical effects of exact, Hermite,
and Lagrange-interpolated Kerr ADM data. It evolves the canonical variables

```text
(x, y, z, px, py, pz)
```

with the coordinate-time Hamilton equations derived in
[`../geodesic_notes/`](../geodesic_notes/). The fixed spacetime has `M=1` and
`a=1/2`. The particle starts on the stable prograde equatorial circular orbit
at Boyer--Lindquist radius `R_BL=5`, well outside both the horizon and the
prograde ISCO.

## Three ADM backends

- `EXACT` evaluates the generated analytic Cartesian ADM fields and gradients.
- `HERMITE` obtains all ten fields and their gradients from one tensor-product
  quintic interpolation call. Its gradients differentiate the same globally
  `C2` piecewise interpolant used for the values.
- `LAGRANGE+FD` uses the copied seven-point-per-axis Lagrange routine for field
  values. Each gradient uses twelve additional Lagrange evaluations and the
  fourth-order centered formula with displacement equal to the grid spacing
  `H`.

The Lagrange construction intentionally matches the existing ADM comparison.
It is not generally the derivative of the center's local value polynomial
because displaced evaluations may select different stencils.

## Build and run

Choose one or more strictly increasing even resolutions of at least 16:

```sh
make -C kerr_geodesic
OMP_NUM_THREADS=8 make -C kerr_geodesic run \
  ARGS='--resolutions 64,96,128'
```

The default evolution covers ten orbital periods with 4096 classical RK4
steps per period. Both values are configurable:

```sh
make -C kerr_geodesic run \
  ARGS='--resolutions 96 --orbits 5 --steps-per-orbit 8192'
```

Grid initialization uses OpenMP; each trajectory is integrated serially.
`OMP_NUM_THREADS`, `OMP_DYNAMIC`, `OMP_PROC_BIND`, and `OMP_PLACES` therefore
control grid-population throughput without changing the ODE calculation.

## Grid and noise

The active cell-centered cube is `[-5,5]^3`, with `H=10/N`. Six analytic ghost
layers on each face supply every Hermite stencil and every Lagrange stencil at
the `+/-2H` finite-difference samples. The allocation is

```text
10 (N + 12)^3 sizeof(double).
```

Even resolutions prevent a stored point from coinciding with the puncture.
All trajectory stages must remain in the active exterior region `r>s`; the
program aborts rather than clamping, extrapolating, or hiding a failed orbit.

Optional deterministic multiplicative noise uses the same model as the ADM
pointwise comparison:

```sh
make -C kerr_geodesic run \
  ARGS='--resolutions 64,96 --noise-epsilon 1e-8'
```

Both interpolation methods receive identical noisy arrays. The exact backend,
closed circular orbit, and invariant references remain noise-free. Select an
independent realization with `--noise-seed UINT64`.

## Diagnostics

The report separates three effects:

1. The exact-field RK4 trajectory is compared with the closed circular orbit,
   directly measuring time-integration error.
2. On 256 fixed points around the closed orbit, interpolated coordinate
   velocities and momentum forces are compared with the exact RHS. This is a
   local force-noise diagnostic without trajectory feedback.
3. Hermite and Lagrange trajectories are compared with the same-step exact
   RK4 trajectory. The report gives RMS and maximum position, momentum,
   radial, vertical, and unwrapped phase errors with pairwise grid orders.

It also prints drift in the exact Kerr energy, each method's own Hamiltonian,
`Lz=x*py-y*px`, and the exact-metric mass-shell residual formed using the fixed
initial `pt=-E0`. Conserving a method Hamiltonian is not the same as following
the exact trajectory: Hermite is internally Hamiltonian, whereas the chosen
Lagrange-plus-FD force is generally not.

If the exact RK4 position or momentum error exceeds one tenth of the smaller
finest-grid interpolation error, the program warns that
`--steps-per-orbit` should be increased.

## CSV and plots

Write sampled long-form data with:

```sh
make -C kerr_geodesic run \
  ARGS='--resolutions 64,96 --trajectory-csv build/trajectory.csv'
```

`--output-every` selects the number of RK4 steps between rows and defaults to
16. The CSV includes all states, the closed and exact-RK4 references, radial
and phase errors, exact and method energies, `Lz`, and the mass-shell residual.

The optional Matplotlib plot shows the equatorial path, radial drift,
unwrapped phase error, and exact-invariant drift. It selects the finest CSV
resolution unless `--resolution` is supplied:

```sh
make -C kerr_geodesic plot \
  ARGS='--resolutions 64 --orbits 2 --steps-per-orbit 4096'

python3 kerr_geodesic/plot_trajectory.py \
  kerr_geodesic/build/trajectory.csv --resolution 64
```

## Verification

```sh
make -C kerr_geodesic check
make -C kerr_geodesic sanitize
```

The checks cover the Hamiltonian RHS, circular data, analytic momentum,
invariants, phase unwrapping, fourth-order RK4 convergence, clean/noisy
determinism, OpenMP independence, CSV structure, invalid inputs, and all
existing interpolation and symbolic dependencies.
