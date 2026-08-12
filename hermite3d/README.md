# hermite3d

`hermite3d` is a small C99 library for value-only, tensor-product quintic
Hermite interpolation on a uniform, logically Cartesian three-dimensional
grid. A single call evaluates any positive number of scalar functions at the
same point, reusing the coordinate-dependent interpolation weights.

The library has no runtime dependencies, performs no allocation, and keeps no
mutable global state. Its public header can be included from both C and C++.

## Grid and storage

A `hermite3d_grid` describes the logical grid dimensions, the coordinate of
sample `(0,0,0)`, the three positive spacings, and the element strides:

```c
hermite3d_grid grid = {
    .nxx0 = nxx0, .nxx1 = nxx1, .nxx2 = nxx2,
    .xx0_start = xx0_start,
    .xx1_start = xx1_start,
    .xx2_start = xx2_start,
    .dxx0 = dxx0, .dxx1 = dxx1, .dxx2 = dxx2,
    .stride0 = stride0,
    .stride1 = stride1,
    .stride2 = stride2
};
```

An input function is a `const double *` pointing to its logical sample
`(0,0,0)`. Sample `(i0,i1,i2)` is read at

```text
i0*stride0 + i1*stride1 + i2*stride2
```

where strides are measured in `double` elements, not bytes. All functions in
one call share the same grid and strides. For an ordinary contiguous layout
with `xx0` varying fastest, use `stride0 = 1`, `stride1 = nxx0`, and
`stride2 = nxx0*nxx1`. Padded and permuted layouts can be represented by other
positive, non-overlapping strides.

## Interpolation domain

For a query coordinate, define

```text
u0 = (xx0 - xx0_start)/dxx0
u1 = (xx1 - xx1_start)/dxx1
u2 = (xx2 - xx2_start)/dxx2
```

Each coordinate must satisfy `2 <= ua <= nxxa - 3`. The library locates the
containing cell and uses the six samples with offsets `-2` through `+3` from
its lower endpoint in each direction, for a total stencil of 216 samples. The
upper boundary is included and is evaluated as the upper endpoint of the
preceding cell. Queries outside this domain are rejected; the library does not
clamp, extrapolate, or choose a boundary closure.

The value-only construction estimates the endpoint derivative data from these
same grid samples and incorporates that reconstruction into six direct weights
per direction. No caller-supplied derivative arrays are needed.

## Calling the library

Include `hermite3d.h` and pass an array of input field pointers together with a
matching output array:

```c
const double *functions[] = {field_a, field_b};
hermite3d_value_gradient result[2];

hermite3d_status status = hermite3d_interpolate_value_gradient(
    &grid, xx0, xx1, xx2, 2, functions, result);
```

For each function, this routine returns its value and
`gradient[0]`, `gradient[1]`, and `gradient[2]`, the derivatives with respect
to `xx0`, `xx1`, and `xx2`.

The testing and analysis interface returns the complete tensor-product jet:

```c
hermite3d_jet jet[2];
status = hermite3d_interpolate_jet(
    &grid, xx0, xx1, xx2, 2, functions, jet);
```

`jet[f].derivative[q0][q1][q2]` is

```text
d^(q0+q1+q2) functions[f]
---------------------------------,  with q0,q1,q2 in {0,1,2}.
 dxx0^q0 dxx1^q1 dxx2^q2
```

Consequently, `[0][0][0]` is the value and entries `[1][0][0]`,
`[0][1][0]`, and `[0][0][1]` are the gradient. The 27 entries include mixed
derivatives whose total order can exceed two; each coordinate's derivative
order is at most two.

All outputs are derivatives of the same piecewise polynomial. In particular,
the interpolated field is continuous through second derivatives at shared grid
planes. Physical spacing factors are already included in returned derivatives.

On failure, a function returns a non-success `hermite3d_status` and leaves the
entire output array unchanged. The caller owns all input and output storage,
and those regions must not overlap.

## Build and test

From this directory:

```sh
make
make check
```

`make` creates `build/libhermite3d.a`. `make check` also builds and runs the
dependency-free C test suite and a C++ header/linkage smoke test. Build
products remain under `build/`; remove them with `make clean`.
