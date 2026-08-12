#ifndef HERMITE3D_H
#define HERMITE3D_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Status returned by every hermite3d interpolation routine. */
typedef enum hermite3d_status {
  HERMITE3D_SUCCESS = 0,
  HERMITE3D_INVALID_ARGUMENT,
  HERMITE3D_INVALID_GRID,
  HERMITE3D_INDEX_OVERFLOW,
  HERMITE3D_STENCIL_UNAVAILABLE
} hermite3d_status;

/**
 * Geometry and memory layout of a uniform, logically Cartesian grid.
 *
 * xx0_start, xx1_start, and xx2_start are the coordinates of logical sample
 * (0, 0, 0).  The strides are measured in double elements, not bytes, so the
 * sample (i0, i1, i2) is found at
 *
 *   i0 * stride0 + i1 * stride1 + i2 * stride2.
 *
 * Every function passed to an interpolation call must use this same layout.
 */
typedef struct hermite3d_grid {
  size_t nxx0;
  size_t nxx1;
  size_t nxx2;

  double xx0_start;
  double xx1_start;
  double xx2_start;

  double dxx0;
  double dxx1;
  double dxx2;

  size_t stride0;
  size_t stride1;
  size_t stride2;
} hermite3d_grid;

/** Value and coordinate gradient of one interpolated function. */
typedef struct hermite3d_value_gradient {
  double value;
  double gradient[3];
} hermite3d_value_gradient;

/**
 * Componentwise-order-two jet of one interpolated function.
 *
 * derivative[q0][q1][q2] is
 *
 *   d^(q0+q1+q2) f / (dxx0^q0 dxx1^q1 dxx2^q2),
 *
 * where each q is independently 0, 1, or 2.  In particular,
 * derivative[0][0][0] is the value.  The array contains 27 entries, including
 * mixed derivatives whose total order is greater than two.
 */
typedef struct hermite3d_jet {
  double derivative[3][3][3];
} hermite3d_jet;

/**
 * Interpolate values and gradients of several functions at one point.
 *
 * functions is an array of num_functions pointers.  Each pointer addresses
 * logical grid sample (0, 0, 0) of one scalar function.  results contains one
 * output structure per input function.  The weights are computed once and
 * reused for the entire function array.
 *
 * On failure, results is left unchanged.
 */
hermite3d_status hermite3d_interpolate_value_gradient(
    const hermite3d_grid *grid, double xx0, double xx1, double xx2,
    size_t num_functions, const double *const functions[],
    hermite3d_value_gradient results[]);

/**
 * Interpolate the complete componentwise-order-two jet of several functions.
 *
 * The input layout and error behavior are the same as for
 * hermite3d_interpolate_value_gradient.  All returned entries are analytic
 * derivatives of the same tensor-product interpolating polynomial.
 */
hermite3d_status hermite3d_interpolate_jet(
    const hermite3d_grid *grid, double xx0, double xx1, double xx2,
    size_t num_functions, const double *const functions[],
    hermite3d_jet results[]);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* HERMITE3D_H */
