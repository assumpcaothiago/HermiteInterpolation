#include "hermite3d.h"

#include <math.h>
#include <stdint.h>

enum {
  HERMITE3D_STENCIL_WIDTH = 6,
  HERMITE3D_DERIVATIVE_COUNT = 3
};

typedef struct hermite3d_prepared_query {
  size_t base[3];
  double weight[3][HERMITE3D_DERIVATIVE_COUNT][HERMITE3D_STENCIL_WIDTH];
} hermite3d_prepared_query;

/*
 * Coefficients of the six one-dimensional value-only quintic weights.
 *
 * Row s corresponds to grid offset s - 2 relative to the lower cell endpoint.
 * Columns contain coefficients of t^5 through t^0.  These weights result from
 * inserting fourth-order centered estimates of the endpoint slopes and
 * curvatures into the endpoint-value/slope/curvature Hermite polynomial.  The
 * union of the two centered endpoint stencils is therefore {-2,...,+3}.
 */
static const double weight_coefficient[HERMITE3D_STENCIL_WIDTH][6] = {
    {-5.0 / 24.0, 13.0 / 24.0, -3.0 / 8.0, -1.0 / 24.0,
     1.0 / 12.0, 0.0},
    {25.0 / 24.0, -8.0 / 3.0, 13.0 / 8.0, 2.0 / 3.0,
     -2.0 / 3.0, 0.0},
    {-25.0 / 12.0, 21.0 / 4.0, -35.0 / 12.0, -5.0 / 4.0, 0.0,
     1.0},
    {25.0 / 12.0, -31.0 / 6.0, 11.0 / 4.0, 2.0 / 3.0, 2.0 / 3.0,
     0.0},
    {-25.0 / 24.0, 61.0 / 24.0, -11.0 / 8.0, -1.0 / 24.0,
     -1.0 / 12.0, 0.0},
    {5.0 / 24.0, -1.0 / 2.0, 7.0 / 24.0, 0.0, 0.0, 0.0}};

static int add_overflows_size_t(size_t left, size_t right) {
  return left > SIZE_MAX - right;
}

static int multiply_overflows_size_t(size_t left, size_t right) {
  return left != 0 && right > SIZE_MAX / left;
}

/*
 * Verify that every logical index can be converted to one element offset.
 * Positive arbitrary strides permit padding and any chosen fastest axis.  The
 * caller remains responsible for providing storage that covers the resulting
 * maximum offset and for choosing a nonoverlapping logical layout.
 */
static hermite3d_status validate_grid(const hermite3d_grid *grid) {
  size_t maximum_offset = 0;
  const size_t dimensions[3] = {grid->nxx0, grid->nxx1, grid->nxx2};
  const size_t strides[3] = {grid->stride0, grid->stride1, grid->stride2};
  const double starts[3] = {grid->xx0_start, grid->xx1_start,
                            grid->xx2_start};
  const double spacings[3] = {grid->dxx0, grid->dxx1, grid->dxx2};

  for (size_t axis = 0; axis < 3; ++axis) {
    size_t axis_offset;

    if (dimensions[axis] < HERMITE3D_STENCIL_WIDTH || strides[axis] == 0 ||
        !isfinite(starts[axis]) || !isfinite(spacings[axis]) ||
        spacings[axis] <= 0.0 || !isfinite(1.0 / spacings[axis])) {
      return HERMITE3D_INVALID_GRID;
    }

    if (multiply_overflows_size_t(dimensions[axis] - 1, strides[axis])) {
      return HERMITE3D_INDEX_OVERFLOW;
    }
    axis_offset = (dimensions[axis] - 1) * strides[axis];
    if (add_overflows_size_t(maximum_offset, axis_offset)) {
      return HERMITE3D_INDEX_OVERFLOW;
    }
    maximum_offset += axis_offset;
  }

  return HERMITE3D_SUCCESS;
}

/*
 * Locate one coordinate in the closed region covered by valid six-point
 * stencils.  At the upper endpoint, the cell on the left must be selected:
 * choosing floor(u) there would request a nonexistent point beyond +3.
 */
static hermite3d_status locate_axis(size_t dimension, double start,
                                    double spacing, double coordinate,
                                    size_t *cell, double *t) {
  const double lower_coordinate = start + 2.0 * spacing;
  const double upper_coordinate =
      start + (double)(dimension - 3) * spacing;
  const double upper = (double)(dimension - 3);
  double u;
  double lower_as_double;

  if (!isfinite(lower_coordinate) || !isfinite(upper_coordinate)) {
    return HERMITE3D_INDEX_OVERFLOW;
  }
  if (coordinate < lower_coordinate || coordinate > upper_coordinate) {
    return HERMITE3D_STENCIL_UNAVAILABLE;
  }

  /*
   * Test coordinate-space endpoints before normalizing.  For decimal starts
   * and spacings, an endpoint constructed as start + k * spacing can otherwise
   * normalize to 1.9999999999999998 or to upper plus one rounding unit.
   */
  if (coordinate == lower_coordinate) {
    *cell = 2;
    *t = 0.0;
    return HERMITE3D_SUCCESS;
  }
  if (coordinate == upper_coordinate) {
    *cell = dimension - 4;
    *t = 1.0;
    return HERMITE3D_SUCCESS;
  }

  u = (coordinate - start) / spacing;
  if (!isfinite(u)) {
    return HERMITE3D_INDEX_OVERFLOW;
  }

  /* The coordinate comparison proved the point is in range.  These two
   * assignments correct only roundoff introduced by the normalization. */
  if (u < 2.0)
    u = 2.0;
  if (u > upper)
    u = upper;
  if (u == upper) {
    *cell = dimension - 4;
    *t = 1.0;
    return HERMITE3D_SUCCESS;
  }

  lower_as_double = floor(u);
  if (lower_as_double < 0.0 || lower_as_double > (double)SIZE_MAX) {
    return HERMITE3D_INDEX_OVERFLOW;
  }
  *cell = (size_t)lower_as_double;
  *t = u - lower_as_double;
  return HERMITE3D_SUCCESS;
}

/*
 * Evaluate every weight and its first two coordinate derivatives by Horner's
 * rule.  Differentiating the coefficient form directly ensures that all
 * three tables describe one polynomial.  The inverse-spacing factors convert
 * derivatives with respect to t into derivatives with respect to xx.
 */
static void evaluate_axis_weights(
    double t, double spacing,
    double weight[HERMITE3D_DERIVATIVE_COUNT][HERMITE3D_STENCIL_WIDTH]) {
  const double inverse_spacing = 1.0 / spacing;
  const double inverse_spacing_squared = inverse_spacing * inverse_spacing;

  for (size_t stencil = 0; stencil < HERMITE3D_STENCIL_WIDTH; ++stencil) {
    const double *coefficient = weight_coefficient[stencil];

    weight[0][stencil] =
        (((((coefficient[0] * t + coefficient[1]) * t + coefficient[2]) * t +
            coefficient[3]) *
               t +
           coefficient[4]) *
              t +
          coefficient[5]);

    weight[1][stencil] =
        ((((5.0 * coefficient[0] * t + 4.0 * coefficient[1]) * t +
           3.0 * coefficient[2]) *
              t +
          2.0 * coefficient[3]) *
             t +
         coefficient[4]) *
        inverse_spacing;

    weight[2][stencil] =
        (((20.0 * coefficient[0] * t + 12.0 * coefficient[1]) * t +
          6.0 * coefficient[2]) *
             t +
         2.0 * coefficient[3]) *
        inverse_spacing_squared;
  }
}

static hermite3d_status prepare_query(const hermite3d_grid *grid, double xx0,
                                      double xx1, double xx2,
                                      hermite3d_prepared_query *query) {
  const double coordinates[3] = {xx0, xx1, xx2};
  const size_t dimensions[3] = {grid->nxx0, grid->nxx1, grid->nxx2};
  const double starts[3] = {grid->xx0_start, grid->xx1_start,
                            grid->xx2_start};
  const double spacings[3] = {grid->dxx0, grid->dxx1, grid->dxx2};
  double t[3];

  for (size_t axis = 0; axis < 3; ++axis) {
    size_t cell;
    hermite3d_status status = locate_axis(dimensions[axis], starts[axis],
                                          spacings[axis], coordinates[axis],
                                          &cell, &t[axis]);
    if (status != HERMITE3D_SUCCESS) {
      return status;
    }

    /* cell >= 2 follows from locate_axis, so this subtraction is safe. */
    query->base[axis] = cell - 2;
    evaluate_axis_weights(t[axis], spacings[axis], query->weight[axis]);
  }

  return HERMITE3D_SUCCESS;
}

static hermite3d_status validate_and_prepare(
    const hermite3d_grid *grid, double xx0, double xx1, double xx2,
    size_t num_functions, const double *const functions[], const void *results,
    hermite3d_prepared_query *query) {
  hermite3d_status status;

  if (grid == NULL || functions == NULL || results == NULL ||
      num_functions == 0 || !isfinite(xx0) || !isfinite(xx1) ||
      !isfinite(xx2)) {
    return HERMITE3D_INVALID_ARGUMENT;
  }
  for (size_t function = 0; function < num_functions; ++function) {
    if (functions[function] == NULL) {
      return HERMITE3D_INVALID_ARGUMENT;
    }
  }

  status = validate_grid(grid);
  if (status != HERMITE3D_SUCCESS) {
    return status;
  }
  return prepare_query(grid, xx0, xx1, xx2, query);
}

hermite3d_status hermite3d_interpolate_value_gradient(
    const hermite3d_grid *grid, double xx0, double xx1, double xx2,
    size_t num_functions, const double *const functions[],
    hermite3d_value_gradient results[]) {
  hermite3d_prepared_query query;
  hermite3d_status status =
      validate_and_prepare(grid, xx0, xx1, xx2, num_functions, functions,
                           results, &query);

  if (status != HERMITE3D_SUCCESS) {
    return status;
  }

  for (size_t function = 0; function < num_functions; ++function) {
    double after_xx0[2][HERMITE3D_STENCIL_WIDTH]
                    [HERMITE3D_STENCIL_WIDTH] = {{{0.0}}};
    double after_xx1[3][HERMITE3D_STENCIL_WIDTH] = {{0.0}};
    hermite3d_value_gradient result = {0.0, {0.0, 0.0, 0.0}};
    const double *data = functions[function];

    /*
     * Contracting one axis at a time is algebraically identical to applying
     * four precomputed 216-point tensor weights.  It needs fewer operations
     * and makes the separability of the interpolant explicit.
     */
    for (size_t s2 = 0; s2 < HERMITE3D_STENCIL_WIDTH; ++s2) {
      const size_t offset2 = (query.base[2] + s2) * grid->stride2;
      for (size_t s1 = 0; s1 < HERMITE3D_STENCIL_WIDTH; ++s1) {
        const size_t offset1 = (query.base[1] + s1) * grid->stride1;
        double value_sum = 0.0;
        double derivative0_sum = 0.0;

        for (size_t s0 = 0; s0 < HERMITE3D_STENCIL_WIDTH; ++s0) {
          const size_t offset = offset2 + offset1 +
                                (query.base[0] + s0) * grid->stride0;
          const double sample = data[offset];
          value_sum += query.weight[0][0][s0] * sample;
          derivative0_sum += query.weight[0][1][s0] * sample;
        }
        after_xx0[0][s1][s2] = value_sum;
        after_xx0[1][s1][s2] = derivative0_sum;
      }
    }

    for (size_t s2 = 0; s2 < HERMITE3D_STENCIL_WIDTH; ++s2) {
      for (size_t s1 = 0; s1 < HERMITE3D_STENCIL_WIDTH; ++s1) {
        after_xx1[0][s2] +=
            query.weight[1][0][s1] * after_xx0[0][s1][s2];
        after_xx1[1][s2] +=
            query.weight[1][0][s1] * after_xx0[1][s1][s2];
        after_xx1[2][s2] +=
            query.weight[1][1][s1] * after_xx0[0][s1][s2];
      }
    }

    for (size_t s2 = 0; s2 < HERMITE3D_STENCIL_WIDTH; ++s2) {
      result.value += query.weight[2][0][s2] * after_xx1[0][s2];
      result.gradient[0] += query.weight[2][0][s2] * after_xx1[1][s2];
      result.gradient[1] += query.weight[2][0][s2] * after_xx1[2][s2];
      result.gradient[2] += query.weight[2][1][s2] * after_xx1[0][s2];
    }

    results[function] = result;
  }

  return HERMITE3D_SUCCESS;
}

hermite3d_status hermite3d_interpolate_jet(
    const hermite3d_grid *grid, double xx0, double xx1, double xx2,
    size_t num_functions, const double *const functions[],
    hermite3d_jet results[]) {
  hermite3d_prepared_query query;
  hermite3d_status status =
      validate_and_prepare(grid, xx0, xx1, xx2, num_functions, functions,
                           results, &query);

  if (status != HERMITE3D_SUCCESS) {
    return status;
  }

  for (size_t function = 0; function < num_functions; ++function) {
    double after_xx0[HERMITE3D_DERIVATIVE_COUNT]
                    [HERMITE3D_STENCIL_WIDTH]
                    [HERMITE3D_STENCIL_WIDTH] = {{{0.0}}};
    double after_xx1[HERMITE3D_DERIVATIVE_COUNT]
                    [HERMITE3D_DERIVATIVE_COUNT]
                    [HERMITE3D_STENCIL_WIDTH] = {{{0.0}}};
    hermite3d_jet result = {{{{0.0}}}};
    const double *data = functions[function];

    /* First contract the 216 values along xx0 for derivative orders 0..2. */
    for (size_t s2 = 0; s2 < HERMITE3D_STENCIL_WIDTH; ++s2) {
      const size_t offset2 = (query.base[2] + s2) * grid->stride2;
      for (size_t s1 = 0; s1 < HERMITE3D_STENCIL_WIDTH; ++s1) {
        const size_t offset1 = (query.base[1] + s1) * grid->stride1;
        for (size_t s0 = 0; s0 < HERMITE3D_STENCIL_WIDTH; ++s0) {
          const size_t offset = offset2 + offset1 +
                                (query.base[0] + s0) * grid->stride0;
          const double sample = data[offset];
          for (size_t q0 = 0; q0 < HERMITE3D_DERIVATIVE_COUNT; ++q0) {
            after_xx0[q0][s1][s2] += query.weight[0][q0][s0] * sample;
          }
        }
      }
    }

    /* Continue along xx1, producing all nine (q0,q1) combinations. */
    for (size_t q0 = 0; q0 < HERMITE3D_DERIVATIVE_COUNT; ++q0) {
      for (size_t q1 = 0; q1 < HERMITE3D_DERIVATIVE_COUNT; ++q1) {
        for (size_t s2 = 0; s2 < HERMITE3D_STENCIL_WIDTH; ++s2) {
          for (size_t s1 = 0; s1 < HERMITE3D_STENCIL_WIDTH; ++s1) {
            after_xx1[q0][q1][s2] +=
                query.weight[1][q1][s1] * after_xx0[q0][s1][s2];
          }
        }
      }
    }

    /* The xx2 contraction completes every entry in the 3 x 3 x 3 jet. */
    for (size_t q0 = 0; q0 < HERMITE3D_DERIVATIVE_COUNT; ++q0) {
      for (size_t q1 = 0; q1 < HERMITE3D_DERIVATIVE_COUNT; ++q1) {
        for (size_t q2 = 0; q2 < HERMITE3D_DERIVATIVE_COUNT; ++q2) {
          for (size_t s2 = 0; s2 < HERMITE3D_STENCIL_WIDTH; ++s2) {
            result.derivative[q0][q1][q2] +=
                query.weight[2][q2][s2] * after_xx1[q0][q1][s2];
          }
        }
      }
    }

    results[function] = result;
  }

  return HERMITE3D_SUCCESS;
}
