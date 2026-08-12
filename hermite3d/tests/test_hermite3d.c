/*
 * Black-box regression tests for the public C interface.
 *
 * The tests deliberately avoid private implementation details.  They check
 * the interpolant through mathematical invariants (polynomial reproduction,
 * symmetry, endpoint data, and C2 matching), through agreement between the
 * two public entry points, and through the API's documented failure behavior.
 */
#include "hermite3d.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ARRAY_COUNT(array) (sizeof(array) / sizeof((array)[0]))

static int failures = 0;

/* Record all failures in one run so a single invocation can expose more than
 * one regression.  Tests that need richer diagnostics use the same counter
 * directly and print the relevant polynomial and derivative multi-indices. */
#define CHECK(condition)                                               \
  do {                                                                 \
    if (!(condition)) {                                                \
      fprintf(stderr, "%s:%d: check failed: %s\n", __FILE__, __LINE__, \
              #condition);                                             \
      ++failures;                                                      \
    }                                                                  \
  } while (0)

typedef struct test_grid {
  hermite3d_grid descriptor;
  /* Number of double elements needed to reach the largest strided offset. */
  size_t storage_size;
} test_grid;

/*
 * Construct one intentionally nontrivial grid shared by most tests.
 *
 * Different dimensions and spacings catch accidental axis permutations.
 * Shifted, decimal-valued coordinates exercise physical-spacing factors and
 * floating-point endpoint location.  The dimensions leave several eligible
 * cells after reserving the two samples below and three samples above each
 * lower cell endpoint.
 */
static test_grid make_test_grid(void) {
  test_grid result;

  result.descriptor.nxx0 = 9;
  result.descriptor.nxx1 = 10;
  result.descriptor.nxx2 = 11;
  result.descriptor.xx0_start = -0.80;
  result.descriptor.xx1_start = -0.70;
  result.descriptor.xx2_start = -0.90;
  result.descriptor.dxx0 = 0.19;
  result.descriptor.dxx1 = 0.13;
  result.descriptor.dxx2 = 0.11;

  /*
   * None of these strides describes a tightly packed array.  The unused
   * elements between logical samples ensure the implementation honors the
   * descriptor instead of silently assuming xx0-contiguous storage.
   */
  result.descriptor.stride0 = 2;
  result.descriptor.stride1 = 21;
  result.descriptor.stride2 = 215;
  result.storage_size =
      (result.descriptor.nxx0 - 1) * result.descriptor.stride0 +
      (result.descriptor.nxx1 - 1) * result.descriptor.stride1 +
      (result.descriptor.nxx2 - 1) * result.descriptor.stride2 + 1;
  return result;
}

static size_t offset(const hermite3d_grid *grid, size_t i0, size_t i1,
                     size_t i2) {
  return i0 * grid->stride0 + i1 * grid->stride1 + i2 * grid->stride2;
}

/* Construct sample coordinates exactly as an ordinary caller would. */
static double coordinate(double start, double spacing, size_t index) {
  return start + spacing * (double)index;
}

static double integer_power(double value, unsigned int power) {
  double result = 1.0;
  unsigned int q;
  for (q = 0; q < power; ++q) {
    result *= value;
  }
  return result;
}

static double monomial_derivative(double value, unsigned int power,
                                  unsigned int derivative_order) {
  double coefficient = 1.0;
  unsigned int q;

  if (derivative_order > power) {
    return 0.0;
  }
  for (q = 0; q < derivative_order; ++q) {
    coefficient *= (double)(power - q);
  }
  return coefficient * integer_power(value, power - derivative_order);
}

static int nearly_equal(double actual, double expected, double tolerance) {
  /* The floor at one makes tolerance meaningful for derivatives that vanish. */
  const double scale = fmax(1.0, fmax(fabs(actual), fabs(expected)));
  return fabs(actual - expected) <= tolerance * scale;
}

static double *allocate_field(const test_grid *grid) {
  /* calloc also gives padding elements a recognizable value that is unrelated
   * to the analytic fields stored at logical samples. */
  double *field = (double *)calloc(grid->storage_size, sizeof(*field));
  CHECK(field != NULL);
  return field;
}

static void fill_monomial(const test_grid *test, double *field,
                          unsigned int p0, unsigned int p1,
                          unsigned int p2) {
  /* Store xx0^p0 xx1^p1 xx2^p2 using the public strided layout. */
  const hermite3d_grid *grid = &test->descriptor;
  size_t i0;
  size_t i1;
  size_t i2;

  for (i2 = 0; i2 < grid->nxx2; ++i2) {
    const double xx2 = coordinate(grid->xx2_start, grid->dxx2, i2);
    for (i1 = 0; i1 < grid->nxx1; ++i1) {
      const double xx1 = coordinate(grid->xx1_start, grid->dxx1, i1);
      for (i0 = 0; i0 < grid->nxx0; ++i0) {
        const double xx0 = coordinate(grid->xx0_start, grid->dxx0, i0);
        field[offset(grid, i0, i1, i2)] =
            integer_power(xx0, p0) * integer_power(xx1, p1) *
            integer_power(xx2, p2);
      }
    }
  }
}

static void fill_smooth_field(const test_grid *test, double *field,
                              double phase) {
  /*
   * A nonpolynomial, nonseparable field is useful when comparing independent
   * API paths: no exactness property can accidentally hide a missing axis or
   * a contraction performed in the wrong order.  phase creates distinct
   * functions while preserving comparable magnitudes.
   */
  const hermite3d_grid *grid = &test->descriptor;
  size_t i0;
  size_t i1;
  size_t i2;

  for (i2 = 0; i2 < grid->nxx2; ++i2) {
    const double xx2 = coordinate(grid->xx2_start, grid->dxx2, i2);
    for (i1 = 0; i1 < grid->nxx1; ++i1) {
      const double xx1 = coordinate(grid->xx1_start, grid->dxx1, i1);
      for (i0 = 0; i0 < grid->nxx0; ++i0) {
        const double xx0 = coordinate(grid->xx0_start, grid->dxx0, i0);
        field[offset(grid, i0, i1, i2)] =
            sin(0.7 * xx0 + phase) + cos(0.4 * xx1 - 0.3 * phase) +
            exp(0.2 * xx2) + 0.15 * xx0 * xx1 * xx2;
      }
    }
  }
}

/*
 * The one-dimensional value-only operator reproduces every polynomial through
 * degree four.  Tensor separability therefore requires exact reproduction of
 * xx0^p0 xx1^p1 xx2^p2 for every p0,p1,p2 <= 4.  Differentiating that identity
 * also fixes all 27 entries of the returned jet, including high-total-order
 * mixed derivatives.  This single sweep is the main mathematical regression.
 */
static void test_tensor_monomials(void) {
  const test_grid test = make_test_grid();
  const hermite3d_grid *grid = &test.descriptor;
  double *field = allocate_field(&test);
  const double *functions[1] = {field};
  const double xx0 = grid->xx0_start + 3.37 * grid->dxx0;
  const double xx1 = grid->xx1_start + 4.23 * grid->dxx1;
  const double xx2 = grid->xx2_start + 5.61 * grid->dxx2;
  unsigned int p0;
  unsigned int p1;
  unsigned int p2;

  if (field == NULL) {
    return;
  }

  for (p2 = 0; p2 <= 4; ++p2) {
    for (p1 = 0; p1 <= 4; ++p1) {
      for (p0 = 0; p0 <= 4; ++p0) {
        hermite3d_jet result;
        unsigned int q0;
        unsigned int q1;
        unsigned int q2;

        fill_monomial(&test, field, p0, p1, p2);
        CHECK(hermite3d_interpolate_jet(grid, xx0, xx1, xx2, 1,
                                        functions, &result) ==
              HERMITE3D_SUCCESS);

        /* q0, q1, and q2 are componentwise derivative orders, not a bound on
         * their sum.  For example, (2,2,2) is intentionally included. */
        for (q2 = 0; q2 <= 2; ++q2) {
          for (q1 = 0; q1 <= 2; ++q1) {
            for (q0 = 0; q0 <= 2; ++q0) {
              const double expected =
                  monomial_derivative(xx0, p0, q0) *
                  monomial_derivative(xx1, p1, q1) *
                  monomial_derivative(xx2, p2, q2);
              /* High mixed derivatives multiply several inverse spacings and
               * amplify ordinary rounding from the three staged contractions.
               * The identity being tested is exact; this tolerance covers only
               * its floating-point evaluation on the anisotropic test grid. */
              if (!nearly_equal(result.derivative[q0][q1][q2], expected,
                                2.0e-7)) {
                fprintf(stderr,
                        "monomial (%u,%u,%u), derivative (%u,%u,%u): "
                        "got %.17g, expected %.17g\n",
                        p0, p1, p2, q0, q1, q2,
                        result.derivative[q0][q1][q2], expected);
                ++failures;
              }
            }
          }
        }
      }
    }
  }

  free(field);
}

/*
 * The compact routine has a specialized contraction that computes only four
 * outputs.  Compare it with the corresponding entries of the general jet so
 * the optimization cannot change its mathematics.  Passing three distinct
 * fields at once also verifies ordering and the advertised weight reuse API.
 */
static void test_value_gradient_matches_jet_and_multiple_functions(void) {
  const test_grid test = make_test_grid();
  const hermite3d_grid *grid = &test.descriptor;
  double *fields[3];
  const double *functions[3];
  hermite3d_value_gradient compact[3];
  hermite3d_jet full[3];
  const double xx0 = grid->xx0_start + 3.19 * grid->dxx0;
  const double xx1 = grid->xx1_start + 4.73 * grid->dxx1;
  const double xx2 = grid->xx2_start + 6.11 * grid->dxx2;
  size_t function;

  for (function = 0; function < ARRAY_COUNT(fields); ++function) {
    fields[function] = allocate_field(&test);
    if (fields[function] == NULL) {
      while (function > 0) {
        free(fields[--function]);
      }
      return;
    }
    fill_smooth_field(&test, fields[function], 0.3 * (double)function);
    functions[function] = fields[function];
  }

  CHECK(hermite3d_interpolate_value_gradient(
            grid, xx0, xx1, xx2, ARRAY_COUNT(functions), functions, compact) ==
        HERMITE3D_SUCCESS);
  CHECK(hermite3d_interpolate_jet(grid, xx0, xx1, xx2,
                                  ARRAY_COUNT(functions), functions, full) ==
        HERMITE3D_SUCCESS);

  for (function = 0; function < ARRAY_COUNT(functions); ++function) {
    CHECK(nearly_equal(compact[function].value,
                       full[function].derivative[0][0][0], 2.0e-13));
    CHECK(nearly_equal(compact[function].gradient[0],
                       full[function].derivative[1][0][0], 2.0e-13));
    CHECK(nearly_equal(compact[function].gradient[1],
                       full[function].derivative[0][1][0], 2.0e-13));
    CHECK(nearly_equal(compact[function].gradient[2],
                       full[function].derivative[0][0][1], 2.0e-13));
    free(fields[function]);
  }
}

/*
 * A constant probes partition of unity: its value must survive and all 26
 * nonzero-order jet entries must vanish.  At a grid node, the endpoint
 * cardinal conditions additionally require the interpolant's first two
 * derivatives to equal the centered finite-difference data from which the
 * value-only Hermite polynomial is constructed.
 */
static void test_constant_and_endpoint_derivatives(void) {
  const test_grid test = make_test_grid();
  const hermite3d_grid *grid = &test.descriptor;
  double *field = allocate_field(&test);
  const double *functions[1] = {field};
  hermite3d_jet result;
  size_t q0;
  size_t q1;
  size_t q2;
  size_t i0;
  size_t i1;
  size_t i2;

  if (field == NULL) {
    return;
  }

  for (i2 = 0; i2 < grid->nxx2; ++i2) {
    for (i1 = 0; i1 < grid->nxx1; ++i1) {
      for (i0 = 0; i0 < grid->nxx0; ++i0) {
        field[offset(grid, i0, i1, i2)] = 2.75;
      }
    }
  }

  CHECK(hermite3d_interpolate_jet(
            grid, grid->xx0_start + 3.4 * grid->dxx0,
            grid->xx1_start + 4.2 * grid->dxx1,
            grid->xx2_start + 5.7 * grid->dxx2, 1, functions, &result) ==
        HERMITE3D_SUCCESS);
  for (q2 = 0; q2 <= 2; ++q2) {
    for (q1 = 0; q1 <= 2; ++q1) {
      for (q0 = 0; q0 <= 2; ++q0) {
        const double expected = (q0 == 0 && q1 == 0 && q2 == 0) ? 2.75 : 0.0;
        CHECK(nearly_equal(result.derivative[q0][q1][q2], expected, 2.0e-10));
      }
    }
  }

  /*
   * Check one direction explicitly against the familiar five-point formulas.
   * The tensor factors in the other directions evaluate to their node values,
   * isolating the xx0 endpoint derivative weights.
   */
  fill_smooth_field(&test, field, 0.17);
  i0 = 4;
  i1 = 5;
  i2 = 6;
  CHECK(hermite3d_interpolate_jet(
            grid, coordinate(grid->xx0_start, grid->dxx0, i0),
            coordinate(grid->xx1_start, grid->dxx1, i1),
            coordinate(grid->xx2_start, grid->dxx2, i2), 1, functions,
            &result) == HERMITE3D_SUCCESS);
  {
    const double fm2 = field[offset(grid, i0 - 2, i1, i2)];
    const double fm1 = field[offset(grid, i0 - 1, i1, i2)];
    const double f0 = field[offset(grid, i0, i1, i2)];
    const double fp1 = field[offset(grid, i0 + 1, i1, i2)];
    const double fp2 = field[offset(grid, i0 + 2, i1, i2)];
    const double expected_first =
        (fm2 - 8.0 * fm1 + 8.0 * fp1 - fp2) / (12.0 * grid->dxx0);
    const double expected_second =
        (-fm2 + 16.0 * fm1 - 30.0 * f0 + 16.0 * fp1 - fp2) /
        (12.0 * grid->dxx0 * grid->dxx0);
    CHECK(nearly_equal(result.derivative[0][0][0], f0, 2.0e-13));
    CHECK(nearly_equal(result.derivative[1][0][0], expected_first, 2.0e-12));
    CHECK(nearly_equal(result.derivative[2][0][0], expected_second, 2.0e-11));
  }

  free(field);
}

/*
 * Reflect both the samples and query across the center of every grid axis.
 * Values and even-total-order derivatives are invariant, while every odd
 * total derivative changes sign by the chain rule.  Testing the complete jet
 * exercises orientation signs in all one-dimensional weight tables.
 */
static void test_reflection_symmetry(void) {
  const test_grid test = make_test_grid();
  const hermite3d_grid *grid = &test.descriptor;
  double *field = allocate_field(&test);
  double *reflected = allocate_field(&test);
  const double *field_pointer[1] = {field};
  const double *reflected_pointer[1] = {reflected};
  hermite3d_jet original_result;
  hermite3d_jet reflected_result;
  const double u0 = 3.27;
  const double u1 = 4.38;
  const double u2 = 5.49;
  size_t i0;
  size_t i1;
  size_t i2;
  size_t q0;
  size_t q1;
  size_t q2;

  if (field == NULL || reflected == NULL) {
    free(field);
    free(reflected);
    return;
  }

  fill_smooth_field(&test, field, 0.41);
  for (i2 = 0; i2 < grid->nxx2; ++i2) {
    for (i1 = 0; i1 < grid->nxx1; ++i1) {
      for (i0 = 0; i0 < grid->nxx0; ++i0) {
        reflected[offset(grid, i0, i1, i2)] =
            field[offset(grid, grid->nxx0 - 1 - i0, grid->nxx1 - 1 - i1,
                         grid->nxx2 - 1 - i2)];
      }
    }
  }

  CHECK(hermite3d_interpolate_jet(
            grid, grid->xx0_start + u0 * grid->dxx0,
            grid->xx1_start + u1 * grid->dxx1,
            grid->xx2_start + u2 * grid->dxx2, 1, field_pointer,
            &original_result) == HERMITE3D_SUCCESS);
  CHECK(hermite3d_interpolate_jet(
            grid,
            grid->xx0_start + ((double)(grid->nxx0 - 1) - u0) * grid->dxx0,
            grid->xx1_start + ((double)(grid->nxx1 - 1) - u1) * grid->dxx1,
            grid->xx2_start + ((double)(grid->nxx2 - 1) - u2) * grid->dxx2,
            1, reflected_pointer, &reflected_result) == HERMITE3D_SUCCESS);

  for (q2 = 0; q2 <= 2; ++q2) {
    for (q1 = 0; q1 <= 2; ++q1) {
      for (q0 = 0; q0 <= 2; ++q0) {
        const double sign = ((q0 + q1 + q2) % 2 == 0) ? 1.0 : -1.0;
        CHECK(nearly_equal(reflected_result.derivative[q0][q1][q2],
                           sign * original_result.derivative[q0][q1][q2],
                           2.0e-9));
      }
    }
  }

  free(field);
  free(reflected);
}

/*
 * Adjacent cell polynomials share value, slope, and curvature data at a grid
 * plane.  Sample the interpolant just to either side of one xx0 plane and at
 * the plane itself to detect a discontinuity in normal derivatives q0=0,1,2.
 * epsilon is measured in normalized grid coordinates, so both nearby queries
 * remain well inside their respective cells.
 */
static void test_c2_continuity(void) {
  const test_grid test = make_test_grid();
  const hermite3d_grid *grid = &test.descriptor;
  double *field = allocate_field(&test);
  const double *functions[1] = {field};
  hermite3d_jet left;
  hermite3d_jet node;
  hermite3d_jet right;
  const double epsilon = 1.0e-7;
  const double node_u0 = 4.0;
  const double u1 = 4.31;
  const double u2 = 5.62;
  size_t q0;

  if (field == NULL) {
    return;
  }
  fill_smooth_field(&test, field, 0.73);

  CHECK(hermite3d_interpolate_jet(
            grid, grid->xx0_start + (node_u0 - epsilon) * grid->dxx0,
            grid->xx1_start + u1 * grid->dxx1,
            grid->xx2_start + u2 * grid->dxx2, 1, functions, &left) ==
        HERMITE3D_SUCCESS);
  CHECK(hermite3d_interpolate_jet(
            grid, grid->xx0_start + node_u0 * grid->dxx0,
            grid->xx1_start + u1 * grid->dxx1,
            grid->xx2_start + u2 * grid->dxx2, 1, functions, &node) ==
        HERMITE3D_SUCCESS);
  CHECK(hermite3d_interpolate_jet(
            grid, grid->xx0_start + (node_u0 + epsilon) * grid->dxx0,
            grid->xx1_start + u1 * grid->dxx1,
            grid->xx2_start + u2 * grid->dxx2, 1, functions, &right) ==
        HERMITE3D_SUCCESS);

  for (q0 = 0; q0 <= 2; ++q0) {
    /* Nearby values need not be bit-identical to the node value because they
     * are evaluated epsilon away; the tolerance scales with that displacement
     * and is tight enough to reveal a finite jump between the two cells. */
    CHECK(nearly_equal(left.derivative[q0][0][0], node.derivative[q0][0][0],
                       2.0e-6));
    CHECK(nearly_equal(right.derivative[q0][0][0], node.derivative[q0][0][0],
                       2.0e-6));
  }

  free(field);
}

/* Distinct sentinels make partial writes visible byte-for-byte after errors. */
static void set_value_gradient_sentinel(hermite3d_value_gradient *result) {
  result->value = 123456.25;
  result->gradient[0] = -31.0;
  result->gradient[1] = -32.0;
  result->gradient[2] = -33.0;
}

static void set_jet_sentinel(hermite3d_jet *result) {
  size_t q0;
  size_t q1;
  size_t q2;
  for (q0 = 0; q0 < 3; ++q0) {
    for (q1 = 0; q1 < 3; ++q1) {
      for (q2 = 0; q2 < 3; ++q2) {
        result->derivative[q0][q1][q2] =
            1000.0 + 100.0 * (double)q0 + 10.0 * (double)q1 + (double)q2;
      }
    }
  }
}

/*
 * Exercise the closed complete-stencil domain and each public error class.
 * The sentinel comparisons enforce transactional behavior: validation must
 * finish before the first result is written, so callers never receive a
 * mixture of old and partially computed outputs.
 */
static void test_valid_domain_and_errors(void) {
  const test_grid test = make_test_grid();
  const hermite3d_grid *grid = &test.descriptor;
  double *field = allocate_field(&test);
  const double *functions[1] = {field};
  const double *null_function[1] = {NULL};
  hermite3d_value_gradient result;
  hermite3d_value_gradient saved;
  hermite3d_jet jet_result;
  hermite3d_jet saved_jet;
  hermite3d_grid invalid;
  const double valid_xx0 = grid->xx0_start + 3.2 * grid->dxx0;
  const double valid_xx1 = grid->xx1_start + 4.2 * grid->dxx1;
  const double valid_xx2 = grid->xx2_start + 5.2 * grid->dxx2;

  if (field == NULL) {
    return;
  }
  fill_smooth_field(&test, field, 0.0);

  /*
   * Both complete-stencil endpoints are inclusive.  Decimal starts and
   * spacings make these calls regress the locator's endpoint-roundoff handling
   * as well as the lower-cell choice at the upper endpoint.
   */
  CHECK(hermite3d_interpolate_value_gradient(
            grid, grid->xx0_start + 2.0 * grid->dxx0,
            grid->xx1_start + 2.0 * grid->dxx1,
            grid->xx2_start + 2.0 * grid->dxx2, 1, functions, &result) ==
        HERMITE3D_SUCCESS);
  CHECK(hermite3d_interpolate_value_gradient(
            grid, grid->xx0_start + (double)(grid->nxx0 - 3) * grid->dxx0,
            grid->xx1_start + (double)(grid->nxx1 - 3) * grid->dxx1,
            grid->xx2_start + (double)(grid->nxx2 - 3) * grid->dxx2, 1,
            functions, &result) == HERMITE3D_SUCCESS);

  set_value_gradient_sentinel(&result);
  saved = result;
  /* Every failing compact call below must preserve the full result object. */
#define CHECK_VALUE_ERROR(expected_status, call)         \
  do {                                                   \
    CHECK((call) == (expected_status));                  \
    CHECK(memcmp(&result, &saved, sizeof(result)) == 0); \
  } while (0)

  CHECK_VALUE_ERROR(
      HERMITE3D_STENCIL_UNAVAILABLE,
      hermite3d_interpolate_value_gradient(
          grid, grid->xx0_start + 1.99 * grid->dxx0, valid_xx1, valid_xx2, 1,
          functions, &result));
  CHECK_VALUE_ERROR(
      HERMITE3D_STENCIL_UNAVAILABLE,
      hermite3d_interpolate_value_gradient(
          grid, grid->xx0_start + ((double)(grid->nxx0 - 3) + 0.01) * grid->dxx0,
          valid_xx1, valid_xx2, 1, functions, &result));
  CHECK_VALUE_ERROR(
      HERMITE3D_INVALID_ARGUMENT,
      hermite3d_interpolate_value_gradient(NULL, valid_xx0, valid_xx1,
                                           valid_xx2, 1, functions, &result));
  CHECK_VALUE_ERROR(
      HERMITE3D_INVALID_ARGUMENT,
      hermite3d_interpolate_value_gradient(grid, valid_xx0, valid_xx1,
                                           valid_xx2, 0, functions, &result));
  CHECK_VALUE_ERROR(
      HERMITE3D_INVALID_ARGUMENT,
      hermite3d_interpolate_value_gradient(grid, valid_xx0, valid_xx1,
                                           valid_xx2, 1, NULL, &result));
  CHECK_VALUE_ERROR(
      HERMITE3D_INVALID_ARGUMENT,
      hermite3d_interpolate_value_gradient(grid, valid_xx0, valid_xx1,
                                           valid_xx2, 1, null_function,
                                           &result));
  CHECK(hermite3d_interpolate_value_gradient(grid, valid_xx0, valid_xx1,
                                             valid_xx2, 1, functions, NULL) ==
        HERMITE3D_INVALID_ARGUMENT);
  CHECK_VALUE_ERROR(
      HERMITE3D_INVALID_ARGUMENT,
      hermite3d_interpolate_value_gradient(grid, NAN, valid_xx1, valid_xx2, 1,
                                           functions, &result));

  /* Distinguish malformed geometry from a valid grid lacking the query's
   * complete six-point stencil. */
  invalid = *grid;
  invalid.nxx0 = 5;
  CHECK_VALUE_ERROR(HERMITE3D_INVALID_GRID,
                    hermite3d_interpolate_value_gradient(
                        &invalid, valid_xx0, valid_xx1, valid_xx2, 1, functions,
                        &result));
  invalid = *grid;
  invalid.dxx1 = 0.0;
  CHECK_VALUE_ERROR(HERMITE3D_INVALID_GRID,
                    hermite3d_interpolate_value_gradient(
                        &invalid, valid_xx0, valid_xx1, valid_xx2, 1, functions,
                        &result));
  invalid = *grid;
  invalid.dxx2 = -1.0;
  CHECK_VALUE_ERROR(HERMITE3D_INVALID_GRID,
                    hermite3d_interpolate_value_gradient(
                        &invalid, valid_xx0, valid_xx1, valid_xx2, 1, functions,
                        &result));
  invalid = *grid;
  invalid.xx0_start = INFINITY;
  CHECK_VALUE_ERROR(HERMITE3D_INVALID_GRID,
                    hermite3d_interpolate_value_gradient(
                        &invalid, valid_xx0, valid_xx1, valid_xx2, 1, functions,
                        &result));
  invalid = *grid;
  invalid.stride0 = 0;
  CHECK_VALUE_ERROR(HERMITE3D_INVALID_GRID,
                    hermite3d_interpolate_value_gradient(
                        &invalid, valid_xx0, valid_xx1, valid_xx2, 1, functions,
                        &result));

  /* A descriptor whose maximum element offset would overflow size_t must fail
   * before any pointer arithmetic is attempted. */
  invalid = *grid;
  invalid.nxx0 = SIZE_MAX;
  invalid.stride0 = 2;
  CHECK_VALUE_ERROR(HERMITE3D_INDEX_OVERFLOW,
                    hermite3d_interpolate_value_gradient(
                        &invalid, valid_xx0, valid_xx1, valid_xx2, 1, functions,
                        &result));

  /* Repeat the unchanged-output contract through the full-jet entry point. */
  set_jet_sentinel(&jet_result);
  saved_jet = jet_result;
  CHECK(hermite3d_interpolate_jet(grid, valid_xx0, valid_xx1, valid_xx2, 1,
                                  null_function, &jet_result) ==
        HERMITE3D_INVALID_ARGUMENT);
  CHECK(memcmp(&jet_result, &saved_jet, sizeof(jet_result)) == 0);

#undef CHECK_VALUE_ERROR
  free(field);
}

int main(void) {
  /* Keep the calls explicit so a failure can be associated with one invariant
   * without introducing a test framework dependency. */
  test_tensor_monomials();
  test_value_gradient_matches_jet_and_multiple_functions();
  test_constant_and_endpoint_derivatives();
  test_reflection_symmetry();
  test_c2_continuity();
  test_valid_domain_and_errors();

  if (failures != 0) {
    fprintf(stderr, "hermite3d tests: %d failure(s)\n", failures);
    return EXIT_FAILURE;
  }

  puts("hermite3d tests: all checks passed");
  return EXIT_SUCCESS;
}
