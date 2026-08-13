#include "lagrange3d_reference.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { DIMENSION = 15 };

static int checks = 0;
static int failures = 0;

static void check(int condition, const char *message) {
  ++checks;
  if (!condition) {
    ++failures;
    fprintf(stderr, "FAIL: %s\n", message);
  }
}

static double integer_power(double value, unsigned exponent) {
  double result = 1.0;
  while (exponent-- > 0) result *= value;
  return result;
}

static void set_monomial(double *storage, const double *xx0, const double *xx1,
                         const double *xx2, unsigned p0, unsigned p1,
                         unsigned p2) {
  for (size_t i2 = 0; i2 < DIMENSION; ++i2)
    for (size_t i1 = 0; i1 < DIMENSION; ++i1)
      for (size_t i0 = 0; i0 < DIMENSION; ++i0) {
        const size_t offset = i0 + DIMENSION * (i1 + DIMENSION * i2);
        storage[offset] = integer_power(xx0[i0], p0) *
                          integer_power(xx1[i1], p1) *
                          integer_power(xx2[i2], p2);
      }
}

static void initialize_grid(lagrange3d_reference_grid *grid,
                            double *storage[KERR_ADM_FIELD_COUNT]) {
  static double coordinate[3][DIMENSION];
  static const double start[3] = {-1.4, -1.7, -2.0};
  static const double spacing[3] = {0.2, 0.25, 0.3};

  for (size_t axis = 0; axis < 3; ++axis) {
    grid->dimension[axis] = DIMENSION;
    grid->spacing[axis] = spacing[axis];
    grid->coordinate[axis] = coordinate[axis];
    for (size_t index = 0; index < DIMENSION; ++index)
      coordinate[axis][index] = start[axis] + (double)index * spacing[axis];
  }
  for (size_t field = 0; field < KERR_ADM_FIELD_COUNT; ++field) {
    storage[field] =
        (double *)malloc(DIMENSION * DIMENSION * DIMENSION * sizeof(double));
    check(storage[field] != NULL, "monomial storage allocation");
    grid->function[field] = storage[field];
  }
}

static void test_tensor_monomials(lagrange3d_reference_grid *grid,
                                  double *storage[KERR_ADM_FIELD_COUNT]) {
  static const double query[3] = {0.17, -0.23, 0.31};

  for (unsigned first = 0; first < 343; first += KERR_ADM_FIELD_COUNT) {
    double values[KERR_ADM_FIELD_COUNT];
    for (size_t field = 0; field < KERR_ADM_FIELD_COUNT; ++field) {
      const unsigned monomial = first + (unsigned)field;
      const unsigned active = monomial < 343 ? monomial : 0;
      const unsigned p0 = active % 7;
      const unsigned p1 = (active / 7) % 7;
      const unsigned p2 = active / 49;
      set_monomial(storage[field], grid->coordinate[0], grid->coordinate[1],
                   grid->coordinate[2], p0, p1, p2);
    }
    check(lagrange3d_reference_interpolate(grid, query[0], query[1], query[2],
                                           values) ==
              LAGRANGE3D_REFERENCE_SUCCESS,
          "tensor monomial interpolation status");
    for (size_t field = 0; field < KERR_ADM_FIELD_COUNT; ++field) {
      const unsigned monomial = first + (unsigned)field;
      if (monomial < 343) {
        const unsigned p0 = monomial % 7;
        const unsigned p1 = (monomial / 7) % 7;
        const unsigned p2 = monomial / 49;
        const double expected = integer_power(query[0], p0) *
                                integer_power(query[1], p1) *
                                integer_power(query[2], p2);
        check(fabs(values[field] - expected) <=
                  2e-10 * fmax(1.0, fabs(expected)),
              "tensor monomial reproduction through degree six");
      }
    }
  }
}

static void test_fd_formula(void) {
  const double x = 0.37;
  const double h = 0.02;

  for (unsigned degree = 0; degree <= 5; ++degree) {
    const double derivative =
        (integer_power(x - 2 * h, degree) -
         8 * integer_power(x - h, degree) +
         8 * integer_power(x + h, degree) -
         integer_power(x + 2 * h, degree)) /
        (12 * h);
    const double exact = degree == 0
                             ? 0.0
                             : (double)degree * integer_power(x, degree - 1);
    if (degree <= 4)
      check(fabs(derivative - exact) <= 2e-13,
            "fourth-order FD polynomial exactness");
    else
      check(fabs((derivative - exact) + 4 * integer_power(h, 4)) <= 2e-13,
            "degree-five FD defect");
  }
}

/*
 * Six stored layers cover the floating-point endpoint case in which a query
 * just inside the active face has its +2H coordinate round to the outward
 * midpoint.  The copied nearest-node locator then chooses the outer node.
 * The current test grid has six layers around active bounds (-0.3,0.3).
 */
static void test_six_layer_padding(lagrange3d_reference_grid *grid) {
  double values[KERR_ADM_FIELD_COUNT];
  const double near_upper_face = nextafter(0.3, 0.0);
  const double displaced = near_upper_face + 2.0 * grid->spacing[0];

  check(lagrange3d_reference_interpolate(grid, displaced, -0.23, 0.31,
                                         values) ==
            LAGRANGE3D_REFERENCE_SUCCESS,
        "six-layer upper-face padding");
  check(lagrange3d_reference_interpolate(grid, -displaced, -0.23, 0.31,
                                         values) ==
            LAGRANGE3D_REFERENCE_SUCCESS,
        "six-layer lower-face padding");
}

static void test_adapter_errors(lagrange3d_reference_grid *grid) {
  double output[KERR_ADM_FIELD_COUNT];
  double original[KERR_ADM_FIELD_COUNT];
  lagrange3d_reference_grid invalid = *grid;

  for (size_t field = 0; field < KERR_ADM_FIELD_COUNT; ++field)
    output[field] = 100.0 + (double)field;
  memcpy(original, output, sizeof(output));
  check(lagrange3d_reference_interpolate(NULL, 0, 0, 0, output) ==
            LAGRANGE3D_REFERENCE_INVALID_ARGUMENT,
        "null grid rejected");
  check(memcmp(output, original, sizeof(output)) == 0,
        "output unchanged after null grid");
  check(lagrange3d_reference_interpolate(grid, NAN, 0, 0, output) ==
            LAGRANGE3D_REFERENCE_INVALID_ARGUMENT,
        "nonfinite query rejected");
  check(memcmp(output, original, sizeof(output)) == 0,
        "output unchanged after nonfinite query");
  invalid.dimension[0] = 6;
  check(lagrange3d_reference_interpolate(&invalid, 0, 0, 0, output) ==
            LAGRANGE3D_REFERENCE_INVALID_GRID,
        "short grid rejected");
  invalid = *grid;
  invalid.spacing[1] = 0.0;
  check(lagrange3d_reference_interpolate(&invalid, 0, 0, 0, output) ==
            LAGRANGE3D_REFERENCE_INVALID_GRID,
        "zero spacing rejected");
  invalid = *grid;
  invalid.function[3] = NULL;
  check(lagrange3d_reference_interpolate(&invalid, 0, 0, 0, output) ==
            LAGRANGE3D_REFERENCE_INVALID_ARGUMENT,
        "null function rejected");
  check(lagrange3d_reference_interpolate(grid, -10, 0, 0, output) ==
            LAGRANGE3D_REFERENCE_STENCIL_UNAVAILABLE,
        "unavailable stencil reported");
  check(memcmp(output, original, sizeof(output)) == 0,
        "output unchanged after unavailable stencil");
}

int main(void) {
  lagrange3d_reference_grid grid;
  double *storage[KERR_ADM_FIELD_COUNT] = {0};

  initialize_grid(&grid, storage);
  if (failures == 0) {
    test_tensor_monomials(&grid, storage);
    test_six_layer_padding(&grid);
    test_adapter_errors(&grid);
  }
  test_fd_formula();
  for (size_t field = 0; field < KERR_ADM_FIELD_COUNT; ++field)
    free(storage[field]);
  if (failures != 0) {
    fprintf(stderr, "%d of %d Lagrange/FD checks failed\n", failures, checks);
    return 1;
  }
  printf("Lagrange reference and FD: %d checks passed\n", checks);
  return 0;
}
