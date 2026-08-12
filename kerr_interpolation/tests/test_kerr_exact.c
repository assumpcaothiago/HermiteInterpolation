#include "kerr_exact.h"
#include "kerr_exact_fixtures.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

static int failures = 0;

#define CHECK(condition, description)                                      \
  do {                                                                     \
    if (!(condition)) {                                                    \
      fprintf(stderr, "[FAIL] %s (line %d)\n", description, __LINE__);     \
      ++failures;                                                          \
    }                                                                      \
  } while (0)

static int close_scaled(double actual, double expected, double tolerance) {
  const double scale = fmax(1.0, fabs(expected));
  return isfinite(actual) && fabs(actual - expected) <= tolerance * scale;
}

static void test_high_precision_fixtures(void) {
  for (size_t fixture_index = 0; fixture_index < KERR_EXACT_FIXTURE_COUNT;
       ++fixture_index) {
    const kerr_exact_fixture *fixture = &kerr_exact_fixtures[fixture_index];
    double values[KERR_EXACT_COMPONENT_COUNT];
    kerr_exact_value_gradient jet[KERR_EXACT_COMPONENT_COUNT];

    CHECK(kerr_exact_metric(fixture->coordinate[0], fixture->coordinate[1],
                            fixture->coordinate[2], values) ==
              KERR_EXACT_SUCCESS,
          "fixture value evaluation succeeds");
    CHECK(kerr_exact_metric_gradient(
              fixture->coordinate[0], fixture->coordinate[1],
              fixture->coordinate[2], jet) == KERR_EXACT_SUCCESS,
          "fixture value-gradient evaluation succeeds");

    for (size_t component = 0; component < KERR_EXACT_COMPONENT_COUNT;
         ++component) {
      CHECK(close_scaled(values[component], fixture->value[component], 3e-12),
            "value-only evaluator agrees with 80-digit fixture");
      CHECK(close_scaled(jet[component].value, fixture->value[component],
                         3e-12),
            "value-gradient evaluator agrees with value fixture");
      CHECK(close_scaled(values[component], jet[component].value, 2e-14),
            "independently optimized value kernels agree within roundoff");
      for (size_t axis = 0; axis < KERR_EXACT_DIMENSION; ++axis) {
        CHECK(close_scaled(jet[component].gradient[axis],
                           fixture->gradient[component][axis], 3e-11),
              "analytic gradient agrees with 80-digit fixture");
      }
    }
  }
}

static int finite_difference_component(const double point[3], size_t axis,
                                       double step, size_t component,
                                       double *derivative) {
  double shifted[3];
  double minus_two[KERR_EXACT_COMPONENT_COUNT];
  double minus_one[KERR_EXACT_COMPONENT_COUNT];
  double plus_one[KERR_EXACT_COMPONENT_COUNT];
  double plus_two[KERR_EXACT_COMPONENT_COUNT];

  memcpy(shifted, point, sizeof(shifted));
  shifted[axis] = point[axis] - 2.0 * step;
  if (kerr_exact_metric(shifted[0], shifted[1], shifted[2], minus_two) !=
      KERR_EXACT_SUCCESS)
    return 0;
  shifted[axis] = point[axis] - step;
  if (kerr_exact_metric(shifted[0], shifted[1], shifted[2], minus_one) !=
      KERR_EXACT_SUCCESS)
    return 0;
  shifted[axis] = point[axis] + step;
  if (kerr_exact_metric(shifted[0], shifted[1], shifted[2], plus_one) !=
      KERR_EXACT_SUCCESS)
    return 0;
  shifted[axis] = point[axis] + 2.0 * step;
  if (kerr_exact_metric(shifted[0], shifted[1], shifted[2], plus_two) !=
      KERR_EXACT_SUCCESS)
    return 0;

  *derivative =
      (minus_two[component] - 8.0 * minus_one[component] +
       8.0 * plus_one[component] - plus_two[component]) /
      (12.0 * step);
  return 1;
}

static void test_gradients_against_finite_differences(void) {
  static const double points[][3] = {
      {0.81, -0.47, 1.13}, {-0.19, 0.31, 0.27}, {1.41, 0.73, -0.88}};
  static const double steps[] = {3e-4, 1e-4, 3e-5, 1e-5};

  for (size_t point_index = 0; point_index < sizeof(points) / sizeof(points[0]);
       ++point_index) {
    kerr_exact_value_gradient exact[KERR_EXACT_COMPONENT_COUNT];
    CHECK(kerr_exact_metric_gradient(
              points[point_index][0], points[point_index][1],
              points[point_index][2], exact) == KERR_EXACT_SUCCESS,
          "finite-difference reference evaluation succeeds");

    for (size_t component = 0; component < KERR_EXACT_COMPONENT_COUNT;
         ++component) {
      for (size_t axis = 0; axis < KERR_EXACT_DIMENSION; ++axis) {
        double best_error = INFINITY;
        for (size_t step_index = 0;
             step_index < sizeof(steps) / sizeof(steps[0]); ++step_index) {
          double approximation = 0.0;
          const int evaluated = finite_difference_component(
              points[point_index], axis, steps[step_index], component,
              &approximation);
          CHECK(evaluated,
                "finite-difference samples are evaluable");
          if (evaluated && isfinite(approximation)) {
            const double error =
                fabs(approximation - exact[component].gradient[axis]);
            if (error < best_error) best_error = error;
          }
        }
        CHECK(best_error <=
                  2e-8 * fmax(1.0, fabs(exact[component].gradient[axis])),
              "analytic derivative agrees with independent five-point difference");
      }
    }
  }
}

static void test_symmetry_and_zero_channel(void) {
  const double point[3] = {0.83, -0.52, 0.91};
  kerr_exact_value_gradient positive[KERR_EXACT_COMPONENT_COUNT];
  kerr_exact_value_gradient negative[KERR_EXACT_COMPONENT_COUNT];

  CHECK(kerr_exact_metric_gradient(point[0], point[1], point[2], positive) ==
            KERR_EXACT_SUCCESS,
        "positive parity point is evaluable");
  CHECK(kerr_exact_metric_gradient(-point[0], -point[1], -point[2], negative) ==
            KERR_EXACT_SUCCESS,
        "inverted parity point is evaluable");

  for (size_t component = 0; component < KERR_EXACT_COMPONENT_COUNT;
       ++component) {
    const int odd = component == KERR_EXACT_TX ||
                    component == KERR_EXACT_TY ||
                    component == KERR_EXACT_TZ;
    const double parity = odd ? -1.0 : 1.0;
    CHECK(close_scaled(negative[component].value,
                       parity * positive[component].value, 2e-13),
          "metric component has inversion parity");
    for (size_t axis = 0; axis < KERR_EXACT_DIMENSION; ++axis) {
      CHECK(close_scaled(negative[component].gradient[axis],
                         -parity * positive[component].gradient[axis], 2e-12),
            "metric derivative has differentiated inversion parity");
    }
  }

  CHECK(positive[KERR_EXACT_TZ].value == 0.0,
        "g_tz is represented as exact zero");
  for (size_t axis = 0; axis < KERR_EXACT_DIMENSION; ++axis) {
    CHECK(positive[KERR_EXACT_TZ].gradient[axis] == 0.0,
          "all g_tz derivatives are represented as exact zero");
  }
}

static void test_api_errors_and_names(void) {
  static const char *const names[KERR_EXACT_COMPONENT_COUNT] = {
      "tt", "tx", "ty", "tz", "xx", "xy", "xz", "yy", "yz", "zz"};
  double values[KERR_EXACT_COMPONENT_COUNT];
  double sentinel_values[KERR_EXACT_COMPONENT_COUNT];
  kerr_exact_value_gradient results[KERR_EXACT_COMPONENT_COUNT];
  kerr_exact_value_gradient sentinel_results[KERR_EXACT_COMPONENT_COUNT];

  for (size_t component = 0; component < KERR_EXACT_COMPONENT_COUNT;
       ++component) {
    values[component] = 1234.5 + (double)component;
    results[component].value = -4321.5 - (double)component;
    for (size_t axis = 0; axis < KERR_EXACT_DIMENSION; ++axis)
      results[component].gradient[axis] = 91.0 + (double)(3 * component + axis);
    CHECK(strcmp(kerr_exact_component_name((kerr_exact_component)component),
                 names[component]) == 0,
          "component name follows public ordering");
  }
  memcpy(sentinel_values, values, sizeof(values));
  memcpy(sentinel_results, results, sizeof(results));

  CHECK(kerr_exact_metric(1.0, 2.0, 3.0, NULL) ==
            KERR_EXACT_INVALID_ARGUMENT,
        "null value output is rejected");
  CHECK(kerr_exact_metric_gradient(1.0, 2.0, 3.0, NULL) ==
            KERR_EXACT_INVALID_ARGUMENT,
        "null gradient output is rejected");
  CHECK(kerr_exact_metric(NAN, 0.0, 1.0, values) ==
            KERR_EXACT_INVALID_ARGUMENT,
        "NaN coordinate is rejected");
  CHECK(kerr_exact_metric_gradient(INFINITY, 0.0, 1.0, results) ==
            KERR_EXACT_INVALID_ARGUMENT,
        "infinite coordinate is rejected");
  CHECK(kerr_exact_metric(0.0, 0.0, 0.0, values) ==
            KERR_EXACT_DOMAIN_ERROR,
        "puncture is rejected");
  CHECK(kerr_exact_metric_gradient(0.0, 0.0, 0.0, results) ==
            KERR_EXACT_DOMAIN_ERROR,
        "puncture gradient is rejected");
  CHECK(memcmp(values, sentinel_values, sizeof(values)) == 0,
        "value buffer remains unchanged on errors");
  CHECK(memcmp(results, sentinel_results, sizeof(results)) == 0,
        "gradient buffer remains unchanged on errors");

  memcpy(values, sentinel_values, sizeof(values));
  CHECK(kerr_exact_metric(1e-150, 0.0, 0.0, values) ==
            KERR_EXACT_NONFINITE_RESULT,
        "finite coordinates producing double overflow are reported");
  CHECK(memcmp(values, sentinel_values, sizeof(values)) == 0,
        "value buffer remains unchanged after nonfinite result");

  CHECK(kerr_exact_component_name((kerr_exact_component)-1) == NULL,
        "negative component name is rejected");
  CHECK(kerr_exact_component_name((kerr_exact_component)KERR_EXACT_COMPONENT_COUNT) ==
            NULL,
        "past-end component name is rejected");
  CHECK(close_scaled(kerr_exact_throat_radius(),
                     (2.0 + sqrt(3.0)) / 8.0, 2e-16),
        "public throat radius has the documented value");
}

int main(void) {
  test_high_precision_fixtures();
  test_gradients_against_finite_differences();
  test_symmetry_and_zero_channel();
  test_api_errors_and_names();

  if (failures != 0) {
    fprintf(stderr, "%d Kerr exact-evaluator test(s) failed\n", failures);
    return 1;
  }
  printf("All Kerr exact-evaluator tests passed.\n");
  return 0;
}
