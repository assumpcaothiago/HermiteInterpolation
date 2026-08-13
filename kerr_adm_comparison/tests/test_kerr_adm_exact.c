#include "kerr_adm_exact.h"
#include "kerr_adm_exact_fixtures.h"
#include "kerr_exact.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static int checks = 0;
static int failures = 0;

static void check(int condition, const char *message) {
  ++checks;
  if (!condition) {
    ++failures;
    fprintf(stderr, "FAIL: %s\n", message);
  }
}

static int close_scaled(double actual, double expected, double tolerance) {
  const double scale = fmax(1.0, fabs(expected));
  return isfinite(actual) && fabs(actual - expected) <= tolerance * scale;
}

static size_t gamma_field(size_t row, size_t column) {
  static const size_t map[3][3] = {
      {KERR_ADM_GAMMAXX, KERR_ADM_GAMMAXY, KERR_ADM_GAMMAXZ},
      {KERR_ADM_GAMMAXY, KERR_ADM_GAMMAYY, KERR_ADM_GAMMAYZ},
      {KERR_ADM_GAMMAXZ, KERR_ADM_GAMMAYZ, KERR_ADM_GAMMAZZ}};
  return map[row][column];
}

static size_t spatial_metric_component(size_t row, size_t column) {
  return gamma_field(row, column) - KERR_ADM_GAMMAXX + KERR_EXACT_XX;
}

/*
 * Reconstruct g_(mu nu) and its gradient by the ADM identities.  Comparing
 * these results with the independently generated four-metric evaluator catches
 * field-ordering, shift-index, lapse-sign, and product-rule mistakes.
 */
static void check_adm_reconstruction(
    const kerr_adm_value_gradient adm[KERR_ADM_FIELD_COUNT],
    const kerr_exact_value_gradient metric[KERR_EXACT_COMPONENT_COUNT],
    const char *name) {
  double beta_down[3] = {0.0, 0.0, 0.0};
  double gtt = -adm[KERR_ADM_ALPHA].value * adm[KERR_ADM_ALPHA].value;
  double d_beta_down[3][3] = {{0.0}};
  double d_gtt[3];

  for (size_t row = 0; row < 3; ++row) {
    for (size_t column = 0; column < 3; ++column) {
      const size_t gamma = gamma_field(row, column);
      const size_t beta = KERR_ADM_BETAX + column;
      beta_down[row] += adm[gamma].value * adm[beta].value;
      for (size_t axis = 0; axis < 3; ++axis) {
        d_beta_down[row][axis] +=
            adm[gamma].gradient[axis] * adm[beta].value +
            adm[gamma].value * adm[beta].gradient[axis];
      }
    }
  }
  for (size_t row = 0; row < 3; ++row) {
    const size_t metric_component = KERR_EXACT_TX + row;
    char message[160];
    snprintf(message, sizeof(message), "%s reconstructed g_t%zu", name, row);
    check(close_scaled(beta_down[row], metric[metric_component].value, 3e-12),
          message);
    for (size_t axis = 0; axis < 3; ++axis) {
      snprintf(message, sizeof(message), "%s reconstructed d%zu g_t%zu", name,
               axis, row);
      check(close_scaled(d_beta_down[row][axis],
                         metric[metric_component].gradient[axis], 2e-11),
            message);
    }
  }

  for (size_t row = 0; row < 3; ++row) {
    const size_t beta = KERR_ADM_BETAX + row;
    gtt += beta_down[row] * adm[beta].value;
  }
  for (size_t axis = 0; axis < 3; ++axis) {
    d_gtt[axis] = -2.0 * adm[KERR_ADM_ALPHA].value *
                   adm[KERR_ADM_ALPHA].gradient[axis];
    for (size_t row = 0; row < 3; ++row) {
      const size_t beta = KERR_ADM_BETAX + row;
      d_gtt[axis] += d_beta_down[row][axis] * adm[beta].value +
                     beta_down[row] * adm[beta].gradient[axis];
    }
  }
  check(close_scaled(gtt, metric[KERR_EXACT_TT].value, 3e-12),
        "reconstructed g_tt");
  for (size_t axis = 0; axis < 3; ++axis)
    check(close_scaled(d_gtt[axis], metric[KERR_EXACT_TT].gradient[axis],
                       2e-11),
          "reconstructed gradient g_tt");

  for (size_t row = 0; row < 3; ++row) {
    for (size_t column = row; column < 3; ++column) {
      const size_t adm_component = gamma_field(row, column);
      const size_t metric_component = spatial_metric_component(row, column);
      check(close_scaled(adm[adm_component].value,
                         metric[metric_component].value, 3e-12),
            "spatial metric reconstruction");
      for (size_t axis = 0; axis < 3; ++axis)
        check(close_scaled(adm[adm_component].gradient[axis],
                           metric[metric_component].gradient[axis], 2e-11),
              "spatial metric gradient reconstruction");
    }
  }
}

static void test_fixtures(void) {
  for (size_t fixture_index = 0; fixture_index < kerr_adm_exact_fixture_count;
       ++fixture_index) {
    const kerr_adm_exact_fixture *fixture =
        &kerr_adm_exact_fixtures[fixture_index];
    double values[KERR_ADM_FIELD_COUNT];
    kerr_adm_value_gradient adm[KERR_ADM_FIELD_COUNT];
    const kerr_adm_exact_status value_status = kerr_adm_exact_values(
        fixture->coordinate[0], fixture->coordinate[1], fixture->coordinate[2],
        values);

    check(value_status == KERR_ADM_EXACT_SUCCESS, "fixture value status");
    for (size_t field = 0; field < KERR_ADM_FIELD_COUNT; ++field) {
      char message[160];
      snprintf(message, sizeof(message), "%s fixture value %s", fixture->name,
               kerr_adm_field_name((kerr_adm_field)field));
      check(close_scaled(values[field], fixture->value[field], 3e-13), message);
    }

    if (fixture->has_gradient) {
      kerr_exact_value_gradient metric[KERR_EXACT_COMPONENT_COUNT];
      check(kerr_adm_exact_value_gradient(
                fixture->coordinate[0], fixture->coordinate[1],
                fixture->coordinate[2], adm) == KERR_ADM_EXACT_SUCCESS,
            "fixture ADM gradient status");
      check(kerr_exact_metric_gradient(
                fixture->coordinate[0], fixture->coordinate[1],
                fixture->coordinate[2], metric) == KERR_EXACT_SUCCESS,
            "fixture four-metric gradient status");
      for (size_t field = 0; field < KERR_ADM_FIELD_COUNT; ++field) {
        check(close_scaled(adm[field].value, fixture->value[field], 3e-13),
              "fixture combined value");
        for (size_t axis = 0; axis < KERR_ADM_DIMENSION; ++axis)
          check(close_scaled(adm[field].gradient[axis],
                             fixture->gradient[field][axis], 2e-12),
                "fixture gradient");
      }
      check_adm_reconstruction(adm, metric, fixture->name);
    }
  }
}

static void test_errors_and_names(void) {
  double values[KERR_ADM_FIELD_COUNT];
  double original_values[KERR_ADM_FIELD_COUNT];
  kerr_adm_value_gradient results[KERR_ADM_FIELD_COUNT];
  kerr_adm_value_gradient original_results[KERR_ADM_FIELD_COUNT];
  const double throat = kerr_adm_exact_throat_radius();

  for (size_t field = 0; field < KERR_ADM_FIELD_COUNT; ++field) {
    values[field] = 1234.5 + (double)field;
    results[field].value = -987.0 - (double)(4 * field);
    for (size_t axis = 0; axis < KERR_ADM_DIMENSION; ++axis)
      results[field].gradient[axis] =
          -988.0 - (double)(4 * field + axis);
    check(kerr_adm_field_name((kerr_adm_field)field) != NULL, "field name");
  }
  memcpy(original_values, values, sizeof(values));
  memcpy(original_results, results, sizeof(results));
  check(kerr_adm_field_name((kerr_adm_field)-1) == NULL, "negative field name");
  check(kerr_adm_field_name((kerr_adm_field)KERR_ADM_FIELD_COUNT) == NULL,
        "large field name");
  check(kerr_adm_exact_values(1.0, 0.0, 0.0, NULL) ==
            KERR_ADM_EXACT_INVALID_ARGUMENT,
        "null value output");
  check(kerr_adm_exact_value_gradient(1.0, 0.0, 0.0, NULL) ==
            KERR_ADM_EXACT_INVALID_ARGUMENT,
        "null gradient output");
  check(kerr_adm_exact_values(NAN, 0.0, 1.0, values) ==
            KERR_ADM_EXACT_INVALID_ARGUMENT,
        "nan coordinate");
  check(memcmp(values, original_values, sizeof(values)) == 0,
        "value output unchanged after nan");
  check(kerr_adm_exact_values(0.0, 0.0, 0.0, values) ==
            KERR_ADM_EXACT_DOMAIN_ERROR,
        "origin rejected");
  check(memcmp(values, original_values, sizeof(values)) == 0,
        "value output unchanged after origin");
  check(kerr_adm_exact_values(throat, 0.0, 0.0, values) ==
            KERR_ADM_EXACT_SUCCESS,
        "throat value accepted");
  check(values[KERR_ADM_ALPHA] == 0.0, "throat lapse is zero");
  check(kerr_adm_exact_value_gradient(throat, 0.0, 0.0, results) ==
            KERR_ADM_EXACT_NONDIFFERENTIABLE,
        "throat gradient rejected");
  check(memcmp(results, original_results, sizeof(results)) == 0,
        "gradient output unchanged at throat");
}

int main(void) {
  test_fixtures();
  test_errors_and_names();
  if (failures != 0) {
    fprintf(stderr, "%d of %d Kerr ADM evaluator checks failed\n", failures,
            checks);
    return 1;
  }
  printf("Kerr ADM evaluator: %d checks passed\n", checks);
  return 0;
}
