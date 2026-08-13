#include "kerr_adm_noise.h"

#include <math.h>
#include <stdint.h>
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

static void test_regression_values(void) {
  const uint64_t seed = UINT64_C(0x4e4f49534541444d);
  const double first = kerr_adm_noise_sample(seed, 200, 0, 1, 2, 3);
  const double different_field =
      kerr_adm_noise_sample(seed, 200, 1, 1, 2, 3);
  const double different_resolution =
      kerr_adm_noise_sample(seed, 400, 0, 1, 2, 3);
  const double different_node =
      kerr_adm_noise_sample(seed, 200, 0, 2, 2, 3);
  const double different_seed =
      kerr_adm_noise_sample(UINT64_C(123), 200, 0, 1, 2, 3);

  check(first == -0x1.8f0c12f33f166p-1, "fixed noise regression value");
  check(different_field == 0x1.961eab44fed4ep-1,
        "field-keyed regression value");
  check(different_resolution == -0x1.516b705bc9b38p-3,
        "resolution-keyed regression value");
  check(different_seed == -0x1.cfb012c4d33aep-1,
        "seed-keyed regression value");
  check(first == kerr_adm_noise_sample(seed, 200, 0, 1, 2, 3),
        "noise sample is repeatable");
  check(first != different_field, "fields receive independent samples");
  check(first != different_node, "grid nodes receive independent samples");
  check(first != different_resolution,
        "resolutions receive distinct realizations");
  check(first != different_seed, "noise seed changes realization");
}

static void test_range_and_distribution(void) {
  const uint64_t seed = UINT64_C(0x4e4f49534541444d);
  long double sum_squares = 0.0L;
  const size_t samples = 100000;

  for (size_t index = 0; index < samples; ++index) {
    const double sample =
        kerr_adm_noise_sample(seed, 64, index % 10, index, index / 7,
                              index / 97);
    check(sample > -1.0 && sample < 1.0, "noise sample lies in (-1,1)");
    sum_squares += (long double)sample * (long double)sample;
  }
  check(fabsl(sqrtl(sum_squares / (long double)samples) -
              1.0L / sqrtl(3.0L)) < 0.005L,
        "noise RMS is consistent with a uniform distribution");
}

static void test_application(void) {
  const uint64_t seed = UINT64_C(0x4e4f49534541444d);
  const double value = 3.25;
  const double epsilon = 1.0e-8;
  const double perturbed =
      kerr_adm_apply_noise(value, epsilon, seed, 200, 4, 5, 6, 7);
  const double relative = perturbed / value - 1.0;
  const double unchanged =
      kerr_adm_apply_noise(value, 0.0, seed, 200, 4, 5, 6, 7);
  const double zero =
      kerr_adm_apply_noise(0.0, epsilon, seed, 200, 3, 5, 6, 7);

  check(fabs(relative) <= epsilon * (1.0 + 8.0e-8),
        "relative perturbation is bounded by epsilon");
  check(memcmp(&unchanged, &value, sizeof(value)) == 0,
        "zero epsilon preserves the value bit-for-bit");
  check(zero == 0.0 && !signbit(zero),
        "multiplicative noise preserves an exact positive zero");
}

int main(void) {
  test_regression_values();
  test_range_and_distribution();
  test_application();
  if (failures != 0) {
    fprintf(stderr, "%d of %d noise checks failed\n", failures, checks);
    return 1;
  }
  printf("Kerr ADM noise: %d checks passed\n", checks);
  return 0;
}
