#ifndef KERR_ADM_NOISE_H
#define KERR_ADM_NOISE_H

#include <stddef.h>
#include <stdint.h>

/* SplitMix64's stateless finalizer provides a reproducible 64-bit hash. */
static inline uint64_t kerr_adm_noise_mix(uint64_t value) {
  value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
  return value ^ (value >> 31);
}

static inline uint64_t kerr_adm_noise_combine(uint64_t state,
                                               uint64_t value) {
  return kerr_adm_noise_mix(
      state ^ (kerr_adm_noise_mix(value + UINT64_C(0x9e3779b97f4a7c15)) +
               UINT64_C(0x9e3779b97f4a7c15)));
}

/*
 * Return a deterministic independent sample in (-1,1).  The key contains the
 * resolution as well as every logical sample index, so refinement levels have
 * distinct realizations and OpenMP scheduling cannot affect the result.
 */
static inline double kerr_adm_noise_sample(uint64_t seed, size_t resolution,
                                           size_t field, size_t index0,
                                           size_t index1, size_t index2) {
  uint64_t state = kerr_adm_noise_mix(seed ^ UINT64_C(0x4b45525241444d4e));
  uint64_t bin;
  double unit;

  state = kerr_adm_noise_combine(state, (uint64_t)resolution);
  state = kerr_adm_noise_combine(state, (uint64_t)field);
  state = kerr_adm_noise_combine(state, (uint64_t)index0);
  state = kerr_adm_noise_combine(state, (uint64_t)index1);
  state = kerr_adm_noise_combine(state, (uint64_t)index2);
  bin = state >> 12;
  unit = (double)(2 * bin + 1) * 0x1p-53;
  return 2.0 * unit - 1.0;
}

static inline double kerr_adm_apply_noise(double value, double epsilon,
                                           uint64_t seed, size_t resolution,
                                           size_t field, size_t index0,
                                           size_t index1, size_t index2) {
  if (epsilon == 0.0) return value;
  return value *
         (1.0 + epsilon *
                    kerr_adm_noise_sample(seed, resolution, field, index0,
                                          index1, index2));
}

#endif /* KERR_ADM_NOISE_H */
