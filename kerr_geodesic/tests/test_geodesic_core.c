#include "kerr_adm_exact.h"
#include "kerr_geodesic_core.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static size_t checks = 0;

#define CHECK(condition)                                                        \
  do {                                                                          \
    ++checks;                                                                    \
    if (!(condition)) {                                                          \
      fprintf(stderr, "check failed at %s:%d: %s\n", __FILE__, __LINE__,       \
              #condition);                                                       \
      return 0;                                                                  \
    }                                                                            \
  } while (0)

static int close_scaled(double actual, double expected, double tolerance) {
  return fabs(actual - expected) <=
         tolerance * fmax(1.0, fmax(fabs(actual), fabs(expected)));
}

static int exact_backend(void *context, const double position[3],
                         kerr_geodesic_adm *adm) {
  kerr_adm_value_gradient fields[KERR_ADM_FIELD_COUNT];
  double values[KERR_ADM_FIELD_COUNT];
  double gradients[KERR_ADM_FIELD_COUNT * 3];
  size_t field;
  size_t axis;
  (void)context;
  if (kerr_adm_exact_value_gradient(position[0], position[1], position[2],
                                    fields) != KERR_ADM_EXACT_SUCCESS)
    return 0;
  for (field = 0; field < KERR_ADM_FIELD_COUNT; ++field) {
    values[field] = fields[field].value;
    for (axis = 0; axis < 3; ++axis)
      gradients[field * 3 + axis] = fields[field].gradient[axis];
  }
  return kerr_geodesic_adm_from_ordered_fields(values, gradients, adm) ==
         KERR_GEODESIC_SUCCESS;
}

static int failing_backend(void *context, const double position[3],
                           kerr_geodesic_adm *adm) {
  (void)context;
  (void)position;
  (void)adm;
  return 0;
}

static int test_field_mapping_and_flat_rhs(void) {
  double values[10] = {1.0, 0.0, 0.0, 0.0, 1.0,
                       0.0, 0.0, 1.0, 0.0, 1.0};
  double gradients[30] = {0.0};
  kerr_geodesic_adm adm;
  kerr_geodesic_state state = {{1.0, 2.0, 3.0}, {0.3, -0.4, 0.5}};
  kerr_geodesic_state rhs;
  const kerr_geodesic_state sentinel = {{-9.0, -8.0, -7.0},
                                        {-6.0, -5.0, -4.0}};
  double hamiltonian;
  double W = sqrt(1.0 + 0.3 * 0.3 + 0.4 * 0.4 + 0.5 * 0.5);
  size_t axis;

  CHECK(kerr_geodesic_adm_from_ordered_fields(values, gradients, &adm) ==
        KERR_GEODESIC_SUCCESS);
  CHECK(kerr_geodesic_rhs_from_adm(&adm, &state, &rhs, &hamiltonian) ==
        KERR_GEODESIC_SUCCESS);
  CHECK(close_scaled(hamiltonian, W, 2e-15));
  for (axis = 0; axis < 3; ++axis) {
    CHECK(close_scaled(rhs.position[axis], state.momentum[axis] / W, 2e-15));
    CHECK(rhs.momentum[axis] == 0.0);
  }
  adm.gamma[2][2] = -1.0;
  rhs = sentinel;
  CHECK(kerr_geodesic_rhs_from_adm(&adm, &state, &rhs, NULL) ==
        KERR_GEODESIC_NON_POSITIVE_METRIC);
  CHECK(memcmp(&rhs, &sentinel, sizeof(rhs)) == 0);
  values[0] = NAN;
  memset(&adm, 0x5a, sizeof(adm));
  {
    const kerr_geodesic_adm adm_sentinel = adm;
    CHECK(kerr_geodesic_adm_from_ordered_fields(values, gradients, &adm) ==
          KERR_GEODESIC_NONFINITE_RESULT);
    CHECK(memcmp(&adm, &adm_sentinel, sizeof(adm)) == 0);
  }
  CHECK(kerr_geodesic_adm_from_ordered_fields(NULL, gradients, &adm) ==
        KERR_GEODESIC_INVALID_ARGUMENT);
  return 1;
}

static int test_circular_data_and_invariants(void) {
  kerr_geodesic_circular_orbit orbit;
  kerr_geodesic_state rhs;
  kerr_geodesic_state closed;
  kerr_geodesic_adm adm;
  double hamiltonian;
  double residual;
  const double radius_bl = 5.0;
  const double spin = 0.5;
  const double sqrt_radius = sqrt(radius_bl);
  const double common =
      pow(radius_bl, 0.75) *
      sqrt(pow(radius_bl, 1.5) - 3.0 * sqrt_radius +
           2.0 * spin);
  const double expected_energy =
      (pow(radius_bl, 1.5) - 2.0 * sqrt_radius + spin) / common;
  const double expected_lz =
      (radius_bl * radius_bl - 2.0 * spin * sqrt_radius +
       spin * spin) /
      common;

  CHECK(kerr_geodesic_initialize_circular_orbit(5.0, exact_backend, NULL,
                                                &orbit) ==
        KERR_GEODESIC_SUCCESS);
  CHECK(close_scaled(orbit.les_radius, 4.012753168487724, 2e-15));
  CHECK(close_scaled(orbit.angular_frequency, 0.08561394699397955, 2e-15));
  CHECK(close_scaled(orbit.energy, expected_energy, 3e-14));
  CHECK(close_scaled(orbit.angular_momentum, expected_lz, 3e-14));
  CHECK(orbit.initial_state.position[0] >
        8.0 * kerr_adm_exact_throat_radius());
  CHECK(kerr_geodesic_rhs(exact_backend, NULL, &orbit.initial_state, &rhs,
                          &hamiltonian) == KERR_GEODESIC_SUCCESS);
  CHECK(close_scaled(rhs.position[0], 0.0, 2e-12));
  CHECK(close_scaled(rhs.position[1],
                     orbit.les_radius * orbit.angular_frequency, 2e-12));
  CHECK(close_scaled(rhs.position[2], 0.0, 2e-12));
  CHECK(close_scaled(rhs.momentum[0],
                     -orbit.angular_frequency *
                         orbit.initial_state.momentum[1],
                     2e-11));
  CHECK(close_scaled(rhs.momentum[1], 0.0, 2e-11));
  CHECK(close_scaled(rhs.momentum[2], 0.0, 2e-11));
  CHECK(close_scaled(hamiltonian, orbit.energy, 2e-14));
  CHECK(close_scaled(kerr_geodesic_angular_momentum(&orbit.initial_state),
                     orbit.angular_momentum, 2e-15));
  CHECK(exact_backend(NULL, orbit.initial_state.position, &adm));
  CHECK(kerr_geodesic_mass_shell_residual(
            &adm, &orbit.initial_state, -orbit.energy, &residual) ==
        KERR_GEODESIC_SUCCESS);
  CHECK(fabs(residual) < 2e-14);
  kerr_geodesic_circular_state(&orbit, 0.25 * orbit.period, &closed);
  CHECK(close_scaled(closed.position[0], 0.0, 2e-14));
  CHECK(close_scaled(closed.position[1], orbit.les_radius, 2e-14));
  CHECK(close_scaled(closed.momentum[0],
                     -orbit.initial_state.momentum[1], 2e-14));
  CHECK(close_scaled(closed.momentum[1], 0.0, 2e-14));
  return 1;
}

static double integrate_one_orbit(size_t steps) {
  kerr_geodesic_circular_orbit orbit;
  kerr_geodesic_state state;
  kerr_geodesic_state exact;
  double error = 0.0;
  size_t step;
  size_t axis;
  if (kerr_geodesic_initialize_circular_orbit(5.0, exact_backend, NULL,
                                              &orbit) !=
      KERR_GEODESIC_SUCCESS)
    return NAN;
  state = orbit.initial_state;
  for (step = 0; step < steps; ++step) {
    if (kerr_geodesic_rk4_step(exact_backend, NULL,
                               orbit.period / (double)steps, &state) !=
        KERR_GEODESIC_SUCCESS)
      return NAN;
  }
  kerr_geodesic_circular_state(&orbit, orbit.period, &exact);
  for (axis = 0; axis < 3; ++axis) {
    const double dx = state.position[axis] - exact.position[axis];
    const double dp = state.momentum[axis] - exact.momentum[axis];
    error += dx * dx + dp * dp;
  }
  return sqrt(error);
}

static int test_rk4_and_phase(void) {
  const double coarse = integrate_one_orbit(64);
  const double medium = integrate_one_orbit(128);
  const double fine = integrate_one_orbit(256);
  kerr_geodesic_phase phase = {0.0, 0.0, 0.0, 0};
  kerr_geodesic_state state = {{1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  kerr_geodesic_state unchanged = state;

  CHECK(isfinite(coarse) && isfinite(medium) && isfinite(fine));
  CHECK(log(coarse / medium) / log(2.0) > 3.8);
  CHECK(log(medium / fine) / log(2.0) > 3.8);
  kerr_geodesic_phase_initialize(&phase, &state);
  state.position[0] = 0.0;
  state.position[1] = 1.0;
  CHECK(close_scaled(kerr_geodesic_phase_update(&phase, &state),
                     0.5 * acos(-1.0), 2e-15));
  state.position[0] = -1.0;
  state.position[1] = 0.0;
  CHECK(close_scaled(kerr_geodesic_phase_update(&phase, &state), acos(-1.0),
                     2e-15));
  state.position[0] = 0.0;
  state.position[1] = -1.0;
  CHECK(close_scaled(kerr_geodesic_phase_update(&phase, &state),
                     1.5 * acos(-1.0), 2e-15));
  CHECK(kerr_geodesic_rk4_step(failing_backend, NULL, 0.1, &unchanged) ==
        KERR_GEODESIC_BACKEND_FAILURE);
  CHECK(memcmp(&unchanged, &(kerr_geodesic_state){{1.0, 0.0, 0.0},
                                                  {0.0, 0.0, 0.0}},
               sizeof(unchanged)) == 0);
  CHECK(kerr_geodesic_rk4_step(exact_backend, NULL, 0.0, &unchanged) ==
        KERR_GEODESIC_INVALID_ARGUMENT);
  CHECK(memcmp(&unchanged, &(kerr_geodesic_state){{1.0, 0.0, 0.0},
                                                  {0.0, 0.0, 0.0}},
               sizeof(unchanged)) == 0);
  {
    kerr_geodesic_circular_orbit orbit;
    CHECK(kerr_geodesic_initialize_circular_orbit(
              1.0, exact_backend, NULL, &orbit) ==
          KERR_GEODESIC_INVALID_ARGUMENT);
  }
  return 1;
}

int main(void) {
  if (!test_field_mapping_and_flat_rhs() ||
      !test_circular_data_and_invariants() || !test_rk4_and_phase())
    return EXIT_FAILURE;
  printf("kerr geodesic core: %zu checks passed\n", checks);
  return EXIT_SUCCESS;
}
