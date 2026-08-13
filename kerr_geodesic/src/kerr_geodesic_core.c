#include "kerr_geodesic_core.h"

#include <math.h>
#include <stddef.h>
#include <string.h>

enum {
  ADM_ALPHA = 0,
  ADM_BETAX,
  ADM_BETAY,
  ADM_BETAZ,
  ADM_GAMMAXX,
  ADM_GAMMAXY,
  ADM_GAMMAXZ,
  ADM_GAMMAYY,
  ADM_GAMMAYZ,
  ADM_GAMMAZZ
};

static int finite_state(const kerr_geodesic_state *state) {
  size_t axis;
  for (axis = 0; axis < KERR_GEODESIC_DIMENSION; ++axis) {
    if (!isfinite(state->position[axis]) ||
        !isfinite(state->momentum[axis]))
      return 0;
  }
  return 1;
}

static int finite_adm(const kerr_geodesic_adm *adm) {
  size_t i;
  size_t j;
  size_t k;
  if (!isfinite(adm->alpha)) return 0;
  for (i = 0; i < KERR_GEODESIC_DIMENSION; ++i) {
    if (!isfinite(adm->beta[i]) || !isfinite(adm->d_alpha[i])) return 0;
    for (j = 0; j < KERR_GEODESIC_DIMENSION; ++j) {
      if (!isfinite(adm->gamma[i][j]) || !isfinite(adm->d_beta[i][j]))
        return 0;
      for (k = 0; k < KERR_GEODESIC_DIMENSION; ++k) {
        if (!isfinite(adm->d_gamma[i][j][k])) return 0;
      }
    }
  }
  return 1;
}

/* Solve gamma q=p by Cholesky factorization, checking positive definiteness. */
static kerr_geodesic_status solve_spatial_metric(
    const double gamma[KERR_GEODESIC_DIMENSION][KERR_GEODESIC_DIMENSION],
    const double momentum[KERR_GEODESIC_DIMENSION],
    double solution[KERR_GEODESIC_DIMENSION]) {
  double lower[KERR_GEODESIC_DIMENSION][KERR_GEODESIC_DIMENSION] = {{0.0}};
  double intermediate[KERR_GEODESIC_DIMENSION];
  size_t i;
  size_t j;
  size_t k;

  for (i = 0; i < KERR_GEODESIC_DIMENSION; ++i) {
    for (j = 0; j <= i; ++j) {
      double value = gamma[i][j];
      if (!isfinite(value) ||
          fabs(gamma[i][j] - gamma[j][i]) >
              64.0 * 0x1p-52 * fmax(1.0, fabs(value)))
        return KERR_GEODESIC_NON_POSITIVE_METRIC;
      for (k = 0; k < j; ++k) value -= lower[i][k] * lower[j][k];
      if (i == j) {
        if (!(value > 0.0) || !isfinite(value))
          return KERR_GEODESIC_NON_POSITIVE_METRIC;
        lower[i][j] = sqrt(value);
      } else {
        lower[i][j] = value / lower[j][j];
      }
    }
  }
  for (i = 0; i < KERR_GEODESIC_DIMENSION; ++i) {
    double value = momentum[i];
    for (j = 0; j < i; ++j) value -= lower[i][j] * intermediate[j];
    intermediate[i] = value / lower[i][i];
  }
  for (i = KERR_GEODESIC_DIMENSION; i-- > 0;) {
    double value = intermediate[i];
    for (j = i + 1; j < KERR_GEODESIC_DIMENSION; ++j)
      value -= lower[j][i] * solution[j];
    solution[i] = value / lower[i][i];
  }
  for (i = 0; i < KERR_GEODESIC_DIMENSION; ++i) {
    if (!isfinite(solution[i])) return KERR_GEODESIC_NONFINITE_RESULT;
  }
  return KERR_GEODESIC_SUCCESS;
}

kerr_geodesic_status kerr_geodesic_adm_from_ordered_fields(
    const double values[KERR_GEODESIC_ADM_FIELD_COUNT],
    const double gradients[KERR_GEODESIC_ADM_FIELD_COUNT *
                           KERR_GEODESIC_DIMENSION],
    kerr_geodesic_adm *adm) {
  static const size_t gamma_field[KERR_GEODESIC_DIMENSION]
                                 [KERR_GEODESIC_DIMENSION] = {
      {ADM_GAMMAXX, ADM_GAMMAXY, ADM_GAMMAXZ},
      {ADM_GAMMAXY, ADM_GAMMAYY, ADM_GAMMAYZ},
      {ADM_GAMMAXZ, ADM_GAMMAYZ, ADM_GAMMAZZ}};
  kerr_geodesic_adm temporary;
  size_t axis;
  size_t i;
  size_t j;

  if (values == NULL || gradients == NULL || adm == NULL)
    return KERR_GEODESIC_INVALID_ARGUMENT;
  temporary.alpha = values[ADM_ALPHA];
  for (i = 0; i < KERR_GEODESIC_DIMENSION; ++i) {
    temporary.beta[i] = values[ADM_BETAX + i];
    temporary.d_alpha[i] =
        gradients[(size_t)ADM_ALPHA * KERR_GEODESIC_DIMENSION + i];
  }
  for (axis = 0; axis < KERR_GEODESIC_DIMENSION; ++axis) {
    for (i = 0; i < KERR_GEODESIC_DIMENSION; ++i) {
      temporary.d_beta[axis][i] =
          gradients[((size_t)ADM_BETAX + i) * KERR_GEODESIC_DIMENSION + axis];
      for (j = 0; j < KERR_GEODESIC_DIMENSION; ++j) {
        const size_t field = gamma_field[i][j];
        temporary.gamma[i][j] = values[field];
        temporary.d_gamma[axis][i][j] =
            gradients[field * KERR_GEODESIC_DIMENSION + axis];
      }
    }
  }
  if (!finite_adm(&temporary)) return KERR_GEODESIC_NONFINITE_RESULT;
  *adm = temporary;
  return KERR_GEODESIC_SUCCESS;
}

kerr_geodesic_status kerr_geodesic_rhs_from_adm(
    const kerr_geodesic_adm *adm, const kerr_geodesic_state *state,
    kerr_geodesic_state *rhs, double *hamiltonian) {
  kerr_geodesic_state temporary;
  double q[KERR_GEODESIC_DIMENSION];
  double momentum_squared = 0.0;
  double shift_momentum = 0.0;
  double W;
  size_t i;
  size_t j;
  size_t k;
  kerr_geodesic_status status;

  if (adm == NULL || state == NULL || rhs == NULL || !finite_adm(adm) ||
      !finite_state(state))
    return KERR_GEODESIC_INVALID_ARGUMENT;
  status = solve_spatial_metric(adm->gamma, state->momentum, q);
  if (status != KERR_GEODESIC_SUCCESS) return status;
  for (i = 0; i < KERR_GEODESIC_DIMENSION; ++i) {
    momentum_squared += state->momentum[i] * q[i];
    shift_momentum += adm->beta[i] * state->momentum[i];
  }
  if (!(momentum_squared > -1.0) || !isfinite(momentum_squared))
    return KERR_GEODESIC_NONFINITE_RESULT;
  W = sqrt(1.0 + momentum_squared);
  if (!(W > 0.0) || !isfinite(W)) return KERR_GEODESIC_NONFINITE_RESULT;

  for (i = 0; i < KERR_GEODESIC_DIMENSION; ++i) {
    double metric_term = 0.0;
    double shift_term = 0.0;
    temporary.position[i] = adm->alpha * q[i] / W - adm->beta[i];
    for (j = 0; j < KERR_GEODESIC_DIMENSION; ++j) {
      shift_term += state->momentum[j] * adm->d_beta[i][j];
      for (k = 0; k < KERR_GEODESIC_DIMENSION; ++k)
        metric_term += q[j] * q[k] * adm->d_gamma[i][j][k];
    }
    temporary.momentum[i] = -W * adm->d_alpha[i] +
                            adm->alpha * metric_term / (2.0 * W) + shift_term;
  }
  if (!finite_state(&temporary)) return KERR_GEODESIC_NONFINITE_RESULT;
  if (hamiltonian != NULL) {
    const double value = adm->alpha * W - shift_momentum;
    if (!isfinite(value)) return KERR_GEODESIC_NONFINITE_RESULT;
    *hamiltonian = value;
  }
  *rhs = temporary;
  return KERR_GEODESIC_SUCCESS;
}

kerr_geodesic_status kerr_geodesic_rhs(
    kerr_geodesic_adm_evaluator evaluator, void *context,
    const kerr_geodesic_state *state, kerr_geodesic_state *rhs,
    double *hamiltonian) {
  kerr_geodesic_adm adm;
  if (evaluator == NULL || state == NULL || rhs == NULL)
    return KERR_GEODESIC_INVALID_ARGUMENT;
  if (!evaluator(context, state->position, &adm))
    return KERR_GEODESIC_BACKEND_FAILURE;
  return kerr_geodesic_rhs_from_adm(&adm, state, rhs, hamiltonian);
}

static kerr_geodesic_state state_plus_scaled_rhs(
    const kerr_geodesic_state *state, const kerr_geodesic_state *rhs,
    double scale) {
  kerr_geodesic_state result;
  size_t axis;
  for (axis = 0; axis < KERR_GEODESIC_DIMENSION; ++axis) {
    result.position[axis] = state->position[axis] + scale * rhs->position[axis];
    result.momentum[axis] = state->momentum[axis] + scale * rhs->momentum[axis];
  }
  return result;
}

kerr_geodesic_status kerr_geodesic_rk4_step(
    kerr_geodesic_adm_evaluator evaluator, void *context, double step,
    kerr_geodesic_state *state) {
  kerr_geodesic_state k1;
  kerr_geodesic_state k2;
  kerr_geodesic_state k3;
  kerr_geodesic_state k4;
  kerr_geodesic_state stage;
  kerr_geodesic_state result;
  kerr_geodesic_status status;
  size_t axis;

  if (evaluator == NULL || state == NULL || !isfinite(step) || step <= 0.0 ||
      !finite_state(state))
    return KERR_GEODESIC_INVALID_ARGUMENT;
  status = kerr_geodesic_rhs(evaluator, context, state, &k1, NULL);
  if (status != KERR_GEODESIC_SUCCESS) return status;
  stage = state_plus_scaled_rhs(state, &k1, 0.5 * step);
  status = kerr_geodesic_rhs(evaluator, context, &stage, &k2, NULL);
  if (status != KERR_GEODESIC_SUCCESS) return status;
  stage = state_plus_scaled_rhs(state, &k2, 0.5 * step);
  status = kerr_geodesic_rhs(evaluator, context, &stage, &k3, NULL);
  if (status != KERR_GEODESIC_SUCCESS) return status;
  stage = state_plus_scaled_rhs(state, &k3, step);
  status = kerr_geodesic_rhs(evaluator, context, &stage, &k4, NULL);
  if (status != KERR_GEODESIC_SUCCESS) return status;
  for (axis = 0; axis < KERR_GEODESIC_DIMENSION; ++axis) {
    result.position[axis] =
        state->position[axis] +
        step * (k1.position[axis] + 2.0 * k2.position[axis] +
                2.0 * k3.position[axis] + k4.position[axis]) /
            6.0;
    result.momentum[axis] =
        state->momentum[axis] +
        step * (k1.momentum[axis] + 2.0 * k2.momentum[axis] +
                2.0 * k3.momentum[axis] + k4.momentum[axis]) /
            6.0;
  }
  if (!finite_state(&result)) return KERR_GEODESIC_NONFINITE_RESULT;
  *state = result;
  return KERR_GEODESIC_SUCCESS;
}

kerr_geodesic_status kerr_geodesic_initialize_circular_orbit(
    double boyer_lindquist_radius, kerr_geodesic_adm_evaluator exact_evaluator,
    void *exact_context, kerr_geodesic_circular_orbit *orbit) {
  const double mass = 1.0;
  const double spin = 0.5;
  const double throat = (2.0 + sqrt(3.0)) / 8.0;
  kerr_geodesic_circular_orbit temporary;
  kerr_geodesic_adm adm;
  kerr_geodesic_state rhs;
  double velocity[KERR_GEODESIC_DIMENSION] = {0.0, 0.0, 0.0};
  double shifted_velocity[KERR_GEODESIC_DIMENSION];
  double normalization = 0.0;
  double ut;
  size_t i;
  size_t j;
  kerr_geodesic_status status;

  if (orbit == NULL || exact_evaluator == NULL ||
      !isfinite(boyer_lindquist_radius) ||
      !(boyer_lindquist_radius > 4.0 * throat))
    return KERR_GEODESIC_INVALID_ARGUMENT;
  temporary.boyer_lindquist_radius = boyer_lindquist_radius;
  temporary.les_radius =
      (boyer_lindquist_radius - 2.0 * throat +
       sqrt(boyer_lindquist_radius *
            (boyer_lindquist_radius - 4.0 * throat))) /
      2.0;
  temporary.angular_frequency =
      sqrt(mass) /
      (pow(boyer_lindquist_radius, 1.5) + spin * sqrt(mass));
  temporary.period = 2.0 * acos(-1.0) / temporary.angular_frequency;
  memset(&temporary.initial_state, 0, sizeof(temporary.initial_state));
  temporary.initial_state.position[0] = temporary.les_radius;
  velocity[1] = temporary.les_radius * temporary.angular_frequency;
  if (!exact_evaluator(exact_context, temporary.initial_state.position, &adm))
    return KERR_GEODESIC_BACKEND_FAILURE;
  for (i = 0; i < KERR_GEODESIC_DIMENSION; ++i)
    shifted_velocity[i] = velocity[i] + adm.beta[i];
  for (i = 0; i < KERR_GEODESIC_DIMENSION; ++i) {
    for (j = 0; j < KERR_GEODESIC_DIMENSION; ++j)
      normalization +=
          adm.gamma[i][j] * shifted_velocity[i] * shifted_velocity[j];
  }
  normalization = adm.alpha * adm.alpha - normalization;
  if (!(normalization > 0.0) || !isfinite(normalization))
    return KERR_GEODESIC_NONFINITE_RESULT;
  ut = 1.0 / sqrt(normalization);
  for (i = 0; i < KERR_GEODESIC_DIMENSION; ++i) {
    for (j = 0; j < KERR_GEODESIC_DIMENSION; ++j)
      temporary.initial_state.momentum[i] +=
          ut * adm.gamma[i][j] * shifted_velocity[j];
  }
  status = kerr_geodesic_rhs_from_adm(
      &adm, &temporary.initial_state, &rhs, &temporary.energy);
  if (status != KERR_GEODESIC_SUCCESS) return status;
  temporary.angular_momentum =
      kerr_geodesic_angular_momentum(&temporary.initial_state);
  if (!isfinite(temporary.les_radius) ||
      !isfinite(temporary.angular_frequency) || !isfinite(temporary.period) ||
      !isfinite(temporary.energy) || !isfinite(temporary.angular_momentum))
    return KERR_GEODESIC_NONFINITE_RESULT;
  *orbit = temporary;
  return KERR_GEODESIC_SUCCESS;
}

void kerr_geodesic_circular_state(
    const kerr_geodesic_circular_orbit *orbit, double time,
    kerr_geodesic_state *state) {
  const double angle = orbit->angular_frequency * time;
  const double cosine = cos(angle);
  const double sine = sin(angle);
  const kerr_geodesic_state *initial = &orbit->initial_state;
  state->position[0] = orbit->les_radius * cosine;
  state->position[1] = orbit->les_radius * sine;
  state->position[2] = 0.0;
  state->momentum[0] =
      initial->momentum[0] * cosine - initial->momentum[1] * sine;
  state->momentum[1] =
      initial->momentum[0] * sine + initial->momentum[1] * cosine;
  state->momentum[2] = initial->momentum[2];
}

double kerr_geodesic_angular_momentum(const kerr_geodesic_state *state) {
  return state->position[0] * state->momentum[1] -
         state->position[1] * state->momentum[0];
}

kerr_geodesic_status kerr_geodesic_mass_shell_residual(
    const kerr_geodesic_adm *adm, const kerr_geodesic_state *state,
    double fixed_pt, double *residual) {
  double q[KERR_GEODESIC_DIMENSION];
  double momentum_squared = 0.0;
  double shifted = fixed_pt;
  double value;
  size_t axis;
  kerr_geodesic_status status;
  if (adm == NULL || state == NULL || residual == NULL || !isfinite(fixed_pt) ||
      !finite_adm(adm) || !finite_state(state) || adm->alpha == 0.0)
    return KERR_GEODESIC_INVALID_ARGUMENT;
  status = solve_spatial_metric(adm->gamma, state->momentum, q);
  if (status != KERR_GEODESIC_SUCCESS) return status;
  for (axis = 0; axis < KERR_GEODESIC_DIMENSION; ++axis) {
    momentum_squared += state->momentum[axis] * q[axis];
    shifted -= adm->beta[axis] * state->momentum[axis];
  }
  value = momentum_squared - shifted * shifted / (adm->alpha * adm->alpha) + 1.0;
  if (!isfinite(value)) return KERR_GEODESIC_NONFINITE_RESULT;
  *residual = value;
  return KERR_GEODESIC_SUCCESS;
}

void kerr_geodesic_phase_initialize(kerr_geodesic_phase *phase,
                                    const kerr_geodesic_state *state) {
  phase->previous_x = state->position[0];
  phase->previous_y = state->position[1];
  phase->angle = atan2(state->position[1], state->position[0]);
  phase->initialized = 1;
}

double kerr_geodesic_phase_update(kerr_geodesic_phase *phase,
                                  const kerr_geodesic_state *state) {
  double dot;
  double cross;
  if (!phase->initialized) {
    kerr_geodesic_phase_initialize(phase, state);
    return phase->angle;
  }
  dot = phase->previous_x * state->position[0] +
        phase->previous_y * state->position[1];
  cross = phase->previous_x * state->position[1] -
          phase->previous_y * state->position[0];
  phase->angle += atan2(cross, dot);
  phase->previous_x = state->position[0];
  phase->previous_y = state->position[1];
  return phase->angle;
}

const char *kerr_geodesic_status_name(kerr_geodesic_status status) {
  switch (status) {
    case KERR_GEODESIC_SUCCESS:
      return "success";
    case KERR_GEODESIC_INVALID_ARGUMENT:
      return "invalid argument";
    case KERR_GEODESIC_BACKEND_FAILURE:
      return "ADM backend failure";
    case KERR_GEODESIC_NON_POSITIVE_METRIC:
      return "spatial metric is not positive definite";
    case KERR_GEODESIC_NONFINITE_RESULT:
      return "nonfinite geodesic result";
  }
  return "unknown status";
}
