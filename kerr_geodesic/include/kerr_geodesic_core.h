#ifndef KERR_GEODESIC_CORE_H
#define KERR_GEODESIC_CORE_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

enum {
  KERR_GEODESIC_DIMENSION = 3,
  KERR_GEODESIC_STATE_SIZE = 6,
  KERR_GEODESIC_ADM_FIELD_COUNT = 10
};

typedef enum kerr_geodesic_status {
  KERR_GEODESIC_SUCCESS = 0,
  KERR_GEODESIC_INVALID_ARGUMENT,
  KERR_GEODESIC_BACKEND_FAILURE,
  KERR_GEODESIC_NON_POSITIVE_METRIC,
  KERR_GEODESIC_NONFINITE_RESULT
} kerr_geodesic_status;

typedef struct kerr_geodesic_state {
  double position[KERR_GEODESIC_DIMENSION];
  double momentum[KERR_GEODESIC_DIMENSION];
} kerr_geodesic_state;

/*
 * ADM data and spatial derivatives at one Cartesian point.  Derivative
 * indices come first: d_beta[k][i] is partial_k beta^i and
 * d_gamma[k][i][j] is partial_k gamma_ij.
 */
typedef struct kerr_geodesic_adm {
  double alpha;
  double beta[KERR_GEODESIC_DIMENSION];
  double gamma[KERR_GEODESIC_DIMENSION][KERR_GEODESIC_DIMENSION];
  double d_alpha[KERR_GEODESIC_DIMENSION];
  double d_beta[KERR_GEODESIC_DIMENSION][KERR_GEODESIC_DIMENSION];
  double d_gamma[KERR_GEODESIC_DIMENSION][KERR_GEODESIC_DIMENSION]
                [KERR_GEODESIC_DIMENSION];
} kerr_geodesic_adm;

/* The callback returns nonzero on success and must not modify global state. */
typedef int (*kerr_geodesic_adm_evaluator)(
    void *context, const double position[KERR_GEODESIC_DIMENSION],
    kerr_geodesic_adm *adm);

typedef struct kerr_geodesic_circular_orbit {
  double boyer_lindquist_radius;
  double les_radius;
  double angular_frequency;
  double period;
  kerr_geodesic_state initial_state;
  double energy;
  double angular_momentum;
} kerr_geodesic_circular_orbit;

typedef struct kerr_geodesic_phase {
  double previous_x;
  double previous_y;
  double angle;
  int initialized;
} kerr_geodesic_phase;

/*
 * Convert the repository's shared ADM field ordering
 *
 *   alpha, betax, betay, betaz,
 *   gammaxx, gammaxy, gammaxz, gammayy, gammayz, gammazz
 *
 * into the symmetric matrix representation used by the geodesic equations.
 */
kerr_geodesic_status kerr_geodesic_adm_from_ordered_fields(
    const double values[KERR_GEODESIC_ADM_FIELD_COUNT],
    const double gradients[KERR_GEODESIC_ADM_FIELD_COUNT *
                           KERR_GEODESIC_DIMENSION],
    kerr_geodesic_adm *adm);

kerr_geodesic_status kerr_geodesic_rhs_from_adm(
    const kerr_geodesic_adm *adm, const kerr_geodesic_state *state,
    kerr_geodesic_state *rhs, double *hamiltonian);

kerr_geodesic_status kerr_geodesic_rhs(
    kerr_geodesic_adm_evaluator evaluator, void *context,
    const kerr_geodesic_state *state, kerr_geodesic_state *rhs,
    double *hamiltonian);

/* One classical fourth-order Runge--Kutta step.  State is unchanged on error. */
kerr_geodesic_status kerr_geodesic_rk4_step(
    kerr_geodesic_adm_evaluator evaluator, void *context, double step,
    kerr_geodesic_state *state);

/* Fixed M=1, a=1/2, prograde circular data at the supplied BL radius. */
kerr_geodesic_status kerr_geodesic_initialize_circular_orbit(
    double boyer_lindquist_radius, kerr_geodesic_adm_evaluator exact_evaluator,
    void *exact_context, kerr_geodesic_circular_orbit *orbit);

void kerr_geodesic_circular_state(
    const kerr_geodesic_circular_orbit *orbit, double time,
    kerr_geodesic_state *state);

double kerr_geodesic_angular_momentum(const kerr_geodesic_state *state);

kerr_geodesic_status kerr_geodesic_mass_shell_residual(
    const kerr_geodesic_adm *adm, const kerr_geodesic_state *state,
    double fixed_pt, double *residual);

void kerr_geodesic_phase_initialize(kerr_geodesic_phase *phase,
                                    const kerr_geodesic_state *state);
double kerr_geodesic_phase_update(kerr_geodesic_phase *phase,
                                  const kerr_geodesic_state *state);

const char *kerr_geodesic_status_name(kerr_geodesic_status status);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* KERR_GEODESIC_CORE_H */
