#ifndef KERR_ADM_EXACT_H
#define KERR_ADM_EXACT_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

enum {
  KERR_ADM_DIMENSION = 3,
  KERR_ADM_FIELD_COUNT = 10
};

/* Ordering shared by the analytic evaluator and both interpolators. */
typedef enum kerr_adm_field {
  KERR_ADM_ALPHA = 0,
  KERR_ADM_BETAX,
  KERR_ADM_BETAY,
  KERR_ADM_BETAZ,
  KERR_ADM_GAMMAXX,
  KERR_ADM_GAMMAXY,
  KERR_ADM_GAMMAXZ,
  KERR_ADM_GAMMAYY,
  KERR_ADM_GAMMAYZ,
  KERR_ADM_GAMMAZZ
} kerr_adm_field;

typedef enum kerr_adm_exact_status {
  KERR_ADM_EXACT_SUCCESS = 0,
  KERR_ADM_EXACT_INVALID_ARGUMENT,
  KERR_ADM_EXACT_DOMAIN_ERROR,
  KERR_ADM_EXACT_NONFINITE_RESULT
} kerr_adm_exact_status;

typedef struct kerr_adm_value_gradient {
  double value;
  double gradient[KERR_ADM_DIMENSION];
} kerr_adm_value_gradient;

/*
 * Evaluate the Cartesian ADM fields for M=1 and a=1/2.  The lapse is signed:
 * positive on the exterior sheet, negative on the inner sheet, and smooth
 * through the throat.  The origin is outside the coordinate domain.  On
 * failure, the output array is left unchanged.
 */
kerr_adm_exact_status kerr_adm_exact_values(
    double x, double y, double z,
    double values[KERR_ADM_FIELD_COUNT]);

/*
 * Evaluate the fields and their analytic Cartesian gradients, including on
 * the throat r=s.  On failure, results is left unchanged.
 */
kerr_adm_exact_status kerr_adm_exact_value_gradient(
    double x, double y, double z,
    kerr_adm_value_gradient results[KERR_ADM_FIELD_COUNT]);

double kerr_adm_exact_throat_radius(void);
const char *kerr_adm_field_name(kerr_adm_field field);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* KERR_ADM_EXACT_H */
