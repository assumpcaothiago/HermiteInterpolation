#ifndef KERR_EXACT_H
#define KERR_EXACT_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

enum {
  KERR_EXACT_DIMENSION = 3,
  KERR_EXACT_COMPONENT_COUNT = 10
};

/* Independent components of the symmetric four-metric. */
typedef enum kerr_exact_component {
  KERR_EXACT_TT = 0,
  KERR_EXACT_TX,
  KERR_EXACT_TY,
  KERR_EXACT_TZ,
  KERR_EXACT_XX,
  KERR_EXACT_XY,
  KERR_EXACT_XZ,
  KERR_EXACT_YY,
  KERR_EXACT_YZ,
  KERR_EXACT_ZZ
} kerr_exact_component;

typedef enum kerr_exact_status {
  KERR_EXACT_SUCCESS = 0,
  KERR_EXACT_INVALID_ARGUMENT,
  KERR_EXACT_DOMAIN_ERROR,
  KERR_EXACT_NONFINITE_RESULT
} kerr_exact_status;

typedef struct kerr_exact_value_gradient {
  double value;
  /* gradient[0], gradient[1], gradient[2] are d/dx, d/dy, and d/dz. */
  double gradient[KERR_EXACT_DIMENSION];
} kerr_exact_value_gradient;

/*
 * Evaluate the non-boosted LES Cartesian Kerr metric for M=1 and a=1/2.
 * Components follow kerr_exact_component order.  The origin r=0 is outside
 * the coordinate domain.  On failure, the output array is left unchanged.
 */
kerr_exact_status kerr_exact_metric(
    double x, double y, double z,
    double values[KERR_EXACT_COMPONENT_COUNT]);

/* Evaluate the same components and their analytic Cartesian gradients. */
kerr_exact_status kerr_exact_metric_gradient(
    double x, double y, double z,
    kerr_exact_value_gradient results[KERR_EXACT_COMPONENT_COUNT]);

/* LES coordinate radius of the event-horizon throat. */
double kerr_exact_throat_radius(void);

/* Return a stable short name, or NULL for an invalid component. */
const char *kerr_exact_component_name(kerr_exact_component component);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* KERR_EXACT_H */
