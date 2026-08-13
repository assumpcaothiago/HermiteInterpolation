#ifndef LAGRANGE3D_REFERENCE_H
#define LAGRANGE3D_REFERENCE_H

#include "kerr_adm_exact.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum lagrange3d_reference_status {
  LAGRANGE3D_REFERENCE_SUCCESS = 0,
  LAGRANGE3D_REFERENCE_INVALID_ARGUMENT,
  LAGRANGE3D_REFERENCE_INVALID_GRID,
  LAGRANGE3D_REFERENCE_STENCIL_UNAVAILABLE,
  LAGRANGE3D_REFERENCE_UPSTREAM_ERROR
} lagrange3d_reference_status;

/*
 * View of the shared ADM grid required by the copied production routine.
 * The copied routine uses int dimensions and three explicit coordinate arrays.
 */
typedef struct lagrange3d_reference_grid {
  int dimension[3];
  double spacing[3];
  double *coordinate[3];
  const double *function[KERR_ADM_FIELD_COUNT];
} lagrange3d_reference_grid;

/* Seven samples per axis (degree-six tensor Lagrange interpolation). */
lagrange3d_reference_status lagrange3d_reference_interpolate(
    const lagrange3d_reference_grid *grid, double xx0, double xx1, double xx2,
    double values[KERR_ADM_FIELD_COUNT]);

const char *lagrange3d_reference_status_name(
    lagrange3d_reference_status status);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* LAGRANGE3D_REFERENCE_H */
