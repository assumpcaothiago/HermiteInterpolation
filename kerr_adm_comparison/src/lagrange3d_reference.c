#include "lagrange3d_reference.h"

#include <math.h>
#include <stddef.h>

/* Public symbol supplied by the unmodified vendored production source. */
int interpolation_3d_general__uniform_src_grid(
    int n_interp_ghosts, double src_dxx0, double src_dxx1, double src_dxx2,
    int src_nxx0, int src_nxx1, int src_nxx2, int num_functions,
    double *restrict src_coordinates[],
    const double *restrict src_functions[], int num_destination_points,
    const double destination_coordinates[][3],
    double *restrict destination_data[]);

lagrange3d_reference_status lagrange3d_reference_interpolate(
    const lagrange3d_reference_grid *grid, double xx0, double xx1, double xx2,
    double values[KERR_ADM_FIELD_COUNT]) {
  double temporary[KERR_ADM_FIELD_COUNT];
  double *destinations[KERR_ADM_FIELD_COUNT];
  double destination_coordinate[1][3];
  double *coordinates[3];
  const double *functions[KERR_ADM_FIELD_COUNT];
  int upstream_status;

  if (grid == NULL || values == NULL || !isfinite(xx0) || !isfinite(xx1) ||
      !isfinite(xx2)) {
    return LAGRANGE3D_REFERENCE_INVALID_ARGUMENT;
  }
  for (size_t axis = 0; axis < KERR_ADM_DIMENSION; ++axis) {
    if (grid->dimension[axis] < 7 || !isfinite(grid->spacing[axis]) ||
        grid->spacing[axis] <= 0.0 || grid->coordinate[axis] == NULL) {
      return LAGRANGE3D_REFERENCE_INVALID_GRID;
    }
    coordinates[axis] = grid->coordinate[axis];
  }
  for (size_t field = 0; field < KERR_ADM_FIELD_COUNT; ++field) {
    if (grid->function[field] == NULL)
      return LAGRANGE3D_REFERENCE_INVALID_ARGUMENT;
    functions[field] = grid->function[field];
    destinations[field] = &temporary[field];
  }

  destination_coordinate[0][0] = xx0;
  destination_coordinate[0][1] = xx1;
  destination_coordinate[0][2] = xx2;
  upstream_status = interpolation_3d_general__uniform_src_grid(
      3, grid->spacing[0], grid->spacing[1], grid->spacing[2],
      grid->dimension[0], grid->dimension[1], grid->dimension[2],
      KERR_ADM_FIELD_COUNT, coordinates, functions, 1,
      (const double(*)[3])destination_coordinate, destinations);
  if (upstream_status == 3)
    return LAGRANGE3D_REFERENCE_STENCIL_UNAVAILABLE;
  if (upstream_status != 0) return LAGRANGE3D_REFERENCE_UPSTREAM_ERROR;

  for (size_t field = 0; field < KERR_ADM_FIELD_COUNT; ++field) {
    if (!isfinite(temporary[field]))
      return LAGRANGE3D_REFERENCE_UPSTREAM_ERROR;
  }
  for (size_t field = 0; field < KERR_ADM_FIELD_COUNT; ++field)
    values[field] = temporary[field];
  return LAGRANGE3D_REFERENCE_SUCCESS;
}

const char *lagrange3d_reference_status_name(
    lagrange3d_reference_status status) {
  switch (status) {
    case LAGRANGE3D_REFERENCE_SUCCESS:
      return "success";
    case LAGRANGE3D_REFERENCE_INVALID_ARGUMENT:
      return "invalid argument";
    case LAGRANGE3D_REFERENCE_INVALID_GRID:
      return "invalid grid";
    case LAGRANGE3D_REFERENCE_STENCIL_UNAVAILABLE:
      return "stencil unavailable";
    case LAGRANGE3D_REFERENCE_UPSTREAM_ERROR:
      return "upstream error";
  }
  return "unknown status";
}
