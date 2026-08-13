#include "hermite3d.h"
#include "kerr_adm_exact.h"
#include "kerr_adm_noise.h"
#include "kerr_geodesic_core.h"
#include "lagrange3d_reference.h"

#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum {
  GHOST_LAYERS = 6,
  RHS_DIAGNOSTIC_SAMPLES = 256,
  METHOD_EXACT = 0,
  METHOD_HERMITE = 1,
  METHOD_LAGRANGE = 2,
  METHOD_COUNT = 3,
  METRIC_POSITION = 0,
  METRIC_MOMENTUM,
  METRIC_RADIAL,
  METRIC_VERTICAL,
  METRIC_PHASE,
  METRIC_EXACT_ENERGY,
  METRIC_METHOD_HAMILTONIAN,
  METRIC_ANGULAR_MOMENTUM,
  METRIC_MASS_SHELL,
  METRIC_COUNT
};

static const double ACTIVE_HALF_WIDTH = 5.0;
static const double ORBIT_BL_RADIUS = 5.0;
static const uint64_t DEFAULT_NOISE_SEED =
    UINT64_C(0x4e4f49534541444d);

typedef struct options {
  size_t *resolutions;
  size_t resolution_count;
  size_t orbits;
  size_t steps_per_orbit;
  size_t output_every;
  double noise_epsilon;
  uint64_t noise_seed;
  char *trajectory_csv;
} options;

typedef struct scalar_norm {
  long double sum_squares;
  long double maximum;
  size_t count;
} scalar_norm;

typedef struct norm_result {
  long double rms;
  long double maximum;
} norm_result;

typedef struct method_summary {
  norm_result metric[METRIC_COUNT];
  double final_position_error;
  double final_momentum_error;
  double minimum_radius;
  double backend_seconds;
  double integration_seconds;
  size_t backend_evaluations;
} method_summary;

typedef struct rhs_defect {
  norm_result velocity;
  norm_result force;
} rhs_defect;

typedef struct level_summary {
  size_t resolution;
  size_t dimension;
  double spacing;
  long double memory_mib;
  int grid_threads;
  double grid_seconds;
  method_summary method[2];
  rhs_defect defect[2];
} level_summary;

typedef struct reference_solution {
  kerr_geodesic_state *state;
  double *phase;
  size_t steps;
  method_summary summary;
} reference_solution;

typedef struct grid_failure {
  size_t offset;
  kerr_adm_exact_status status;
} grid_failure;

typedef enum backend_kind {
  BACKEND_EXACT = 0,
  BACKEND_HERMITE,
  BACKEND_LAGRANGE
} backend_kind;

typedef struct backend_context {
  backend_kind kind;
  double throat;
  double half_width;
  double spacing;
  const hermite3d_grid *hermite_grid;
  const double *const *functions;
  const lagrange3d_reference_grid *lagrange_grid;
  char failure[192];
  size_t evaluations;
  double seconds;
} backend_context;

static int add_overflows_size_t(size_t left, size_t right) {
  return left > SIZE_MAX - right;
}

static int multiply_overflows_size_t(size_t left, size_t right) {
  return left != 0 && right > SIZE_MAX / left;
}

static char *copy_string(const char *text) {
  const size_t length = strlen(text) + 1;
  char *result;
  if (length == 0 || multiply_overflows_size_t(length, sizeof(*result)))
    return NULL;
  result = (char *)malloc(length);
  if (result != NULL) memcpy(result, text, length);
  return result;
}

static void norm_add(scalar_norm *norm, long double value) {
  const long double magnitude = fabsl(value);
  norm->sum_squares += magnitude * magnitude;
  if (magnitude > norm->maximum) norm->maximum = magnitude;
  ++norm->count;
}

static norm_result norm_finish(const scalar_norm *norm) {
  norm_result result;
  result.rms = norm->count == 0
                   ? 0.0L
                   : sqrtl(norm->sum_squares / (long double)norm->count);
  result.maximum = norm->maximum;
  return result;
}

static double vector_difference(const double left[3], const double right[3]) {
  double square = 0.0;
  size_t axis;
  for (axis = 0; axis < 3; ++axis) {
    const double difference = left[axis] - right[axis];
    square += difference * difference;
  }
  return sqrt(square);
}

static double radius_of(const double position[3]) {
  return sqrt(position[0] * position[0] + position[1] * position[1] +
              position[2] * position[2]);
}

static const char *method_name(size_t method) {
  static const char *const names[METHOD_COUNT] = {"EXACT", "HERMITE",
                                                  "LAGRANGE+FD"};
  return method < METHOD_COUNT ? names[method] : "UNKNOWN";
}

static const char *exact_status_name(kerr_adm_exact_status status) {
  switch (status) {
    case KERR_ADM_EXACT_SUCCESS:
      return "success";
    case KERR_ADM_EXACT_INVALID_ARGUMENT:
      return "invalid argument";
    case KERR_ADM_EXACT_DOMAIN_ERROR:
      return "domain error";
    case KERR_ADM_EXACT_NONFINITE_RESULT:
      return "nonfinite result";
  }
  return "unknown status";
}

static const char *hermite_status_name(hermite3d_status status) {
  switch (status) {
    case HERMITE3D_SUCCESS:
      return "success";
    case HERMITE3D_INVALID_ARGUMENT:
      return "invalid argument";
    case HERMITE3D_INVALID_GRID:
      return "invalid grid";
    case HERMITE3D_INDEX_OVERFLOW:
      return "index overflow";
    case HERMITE3D_STENCIL_UNAVAILABLE:
      return "stencil unavailable";
  }
  return "unknown status";
}

static void print_usage(FILE *stream, const char *program) {
  fprintf(stream,
          "Usage: %s --resolutions N1,N2,... [options]\n"
          "\n"
          "Required:\n"
          "  --resolutions LIST       strictly increasing even N >= 16\n"
          "\n"
          "Options:\n"
          "  --orbits COUNT           orbital periods to integrate (default 10)\n"
          "  --steps-per-orbit COUNT  RK4 steps per period (default 4096)\n"
          "  --noise-epsilon EPS      relative grid noise in [0,1] (default 0)\n"
          "  --noise-seed UINT64      independent noise seed\n"
          "  --trajectory-csv PATH    write sampled long-form trajectory data\n"
          "  --output-every COUNT     RK4 steps between CSV rows (default 16)\n"
          "  --help                    show this message\n",
          program);
}

static const char *attached_value(const char *argument, const char *name) {
  const size_t length = strlen(name);
  if (strncmp(argument, name, length) == 0 && argument[length] == '=')
    return argument + length + 1;
  return NULL;
}

static int parse_size(const char *text, size_t *value) {
  char *end;
  uintmax_t parsed;
  if (text == NULL || *text == '\0' || *text == '-') return 0;
  errno = 0;
  parsed = strtoumax(text, &end, 10);
  if (errno != 0 || *end != '\0' || parsed == 0 || parsed > SIZE_MAX) return 0;
  *value = (size_t)parsed;
  return 1;
}

static int parse_seed(const char *text, uint64_t *value) {
  char *end;
  uintmax_t parsed;
  if (text == NULL || *text == '\0' || *text == '-') return 0;
  errno = 0;
  parsed = strtoumax(text, &end, 0);
  if (errno != 0 || *end != '\0' || parsed > UINT64_MAX) return 0;
  *value = (uint64_t)parsed;
  return 1;
}

static int parse_epsilon(const char *text, double *value) {
  char *end;
  double parsed;
  if (text == NULL || *text == '\0') return 0;
  errno = 0;
  parsed = strtod(text, &end);
  if (errno != 0 || *end != '\0' || !isfinite(parsed) || parsed < 0.0 ||
      parsed > 1.0)
    return 0;
  *value = parsed;
  return 1;
}

static int parse_resolutions(const char *text, size_t **values,
                             size_t *count) {
  size_t entries = 1;
  size_t index = 0;
  const char *cursor;
  size_t *result;
  char *copy;
  char *token;
  if (text == NULL || *text == '\0') return 0;
  for (cursor = text; *cursor != '\0'; ++cursor)
    if (*cursor == ',') ++entries;
  if (multiply_overflows_size_t(entries, sizeof(*result))) return 0;
  result = (size_t *)malloc(entries * sizeof(*result));
  copy = copy_string(text);
  if (result == NULL || copy == NULL) {
    free(result);
    free(copy);
    return 0;
  }
  token = strtok(copy, ",");
  while (token != NULL && index < entries) {
    if (!parse_size(token, &result[index]) || result[index] < 16 ||
        result[index] % 2 != 0 ||
        (index > 0 && result[index] <= result[index - 1])) {
      free(copy);
      free(result);
      return 0;
    }
    ++index;
    token = strtok(NULL, ",");
  }
  free(copy);
  if (index != entries || token != NULL) {
    free(result);
    return 0;
  }
  *values = result;
  *count = entries;
  return 1;
}

static int parse_options(int argc, char **argv, options *result) {
  const char *resolution_text = NULL;
  const char *orbits_text = NULL;
  const char *steps_text = NULL;
  const char *epsilon_text = NULL;
  const char *noise_seed_text = NULL;
  const char *csv_text = NULL;
  const char *output_every_text = NULL;
  int argument;
  memset(result, 0, sizeof(*result));
  result->orbits = 10;
  result->steps_per_orbit = 4096;
  result->output_every = 16;
  result->noise_seed = DEFAULT_NOISE_SEED;
  for (argument = 1; argument < argc; ++argument) {
    const char *value = NULL;
    const char *name = NULL;
    if (strcmp(argv[argument], "--help") == 0) return 2;
#define OPTION_VALUE(option_name, destination)                                  \
  do {                                                                           \
    value = attached_value(argv[argument], option_name);                          \
    if (strcmp(argv[argument], option_name) == 0) {                               \
      if (++argument >= argc) return 0;                                           \
      value = argv[argument];                                                     \
    }                                                                             \
    if (value != NULL) {                                                          \
      name = option_name;                                                         \
      if ((destination) != NULL) return 0;                                        \
      (destination) = value;                                                      \
    }                                                                             \
  } while (0)
    OPTION_VALUE("--resolutions", resolution_text);
    if (name == NULL) OPTION_VALUE("--orbits", orbits_text);
    if (name == NULL) OPTION_VALUE("--steps-per-orbit", steps_text);
    if (name == NULL) OPTION_VALUE("--noise-epsilon", epsilon_text);
    if (name == NULL) OPTION_VALUE("--noise-seed", noise_seed_text);
    if (name == NULL) OPTION_VALUE("--trajectory-csv", csv_text);
    if (name == NULL) OPTION_VALUE("--output-every", output_every_text);
#undef OPTION_VALUE
    if (name == NULL) return 0;
  }
  if (resolution_text == NULL ||
      !parse_resolutions(resolution_text, &result->resolutions,
                         &result->resolution_count) ||
      (orbits_text != NULL && !parse_size(orbits_text, &result->orbits)) ||
      (steps_text != NULL &&
       !parse_size(steps_text, &result->steps_per_orbit)) ||
      (output_every_text != NULL &&
       !parse_size(output_every_text, &result->output_every)) ||
      (epsilon_text != NULL &&
       !parse_epsilon(epsilon_text, &result->noise_epsilon)) ||
      (noise_seed_text != NULL &&
       !parse_seed(noise_seed_text, &result->noise_seed))) {
    free(result->resolutions);
    result->resolutions = NULL;
    return 0;
  }
  if (csv_text != NULL) {
    result->trajectory_csv = copy_string(csv_text);
    if (result->trajectory_csv == NULL || *result->trajectory_csv == '\0') {
      free(result->trajectory_csv);
      free(result->resolutions);
      return 0;
    }
  }
  if (multiply_overflows_size_t(result->orbits, result->steps_per_orbit)) {
    free(result->trajectory_csv);
    free(result->resolutions);
    return 0;
  }
  return 1;
}

static int exact_adm(const double position[3], kerr_geodesic_adm *adm) {
  kerr_adm_value_gradient fields[KERR_ADM_FIELD_COUNT];
  double values[KERR_ADM_FIELD_COUNT];
  double gradients[KERR_ADM_FIELD_COUNT * 3];
  size_t field;
  size_t axis;
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

static int valid_trajectory_position(backend_context *context,
                                     const double position[3]) {
  const double radius = radius_of(position);
  size_t axis;
  if (!isfinite(radius) || !(radius > context->throat)) {
    snprintf(context->failure, sizeof(context->failure),
             "trajectory reached r=%.17g, outside the regular exterior r>s",
             radius);
    return 0;
  }
  for (axis = 0; axis < 3; ++axis) {
    if (!isfinite(position[axis]) || fabs(position[axis]) > context->half_width) {
      snprintf(context->failure, sizeof(context->failure),
               "trajectory coordinate %zu=% .17g left [-%.17g,%.17g]", axis,
               position[axis], context->half_width, context->half_width);
      return 0;
    }
  }
  return 1;
}

static int backend_evaluate(void *opaque, const double position[3],
                            kerr_geodesic_adm *adm) {
  backend_context *context = (backend_context *)opaque;
  double values[KERR_ADM_FIELD_COUNT];
  double gradients[KERR_ADM_FIELD_COUNT * 3];
  const double started = omp_get_wtime();
  int success = 0;
  size_t field;
  size_t axis;
  context->failure[0] = '\0';
  if (!valid_trajectory_position(context, position)) goto finished;
  if (context->kind == BACKEND_EXACT) {
    success = exact_adm(position, adm);
    if (!success)
      snprintf(context->failure, sizeof(context->failure),
               "exact ADM evaluation failed at (% .17g,% .17g,% .17g)",
               position[0], position[1], position[2]);
    goto finished;
  }
  if (context->kind == BACKEND_HERMITE) {
    hermite3d_value_gradient result[KERR_ADM_FIELD_COUNT];
    const hermite3d_status status = hermite3d_interpolate_value_gradient(
        context->hermite_grid, position[0], position[1], position[2],
        KERR_ADM_FIELD_COUNT, context->functions, result);
    if (status != HERMITE3D_SUCCESS) {
      snprintf(context->failure, sizeof(context->failure),
               "Hermite interpolation failed at (% .17g,% .17g,% .17g): %s",
               position[0], position[1], position[2],
               hermite_status_name(status));
      goto finished;
    }
    for (field = 0; field < KERR_ADM_FIELD_COUNT; ++field) {
      values[field] = result[field].value;
      for (axis = 0; axis < 3; ++axis)
        gradients[field * 3 + axis] = result[field].gradient[axis];
    }
  } else {
    static const int displacement[4] = {-2, -1, 1, 2};
    double samples[4][KERR_ADM_FIELD_COUNT];
    lagrange3d_reference_status status = lagrange3d_reference_interpolate(
        context->lagrange_grid, position[0], position[1], position[2], values);
    if (status != LAGRANGE3D_REFERENCE_SUCCESS) {
      snprintf(context->failure, sizeof(context->failure),
               "Lagrange center interpolation failed at (% .17g,% .17g,% .17g): %s",
               position[0], position[1], position[2],
               lagrange3d_reference_status_name(status));
      goto finished;
    }
    for (axis = 0; axis < 3; ++axis) {
      size_t sample;
      for (sample = 0; sample < 4; ++sample) {
        double displaced[3] = {position[0], position[1], position[2]};
        displaced[axis] += (double)displacement[sample] * context->spacing;
        status = lagrange3d_reference_interpolate(
            context->lagrange_grid, displaced[0], displaced[1], displaced[2],
            samples[sample]);
        if (status != LAGRANGE3D_REFERENCE_SUCCESS) {
          snprintf(context->failure, sizeof(context->failure),
                   "Lagrange offset %dH on axis %zu failed: %s",
                   displacement[sample], axis,
                   lagrange3d_reference_status_name(status));
          goto finished;
        }
      }
      for (field = 0; field < KERR_ADM_FIELD_COUNT; ++field) {
        gradients[field * 3 + axis] =
            (samples[0][field] - 8.0 * samples[1][field] +
             8.0 * samples[2][field] - samples[3][field]) /
            (12.0 * context->spacing);
      }
    }
  }
  success = kerr_geodesic_adm_from_ordered_fields(values, gradients, adm) ==
            KERR_GEODESIC_SUCCESS;
  if (!success)
    snprintf(context->failure, sizeof(context->failure),
             "interpolated ADM data were nonfinite");
finished:
  context->seconds += omp_get_wtime() - started;
  ++context->evaluations;
  return success;
}

static int compute_grid_sizes(size_t resolution, size_t *dimension,
                              size_t *grid_points, size_t *total_values) {
  size_t square;
  if (add_overflows_size_t(resolution, 2U * (size_t)GHOST_LAYERS)) return 0;
  *dimension = resolution + 2U * (size_t)GHOST_LAYERS;
  if (multiply_overflows_size_t(*dimension, *dimension)) return 0;
  square = *dimension * *dimension;
  if (multiply_overflows_size_t(square, *dimension)) return 0;
  *grid_points = square * *dimension;
  if (*grid_points > (size_t)INT_MAX ||
      multiply_overflows_size_t(*grid_points, KERR_ADM_FIELD_COUNT))
    return 0;
  *total_values = *grid_points * KERR_ADM_FIELD_COUNT;
  return !multiply_overflows_size_t(*total_values, sizeof(double));
}

static int initialize_coordinates(double *coordinate[3], size_t dimension,
                                  double start, double spacing) {
  size_t axis;
  size_t index;
  for (axis = 0; axis < 3; ++axis) coordinate[axis] = NULL;
  if (multiply_overflows_size_t(dimension, sizeof(double))) return 0;
  for (axis = 0; axis < 3; ++axis) {
    coordinate[axis] = (double *)malloc(dimension * sizeof(double));
    if (coordinate[axis] == NULL) return 0;
    for (index = 0; index < dimension; ++index)
      coordinate[axis][index] = start + (double)index * spacing;
  }
  return 1;
}

static void free_coordinates(double *coordinate[3]) {
  size_t axis;
  for (axis = 0; axis < 3; ++axis) free(coordinate[axis]);
}

static int fill_grid(double *storage, size_t dimension, size_t grid_points,
                     double *coordinate[3], size_t resolution,
                     double noise_epsilon, uint64_t noise_seed,
                     int *team_size) {
  const int maximum_threads = omp_get_max_threads();
  grid_failure *failures;
  const grid_failure *first_failure = NULL;
  int thread;
  if (maximum_threads < 1 ||
      multiply_overflows_size_t((size_t)maximum_threads, sizeof(*failures)))
    return 0;
  failures =
      (grid_failure *)malloc((size_t)maximum_threads * sizeof(*failures));
  if (failures == NULL) return 0;
  for (thread = 0; thread < maximum_threads; ++thread) {
    failures[thread].offset = SIZE_MAX;
    failures[thread].status = KERR_ADM_EXACT_SUCCESS;
  }
  *team_size = 1;
#pragma omp parallel default(none)                                                \
    shared(storage, dimension, grid_points, coordinate, resolution,             \
           noise_epsilon, noise_seed, failures, team_size)
  {
    const int local_thread = omp_get_thread_num();
#pragma omp single
    *team_size = omp_get_num_threads();
#pragma omp for collapse(2) schedule(static)
    for (size_t index2 = 0; index2 < dimension; ++index2) {
      for (size_t index1 = 0; index1 < dimension; ++index1) {
        for (size_t index0 = 0; index0 < dimension; ++index0) {
          const size_t offset =
              index0 + dimension * (index1 + dimension * index2);
          double values[KERR_ADM_FIELD_COUNT];
          const kerr_adm_exact_status status = kerr_adm_exact_values(
              coordinate[0][index0], coordinate[1][index1],
              coordinate[2][index2], values);
          if (status != KERR_ADM_EXACT_SUCCESS) {
            if (offset < failures[local_thread].offset) {
              failures[local_thread].offset = offset;
              failures[local_thread].status = status;
            }
            continue;
          }
          for (size_t field = 0; field < KERR_ADM_FIELD_COUNT; ++field) {
            const double value = kerr_adm_apply_noise(
                values[field], noise_epsilon, noise_seed, resolution, field,
                index0, index1, index2);
            if (!isfinite(value)) {
              if (offset < failures[local_thread].offset) {
                failures[local_thread].offset = offset;
                failures[local_thread].status = KERR_ADM_EXACT_NONFINITE_RESULT;
              }
              continue;
            }
            storage[field * grid_points + offset] = value;
          }
        }
      }
    }
  }
  for (thread = 0; thread < maximum_threads; ++thread) {
    if (failures[thread].offset != SIZE_MAX &&
        (first_failure == NULL ||
         failures[thread].offset < first_failure->offset))
      first_failure = &failures[thread];
  }
  if (first_failure != NULL) {
    const size_t i0 = first_failure->offset % dimension;
    const size_t i1 = (first_failure->offset / dimension) % dimension;
    const size_t i2 = first_failure->offset / (dimension * dimension);
    fprintf(stderr, "grid evaluation failed at (%a,%a,%a): %s\n",
            coordinate[0][i0], coordinate[1][i1], coordinate[2][i2],
            exact_status_name(first_failure->status));
    free(failures);
    return 0;
  }
  free(failures);
  return 1;
}

static int state_diagnostics(
    const kerr_geodesic_circular_orbit *orbit,
    const kerr_geodesic_state *state, double phase, double time,
    backend_context *method_backend, double initial_method_hamiltonian,
    const kerr_geodesic_state *reference, scalar_norm metric[METRIC_COUNT],
    double *method_hamiltonian, double *exact_energy, double *mass_shell) {
  kerr_geodesic_adm method_adm;
  kerr_geodesic_adm analytic_adm;
  kerr_geodesic_state ignored_rhs;
  double local_method_hamiltonian;
  double local_exact_energy;
  double local_mass_shell;
  const double radius = radius_of(state->position);
  const double position_error =
      vector_difference(state->position, reference->position);
  const double momentum_error =
      vector_difference(state->momentum, reference->momentum);
  const double lz = kerr_geodesic_angular_momentum(state);
  kerr_geodesic_status status;
  if (!backend_evaluate(method_backend, state->position, &method_adm)) return 0;
  status = kerr_geodesic_rhs_from_adm(&method_adm, state, &ignored_rhs,
                                      &local_method_hamiltonian);
  if (status != KERR_GEODESIC_SUCCESS) {
    snprintf(method_backend->failure, sizeof(method_backend->failure),
             "geodesic RHS rejected %s ADM data: %s",
             method_backend->kind == BACKEND_HERMITE ? "Hermite" :
             method_backend->kind == BACKEND_LAGRANGE ? "Lagrange" : "exact",
             kerr_geodesic_status_name(status));
    return 0;
  }
  if (!exact_adm(state->position, &analytic_adm))
    return 0;
  status = kerr_geodesic_rhs_from_adm(&analytic_adm, state, &ignored_rhs,
                                      &local_exact_energy);
  if (status != KERR_GEODESIC_SUCCESS ||
      kerr_geodesic_mass_shell_residual(&analytic_adm, state, -orbit->energy,
                                        &local_mass_shell) !=
          KERR_GEODESIC_SUCCESS)
    return 0;
  norm_add(&metric[METRIC_POSITION], position_error);
  norm_add(&metric[METRIC_MOMENTUM], momentum_error);
  norm_add(&metric[METRIC_RADIAL], radius - orbit->les_radius);
  norm_add(&metric[METRIC_VERTICAL], state->position[2]);
  norm_add(&metric[METRIC_PHASE], phase - orbit->angular_frequency * time);
  norm_add(&metric[METRIC_EXACT_ENERGY], local_exact_energy - orbit->energy);
  norm_add(&metric[METRIC_METHOD_HAMILTONIAN],
           local_method_hamiltonian - initial_method_hamiltonian);
  norm_add(&metric[METRIC_ANGULAR_MOMENTUM],
           lz - orbit->angular_momentum);
  norm_add(&metric[METRIC_MASS_SHELL], local_mass_shell);
  *method_hamiltonian = local_method_hamiltonian;
  *exact_energy = local_exact_energy;
  *mass_shell = local_mass_shell;
  return 1;
}

static int write_csv_header(FILE *stream) {
  return fprintf(
             stream,
             "resolution,method,t,x,y,z,px,py,pz,closed_x,closed_y,closed_z,"
             "reference_x,reference_y,reference_z,r,phase,position_error,"
             "momentum_error,radial_error,phase_error,exact_energy,"
             "exact_energy_drift,method_hamiltonian,method_hamiltonian_drift,"
             "lz,lz_drift,mass_shell_residual\n") >= 0;
}

static int write_csv_row(
    FILE *stream, size_t resolution, const char *method, double time,
    const kerr_geodesic_state *state, const kerr_geodesic_state *closed,
    const kerr_geodesic_state *reference, double phase,
    const kerr_geodesic_circular_orbit *orbit, double exact_energy,
    double method_hamiltonian, double initial_method_hamiltonian,
    double mass_shell) {
  const double radius = radius_of(state->position);
  const double position_error =
      vector_difference(state->position, reference->position);
  const double momentum_error =
      vector_difference(state->momentum, reference->momentum);
  const double lz = kerr_geodesic_angular_momentum(state);
  return fprintf(
             stream,
             "%zu,%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
             "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
             "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
             resolution, method, time, state->position[0], state->position[1],
             state->position[2], state->momentum[0], state->momentum[1],
             state->momentum[2], closed->position[0], closed->position[1],
             closed->position[2], reference->position[0],
             reference->position[1], reference->position[2], radius, phase,
             position_error, momentum_error, radius - orbit->les_radius,
             phase - orbit->angular_frequency * time, exact_energy,
             exact_energy - orbit->energy, method_hamiltonian,
             method_hamiltonian - initial_method_hamiltonian, lz,
             lz - orbit->angular_momentum, mass_shell) >= 0;
}

static int build_exact_reference(
    const options *arguments, const kerr_geodesic_circular_orbit *orbit,
    double throat, reference_solution *reference) {
  backend_context backend;
  scalar_norm metrics[METRIC_COUNT] = {{0.0L, 0.0L, 0}};
  kerr_geodesic_phase phase;
  kerr_geodesic_state state = orbit->initial_state;
  size_t step;
  size_t metric;
  double initial_hamiltonian = orbit->energy;
  double method_hamiltonian = 0.0;
  double exact_energy = 0.0;
  double mass_shell = 0.0;
  const double dt = orbit->period / (double)arguments->steps_per_orbit;
  const size_t total_steps = arguments->orbits * arguments->steps_per_orbit;
  const double started = omp_get_wtime();
  memset(reference, 0, sizeof(*reference));
  if (add_overflows_size_t(total_steps, 1) ||
      multiply_overflows_size_t(total_steps + 1, sizeof(*reference->state)) ||
      multiply_overflows_size_t(total_steps + 1, sizeof(*reference->phase)))
    return 0;
  reference->state = (kerr_geodesic_state *)malloc(
      (total_steps + 1) * sizeof(*reference->state));
  reference->phase =
      (double *)malloc((total_steps + 1) * sizeof(*reference->phase));
  if (reference->state == NULL || reference->phase == NULL) return 0;
  memset(&backend, 0, sizeof(backend));
  backend.kind = BACKEND_EXACT;
  backend.throat = throat;
  backend.half_width = ACTIVE_HALF_WIDTH;
  kerr_geodesic_phase_initialize(&phase, &state);
  for (step = 0; step <= total_steps; ++step) {
    kerr_geodesic_state closed;
    const double time = (double)step * dt;
    reference->state[step] = state;
    reference->phase[step] = phase.angle;
    kerr_geodesic_circular_state(orbit, time, &closed);
    if (!state_diagnostics(orbit, &state, phase.angle, time, &backend,
                           initial_hamiltonian, &closed, metrics,
                           &method_hamiltonian, &exact_energy, &mass_shell)) {
      fprintf(stderr, "exact trajectory diagnostic failed at step %zu: %s\n",
              step, backend.failure);
      return 0;
    }
    if (step == total_steps) break;
    {
      const kerr_geodesic_status status =
          kerr_geodesic_rk4_step(backend_evaluate, &backend, dt, &state);
      if (status != KERR_GEODESIC_SUCCESS) {
        fprintf(stderr, "exact RK4 step %zu failed: %s%s%s\n", step,
                kerr_geodesic_status_name(status),
                backend.failure[0] == '\0' ? "" : ": ", backend.failure);
        return 0;
      }
    }
    kerr_geodesic_phase_update(&phase, &state);
  }
  reference->steps = total_steps;
  reference->summary.minimum_radius = HUGE_VAL;
  for (step = 0; step <= total_steps; ++step) {
    const double radius = radius_of(reference->state[step].position);
    if (radius < reference->summary.minimum_radius)
      reference->summary.minimum_radius = radius;
  }
  for (metric = 0; metric < METRIC_COUNT; ++metric)
    reference->summary.metric[metric] = norm_finish(&metrics[metric]);
  reference->summary.backend_seconds = backend.seconds;
  reference->summary.backend_evaluations = backend.evaluations;
  reference->summary.integration_seconds = omp_get_wtime() - started;
  return 1;
}

static void free_reference(reference_solution *reference) {
  free(reference->phase);
  free(reference->state);
  memset(reference, 0, sizeof(*reference));
}

static int initial_backend_hamiltonian(backend_context *backend,
                                       const kerr_geodesic_state *state,
                                       double *hamiltonian) {
  kerr_geodesic_adm adm;
  kerr_geodesic_state rhs;
  kerr_geodesic_status status;
  if (!backend_evaluate(backend, state->position, &adm)) return 0;
  status = kerr_geodesic_rhs_from_adm(&adm, state, &rhs, hamiltonian);
  if (status != KERR_GEODESIC_SUCCESS) {
    snprintf(backend->failure, sizeof(backend->failure),
             "initial ADM data were rejected: %s",
             kerr_geodesic_status_name(status));
    return 0;
  }
  return 1;
}

static int integrate_method(
    size_t resolution, size_t method, const options *arguments,
    const kerr_geodesic_circular_orbit *orbit,
    const reference_solution *reference, backend_context *backend, FILE *csv,
    method_summary *summary) {
  scalar_norm metrics[METRIC_COUNT] = {{0.0L, 0.0L, 0}};
  kerr_geodesic_state state = orbit->initial_state;
  kerr_geodesic_phase phase;
  double initial_method_hamiltonian;
  double method_hamiltonian = 0.0;
  double exact_energy = 0.0;
  double mass_shell = 0.0;
  const double dt = orbit->period / (double)arguments->steps_per_orbit;
  const double started = omp_get_wtime();
  size_t step;
  size_t metric;
  memset(summary, 0, sizeof(*summary));
  summary->minimum_radius = HUGE_VAL;
  backend->seconds = 0.0;
  backend->evaluations = 0;
  if (!initial_backend_hamiltonian(backend, &state,
                                   &initial_method_hamiltonian)) {
    fprintf(stderr, "%s initial Hamiltonian failed: %s\n", method_name(method),
            backend->failure);
    return 0;
  }
  kerr_geodesic_phase_initialize(&phase, &state);
  for (step = 0; step <= reference->steps; ++step) {
    kerr_geodesic_state closed;
    const double time = (double)step * dt;
    const double radius = radius_of(state.position);
    if (radius < summary->minimum_radius) summary->minimum_radius = radius;
    kerr_geodesic_circular_state(orbit, time, &closed);
    if (!state_diagnostics(orbit, &state, phase.angle, time, backend,
                           initial_method_hamiltonian,
                           &reference->state[step], metrics,
                           &method_hamiltonian, &exact_energy, &mass_shell)) {
      fprintf(stderr, "%s diagnostic failed at step %zu: %s\n",
              method_name(method), step, backend->failure);
      return 0;
    }
    if (csv != NULL &&
        (step % arguments->output_every == 0 || step == reference->steps) &&
        !write_csv_row(csv, resolution, method_name(method), time, &state,
                       &closed, &reference->state[step], phase.angle, orbit,
                       exact_energy, method_hamiltonian,
                       initial_method_hamiltonian, mass_shell)) {
      fprintf(stderr, "could not write trajectory CSV\n");
      return 0;
    }
    if (step == reference->steps) break;
    {
      const kerr_geodesic_status status =
          kerr_geodesic_rk4_step(backend_evaluate, backend, dt, &state);
      if (status != KERR_GEODESIC_SUCCESS) {
        fprintf(stderr, "%s RK4 step %zu failed: %s%s%s\n",
                method_name(method), step, kerr_geodesic_status_name(status),
                backend->failure[0] == '\0' ? "" : ": ", backend->failure);
        return 0;
      }
    }
    kerr_geodesic_phase_update(&phase, &state);
  }
  for (metric = 0; metric < METRIC_COUNT; ++metric)
    summary->metric[metric] = norm_finish(&metrics[metric]);
  summary->final_position_error =
      vector_difference(state.position,
                        reference->state[reference->steps].position);
  summary->final_momentum_error =
      vector_difference(state.momentum,
                        reference->state[reference->steps].momentum);
  summary->backend_seconds = backend->seconds;
  summary->backend_evaluations = backend->evaluations;
  summary->integration_seconds = omp_get_wtime() - started;
  return 1;
}

static int write_exact_csv(size_t resolution, const options *arguments,
                           const kerr_geodesic_circular_orbit *orbit,
                           const reference_solution *reference, FILE *csv) {
  const double dt = orbit->period / (double)arguments->steps_per_orbit;
  size_t step;
  for (step = 0; step <= reference->steps; ++step) {
    kerr_geodesic_state closed;
    kerr_geodesic_adm adm;
    kerr_geodesic_state rhs;
    double energy;
    double residual;
    const double time = (double)step * dt;
    if (step % arguments->output_every != 0 && step != reference->steps)
      continue;
    kerr_geodesic_circular_state(orbit, time, &closed);
    if (!exact_adm(reference->state[step].position, &adm) ||
        kerr_geodesic_rhs_from_adm(&adm, &reference->state[step], &rhs,
                                   &energy) != KERR_GEODESIC_SUCCESS ||
        kerr_geodesic_mass_shell_residual(
            &adm, &reference->state[step], -orbit->energy, &residual) !=
            KERR_GEODESIC_SUCCESS ||
        !write_csv_row(csv, resolution, method_name(METHOD_EXACT), time,
                       &reference->state[step], &closed,
                       &reference->state[step], reference->phase[step], orbit,
                       energy, energy, orbit->energy, residual))
      return 0;
  }
  return 1;
}

static int compute_rhs_defect(const kerr_geodesic_circular_orbit *orbit,
                              backend_context *backend, rhs_defect *result) {
  scalar_norm velocity = {0.0L, 0.0L, 0};
  scalar_norm force = {0.0L, 0.0L, 0};
  size_t sample;
  for (sample = 0; sample < RHS_DIAGNOSTIC_SAMPLES; ++sample) {
    const double time = orbit->period * (double)sample /
                        (double)RHS_DIAGNOSTIC_SAMPLES;
    kerr_geodesic_state state;
    kerr_geodesic_state exact_rhs;
    kerr_geodesic_state approximate_rhs;
    size_t axis;
    kerr_geodesic_status status;
    kerr_geodesic_circular_state(orbit, time, &state);
    status = kerr_geodesic_rhs(backend_evaluate, backend, &state,
                               &approximate_rhs, NULL);
    if (status != KERR_GEODESIC_SUCCESS) {
      if (backend->failure[0] == '\0')
        snprintf(backend->failure, sizeof(backend->failure),
                 "RHS defect sample %zu failed: %s", sample,
                 kerr_geodesic_status_name(status));
      return 0;
    }
    {
      kerr_geodesic_adm adm;
      if (!exact_adm(state.position, &adm) ||
          kerr_geodesic_rhs_from_adm(&adm, &state, &exact_rhs, NULL) !=
              KERR_GEODESIC_SUCCESS)
        return 0;
    }
    for (axis = 0; axis < 3; ++axis) {
      norm_add(&velocity,
               approximate_rhs.position[axis] - exact_rhs.position[axis]);
      norm_add(&force,
               approximate_rhs.momentum[axis] - exact_rhs.momentum[axis]);
    }
  }
  result->velocity = norm_finish(&velocity);
  result->force = norm_finish(&force);
  return 1;
}

static int evaluate_level(
    const options *arguments, size_t resolution,
    const kerr_geodesic_circular_orbit *orbit,
    const reference_solution *reference, double throat, FILE *csv,
    level_summary *summary) {
  size_t dimension;
  size_t grid_points;
  size_t total_values;
  double *storage = NULL;
  double *coordinate[3] = {NULL, NULL, NULL};
  const double *functions[KERR_ADM_FIELD_COUNT];
  hermite3d_grid hermite_grid;
  lagrange3d_reference_grid lagrange_grid;
  backend_context hermite_backend;
  backend_context lagrange_backend;
  const double spacing = 2.0 * ACTIVE_HALF_WIDTH / (double)resolution;
  const double start = ACTIVE_HALF_WIDTH * -1.0 - 5.5 * spacing;
  double started;
  size_t field;
  size_t axis;
  memset(summary, 0, sizeof(*summary));
  if (!isfinite(spacing) || !compute_grid_sizes(resolution, &dimension,
                                                &grid_points, &total_values)) {
    fprintf(stderr, "resolution %zu overflows supported grid geometry\n",
            resolution);
    return 0;
  }
  storage = (double *)malloc(total_values * sizeof(*storage));
  if (storage == NULL) {
    fprintf(stderr, "could not allocate %.3Lf MiB for N=%zu\n",
            (long double)total_values * sizeof(*storage) /
                (1024.0L * 1024.0L),
            resolution);
    return 0;
  }
  if (!initialize_coordinates(coordinate, dimension, start, spacing)) {
    fprintf(stderr, "could not allocate coordinate arrays for N=%zu\n",
            resolution);
    free_coordinates(coordinate);
    free(storage);
    return 0;
  }
  started = omp_get_wtime();
  if (!fill_grid(storage, dimension, grid_points, coordinate, resolution,
                 arguments->noise_epsilon, arguments->noise_seed,
                 &summary->grid_threads)) {
    free_coordinates(coordinate);
    free(storage);
    return 0;
  }
  summary->grid_seconds = omp_get_wtime() - started;
  for (field = 0; field < KERR_ADM_FIELD_COUNT; ++field)
    functions[field] = storage + field * grid_points;
  hermite_grid.nxx0 = dimension;
  hermite_grid.nxx1 = dimension;
  hermite_grid.nxx2 = dimension;
  hermite_grid.xx0_start = start;
  hermite_grid.xx1_start = start;
  hermite_grid.xx2_start = start;
  hermite_grid.dxx0 = spacing;
  hermite_grid.dxx1 = spacing;
  hermite_grid.dxx2 = spacing;
  hermite_grid.stride0 = 1;
  hermite_grid.stride1 = dimension;
  hermite_grid.stride2 = dimension * dimension;
  for (axis = 0; axis < 3; ++axis) {
    lagrange_grid.dimension[axis] = (int)dimension;
    lagrange_grid.spacing[axis] = spacing;
    lagrange_grid.coordinate[axis] = coordinate[axis];
  }
  for (field = 0; field < KERR_ADM_FIELD_COUNT; ++field)
    lagrange_grid.function[field] = functions[field];
  memset(&hermite_backend, 0, sizeof(hermite_backend));
  hermite_backend.kind = BACKEND_HERMITE;
  hermite_backend.throat = throat;
  hermite_backend.half_width = ACTIVE_HALF_WIDTH;
  hermite_backend.spacing = spacing;
  hermite_backend.hermite_grid = &hermite_grid;
  hermite_backend.functions = functions;
  memset(&lagrange_backend, 0, sizeof(lagrange_backend));
  lagrange_backend.kind = BACKEND_LAGRANGE;
  lagrange_backend.throat = throat;
  lagrange_backend.half_width = ACTIVE_HALF_WIDTH;
  lagrange_backend.spacing = spacing;
  lagrange_backend.lagrange_grid = &lagrange_grid;
  if (!compute_rhs_defect(orbit, &hermite_backend, &summary->defect[0]) ||
      !compute_rhs_defect(orbit, &lagrange_backend, &summary->defect[1])) {
    fprintf(stderr, "RHS defect sampling failed: %s%s\n",
            hermite_backend.failure,
            lagrange_backend.failure);
    free_coordinates(coordinate);
    free(storage);
    return 0;
  }
  if (csv != NULL &&
      !write_exact_csv(resolution, arguments, orbit, reference, csv)) {
    fprintf(stderr, "could not write exact trajectory CSV rows\n");
    free_coordinates(coordinate);
    free(storage);
    return 0;
  }
  if (!integrate_method(resolution, METHOD_HERMITE, arguments, orbit, reference,
                        &hermite_backend, csv, &summary->method[0]) ||
      !integrate_method(resolution, METHOD_LAGRANGE, arguments, orbit, reference,
                        &lagrange_backend, csv, &summary->method[1])) {
    free_coordinates(coordinate);
    free(storage);
    return 0;
  }
  summary->resolution = resolution;
  summary->dimension = dimension;
  summary->spacing = spacing;
  summary->memory_mib =
      (long double)total_values * sizeof(*storage) / (1024.0L * 1024.0L);
  free_coordinates(coordinate);
  free(storage);
  return 1;
}

static int measured_order(long double coarse, long double fine,
                          double coarse_spacing, double fine_spacing,
                          long double *order) {
  if (!(coarse > 0.0L) || !(fine > 0.0L) || !isfinite(coarse) ||
      !isfinite(fine))
    return 0;
  *order = logl(coarse / fine) /
           logl((long double)coarse_spacing / (long double)fine_spacing);
  return isfinite(*order);
}

static void print_order(const level_summary *levels, size_t level,
                        size_t method, size_t metric, int maximum) {
  long double order;
  long double coarse;
  long double fine;
  if (level == 0) {
    printf(" %8s", "--");
    return;
  }
  coarse = maximum ? levels[level - 1].method[method].metric[metric].maximum
                   : levels[level - 1].method[method].metric[metric].rms;
  fine = maximum ? levels[level].method[method].metric[metric].maximum
                 : levels[level].method[method].metric[metric].rms;
  if (!measured_order(coarse, fine, levels[level - 1].spacing,
                      levels[level].spacing, &order))
    printf(" %8s", "--");
  else
    printf(" %8.3Lf", order);
}

static void print_results(const options *arguments,
                          const kerr_geodesic_circular_orbit *orbit,
                          const reference_solution *reference,
                          const level_summary *levels) {
  static const char *const convergence_name[5] = {
      "position", "momentum", "radial", "vertical", "phase"};
  size_t level;
  size_t method;
  printf("Kerr geodesic interpolation comparison\n");
  printf("  M=1, a=0.5, R_BL=%.17g, r0=%.17g, s=%.17g\n",
         orbit->boyer_lindquist_radius, orbit->les_radius,
         kerr_adm_exact_throat_radius());
  printf("  Omega=%.17g, period=%.17g, orbits=%zu, steps/orbit=%zu\n",
         orbit->angular_frequency, orbit->period, arguments->orbits,
         arguments->steps_per_orbit);
  printf("  active cube=[-5,5]^3; six analytic ghost layers\n");
  printf("  noise epsilon=%.17g, seed=0x%016" PRIx64
         ", relative RMS=%.6g\n",
         arguments->noise_epsilon, arguments->noise_seed,
         arguments->noise_epsilon / sqrt(3.0));
  printf("  exact RK4 audit: position L2=%.6Le Linf=%.6Le; "
         "momentum L2=%.6Le Linf=%.6Le\n",
         reference->summary.metric[METRIC_POSITION].rms,
         reference->summary.metric[METRIC_POSITION].maximum,
         reference->summary.metric[METRIC_MOMENTUM].rms,
         reference->summary.metric[METRIC_MOMENTUM].maximum);
  printf("  exact RHS evaluations=%zu, %.3f us/eval, integration=%.3f s\n",
         reference->summary.backend_evaluations,
         reference->summary.backend_evaluations == 0
             ? 0.0
             : 1.0e6 * reference->summary.backend_seconds /
                   (double)reference->summary.backend_evaluations,
         reference->summary.integration_seconds);
  printf("\nGrid and timing diagnostics\n");
  printf("      N    stored            H    memory[MiB] threads grid[s]   "
         "method       min(r)      evals  us/eval integrate[s]\n");
  for (level = 0; level < arguments->resolution_count; ++level) {
    for (method = 0; method < 2; ++method) {
      const method_summary *data = &levels[level].method[method];
      printf("%7zu %9zu %12.5e %14.3Lf %7d %7.3f   %-11s %10.6f %10zu "
             "%8.3f %12.3f\n",
             levels[level].resolution, levels[level].dimension,
             levels[level].spacing, levels[level].memory_mib,
             levels[level].grid_threads, levels[level].grid_seconds,
             method_name(method + 1), data->minimum_radius,
             data->backend_evaluations,
             data->backend_evaluations == 0
                 ? 0.0
                 : 1.0e6 * data->backend_seconds /
                       (double)data->backend_evaluations,
             data->integration_seconds);
    }
  }
  printf("\nLocal RHS defects on 256 points of the exact circular orbit\n");
  printf("      N method          velocity L2   velocity Linf      force L2   "
         "force Linf\n");
  for (level = 0; level < arguments->resolution_count; ++level) {
    for (method = 0; method < 2; ++method) {
      printf("%7zu %-11s %13.6Le %15.6Le %13.6Le %13.6Le\n",
             levels[level].resolution, method_name(method + 1),
             levels[level].defect[method].velocity.rms,
             levels[level].defect[method].velocity.maximum,
             levels[level].defect[method].force.rms,
             levels[level].defect[method].force.maximum);
    }
  }
  printf("\nTrajectory errors relative to the same-step exact-ADM RK4 solution\n");
  printf("      N method      quantity               L2         Linf    p(L2) "
         " p(Linf)\n");
  for (level = 0; level < arguments->resolution_count; ++level) {
    for (method = 0; method < 2; ++method) {
      size_t metric;
      for (metric = 0; metric < 5; ++metric) {
        const norm_result norm = levels[level].method[method].metric[metric];
        printf("%7zu %-11s %-12s %12.5Le %12.5Le", levels[level].resolution,
               method_name(method + 1), convergence_name[metric], norm.rms,
               norm.maximum);
        print_order(levels, level, method, metric, 0);
        print_order(levels, level, method, metric, 1);
        putchar('\n');
      }
    }
  }
  printf("\nInvariant and Hamiltonian drift (RMS / maximum)\n");
  printf("      N method          exact E             method H            Lz  "
         "            mass shell\n");
  printf("%7s %-11s %.3Le/%.3Le  %.3Le/%.3Le  %.3Le/%.3Le  "
         "%.3Le/%.3Le\n",
         "--", method_name(METHOD_EXACT),
         reference->summary.metric[METRIC_EXACT_ENERGY].rms,
         reference->summary.metric[METRIC_EXACT_ENERGY].maximum,
         reference->summary.metric[METRIC_METHOD_HAMILTONIAN].rms,
         reference->summary.metric[METRIC_METHOD_HAMILTONIAN].maximum,
         reference->summary.metric[METRIC_ANGULAR_MOMENTUM].rms,
         reference->summary.metric[METRIC_ANGULAR_MOMENTUM].maximum,
         reference->summary.metric[METRIC_MASS_SHELL].rms,
         reference->summary.metric[METRIC_MASS_SHELL].maximum);
  for (level = 0; level < arguments->resolution_count; ++level) {
    for (method = 0; method < 2; ++method) {
      const method_summary *data = &levels[level].method[method];
      printf("%7zu %-11s %.3Le/%.3Le  %.3Le/%.3Le  %.3Le/%.3Le  "
             "%.3Le/%.3Le\n",
             levels[level].resolution, method_name(method + 1),
             data->metric[METRIC_EXACT_ENERGY].rms,
             data->metric[METRIC_EXACT_ENERGY].maximum,
             data->metric[METRIC_METHOD_HAMILTONIAN].rms,
             data->metric[METRIC_METHOD_HAMILTONIAN].maximum,
             data->metric[METRIC_ANGULAR_MOMENTUM].rms,
             data->metric[METRIC_ANGULAR_MOMENTUM].maximum,
             data->metric[METRIC_MASS_SHELL].rms,
             data->metric[METRIC_MASS_SHELL].maximum);
    }
  }
  printf("\nFinal state error relative to the exact-ADM RK4 trajectory\n");
  printf("      N method          |delta x|      |delta p|\n");
  for (level = 0; level < arguments->resolution_count; ++level) {
    for (method = 0; method < 2; ++method) {
      printf("%7zu %-11s %14.6e %14.6e\n", levels[level].resolution,
             method_name(method + 1),
             levels[level].method[method].final_position_error,
             levels[level].method[method].final_momentum_error);
    }
  }
  {
    const level_summary *fine = &levels[arguments->resolution_count - 1];
    const long double minimum_position =
        fminl(fine->method[0].metric[METRIC_POSITION].rms,
              fine->method[1].metric[METRIC_POSITION].rms);
    const long double minimum_momentum =
        fminl(fine->method[0].metric[METRIC_MOMENTUM].rms,
              fine->method[1].metric[METRIC_MOMENTUM].rms);
    if ((minimum_position > 0.0L &&
         reference->summary.metric[METRIC_POSITION].rms >
             0.1L * minimum_position) ||
        (minimum_momentum > 0.0L &&
         reference->summary.metric[METRIC_MOMENTUM].rms >
             0.1L * minimum_momentum)) {
      printf("\nWARNING: exact-field RK4 error is not at least ten times smaller "
             "than\n");
      printf("the best finest-grid interpolation trajectory error. Increase "
             "--steps-per-orbit.\n");
    }
  }
  printf("\nHermite differentiates one C2 interpolating polynomial.  The "
         "Lagrange force\n");
  printf("uses fourth-order finite differences of 13 value interpolations with "
         "step H.\n");
}

int main(int argc, char **argv) {
  options arguments;
  kerr_geodesic_circular_orbit orbit;
  reference_solution reference;
  level_summary *levels = NULL;
  backend_context exact_context;
  FILE *csv = NULL;
  size_t level;
  const int parsed = parse_options(argc, argv, &arguments);
  const double throat = kerr_adm_exact_throat_radius();
  int exit_status = EXIT_FAILURE;
  memset(&reference, 0, sizeof(reference));
  if (parsed == 2) {
    print_usage(stdout, argv[0]);
    return EXIT_SUCCESS;
  }
  if (parsed == 0) {
    print_usage(stderr, argv[0]);
    return EXIT_FAILURE;
  }
  memset(&exact_context, 0, sizeof(exact_context));
  exact_context.kind = BACKEND_EXACT;
  exact_context.throat = throat;
  exact_context.half_width = ACTIVE_HALF_WIDTH;
  if (kerr_geodesic_initialize_circular_orbit(
          ORBIT_BL_RADIUS, backend_evaluate, &exact_context, &orbit) !=
      KERR_GEODESIC_SUCCESS) {
    fprintf(stderr, "could not construct circular orbit: %s\n",
            exact_context.failure);
    goto cleanup;
  }
  if (multiply_overflows_size_t(arguments.resolution_count, sizeof(*levels))) {
    fprintf(stderr, "too many resolution levels\n");
    goto cleanup;
  }
  levels =
      (level_summary *)calloc(arguments.resolution_count, sizeof(*levels));
  if (levels == NULL) {
    fprintf(stderr, "could not allocate level summaries\n");
    goto cleanup;
  }
  if (!build_exact_reference(&arguments, &orbit, throat, &reference)) {
    fprintf(stderr, "could not integrate exact reference trajectory\n");
    goto cleanup_reference;
  }
  if (arguments.trajectory_csv != NULL) {
    csv = fopen(arguments.trajectory_csv, "w");
    if (csv == NULL || !write_csv_header(csv)) {
      fprintf(stderr, "could not create trajectory CSV '%s'\n",
              arguments.trajectory_csv);
      goto cleanup_reference;
    }
  }
  for (level = 0; level < arguments.resolution_count; ++level) {
    fprintf(stderr, "evaluating geodesic grid N=%zu ...\n",
            arguments.resolutions[level]);
    if (!evaluate_level(&arguments, arguments.resolutions[level], &orbit,
                        &reference, throat, csv, &levels[level]))
      goto cleanup_reference;
  }
  if (csv != NULL && fclose(csv) != 0) {
    csv = NULL;
    fprintf(stderr, "could not finalize trajectory CSV\n");
    goto cleanup_reference;
  }
  csv = NULL;
  print_results(&arguments, &orbit, &reference, levels);
  if (arguments.trajectory_csv != NULL)
    printf("\nwrote trajectory samples to %s\n", arguments.trajectory_csv);
  exit_status = EXIT_SUCCESS;
cleanup_reference:
  if (csv != NULL) fclose(csv);
  free_reference(&reference);
cleanup:
  free(levels);
  free(arguments.trajectory_csv);
  free(arguments.resolutions);
  return exit_status;
}
