#include "hermite3d.h"
#include "kerr_adm_exact.h"
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
  QUANTITY_COUNT = 4,
  METHOD_COUNT = 3,
  GHOST_LAYERS = 6
};

typedef enum comparison_method {
  METHOD_HERMITE = 0,
  METHOD_LAGRANGE_FD,
  METHOD_ANALYTIC_FD
} comparison_method;

static const uint64_t DEFAULT_SEED = UINT64_C(0x4b4552524845524d);

typedef struct options {
  size_t *resolutions;
  size_t resolution_count;
  size_t point_count;
  uint64_t seed;
} options;

typedef struct query_point {
  double coordinate[3];
  double radius;
  kerr_adm_value_gradient exact[KERR_ADM_FIELD_COUNT];
} query_point;

typedef struct scaled_sum_squares {
  long double scale;
  long double scaled_sum;
} scaled_sum_squares;

typedef struct norm_accumulator {
  scaled_sum_squares error_squares;
  scaled_sum_squares reference_squares;
  long double maximum_error;
  long double maximum_reference;
} norm_accumulator;

typedef struct error_norm {
  long double rms;
  long double maximum;
  long double relative_rms;
  long double relative_maximum;
  int relative_rms_defined;
  int relative_maximum_defined;
} error_norm;

typedef struct level_summary {
  size_t resolution;
  size_t stored_dimension;
  size_t throat_crossing_points;
  double spacing;
  double minimum_radius;
  long double memory_mib;
  int grid_threads;
  double hermite_seconds_per_point;
  double lagrange_value_seconds_per_point;
  double lagrange_full_seconds_per_point;
  error_norm field[METHOD_COUNT][QUANTITY_COUNT][KERR_ADM_FIELD_COUNT];
  error_norm aggregate[METHOD_COUNT][QUANTITY_COUNT];
} level_summary;

typedef struct timed_results {
  hermite3d_value_gradient *hermite;
  double *lagrange_value;
  kerr_adm_value_gradient *lagrange;
  double hermite_seconds;
  double lagrange_value_seconds;
  double lagrange_full_seconds;
} timed_results;

typedef struct grid_failure {
  size_t index;
  kerr_adm_exact_status status;
} grid_failure;

static int add_overflows_size_t(size_t left, size_t right) {
  return left > SIZE_MAX - right;
}

static int multiply_overflows_size_t(size_t left, size_t right) {
  return left != 0 && right > SIZE_MAX / left;
}

static uint64_t splitmix64_next(uint64_t *state) {
  uint64_t value;

  *state += UINT64_C(0x9e3779b97f4a7c15);
  value = *state;
  value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
  return value ^ (value >> 31);
}

/* Midpoints of 2^52 bins exclude 0, 1, and exactly 1/2. */
static double uniform_open(uint64_t *state) {
  const uint64_t bin = splitmix64_next(state) >> 12;
  return (double)(2 * bin + 1) * 0x1p-53;
}

static void scaled_sum_add(scaled_sum_squares *sum, long double value) {
  const long double magnitude = fabsl(value);

  if (magnitude == 0.0L) return;
  if (sum->scale < magnitude) {
    const long double ratio = sum->scale / magnitude;
    sum->scaled_sum = 1.0L + sum->scaled_sum * ratio * ratio;
    sum->scale = magnitude;
  } else {
    const long double ratio = magnitude / sum->scale;
    sum->scaled_sum += ratio * ratio;
  }
}

static long double scaled_sum_rms(const scaled_sum_squares *sum,
                                  size_t count) {
  if (sum->scale == 0.0L) return 0.0L;
  return sum->scale * sqrtl(sum->scaled_sum / (long double)count);
}

static void add_norm_sample(norm_accumulator *accumulator, long double error,
                            long double reference) {
  const long double absolute_error = fabsl(error);
  const long double absolute_reference = fabsl(reference);

  scaled_sum_add(&accumulator->error_squares, error);
  scaled_sum_add(&accumulator->reference_squares, reference);
  if (absolute_error > accumulator->maximum_error)
    accumulator->maximum_error = absolute_error;
  if (absolute_reference > accumulator->maximum_reference)
    accumulator->maximum_reference = absolute_reference;
}

static error_norm finish_norm(const norm_accumulator *accumulator,
                              size_t count) {
  error_norm norm;
  const long double reference_rms =
      scaled_sum_rms(&accumulator->reference_squares, count);

  norm.rms = scaled_sum_rms(&accumulator->error_squares, count);
  norm.maximum = accumulator->maximum_error;
  norm.relative_rms_defined = reference_rms > 0.0L;
  norm.relative_maximum_defined = accumulator->maximum_reference > 0.0L;
  norm.relative_rms = norm.relative_rms_defined ? norm.rms / reference_rms : 0.0L;
  norm.relative_maximum = norm.relative_maximum_defined
                              ? norm.maximum / accumulator->maximum_reference
                              : 0.0L;
  return norm;
}

static const char *exact_status_name(kerr_adm_exact_status status) {
  switch (status) {
    case KERR_ADM_EXACT_SUCCESS:
      return "success";
    case KERR_ADM_EXACT_INVALID_ARGUMENT:
      return "invalid argument";
    case KERR_ADM_EXACT_DOMAIN_ERROR:
      return "coordinate-domain error";
    case KERR_ADM_EXACT_NONDIFFERENTIABLE:
      return "nondifferentiable throat";
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

static const char *method_name(comparison_method method) {
  static const char *const names[METHOD_COUNT] = {
      "HERMITE", "LAGRANGE+FD", "ANALYTIC-FD"};
  return (int)method >= 0 && (int)method < METHOD_COUNT
             ? names[(size_t)method]
             : NULL;
}

static const char *quantity_name(size_t quantity) {
  static const char *const names[QUANTITY_COUNT] = {
      "VALUE", "D/DX", "D/DY", "D/DZ"};
  return quantity < QUANTITY_COUNT ? names[quantity] : NULL;
}

static void print_usage(FILE *stream, const char *program) {
  fprintf(stream,
          "Usage: %s --resolutions N1,N2,... --points COUNT [--seed UINT64]\n"
          "\n"
          "N values must be positive, even, distinct, and strictly increasing.\n"
          "UINT64 accepts decimal or a 0x-prefixed hexadecimal value.\n",
          program);
}

static int parse_size(const char *text, size_t *value) {
  char *end = NULL;
  uintmax_t parsed;

  if (text == NULL || *text == '\0' || *text == '-') return 0;
  errno = 0;
  parsed = strtoumax(text, &end, 10);
  if (errno != 0 || end == text || *end != '\0' || parsed == 0 ||
      parsed > SIZE_MAX) {
    return 0;
  }
  *value = (size_t)parsed;
  return 1;
}

static int parse_seed(const char *text, uint64_t *value) {
  char *end = NULL;
  uintmax_t parsed;

  if (text == NULL || *text == '\0' || *text == '-') return 0;
  errno = 0;
  parsed = strtoumax(text, &end, 0);
  if (errno != 0 || end == text || *end != '\0' || parsed > UINT64_MAX)
    return 0;
  *value = (uint64_t)parsed;
  return 1;
}

static int parse_resolutions(const char *text, size_t **values,
                             size_t *count) {
  size_t entries = 1;
  size_t length;
  size_t *result;
  char *copy;
  char *cursor;

  if (text == NULL || *text == '\0') return 0;
  for (const char *character = text; *character != '\0'; ++character)
    if (*character == ',') ++entries;
  if (multiply_overflows_size_t(entries, sizeof(*result))) return 0;
  length = strlen(text);
  if (add_overflows_size_t(length, 1)) return 0;
  result = (size_t *)malloc(entries * sizeof(*result));
  copy = (char *)malloc(length + 1);
  if (result == NULL || copy == NULL) {
    free(result);
    free(copy);
    return 0;
  }
  memcpy(copy, text, length + 1);
  cursor = copy;
  for (size_t entry = 0; entry < entries; ++entry) {
    char *comma = strchr(cursor, ',');
    if (comma != NULL) *comma = '\0';
    if (!parse_size(cursor, &result[entry]) ||
        (result[entry] & (size_t)1) != 0 ||
        (entry > 0 && result[entry] <= result[entry - 1]) ||
        (entry + 1 < entries && comma == NULL)) {
      free(result);
      free(copy);
      return 0;
    }
    cursor = comma == NULL ? cursor + strlen(cursor) : comma + 1;
  }
  free(copy);
  *values = result;
  *count = entries;
  return 1;
}

static const char *attached_value(const char *argument, const char *name) {
  const size_t length = strlen(name);
  return strncmp(argument, name, length) == 0 && argument[length] == '='
             ? argument + length + 1
             : NULL;
}

static int parse_options(int argc, char **argv, options *result) {
  const char *resolution_text = NULL;
  const char *point_text = NULL;
  const char *seed_text = NULL;

  result->resolutions = NULL;
  result->resolution_count = 0;
  result->point_count = 0;
  result->seed = DEFAULT_SEED;
  for (int argument = 1; argument < argc; ++argument) {
    const char *value;
    if (strcmp(argv[argument], "--help") == 0 ||
        strcmp(argv[argument], "-h") == 0) {
      print_usage(stdout, argv[0]);
      return 1;
    }
    value = attached_value(argv[argument], "--resolutions");
    if (strcmp(argv[argument], "--resolutions") == 0) {
      if (++argument >= argc) return 0;
      value = argv[argument];
    }
    if (value != NULL) {
      if (resolution_text != NULL) return 0;
      resolution_text = value;
      continue;
    }
    value = attached_value(argv[argument], "--points");
    if (strcmp(argv[argument], "--points") == 0) {
      if (++argument >= argc) return 0;
      value = argv[argument];
    }
    if (value != NULL) {
      if (point_text != NULL) return 0;
      point_text = value;
      continue;
    }
    value = attached_value(argv[argument], "--seed");
    if (strcmp(argv[argument], "--seed") == 0) {
      if (++argument >= argc) return 0;
      value = argv[argument];
    }
    if (value != NULL) {
      if (seed_text != NULL) return 0;
      seed_text = value;
      continue;
    }
    return 0;
  }
  if (resolution_text == NULL || point_text == NULL ||
      !parse_resolutions(resolution_text, &result->resolutions,
                         &result->resolution_count) ||
      !parse_size(point_text, &result->point_count) ||
      (seed_text != NULL && !parse_seed(seed_text, &result->seed))) {
    free(result->resolutions);
    result->resolutions = NULL;
    return 0;
  }
  return 2;
}

static query_point *generate_queries(size_t count, uint64_t seed,
                                     double half_width, double throat,
                                     double *minimum_radius) {
  query_point *points;
  uint64_t state = seed;

  if (multiply_overflows_size_t(count, sizeof(*points))) return NULL;
  points = (query_point *)malloc(count * sizeof(*points));
  if (points == NULL) return NULL;
  *minimum_radius = INFINITY;
  for (size_t point = 0; point < count; ++point) {
    kerr_adm_exact_status status;
    do {
      for (size_t axis = 0; axis < KERR_ADM_DIMENSION; ++axis) {
        points[point].coordinate[axis] =
            half_width * (2.0 * uniform_open(&state) - 1.0);
      }
      points[point].radius =
          sqrt(points[point].coordinate[0] * points[point].coordinate[0] +
               points[point].coordinate[1] * points[point].coordinate[1] +
               points[point].coordinate[2] * points[point].coordinate[2]);
    } while (points[point].radius < throat);
    status = kerr_adm_exact_value_gradient(
        points[point].coordinate[0], points[point].coordinate[1],
        points[point].coordinate[2], points[point].exact);
    if (status != KERR_ADM_EXACT_SUCCESS) {
      fprintf(stderr, "exact reference failed at query %zu: %s\n", point,
              exact_status_name(status));
      free(points);
      return NULL;
    }
    if (points[point].radius < *minimum_radius)
      *minimum_radius = points[point].radius;
  }
  return points;
}

static int compute_grid_sizes(size_t resolution, size_t *dimension,
                              size_t *grid_points, size_t *total_values) {
  size_t square;

  if (add_overflows_size_t(resolution, 2 * (size_t)GHOST_LAYERS)) return 0;
  *dimension = resolution + 2 * (size_t)GHOST_LAYERS;
  if (*dimension > (size_t)INT_MAX) return 0;
  if (multiply_overflows_size_t(*dimension, *dimension)) return 0;
  square = *dimension * *dimension;
  if (multiply_overflows_size_t(square, *dimension)) return 0;
  *grid_points = square * *dimension;
  /* The copied routine performs flattened indexing with int arithmetic. */
  if (*grid_points > (size_t)INT_MAX) return 0;
  if (multiply_overflows_size_t(*grid_points,
                                (size_t)KERR_ADM_FIELD_COUNT))
    return 0;
  *total_values = *grid_points * (size_t)KERR_ADM_FIELD_COUNT;
  return !multiply_overflows_size_t(*total_values, sizeof(double));
}

static int initialize_coordinates(double *coordinate[3], size_t dimension,
                                  double start, double spacing) {
  for (size_t axis = 0; axis < KERR_ADM_DIMENSION; ++axis) {
    if (multiply_overflows_size_t(dimension, sizeof(*coordinate[axis])))
      return 0;
    coordinate[axis] = (double *)malloc(dimension * sizeof(*coordinate[axis]));
    if (coordinate[axis] == NULL) {
      for (size_t previous = 0; previous < axis; ++previous)
        free(coordinate[previous]);
      return 0;
    }
    for (size_t index = 0; index < dimension; ++index)
      coordinate[axis][index] = start + (double)index * spacing;
  }
  return 1;
}

static void free_coordinates(double *coordinate[3]) {
  for (size_t axis = 0; axis < KERR_ADM_DIMENSION; ++axis)
    free(coordinate[axis]);
}

static int fill_grid(double *storage, size_t dimension, size_t grid_points,
                     double *coordinate[3], int *team_size) {
  const int max_threads = omp_get_max_threads();
  grid_failure *failures;
  const grid_failure *first_failure = NULL;

  if (max_threads < 1 ||
      multiply_overflows_size_t((size_t)max_threads, sizeof(*failures)))
    return 0;
  failures = (grid_failure *)malloc((size_t)max_threads * sizeof(*failures));
  if (failures == NULL) return 0;
  for (int thread = 0; thread < max_threads; ++thread) {
    failures[thread].index = SIZE_MAX;
    failures[thread].status = KERR_ADM_EXACT_SUCCESS;
  }
  *team_size = 1;
#pragma omp parallel default(none)                                                \
    shared(storage, dimension, grid_points, coordinate, failures, team_size)
  {
    const int thread = omp_get_thread_num();
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
            if (offset < failures[thread].index) {
              failures[thread].index = offset;
              failures[thread].status = status;
            }
            continue;
          }
          for (size_t field = 0; field < KERR_ADM_FIELD_COUNT; ++field)
            storage[field * grid_points + offset] = values[field];
        }
      }
    }
  }
  for (int thread = 0; thread < max_threads; ++thread) {
    if (failures[thread].index != SIZE_MAX &&
        (first_failure == NULL ||
         failures[thread].index < first_failure->index))
      first_failure = &failures[thread];
  }
  if (first_failure != NULL) {
    const size_t index0 = first_failure->index % dimension;
    const size_t index1 = (first_failure->index / dimension) % dimension;
    const size_t index2 = first_failure->index / (dimension * dimension);
    fprintf(stderr, "grid evaluation failed at (%a,%a,%a): %s\n",
            coordinate[0][index0], coordinate[1][index1],
            coordinate[2][index2], exact_status_name(first_failure->status));
    free(failures);
    return 0;
  }
  free(failures);
  return 1;
}

static int lagrange_full_query(const lagrange3d_reference_grid *grid,
                               const double coordinate[3], double spacing,
                               kerr_adm_value_gradient result[KERR_ADM_FIELD_COUNT],
                               size_t point) {
  static const int displacement[4] = {-2, -1, 1, 2};
  double samples[4][KERR_ADM_FIELD_COUNT];
  double center[KERR_ADM_FIELD_COUNT];
  lagrange3d_reference_status status = lagrange3d_reference_interpolate(
      grid, coordinate[0], coordinate[1], coordinate[2], center);

  if (status != LAGRANGE3D_REFERENCE_SUCCESS) {
    fprintf(stderr, "Lagrange center interpolation failed at point %zu: %s\n",
            point, lagrange3d_reference_status_name(status));
    return 0;
  }
  for (size_t field = 0; field < KERR_ADM_FIELD_COUNT; ++field)
    result[field].value = center[field];

  for (size_t axis = 0; axis < KERR_ADM_DIMENSION; ++axis) {
    for (size_t sample = 0; sample < 4; ++sample) {
      double displaced[3] = {coordinate[0], coordinate[1], coordinate[2]};
      displaced[axis] += (double)displacement[sample] * spacing;
      status = lagrange3d_reference_interpolate(
          grid, displaced[0], displaced[1], displaced[2], samples[sample]);
      if (status != LAGRANGE3D_REFERENCE_SUCCESS) {
        fprintf(stderr,
                "Lagrange offset interpolation failed at point %zu, axis %zu, "
                "offset %dH: %s\n",
                point, axis, displacement[sample],
                lagrange3d_reference_status_name(status));
        return 0;
      }
    }
    for (size_t field = 0; field < KERR_ADM_FIELD_COUNT; ++field) {
      result[field].gradient[axis] =
          (samples[0][field] - 8.0 * samples[1][field] +
           8.0 * samples[2][field] - samples[3][field]) /
          (12.0 * spacing);
    }
  }
  return 1;
}

static int analytic_fd_query(const double coordinate[3], double spacing,
                             double gradient[KERR_ADM_FIELD_COUNT][3],
                             size_t point) {
  static const int displacement[4] = {-2, -1, 1, 2};
  double samples[4][KERR_ADM_FIELD_COUNT];

  for (size_t axis = 0; axis < KERR_ADM_DIMENSION; ++axis) {
    for (size_t sample = 0; sample < 4; ++sample) {
      double displaced[3] = {coordinate[0], coordinate[1], coordinate[2]};
      displaced[axis] += (double)displacement[sample] * spacing;
      const kerr_adm_exact_status status = kerr_adm_exact_values(
          displaced[0], displaced[1], displaced[2], samples[sample]);
      if (status != KERR_ADM_EXACT_SUCCESS) {
        fprintf(stderr,
                "analytic FD sample failed at point %zu, axis %zu, offset "
                "%dH: %s\n",
                point, axis, displacement[sample], exact_status_name(status));
        return 0;
      }
    }
    for (size_t field = 0; field < KERR_ADM_FIELD_COUNT; ++field) {
      gradient[field][axis] =
          (samples[0][field] - 8.0 * samples[1][field] +
           8.0 * samples[2][field] - samples[3][field]) /
          (12.0 * spacing);
    }
  }
  return 1;
}

static int point_crosses_throat(const query_point *point, double spacing,
                                double throat) {
  for (size_t axis = 0; axis < KERR_ADM_DIMENSION; ++axis) {
    for (int multiplier = -2; multiplier <= 2; ++multiplier) {
      double displaced[3];
      double radius;
      if (multiplier == 0) continue;
      memcpy(displaced, point->coordinate, sizeof(displaced));
      displaced[axis] += (double)multiplier * spacing;
      radius = sqrt(displaced[0] * displaced[0] +
                    displaced[1] * displaced[1] +
                    displaced[2] * displaced[2]);
      if (radius < throat) return 1;
    }
  }
  return 0;
}

static int allocate_timed_results(size_t point_count, timed_results *results) {
  size_t entries;

  memset(results, 0, sizeof(*results));
  if (multiply_overflows_size_t(point_count,
                                (size_t)KERR_ADM_FIELD_COUNT))
    return 0;
  entries = point_count * (size_t)KERR_ADM_FIELD_COUNT;
  if (multiply_overflows_size_t(entries, sizeof(*results->hermite)) ||
      multiply_overflows_size_t(entries, sizeof(*results->lagrange_value)) ||
      multiply_overflows_size_t(entries, sizeof(*results->lagrange)))
    return 0;
  results->hermite =
      (hermite3d_value_gradient *)malloc(entries * sizeof(*results->hermite));
  results->lagrange_value =
      (double *)malloc(entries * sizeof(*results->lagrange_value));
  results->lagrange =
      (kerr_adm_value_gradient *)malloc(entries * sizeof(*results->lagrange));
  if (results->hermite == NULL || results->lagrange_value == NULL ||
      results->lagrange == NULL) {
    free(results->hermite);
    free(results->lagrange_value);
    free(results->lagrange);
    memset(results, 0, sizeof(*results));
    return 0;
  }
  return 1;
}

static void free_timed_results(timed_results *results) {
  free(results->hermite);
  free(results->lagrange_value);
  free(results->lagrange);
}

static int run_timed_interpolation(
    const query_point *points, size_t point_count, double spacing,
    const hermite3d_grid *hermite_grid,
    const double *const functions[KERR_ADM_FIELD_COUNT],
    const lagrange3d_reference_grid *lagrange_grid, timed_results *results) {
  double warm_values[KERR_ADM_FIELD_COUNT];
  kerr_adm_value_gradient warm_lagrange[KERR_ADM_FIELD_COUNT];
  hermite3d_value_gradient warm_hermite[KERR_ADM_FIELD_COUNT];
  double started;

  if (hermite3d_interpolate_value_gradient(
          hermite_grid, points[0].coordinate[0], points[0].coordinate[1],
          points[0].coordinate[2], KERR_ADM_FIELD_COUNT, functions,
          warm_hermite) != HERMITE3D_SUCCESS ||
      lagrange3d_reference_interpolate(
          lagrange_grid, points[0].coordinate[0], points[0].coordinate[1],
          points[0].coordinate[2], warm_values) !=
          LAGRANGE3D_REFERENCE_SUCCESS ||
      !lagrange_full_query(lagrange_grid, points[0].coordinate, spacing,
                           warm_lagrange, 0)) {
    fprintf(stderr, "interpolation warm-up failed\n");
    return 0;
  }

  started = omp_get_wtime();
  for (size_t point = 0; point < point_count; ++point) {
    const hermite3d_status status = hermite3d_interpolate_value_gradient(
        hermite_grid, points[point].coordinate[0], points[point].coordinate[1],
        points[point].coordinate[2], KERR_ADM_FIELD_COUNT, functions,
        &results->hermite[point * KERR_ADM_FIELD_COUNT]);
    if (status != HERMITE3D_SUCCESS) {
      fprintf(stderr, "Hermite interpolation failed at point %zu: %s\n",
              point, hermite_status_name(status));
      return 0;
    }
  }
  results->hermite_seconds = omp_get_wtime() - started;

  started = omp_get_wtime();
  for (size_t point = 0; point < point_count; ++point) {
    const lagrange3d_reference_status status =
        lagrange3d_reference_interpolate(
            lagrange_grid, points[point].coordinate[0],
            points[point].coordinate[1], points[point].coordinate[2],
            &results->lagrange_value[point * KERR_ADM_FIELD_COUNT]);
    if (status != LAGRANGE3D_REFERENCE_SUCCESS) {
      fprintf(stderr, "Lagrange value interpolation failed at point %zu: %s\n",
              point, lagrange3d_reference_status_name(status));
      return 0;
    }
  }
  results->lagrange_value_seconds = omp_get_wtime() - started;

  started = omp_get_wtime();
  for (size_t point = 0; point < point_count; ++point) {
    if (!lagrange_full_query(
            lagrange_grid, points[point].coordinate, spacing,
            &results->lagrange[point * KERR_ADM_FIELD_COUNT], point))
      return 0;
  }
  results->lagrange_full_seconds = omp_get_wtime() - started;
  return 1;
}

static int accumulate_errors(const query_point *points, size_t point_count,
                             double spacing, const timed_results *results,
                             level_summary *summary) {
  norm_accumulator field[METHOD_COUNT][QUANTITY_COUNT][KERR_ADM_FIELD_COUNT] =
      {0};
  norm_accumulator aggregate[METHOD_COUNT][QUANTITY_COUNT] = {0};

  for (size_t point = 0; point < point_count; ++point) {
    double analytic_fd[KERR_ADM_FIELD_COUNT][3];
    if (!analytic_fd_query(points[point].coordinate, spacing, analytic_fd,
                           point))
      return 0;
    for (size_t field_index = 0; field_index < KERR_ADM_FIELD_COUNT;
         ++field_index) {
      const kerr_adm_value_gradient *const exact =
          &points[point].exact[field_index];
      const hermite3d_value_gradient *const hermite =
          &results->hermite[point * KERR_ADM_FIELD_COUNT + field_index];
      const kerr_adm_value_gradient *const lagrange =
          &results->lagrange[point * KERR_ADM_FIELD_COUNT + field_index];
      const double lagrange_value =
          results->lagrange_value[point * KERR_ADM_FIELD_COUNT + field_index];
      const long double hermite_value_error =
          (long double)hermite->value - (long double)exact->value;
      const long double lagrange_value_error =
          (long double)lagrange_value - (long double)exact->value;

      if (!isfinite(hermite->value) || !isfinite(lagrange_value) ||
          !isfinite(hermite_value_error) || !isfinite(lagrange_value_error)) {
        fprintf(stderr, "nonfinite value result at point %zu, field %s\n",
                point, kerr_adm_field_name((kerr_adm_field)field_index));
        return 0;
      }
      add_norm_sample(&field[METHOD_HERMITE][0][field_index],
                      hermite_value_error, exact->value);
      add_norm_sample(&aggregate[METHOD_HERMITE][0], hermite_value_error,
                      exact->value);
      add_norm_sample(&field[METHOD_LAGRANGE_FD][0][field_index],
                      lagrange_value_error, exact->value);
      add_norm_sample(&aggregate[METHOD_LAGRANGE_FD][0], lagrange_value_error,
                      exact->value);

      for (size_t axis = 0; axis < KERR_ADM_DIMENSION; ++axis) {
        const long double hermite_error =
            (long double)hermite->gradient[axis] - exact->gradient[axis];
        const long double lagrange_error =
            (long double)lagrange->gradient[axis] - exact->gradient[axis];
        const long double baseline_error =
            (long double)analytic_fd[field_index][axis] - exact->gradient[axis];
        if (!isfinite(hermite->gradient[axis]) ||
            !isfinite(lagrange->gradient[axis]) ||
            !isfinite(analytic_fd[field_index][axis]) ||
            !isfinite(hermite_error) || !isfinite(lagrange_error) ||
            !isfinite(baseline_error)) {
          fprintf(stderr,
                  "nonfinite derivative at point %zu, field %s, axis %zu\n",
                  point, kerr_adm_field_name((kerr_adm_field)field_index),
                  axis);
          return 0;
        }
        add_norm_sample(&field[METHOD_HERMITE][axis + 1][field_index],
                        hermite_error, exact->gradient[axis]);
        add_norm_sample(&aggregate[METHOD_HERMITE][axis + 1], hermite_error,
                        exact->gradient[axis]);
        add_norm_sample(&field[METHOD_LAGRANGE_FD][axis + 1][field_index],
                        lagrange_error, exact->gradient[axis]);
        add_norm_sample(&aggregate[METHOD_LAGRANGE_FD][axis + 1],
                        lagrange_error, exact->gradient[axis]);
        add_norm_sample(&field[METHOD_ANALYTIC_FD][axis + 1][field_index],
                        baseline_error, exact->gradient[axis]);
        add_norm_sample(&aggregate[METHOD_ANALYTIC_FD][axis + 1],
                        baseline_error, exact->gradient[axis]);
      }
    }
  }

  for (size_t method = 0; method < METHOD_COUNT; ++method) {
    for (size_t quantity = 0; quantity < QUANTITY_COUNT; ++quantity) {
      if (method == METHOD_ANALYTIC_FD && quantity == 0) continue;
      for (size_t field_index = 0; field_index < KERR_ADM_FIELD_COUNT;
           ++field_index) {
        summary->field[method][quantity][field_index] =
            finish_norm(&field[method][quantity][field_index], point_count);
      }
      summary->aggregate[method][quantity] = finish_norm(
          &aggregate[method][quantity],
          point_count * (size_t)KERR_ADM_FIELD_COUNT);
    }
  }
  return 1;
}

static int evaluate_level(const options *arguments, const query_point *points,
                          double minimum_radius, size_t resolution,
                          double half_width, double throat,
                          level_summary *summary) {
  size_t dimension;
  size_t grid_points;
  size_t total_values;
  double *storage = NULL;
  double *coordinate[3] = {NULL, NULL, NULL};
  const double *functions[KERR_ADM_FIELD_COUNT];
  hermite3d_grid hermite_grid;
  lagrange3d_reference_grid lagrange_grid;
  timed_results results;
  const double spacing = 2.0 * half_width / (double)resolution;
  const double start = -half_width - 5.5 * spacing;

  memset(summary, 0, sizeof(*summary));
  if (!isfinite(spacing) || spacing <= 0.0 || !isfinite(start) ||
      !compute_grid_sizes(resolution, &dimension, &grid_points,
                          &total_values)) {
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
    free(storage);
    return 0;
  }
  if (!fill_grid(storage, dimension, grid_points, coordinate,
                 &summary->grid_threads)) {
    free_coordinates(coordinate);
    free(storage);
    return 0;
  }
  for (size_t field = 0; field < KERR_ADM_FIELD_COUNT; ++field)
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

  for (size_t axis = 0; axis < KERR_ADM_DIMENSION; ++axis) {
    lagrange_grid.dimension[axis] = (int)dimension;
    lagrange_grid.spacing[axis] = spacing;
    lagrange_grid.coordinate[axis] = coordinate[axis];
  }
  for (size_t field = 0; field < KERR_ADM_FIELD_COUNT; ++field)
    lagrange_grid.function[field] = functions[field];

  if (!allocate_timed_results(arguments->point_count, &results)) {
    fprintf(stderr, "could not allocate interpolation result buffers\n");
    free_coordinates(coordinate);
    free(storage);
    return 0;
  }
  if (!run_timed_interpolation(points, arguments->point_count, spacing,
                               &hermite_grid, functions, &lagrange_grid,
                               &results) ||
      !accumulate_errors(points, arguments->point_count, spacing, &results,
                         summary)) {
    free_timed_results(&results);
    free_coordinates(coordinate);
    free(storage);
    return 0;
  }

  summary->resolution = resolution;
  summary->stored_dimension = dimension;
  summary->spacing = spacing;
  summary->minimum_radius = minimum_radius;
  summary->memory_mib =
      (long double)total_values * sizeof(*storage) / (1024.0L * 1024.0L);
  summary->hermite_seconds_per_point =
      results.hermite_seconds / (double)arguments->point_count;
  summary->lagrange_value_seconds_per_point =
      results.lagrange_value_seconds / (double)arguments->point_count;
  summary->lagrange_full_seconds_per_point =
      results.lagrange_full_seconds / (double)arguments->point_count;
  for (size_t point = 0; point < arguments->point_count; ++point)
    if (point_crosses_throat(&points[point], spacing, throat))
      ++summary->throat_crossing_points;

  free_timed_results(&results);
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
           logl((long double)coarse_spacing / fine_spacing);
  return isfinite(*order);
}

static void print_order(long double coarse, long double fine,
                        double coarse_spacing, double fine_spacing) {
  long double order;
  if (measured_order(coarse, fine, coarse_spacing, fine_spacing, &order))
    printf(" %8.3Lf", order);
  else
    printf(" %8s", "--");
}

static void print_relative(long double value, int defined) {
  if (defined)
    printf(" %12.4Le", value);
  else
    printf(" %12s", "--");
}

static void print_norm_rows(const level_summary *levels, size_t level_count,
                            comparison_method method, size_t quantity,
                            size_t field, int aggregate) {
  const char *name = aggregate ? "ALL" : kerr_adm_field_name((kerr_adm_field)field);
  for (size_t level = 0; level < level_count; ++level) {
    const error_norm norm = aggregate
                                ? levels[level].aggregate[method][quantity]
                                : levels[level].field[method][quantity][field];
    printf("%-9s %6zu %11.4e %12.4Le %12.4Le", name,
           levels[level].resolution, levels[level].spacing, norm.rms,
           norm.maximum);
    print_relative(norm.relative_rms, norm.relative_rms_defined);
    print_relative(norm.relative_maximum, norm.relative_maximum_defined);
    if (level == 0) {
      printf(" %8s %8s\n", "--", "--");
    } else {
      const error_norm previous =
          aggregate ? levels[level - 1].aggregate[method][quantity]
                    : levels[level - 1].field[method][quantity][field];
      print_order(previous.rms, norm.rms, levels[level - 1].spacing,
                  levels[level].spacing);
      print_order(previous.maximum, norm.maximum,
                  levels[level - 1].spacing, levels[level].spacing);
      putchar('\n');
    }
  }
}

static void print_results(const options *arguments, const level_summary *levels,
                          double throat, double half_width) {
  printf("Kerr ADM Hermite-Lagrange interpolation comparison\n");
  printf("  M = 1, a = 0.5, throat radius s = %.17g\n", throat);
  printf("  active cube = [%.17g, %.17g]^3; six analytic ghost layers\n",
         -half_width, half_width);
  printf("  random exterior centers = %zu, seed = 0x%016" PRIx64 "\n",
         arguments->point_count, arguments->seed);
  printf("  lapse convention = nonnegative LES lapse on both sheets\n");
  printf("  Lagrange FD spacing equals the local grid spacing H\n");

  printf("\nGrid diagnostics\n");
  printf("     N stored           H  memory(MiB)       min r    min r/H"
         "  FD crosses  grid thr\n");
  for (size_t level = 0; level < arguments->resolution_count; ++level) {
    printf("%6zu %6zu %11.4e %12.3Lf %11.4e %10.3e %11zu %9d\n",
           levels[level].resolution, levels[level].stored_dimension,
           levels[level].spacing, levels[level].memory_mib,
           levels[level].minimum_radius,
           levels[level].minimum_radius / levels[level].spacing,
           levels[level].throat_crossing_points,
           levels[level].grid_threads);
  }

  printf("\nSerial interpolation latency\n");
  printf("     N  Hermite value+grad [us]  Lagrange value [us]"
         "  Lagrange value+FD [us]  ratio\n");
  for (size_t level = 0; level < arguments->resolution_count; ++level) {
    const double ratio = levels[level].lagrange_full_seconds_per_point /
                         levels[level].hermite_seconds_per_point;
    printf("%6zu %24.6f %20.6f %25.6f %8.3f\n",
           levels[level].resolution,
           1.0e6 * levels[level].hermite_seconds_per_point,
           1.0e6 * levels[level].lagrange_value_seconds_per_point,
           1.0e6 * levels[level].lagrange_full_seconds_per_point, ratio);
  }

  printf("\nErrors are finite-cloud sampled absolute and scale-normalized norms.\n");
  printf("Smooth-region expectations: Hermite value p=5, Hermite gradient p=4,\n");
  printf("Lagrange value p=7, and Lagrange+FD gradient p=4.  The lapse cusp\n");
  printf("can disrupt these rates when a stencil crosses r=s; no threshold is enforced.\n");

  for (size_t quantity = 0; quantity < QUANTITY_COUNT; ++quantity) {
    for (size_t method = 0; method < METHOD_COUNT; ++method) {
      if (method == METHOD_ANALYTIC_FD && quantity == 0) continue;
      printf("\n%s / %s\n", method_name((comparison_method)method),
             quantity_name(quantity));
      printf("field          N           H       abs L2     abs Linf"
             "       rel L2     rel Linf    p(L2)  p(Linf)\n");
      for (size_t field = 0; field < KERR_ADM_FIELD_COUNT; ++field)
        print_norm_rows(levels, arguments->resolution_count,
                        (comparison_method)method, quantity, field, 0);
      print_norm_rows(levels, arguments->resolution_count,
                      (comparison_method)method, quantity, 0, 1);
    }
  }
}

int main(int argc, char **argv) {
  options arguments;
  const int parsed = parse_options(argc, argv, &arguments);
  const double throat = kerr_adm_exact_throat_radius();
  const double half_width = 5.0 * throat;
  query_point *points;
  level_summary *levels;
  double minimum_radius;

  if (parsed == 1) return EXIT_SUCCESS;
  if (parsed != 2) {
    print_usage(stderr, argv[0]);
    return EXIT_FAILURE;
  }
  if (multiply_overflows_size_t(arguments.resolution_count, sizeof(*levels))) {
    fprintf(stderr, "too many resolution levels\n");
    free(arguments.resolutions);
    return EXIT_FAILURE;
  }
  levels = (level_summary *)calloc(arguments.resolution_count, sizeof(*levels));
  if (levels == NULL) {
    fprintf(stderr, "could not allocate level summaries\n");
    free(arguments.resolutions);
    return EXIT_FAILURE;
  }
  points = generate_queries(arguments.point_count, arguments.seed, half_width,
                            throat, &minimum_radius);
  if (points == NULL) {
    fprintf(stderr, "could not construct the query cloud\n");
    free(levels);
    free(arguments.resolutions);
    return EXIT_FAILURE;
  }

  for (size_t level = 0; level < arguments.resolution_count; ++level) {
    fprintf(stderr, "evaluating N=%zu ...\n", arguments.resolutions[level]);
    if (!evaluate_level(&arguments, points, minimum_radius,
                        arguments.resolutions[level], half_width, throat,
                        &levels[level])) {
      free(points);
      free(levels);
      free(arguments.resolutions);
      return EXIT_FAILURE;
    }
  }
  print_results(&arguments, levels, throat, half_width);
  free(points);
  free(levels);
  free(arguments.resolutions);
  return EXIT_SUCCESS;
}
