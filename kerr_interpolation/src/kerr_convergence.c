#include "hermite3d.h"
#include "kerr_exact.h"

#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum {
  QUANTITY_COUNT = 4,
  GHOST_LAYERS = 3
};

static const uint64_t DEFAULT_SEED = UINT64_C(0x4b4552524845524d);
static const size_t DEFAULT_Z_PROFILE_SAMPLES = 1000;
static const char DEFAULT_Z_PROFILE_OUTPUT[] = "build/z_profile.csv";

typedef struct options {
  size_t *resolutions;
  size_t resolution_count;
  size_t point_count;
  uint64_t seed;
  int z_profile_enabled;
  size_t z_profile_component;
  size_t z_profile_quantity;
  size_t z_profile_samples;
  const char *z_profile_output;
} options;

/*
 * The scaled representation scale*sqrt(sum_squares) is the same strategy
 * used by robust BLAS norm implementations.  It avoids forming error^2 when
 * a finite interpolation error is large enough that its square would overflow
 * even though the requested norm is still representable.
 */
typedef struct scaled_sum_squares {
  long double scale;
  long double sum_squares;
} scaled_sum_squares;

typedef struct error_accumulator {
  scaled_sum_squares squares;
  long double maximum;
} error_accumulator;

typedef struct error_norm {
  long double rms;
  long double maximum;
} error_norm;

typedef struct level_summary {
  size_t resolution;
  size_t stored_dimension;
  double spacing;
  double minimum_sample_radius;
  double minimum_grid_radius;
  error_norm component[QUANTITY_COUNT][KERR_EXACT_COMPONENT_COUNT];
  error_norm aggregate[QUANTITY_COUNT];
} level_summary;

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

/*
 * Return the midpoint of one of 2^52 equal bins in (0,1).  The numerator is
 * odd, so the binary64 result can never equal 1/2.  Applying this independently
 * to x, y, and z therefore cannot produce the excluded puncture (0,0,0), and
 * it needs no rejection or resampling that would bias the cloud.
 */
static double uniform_open(uint64_t *state) {
  const uint64_t bin = splitmix64_next(state) >> 12;
  return (double)(2 * bin + 1) * 0x1p-53;
}

static void add_squared_error(scaled_sum_squares *sum, long double value) {
  const long double magnitude = fabsl(value);

  if (magnitude == 0.0L) return;
  if (sum->scale < magnitude) {
    const long double ratio = sum->scale / magnitude;
    sum->sum_squares = 1.0L + sum->sum_squares * ratio * ratio;
    sum->scale = magnitude;
  } else {
    const long double ratio = magnitude / sum->scale;
    sum->sum_squares += ratio * ratio;
  }
}

static void add_error(error_accumulator *accumulator, long double error) {
  const long double magnitude = fabsl(error);

  add_squared_error(&accumulator->squares, error);
  if (magnitude > accumulator->maximum) accumulator->maximum = magnitude;
}

static error_norm finish_norm(const error_accumulator *accumulator,
                              size_t count) {
  error_norm result;

  if (accumulator->squares.scale == 0.0L) {
    result.rms = 0.0L;
  } else {
    result.rms = accumulator->squares.scale *
                 sqrtl(accumulator->squares.sum_squares / (long double)count);
  }
  result.maximum = accumulator->maximum;
  return result;
}

static const char *exact_status_name(kerr_exact_status status) {
  switch (status) {
    case KERR_EXACT_SUCCESS:
      return "success";
    case KERR_EXACT_INVALID_ARGUMENT:
      return "invalid argument";
    case KERR_EXACT_DOMAIN_ERROR:
      return "coordinate-domain error";
    case KERR_EXACT_NONFINITE_RESULT:
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
          "Usage: %s --resolutions N1,N2,... --points COUNT [options]\n"
          "\n"
          "N values must be positive, even, distinct, and strictly increasing.\n"
          "One resolution reports errors; later resolutions also report pairwise\n"
          "orders. UINT64 accepts decimal or a 0x-prefixed hexadecimal value.\n"
          "\n"
          "Options:\n"
          "  --seed UINT64              random seed\n"
          "  --z-profile COMP:QUANTITY export a z-axis profile, where COMP is\n"
          "                             tt,tx,ty,tz,xx,xy,xz,yy,yz,zz and\n"
          "                             QUANTITY is value,dx,dy,dz\n"
          "  --z-samples COUNT          even profile sample count (default 1000)\n"
          "  --z-output PATH            profile CSV (default build/z_profile.csv)\n",
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
  if (errno != 0 || end == text || *end != '\0' || parsed > UINT64_MAX) {
    return 0;
  }
  *value = (uint64_t)parsed;
  return 1;
}

static int parse_resolutions(const char *text, size_t **values,
                             size_t *count) {
  size_t entries = 1;
  size_t length;
  size_t *parsed_values;
  char *copy;
  char *cursor;

  if (text == NULL || *text == '\0') return 0;
  for (const char *character = text; *character != '\0'; ++character) {
    if (*character == ',') ++entries;
  }
  if (multiply_overflows_size_t(entries, sizeof(*parsed_values))) return 0;

  length = strlen(text);
  if (add_overflows_size_t(length, 1)) return 0;
  copy = (char *)malloc(length + 1);
  parsed_values = (size_t *)malloc(entries * sizeof(*parsed_values));
  if (copy == NULL || parsed_values == NULL) {
    free(copy);
    free(parsed_values);
    return 0;
  }
  memcpy(copy, text, length + 1);

  cursor = copy;
  for (size_t entry = 0; entry < entries; ++entry) {
    char *comma = strchr(cursor, ',');
    if (comma != NULL) *comma = '\0';
    if (!parse_size(cursor, &parsed_values[entry]) ||
        (parsed_values[entry] & (size_t)1) != 0 ||
        (entry > 0 && parsed_values[entry] <= parsed_values[entry - 1]) ||
        (entry + 1 < entries && comma == NULL)) {
      free(copy);
      free(parsed_values);
      return 0;
    }
    cursor = comma == NULL ? cursor + strlen(cursor) : comma + 1;
  }

  free(copy);
  *values = parsed_values;
  *count = entries;
  return 1;
}

static const char *attached_value(const char *argument, const char *name) {
  const size_t length = strlen(name);
  if (strncmp(argument, name, length) == 0 && argument[length] == '=')
    return argument + length + 1;
  return NULL;
}

static const char *quantity_name(size_t quantity) {
  static const char *const names[QUANTITY_COUNT] = {"value", "dx", "dy",
                                                    "dz"};
  return quantity < QUANTITY_COUNT ? names[quantity] : NULL;
}

static int parse_z_profile(const char *text, size_t *component,
                           size_t *quantity) {
  const char *separator;
  size_t component_length;

  if (text == NULL) return 0;
  separator = strchr(text, ':');
  if (separator == NULL || separator == text || separator[1] == '\0' ||
      strchr(separator + 1, ':') != NULL) {
    return 0;
  }
  component_length = (size_t)(separator - text);
  for (size_t candidate = 0; candidate < KERR_EXACT_COMPONENT_COUNT;
       ++candidate) {
    const char *name =
        kerr_exact_component_name((kerr_exact_component)candidate);
    if (strlen(name) == component_length &&
        strncmp(text, name, component_length) == 0) {
      *component = candidate;
      for (size_t entry = 0; entry < QUANTITY_COUNT; ++entry) {
        if (strcmp(separator + 1, quantity_name(entry)) == 0) {
          *quantity = entry;
          return 1;
        }
      }
      return 0;
    }
  }
  return 0;
}

static int parse_options(int argc, char **argv, options *result) {
  const char *resolution_text = NULL;
  const char *point_text = NULL;
  const char *seed_text = NULL;
  const char *z_profile_text = NULL;
  const char *z_samples_text = NULL;
  const char *z_output_text = NULL;

  result->resolutions = NULL;
  result->resolution_count = 0;
  result->point_count = 0;
  result->seed = DEFAULT_SEED;
  result->z_profile_enabled = 0;
  result->z_profile_component = 0;
  result->z_profile_quantity = 0;
  result->z_profile_samples = DEFAULT_Z_PROFILE_SAMPLES;
  result->z_profile_output = DEFAULT_Z_PROFILE_OUTPUT;

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

    value = attached_value(argv[argument], "--z-profile");
    if (strcmp(argv[argument], "--z-profile") == 0) {
      if (++argument >= argc) return 0;
      value = argv[argument];
    }
    if (value != NULL) {
      if (z_profile_text != NULL) return 0;
      z_profile_text = value;
      continue;
    }

    value = attached_value(argv[argument], "--z-samples");
    if (strcmp(argv[argument], "--z-samples") == 0) {
      if (++argument >= argc) return 0;
      value = argv[argument];
    }
    if (value != NULL) {
      if (z_samples_text != NULL) return 0;
      z_samples_text = value;
      continue;
    }

    value = attached_value(argv[argument], "--z-output");
    if (strcmp(argv[argument], "--z-output") == 0) {
      if (++argument >= argc) return 0;
      value = argv[argument];
    }
    if (value != NULL) {
      if (z_output_text != NULL || *value == '\0') return 0;
      z_output_text = value;
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

  if (z_profile_text == NULL) {
    if (z_samples_text != NULL || z_output_text != NULL) {
      free(result->resolutions);
      result->resolutions = NULL;
      return 0;
    }
  } else {
    result->z_profile_enabled = 1;
    if (!parse_z_profile(z_profile_text, &result->z_profile_component,
                         &result->z_profile_quantity) ||
        (z_samples_text != NULL &&
         !parse_size(z_samples_text, &result->z_profile_samples)) ||
        (result->z_profile_samples & (size_t)1) != 0 ||
        (uintmax_t)result->z_profile_samples > UINT64_C(4503599627370496)) {
      free(result->resolutions);
      result->resolutions = NULL;
      return 0;
    }
    if (z_output_text != NULL) result->z_profile_output = z_output_text;
  }
  return 2;
}

static int compute_grid_sizes(size_t resolution, size_t *stored_dimension,
                              size_t *grid_points, size_t *total_values) {
  size_t square;

  if (add_overflows_size_t(resolution, 2 * (size_t)GHOST_LAYERS)) return 0;
  *stored_dimension = resolution + 2 * (size_t)GHOST_LAYERS;
  if (multiply_overflows_size_t(*stored_dimension, *stored_dimension)) return 0;
  square = *stored_dimension * *stored_dimension;
  if (multiply_overflows_size_t(square, *stored_dimension)) return 0;
  *grid_points = square * *stored_dimension;
  if (multiply_overflows_size_t(*grid_points,
                                (size_t)KERR_EXACT_COMPONENT_COUNT)) {
    return 0;
  }
  *total_values = *grid_points * (size_t)KERR_EXACT_COMPONENT_COUNT;
  if (multiply_overflows_size_t(*total_values, sizeof(double))) return 0;
  return 1;
}

static int fill_grid(double *storage, size_t dimension, size_t grid_points,
                     double start, double spacing) {
  for (size_t index2 = 0; index2 < dimension; ++index2) {
    const double z = start + (double)index2 * spacing;
    for (size_t index1 = 0; index1 < dimension; ++index1) {
      const double y = start + (double)index1 * spacing;
      for (size_t index0 = 0; index0 < dimension; ++index0) {
        const double x = start + (double)index0 * spacing;
        const size_t offset = index0 + dimension * (index1 + dimension * index2);
        double values[KERR_EXACT_COMPONENT_COUNT];
        const kerr_exact_status status = kerr_exact_metric(x, y, z, values);

        if (status != KERR_EXACT_SUCCESS) {
          fprintf(stderr,
                  "exact grid evaluation failed at (%a, %a, %a): %s\n",
                  x, y, z, exact_status_name(status));
          return 0;
        }
        for (size_t component = 0; component < KERR_EXACT_COMPONENT_COUNT;
             ++component) {
          storage[component * grid_points + offset] = values[component];
        }
      }
    }
  }
  return 1;
}

/*
 * Export a diagnostic through the rotation axis without changing the random
 * convergence sample.  Half of the midpoint samples lie on each exterior
 * branch, z <= -s and z >= s.  No profile point lies in the lower LES sheet
 * |z| < s, and the plotting script does not connect the separated branches.
 */
static int write_z_profile(
    FILE *stream, const options *arguments, size_t resolution,
    double half_width, double throat, const hermite3d_grid *grid,
    const double *const functions[KERR_EXACT_COMPONENT_COUNT]) {
  const char *component_name = kerr_exact_component_name(
      (kerr_exact_component)arguments->z_profile_component);
  const char *selected_quantity =
      quantity_name(arguments->z_profile_quantity);
  const size_t branch_samples = arguments->z_profile_samples / 2;
  const long double branch_width =
      (long double)half_width - (long double)throat;

  for (size_t sample = 0; sample < arguments->z_profile_samples; ++sample) {
    const int negative_branch = sample < branch_samples;
    const size_t branch_index =
        negative_branch ? sample : sample - branch_samples;
    const long double fraction =
        ((long double)branch_index + 0.5L) / (long double)branch_samples;
    const double z = negative_branch
                         ? (double)(-(long double)half_width +
                                    branch_width * fraction)
                         : (double)((long double)throat +
                                    branch_width * fraction);
    kerr_exact_value_gradient exact[KERR_EXACT_COMPONENT_COUNT];
    hermite3d_value_gradient interpolated[KERR_EXACT_COMPONENT_COUNT];
    double exact_value;
    double interpolated_value;
    long double error;
    const kerr_exact_status exact_status =
        kerr_exact_metric_gradient(0.0, 0.0, z, exact);
    const hermite3d_status interpolation_status =
        hermite3d_interpolate_value_gradient(
            grid, 0.0, 0.0, z, KERR_EXACT_COMPONENT_COUNT, functions,
            interpolated);

    if (fabs(z) < throat || exact_status != KERR_EXACT_SUCCESS) {
      fprintf(stderr,
              "exact z-profile evaluation failed at sample %zu, z=%a: %s\n",
              sample, z, exact_status_name(exact_status));
      return 0;
    }
    if (interpolation_status != HERMITE3D_SUCCESS) {
      fprintf(stderr,
              "z-profile interpolation failed at sample %zu, z=%a: %s\n",
              sample, z, hermite_status_name(interpolation_status));
      return 0;
    }

    if (arguments->z_profile_quantity == 0) {
      exact_value = exact[arguments->z_profile_component].value;
      interpolated_value =
          interpolated[arguments->z_profile_component].value;
    } else {
      const size_t axis = arguments->z_profile_quantity - 1;
      exact_value = exact[arguments->z_profile_component].gradient[axis];
      interpolated_value =
          interpolated[arguments->z_profile_component].gradient[axis];
    }
    error = (long double)interpolated_value - (long double)exact_value;
    if (!isfinite(exact_value) || !isfinite(interpolated_value) ||
        !isfinite(error)) {
      fprintf(stderr,
              "nonfinite z-profile value for %s:%s at sample %zu, z=%a\n",
              component_name, selected_quantity, sample, z);
      return 0;
    }
    if (fprintf(stream, "%s,%s,%zu,%.17g,%.17g,%.17g,%.21Lg\n",
                component_name, selected_quantity, resolution, z, exact_value,
                interpolated_value, error) < 0) {
      fprintf(stderr, "could not write z-profile CSV\n");
      return 0;
    }
  }
  return 1;
}

static int evaluate_level(const options *arguments, size_t resolution,
                          double half_width, double throat,
                          FILE *z_profile_stream,
                          level_summary *summary) {
  size_t dimension;
  size_t grid_points;
  size_t total_values;
  double *storage;
  const double spacing = 2.0 * half_width / (double)resolution;
  const double start = -half_width - 2.5 * spacing;
  const double minimum_grid_radius = sqrt(3.0) * spacing / 2.0;
  const double *functions[KERR_EXACT_COMPONENT_COUNT];
  hermite3d_grid grid;
  error_accumulator accumulators[QUANTITY_COUNT][KERR_EXACT_COMPONENT_COUNT] =
      {0};
  error_accumulator aggregate[QUANTITY_COUNT] = {0};
  uint64_t random_state = arguments->seed;
  double minimum_sample_radius = INFINITY;

  if (!isfinite(spacing) || spacing <= 0.0 || !isfinite(start) ||
      !compute_grid_sizes(resolution, &dimension, &grid_points, &total_values)) {
    fprintf(stderr, "resolution %zu overflows grid geometry or allocation\n",
            resolution);
    return 0;
  }

  storage = (double *)malloc(total_values * sizeof(*storage));
  if (storage == NULL) {
    fprintf(stderr, "could not allocate %.3Lf MiB for resolution %zu\n",
            (long double)total_values * (long double)sizeof(*storage) /
                (1024.0L * 1024.0L),
            resolution);
    return 0;
  }
  if (!fill_grid(storage, dimension, grid_points, start, spacing)) {
    free(storage);
    return 0;
  }
  for (size_t component = 0; component < KERR_EXACT_COMPONENT_COUNT;
       ++component) {
    functions[component] = storage + component * grid_points;
  }

  grid.nxx0 = dimension;
  grid.nxx1 = dimension;
  grid.nxx2 = dimension;
  grid.xx0_start = start;
  grid.xx1_start = start;
  grid.xx2_start = start;
  grid.dxx0 = spacing;
  grid.dxx1 = spacing;
  grid.dxx2 = spacing;
  grid.stride0 = 1;
  grid.stride1 = dimension;
  grid.stride2 = dimension * dimension;

  if (arguments->z_profile_enabled &&
      !write_z_profile(z_profile_stream, arguments, resolution, half_width,
                       throat, &grid, functions)) {
    free(storage);
    return 0;
  }

  for (size_t point = 0; point < arguments->point_count; ++point) {
    const double x = half_width * (2.0 * uniform_open(&random_state) - 1.0);
    const double y = half_width * (2.0 * uniform_open(&random_state) - 1.0);
    const double z = half_width * (2.0 * uniform_open(&random_state) - 1.0);
    const double radius = sqrt(x * x + y * y + z * z);
    kerr_exact_value_gradient exact[KERR_EXACT_COMPONENT_COUNT];
    hermite3d_value_gradient interpolated[KERR_EXACT_COMPONENT_COUNT];
    const kerr_exact_status exact_status =
        kerr_exact_metric_gradient(x, y, z, exact);
    const hermite3d_status interpolation_status =
        hermite3d_interpolate_value_gradient(
            &grid, x, y, z, KERR_EXACT_COMPONENT_COUNT, functions,
            interpolated);

    if (radius < minimum_sample_radius) minimum_sample_radius = radius;
    if (exact_status != KERR_EXACT_SUCCESS) {
      fprintf(stderr,
              "exact query evaluation failed at point %zu (%a, %a, %a): %s\n",
              point, x, y, z, exact_status_name(exact_status));
      free(storage);
      return 0;
    }
    if (interpolation_status != HERMITE3D_SUCCESS) {
      fprintf(stderr,
              "interpolation failed at point %zu (%a, %a, %a): %s\n",
              point, x, y, z, hermite_status_name(interpolation_status));
      free(storage);
      return 0;
    }

    for (size_t component = 0; component < KERR_EXACT_COMPONENT_COUNT;
         ++component) {
      const long double value_error =
          (long double)interpolated[component].value -
          (long double)exact[component].value;
      if (!isfinite(value_error)) {
        fprintf(stderr,
                "nonfinite value error at point %zu, component %s, r=%a\n",
                point,
                kerr_exact_component_name((kerr_exact_component)component),
                radius);
        free(storage);
        return 0;
      }
      add_error(&accumulators[0][component], value_error);
      add_error(&aggregate[0], value_error);

      for (size_t axis = 0; axis < KERR_EXACT_DIMENSION; ++axis) {
        const long double derivative_error =
            (long double)interpolated[component].gradient[axis] -
            (long double)exact[component].gradient[axis];
        if (!isfinite(derivative_error)) {
          fprintf(stderr,
                  "nonfinite derivative error at point %zu, component %s, "
                  "axis %zu, r=%a\n",
                  point,
                  kerr_exact_component_name((kerr_exact_component)component),
                  axis, radius);
          free(storage);
          return 0;
        }
        add_error(&accumulators[axis + 1][component], derivative_error);
        add_error(&aggregate[axis + 1], derivative_error);
      }
    }
  }

  summary->resolution = resolution;
  summary->stored_dimension = dimension;
  summary->spacing = spacing;
  summary->minimum_sample_radius = minimum_sample_radius;
  summary->minimum_grid_radius = minimum_grid_radius;
  for (size_t quantity = 0; quantity < QUANTITY_COUNT; ++quantity) {
    for (size_t component = 0; component < KERR_EXACT_COMPONENT_COUNT;
         ++component) {
      summary->component[quantity][component] =
          finish_norm(&accumulators[quantity][component],
                      arguments->point_count);
    }
    summary->aggregate[quantity] =
        finish_norm(&aggregate[quantity],
                    arguments->point_count *
                        (size_t)KERR_EXACT_COMPONENT_COUNT);
  }

  free(storage);
  return 1;
}

static int measured_order(long double coarse, long double fine,
                          double coarse_spacing, double fine_spacing,
                          long double *order) {
  if (!(coarse > 0.0L) || !(fine > 0.0L) || !isfinite(coarse) ||
      !isfinite(fine)) {
    return 0;
  }
  *order = logl(coarse / fine) /
           logl((long double)coarse_spacing / (long double)fine_spacing);
  return isfinite(*order);
}

static void print_order(long double coarse, long double fine,
                        double coarse_spacing, double fine_spacing) {
  long double order;
  if (measured_order(coarse, fine, coarse_spacing, fine_spacing, &order))
    printf(" %9.4Lf", order);
  else
    printf(" %9s", "--");
}

static void print_norm_row(const char *name, const level_summary *levels,
                           size_t level_count, size_t quantity,
                           size_t component, int aggregate_row) {
  for (size_t level = 0; level < level_count; ++level) {
    const error_norm norm = aggregate_row
                                ? levels[level].aggregate[quantity]
                                : levels[level].component[quantity][component];
    printf("%-4s %8zu %12.5e %14.6Le %14.6Le", name,
           levels[level].resolution, levels[level].spacing, norm.rms,
           norm.maximum);
    if (level == 0) {
      printf(" %9s %9s\n", "--", "--");
    } else {
      const error_norm previous =
          aggregate_row
              ? levels[level - 1].aggregate[quantity]
              : levels[level - 1].component[quantity][component];
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
  static const char *const headings[QUANTITY_COUNT] = {
      "VALUE", "D/DX", "D/DY", "D/DZ"};

  printf("Kerr metric / quintic Hermite interpolation convergence\n");
  printf("  M = 1, a = 0.5, throat radius s = %.17g\n", throat);
  printf("  active cube = [%.17g, %.17g]^3\n", -half_width,
         half_width);
  printf("  side length = 10s = five horizon diameters\n");
  printf("  random points = %zu, seed = 0x%016" PRIx64 "\n",
         arguments->point_count, arguments->seed);
  printf("  random cloud is uniform, unfiltered, and identical at every level\n");
  if (arguments->z_profile_enabled) {
    printf("  z profile = %s:%s, %zu midpoint samples, %s\n",
           kerr_exact_component_name(
               (kerr_exact_component)arguments->z_profile_component),
           quantity_name(arguments->z_profile_quantity),
           arguments->z_profile_samples, arguments->z_profile_output);
  }
  printf("\nGrid diagnostics\n");
  printf("       N   stored              h          min r      min r/h  min grid r\n");
  for (size_t level = 0; level < arguments->resolution_count; ++level) {
    printf("%8zu %8zu %14.6e %14.6e %12.5e %12.5e\n",
           levels[level].resolution, levels[level].stored_dimension,
           levels[level].spacing, levels[level].minimum_sample_radius,
           levels[level].minimum_sample_radius / levels[level].spacing,
           levels[level].minimum_grid_radius);
  }

  printf("\nThe norms below are finite-cloud sampled RMS (L2) and sampled maxima\n");
  printf("(Linf), not continuum norms over the punctured cube. Expected local\n");
  printf("smooth-function orders are 5 for values and 4 for gradients; no order\n");
  printf("threshold is enforced, and negative or erratic orders remain visible.\n");

  for (size_t quantity = 0; quantity < QUANTITY_COUNT; ++quantity) {
    printf("\n%s\n", headings[quantity]);
    printf("comp        N            h     sampled L2   sampled Linf      p(L2)   p(Linf)\n");
    for (size_t component = 0; component < KERR_EXACT_COMPONENT_COUNT;
         ++component) {
      print_norm_row(
          kerr_exact_component_name((kerr_exact_component)component), levels,
          arguments->resolution_count, quantity, component, 0);
    }
    print_norm_row("ALL", levels, arguments->resolution_count, quantity, 0, 1);
  }
}

int main(int argc, char **argv) {
  options arguments;
  const int parsed = parse_options(argc, argv, &arguments);
  const double throat = kerr_exact_throat_radius();
  const double half_width = 5.0 * throat;
  level_summary *levels;
  FILE *z_profile_stream = NULL;

  if (parsed == 1) return EXIT_SUCCESS;
  if (parsed != 2) {
    print_usage(stderr, argv[0]);
    return EXIT_FAILURE;
  }
  if (arguments.point_count >
      SIZE_MAX / (size_t)KERR_EXACT_COMPONENT_COUNT) {
    fprintf(stderr, "point count overflows aggregate norm bookkeeping\n");
    free(arguments.resolutions);
    return EXIT_FAILURE;
  }
  if (multiply_overflows_size_t(arguments.resolution_count, sizeof(*levels))) {
    fprintf(stderr, "too many resolution levels\n");
    free(arguments.resolutions);
    return EXIT_FAILURE;
  }
  levels = (level_summary *)calloc(arguments.resolution_count, sizeof(*levels));
  if (levels == NULL) {
    fprintf(stderr, "could not allocate convergence summaries\n");
    free(arguments.resolutions);
    return EXIT_FAILURE;
  }

  if (arguments.z_profile_enabled) {
    z_profile_stream = fopen(arguments.z_profile_output, "w");
    if (z_profile_stream == NULL) {
      fprintf(stderr, "could not open z-profile output '%s': %s\n",
              arguments.z_profile_output, strerror(errno));
      free(levels);
      free(arguments.resolutions);
      return EXIT_FAILURE;
    }
    if (fprintf(z_profile_stream,
                "component,quantity,resolution,z,analytic,interpolated,error\n") <
        0) {
      fprintf(stderr, "could not write z-profile header to '%s'\n",
              arguments.z_profile_output);
      fclose(z_profile_stream);
      free(levels);
      free(arguments.resolutions);
      return EXIT_FAILURE;
    }
  }

  for (size_t level = 0; level < arguments.resolution_count; ++level) {
    fprintf(stderr, "evaluating N=%zu ...\n", arguments.resolutions[level]);
    if (!evaluate_level(&arguments, arguments.resolutions[level], half_width,
                        throat, z_profile_stream, &levels[level])) {
      if (z_profile_stream != NULL) fclose(z_profile_stream);
      free(levels);
      free(arguments.resolutions);
      return EXIT_FAILURE;
    }
  }

  if (z_profile_stream != NULL) {
    if (fclose(z_profile_stream) != 0) {
      fprintf(stderr, "could not finish z-profile output '%s': %s\n",
              arguments.z_profile_output, strerror(errno));
      free(levels);
      free(arguments.resolutions);
      return EXIT_FAILURE;
    }
    fprintf(stderr, "wrote z profile to %s\n", arguments.z_profile_output);
  }

  print_results(&arguments, levels, throat, half_width);
  free(levels);
  free(arguments.resolutions);
  return EXIT_SUCCESS;
}
