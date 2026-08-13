#ifndef KERR_PARALLEL_HELPERS_H
#define KERR_PARALLEL_HELPERS_H

#include <math.h>
#include <stddef.h>
#include <stdint.h>

/*
 * A scaled sum of squares represents
 *
 *   scale^2 * sum_squares.
 *
 * Keeping the common scale outside the sum avoids overflow when individual
 * errors are finite but their ordinary squares are not representable.
 */
typedef struct kerr_scaled_sum_squares {
  long double scale;
  long double sum_squares;
} kerr_scaled_sum_squares;

static void kerr_scaled_sum_add(kerr_scaled_sum_squares *sum,
                                long double value) {
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

/* Merge scaled sums without first converting either to an unscaled sum. */
static void kerr_scaled_sum_merge(kerr_scaled_sum_squares *destination,
                                  const kerr_scaled_sum_squares *source) {
  if (source->scale == 0.0L) return;
  if (destination->scale == 0.0L) {
    *destination = *source;
  } else if (destination->scale < source->scale) {
    const long double ratio = destination->scale / source->scale;
    destination->sum_squares =
        source->sum_squares + destination->sum_squares * ratio * ratio;
    destination->scale = source->scale;
  } else {
    const long double ratio = source->scale / destination->scale;
    destination->sum_squares += source->sum_squares * ratio * ratio;
  }
}

/*
 * Parallel workers retain diagnostics privately.  The serial caller selects
 * the lowest failed logical index after the OpenMP region, making diagnostics
 * deterministic and avoiding output or locks in the successful hot path.
 */
typedef struct kerr_indexed_failure {
  size_t index;
  int kind;
  size_t detail0;
  size_t detail1;
  int status;
} kerr_indexed_failure;

static void kerr_failure_clear(kerr_indexed_failure *failure) {
  failure->index = SIZE_MAX;
  failure->kind = 0;
  failure->detail0 = 0;
  failure->detail1 = 0;
  failure->status = 0;
}

static void kerr_failure_record(kerr_indexed_failure *failure, size_t index,
                                int kind, size_t detail0, size_t detail1,
                                int status) {
  if (index < failure->index) {
    failure->index = index;
    failure->kind = kind;
    failure->detail0 = detail0;
    failure->detail1 = detail1;
    failure->status = status;
  }
}

static const kerr_indexed_failure *kerr_first_failure(
    const kerr_indexed_failure *failures, size_t count) {
  const kerr_indexed_failure *first = NULL;

  for (size_t entry = 0; entry < count; ++entry) {
    if (failures[entry].index != SIZE_MAX &&
        (first == NULL || failures[entry].index < first->index)) {
      first = &failures[entry];
    }
  }
  return first;
}

#endif
