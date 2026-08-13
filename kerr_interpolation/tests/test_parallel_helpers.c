#include "kerr_parallel_helpers.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static void require(int condition, const char *message) {
  if (!condition) {
    fprintf(stderr, "parallel-helper test failed: %s\n", message);
    exit(EXIT_FAILURE);
  }
}

int main(void) {
  kerr_scaled_sum_squares left = {0};
  kerr_scaled_sum_squares right = {0};
  kerr_scaled_sum_squares empty = {0};
  kerr_scaled_sum_squares large_left = {0};
  kerr_scaled_sum_squares large_right = {0};
  kerr_indexed_failure failures[4];
  const kerr_indexed_failure *first;
  long double norm;

  kerr_scaled_sum_add(&left, 3.0L);
  kerr_scaled_sum_add(&left, 4.0L);
  kerr_scaled_sum_add(&right, 12.0L);
  kerr_scaled_sum_merge(&left, &right);
  norm = left.scale * sqrtl(left.sum_squares);
  require(fabsl(norm - 13.0L) <= 16.0L * LDBL_EPSILON,
          "merged scaled sum must retain the Euclidean norm");

  kerr_scaled_sum_merge(&left, &empty);
  require(fabsl(left.scale * sqrtl(left.sum_squares) - 13.0L) <=
              16.0L * LDBL_EPSILON,
          "merging an empty sum must have no effect");
  kerr_scaled_sum_merge(&empty, &right);
  require(empty.scale == right.scale &&
              empty.sum_squares == right.sum_squares,
          "an empty destination must copy a nonempty sum");

  kerr_scaled_sum_add(&large_left, LDBL_MAX / 4.0L);
  kerr_scaled_sum_add(&large_right, LDBL_MAX / 4.0L);
  kerr_scaled_sum_merge(&large_left, &large_right);
  require(isfinite(large_left.scale) && large_left.sum_squares == 2.0L,
          "merging large finite inputs must not form their squares");

  for (size_t entry = 0; entry < 4; ++entry)
    kerr_failure_clear(&failures[entry]);
  require(kerr_first_failure(failures, 4) == NULL,
          "an empty failure set must have no first entry");

  kerr_failure_record(&failures[0], 19, 1, 2, 3, 4);
  kerr_failure_record(&failures[1], 15, 13, 14, 15, 16);
  kerr_failure_record(&failures[1], 7, 5, 6, 7, 8);
  kerr_failure_record(&failures[2], 12, 9, 10, 11, 12);
  first = kerr_first_failure(failures, 4);
  require(first == &failures[1] && first->index == 7 && first->kind == 5 &&
              first->detail0 == 6 && first->detail1 == 7 &&
              first->status == 8,
          "selection must retain the lowest failure and its diagnostics");

  puts("OpenMP reduction and failure-selection helpers passed.");
  return EXIT_SUCCESS;
}
