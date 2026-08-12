/*
 * Compile and link this translation unit as C++ while the implementation is
 * compiled as C.  Successful linkage verifies the public header's extern "C"
 * boundary, and the value initialization below confirms that the public
 * aggregate types are straightforward to use from C++.
 */
#include "hermite3d.h"

#include <cstddef>

int main() {
  hermite3d_grid grid{};

  /* num_functions == 0 is rejected before the deliberately empty grid or null
   * arrays are examined further.  This gives the smoke test a defined runtime
   * assertion without allocating a real interpolation grid. */
  const hermite3d_status status = hermite3d_interpolate_value_gradient(
      &grid, 0.0, 0.0, 0.0, std::size_t{0}, nullptr, nullptr);
  return status == HERMITE3D_INVALID_ARGUMENT ? 0 : 1;
}
