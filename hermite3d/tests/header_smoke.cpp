#include "hermite3d.h"

#include <cstddef>

int main() {
  hermite3d_grid grid{};
  const hermite3d_status status = hermite3d_interpolate_value_gradient(
      &grid, 0.0, 0.0, 0.0, std::size_t{0}, nullptr, nullptr);
  return status == HERMITE3D_INVALID_ARGUMENT ? 0 : 1;
}
