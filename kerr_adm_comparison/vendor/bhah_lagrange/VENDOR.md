# Vendored BHaH Lagrange interpolator

These files are an unchanged snapshot from the read-only reference tree at
Git commit `00b4ce5d803aacf2ea8f3876ff358316c15ab498`:

| File | SHA-256 |
|---|---|
| `interpolation_3d_general__uniform_src_grid.c` | `bc4fc55764c5f0d55f1905ba43ef74a59ff27b25eb9f810b16e90b6d6045d69e` |
| `interpolation_lagrange_uniform.h` | `961125bd1ea96ac6266d5f10434ca059955621da3d238b09f3d978d49b577112` |
| `intrinsics/simd_intrinsics.h` | `00ec0b57e5cc36b4038735e5f4a4e7c1002f2d9172efc2202ea4379c30ebc422` |

The local checked adapter is deliberately kept outside this directory.  The
vendored source is compiled as GNU99 because its SIMD reduction macro uses a
GNU expression block.
