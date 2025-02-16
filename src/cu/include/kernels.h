/* File:   eulernv.cu
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 * USA.
 * */
#pragma once
#include "constants.h"
#include "grid.hpp"
#include "matrix3d.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cuda.h>
#include <type_traits>


template <typename T, T Volume, std::size_t N>
__global__ void kernel_calc_conserved(std::array<T *, N> primitives,
                                      std::array<T *, N> conserved,
                                      std::size_t len) {

  const std::size_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid >= len) {
    return;
  }
  T store = primitives[0][tid];
  conserved[0][tid] = store * Volume;
  conserved[1][tid] = store * primitives[1][tid] * Volume;
  conserved[2][tid] = store * primitives[2][tid] * Volume;
  conserved[3][tid] = store * primitives[3][tid] * Volume;
  conserved[4][tid] =
      Volume * (store *
                    (primitives[1][tid] * primitives[1][tid] +
                     primitives[2][tid] * primitives[2][tid] +
                     primitives[3][tid] * primitives[3][tid]) /
                    T(2.0f) +
                primitives[4][tid] / (EULERCFD::CONSTS::GAMMA - T(1.0f)));
}

template <typename T, T Volume, std::size_t N>
__global__ void kernel_calc_primitives(std::array<T *, N> conserved,
                                       std::array<T *, N> primitives,
                                       std::size_t len) {

  const std::size_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid >= len) {
    return;
  }
  primitives[0][tid] = conserved[0][tid] / Volume;
  primitives[1][tid] = conserved[1][tid] / (primitives[0][tid] * Volume);
  primitives[2][tid] = conserved[2][tid] / (primitives[0][tid] * Volume);
  primitives[3][tid] = conserved[3][tid] / (primitives[0][tid] * Volume);
  primitives[4][tid] = (conserved[4][tid] / Volume -
                        T(0.5) * primitives[0][tid] *
                            (primitives[1][tid] * primitives[1][tid] +
                             primitives[2][tid] * primitives[2][tid] +
                             primitives[3][tid] * primitives[3][tid])) *
                       (EULERCFD::CONSTS::GAMMA - T(1.0));
}

template <typename T> __inline__ __device__ T warp_reduce_min(T val) {
  for (int offset = EULERCFD::DEVICE_PARAMETERS::WARPSIZE / 2; offset > 0;
       offset /= 2) {
    val = std::min(val, __shfl_down_sync(0xffffffff, val, offset));
  }
  return val;
}

template <typename T, std::size_t N>
__global__ void reduce_min_kernel(std::array<T *, N> primitives, T *block_mins,
                                  std::size_t len) {
  __shared__ T shared_min[EULERCFD::DEVICE_PARAMETERS::WARPSIZE];
  std::size_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  std::size_t lane = threadIdx.x % EULERCFD::DEVICE_PARAMETERS::WARPSIZE;
  std::size_t warp_id = threadIdx.x / EULERCFD::DEVICE_PARAMETERS::WARPSIZE;

  if (tid >= len) {
    return;
  }

  T cand = EULERCFD::CONSTS::CFL * EULERCFD::CONSTS::DELTA /
           (std::sqrt(EULERCFD::CONSTS::GAMMA * (primitives[4][tid]) /
                      primitives[0][tid]) +
            std::sqrt(primitives[1][tid] * primitives[1][tid] +
                      primitives[2][tid] * primitives[2][tid] +
                      primitives[3][tid] * primitives[3][tid]));

  T min_val = warp_reduce_min(cand);

  if (lane == 0) {
    shared_min[warp_id] = min_val;
  }
  __syncthreads();

  if (warp_id == 0) {
    min_val = (lane < blockDim.x / EULERCFD::DEVICE_PARAMETERS::WARPSIZE)
                  ? shared_min[lane]
                  : std::numeric_limits<T>::max();
    min_val = warp_reduce_min(min_val);
  }

  if (threadIdx.x == 0) {
    block_mins[blockIdx.x] = min_val;
  }
}

template <typename T, std::size_t N>
T calc_timestep(std::array<T *, N> primitives, std::size_t len,
                std::array<std::size_t, 2> lp) {
  const auto threads_per_block = lp[1];
  const auto num_blocks = lp[0];
  assert(num_blocks > 1);

  T *d_block_mins;
  DEVICE_MATRIX_MALLOC(&d_block_mins, num_blocks * sizeof(T));

  reduce_min_kernel<<<num_blocks, threads_per_block>>>(primitives, d_block_mins,
                                                       len);

  std::vector<T> h_block_mins(num_blocks, 0);
  cudaMemcpy(h_block_mins.data(), d_block_mins, num_blocks * sizeof(T),
             cudaMemcpyDeviceToHost);

  auto global_min_ptr =
      std::min_element(h_block_mins.cbegin(), h_block_mins.cend());
  if (global_min_ptr == h_block_mins.end()) {
    printf("ERROR: global min is bad");
    abort();
  }
  T global_min = *global_min_ptr;
  DEVICE_MATRIX_FREE(d_block_mins);
  return global_min;
}

template <typename T, T DS, std::size_t N, std::size_t N2>
__global__ void kernel_calc_gradients(std::array<T *, N> src,
                                      std::array<T *, N2> gradients,
                                      std::size_t len) {

  const std::size_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid >= len) {
    return;
  }

  std::size_t i, j, k;
  _1d23dindex_(tid, i, j, k);

  auto id = [](std::size_t i, std::size_t j, std::size_t k) -> std::size_t {
    return id_f(i, j, k);
  };

  if (!(i >= 1 && i < EULERCFD::CONSTS::NX - 1 && j >= 1 &&
        j < EULERCFD::CONSTS::NY - 1 && k >= 1 &&
        k < EULERCFD::CONSTS::NZ - 1)) {
    return;
  }

  constexpr T iscale = 1.0f / (2.0f * DS);
  gradients[0][id(i, j, k)] =
      (src[0][id(i + 1, j, k)] - src[0][id(i - 1, j, k)]) * (iscale);
  gradients[1][id(i, j, k)] =
      (src[0][id(i, j + 1, k)] - src[0][id(i, j - 1, k)]) * (iscale);
  gradients[2][id(i, j, k)] =
      (src[0][id(i, j, k + 1)] - src[0][id(i, j, k - 1)]) * (iscale);

  gradients[3][id(i, j, k)] =
      (src[1][id(i + 1, j, k)] - src[1][id(i - 1, j, k)]) * (iscale);
  gradients[4][id(i, j, k)] =
      (src[1][id(i, j + 1, k)] - src[1][id(i, j - 1, k)]) * (iscale);
  gradients[5][id(i, j, k)] =
      (src[1][id(i, j, k + 1)] - src[1][id(i, j, k - 1)]) * (iscale);

  gradients[6][id(i, j, k)] =
      (src[2][id(i + 1, j, k)] - src[2][id(i - 1, j, k)]) * (iscale);
  gradients[7][id(i, j, k)] =
      (src[2][id(i, j + 1, k)] - src[2][id(i, j - 1, k)]) * (iscale);
  gradients[8][id(i, j, k)] =
      (src[2][id(i, j, k + 1)] - src[2][id(i, j, k - 1)]) * (iscale);

  gradients[9][id(i, j, k)] =
      (src[3][id(i + 1, j, k)] - src[3][id(i - 1, j, k)]) * (iscale);
  gradients[10][id(i, j, k)] =
      (src[3][id(i, j + 1, k)] - src[3][id(i, j - 1, k)]) * (iscale);
  gradients[11][id(i, j, k)] =
      (src[3][id(i, j, k + 1)] - src[3][id(i, j, k - 1)]) * (iscale);

  gradients[12][id(i, j, k)] =
      (src[4][id(i + 1, j, k)] - src[4][id(i - 1, j, k)]) * (iscale);
  gradients[13][id(i, j, k)] =
      (src[4][id(i, j + 1, k)] - src[4][id(i, j - 1, k)]) * (iscale);
  gradients[14][id(i, j, k)] =
      (src[4][id(i, j, k + 1)] - src[4][id(i, j, k - 1)]) * (iscale);
}

template <typename T, int D, EULERCFD::BC B, int NG, std::size_t N>
__global__ void kernel_update_ghosts(std::array<T *, N> src, std::size_t len) {
  const std::size_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid >= len) {
    return;
  }
  std::size_t i, j, k;
  _1d23dindex_(tid, i, j, k);

  auto id = [](std::size_t i, std::size_t j, std::size_t k) -> std::size_t {
    return id_f(i, j, k);
  };

  if constexpr (B == EULERCFD::BC::PERIODIC) {
    std::size_t jump_index;
    if constexpr (D == 0) {
      if ((i >= NG)) {
        return;
      }
      jump_index = EULERCFD::CONSTS::NX - NG - 1 - (i == 0);
      for (std::size_t q = 0; q < N; ++q) {
        src[q][id(i, j, k)] = src[q][id(jump_index, j, k)];
      }
    } else if constexpr (D == 1) {
      if (((j >= NG))) {
        return;
      }
      jump_index = EULERCFD::CONSTS::NY - NG - 1 - (j == 0);
      for (std::size_t q = 0; q < N; ++q) {
        src[q][id(i, j, k)] = src[q][id(i, jump_index, k)];
      }
    } else if constexpr (D == 2) {
      if ((((k >= NG)))) {
        return;
      }
      jump_index = EULERCFD::CONSTS::NZ - NG - 1 - (k == 0);
      for (std::size_t q = 0; q < N; ++q) {
        src[q][id(i, j, k)] = src[q][id(i, j, jump_index)];
      }
    } else if constexpr (D == 3) {
      if (((i < EULERCFD::CONSTS::NX - NG))) {
        return;
      }
      jump_index = 2 + !(i == EULERCFD::CONSTS::NX - NG);
      for (std::size_t q = 0; q < N; ++q) {
        src[q][id(i, j, k)] = src[q][id(jump_index, j, k)];
      }
    } else if constexpr (D == 4) {
      if (((j < EULERCFD::CONSTS::NY - NG))) {
        return;
      }
      jump_index = 2 + !(j == EULERCFD::CONSTS::NY - NG);
      for (std::size_t q = 0; q < N; ++q) {
        src[q][id(i, j, k)] = src[q][id(i, jump_index, k)];
      }
    } else if constexpr (D == 5) {
      if ((((k < EULERCFD::CONSTS::NZ - NG)))) {
        return;
      }
      jump_index = 2 + !(k == EULERCFD::CONSTS::NZ - NG);

      for (std::size_t q = 0; q < N; ++q) {
        src[q][id(i, j, k)] = src[q][id(i, j, jump_index)];
      }
    }
  }

  if constexpr (B == EULERCFD::BC::OUTFLOW) {
    std::size_t jump_index;
    if constexpr (D == 0) {
      if ((i >= NG)) {
        return;
      }
      jump_index = 2;
      for (std::size_t q = 0; q < N; ++q) {
        src[q][id(i, j, k)] = src[q][id(jump_index, j, k)];
      }
    } else if constexpr (D == 1) {
      if (((j >= NG))) {
        return;
      }
      jump_index = 2; // Lower boundary (j == 0)
      for (std::size_t q = 0; q < N; ++q) {
        src[q][id(i, j, k)] = src[q][id(i, jump_index, k)];
      }
    } else if constexpr (D == 2) {
      if ((((k >= NG)))) {
        return;
      }
      jump_index = 2; // Lower boundary (k == 0)
      for (std::size_t q = 0; q < N; ++q) {
        src[q][id(i, j, k)] = src[q][id(i, j, jump_index)];
      }
    } else if constexpr (D == 3) {
      if (((i < EULERCFD::CONSTS::NX - NG))) {
        return;
      }
      jump_index = EULERCFD::CONSTS::NX - NG -
                   1; // Upper boundary (i == EULERCFD::CONSTS::NX - 1)
      for (std::size_t q = 0; q < N; ++q) {
        src[q][id(i, j, k)] = src[q][id(jump_index, j, k)];
      }
    } else if constexpr (D == 4) {
      if (((j < EULERCFD::CONSTS::NY - NG))) {
        return;
      }
      jump_index = EULERCFD::CONSTS::NY - NG - 1;
      for (std::size_t q = 0; q < N; ++q) {
        src[q][id(i, j, k)] = src[q][id(i, jump_index, k)];
      }
    } else if constexpr (D == 5) {
      if ((((k < EULERCFD::CONSTS::NZ - NG)))) {
        return;
      }
      jump_index = EULERCFD::CONSTS::NZ - NG - 1;
      for (std::size_t q = 0; q < N; ++q) {
        src[q][id(i, j, k)] = src[q][id(i, j, jump_index)];
      }
    }
  }

  if constexpr (B == EULERCFD::BC::WALL || B==EULERCFD::BC::CONDUCTING_WALL) {
    std::size_t jump_index;
    if constexpr (D == 0) {
      if ((i >= NG)) {
        return;
      }
      jump_index = 2;
      for (std::size_t q = 0; q < N; ++q) {
        T scale = (q == 1) ? -1.0f : 1.0f;
        src[q][id(i, j, k)] = scale * src[q][id(jump_index, j, k)];
      }
    } else if constexpr (D == 1) {
      if (((j >= NG))) {
        return;
      }
      jump_index = 2; // Lower boundary (j == 0)
      for (std::size_t q = 0; q < N; ++q) {
        T scale = (q == 2) ? -1.0f : 1.0f;
        src[q][id(i, j, k)] = scale * src[q][id(i, jump_index, k)];
      }
    } else if constexpr (D == 2) {
      if ((((k >= NG)))) {
        return;
      }
      jump_index = 2; // Lower boundary (k == 0)
      for (std::size_t q = 0; q < N; ++q) {
        T scale = (q == 3) ? -1.0f : 1.0f;
        src[q][id(i, j, k)] = scale * src[q][id(i, j, jump_index)];
      }
    } else if constexpr (D == 3) {
      if (((i < EULERCFD::CONSTS::NX - NG))) {
        return;
      }
      jump_index = EULERCFD::CONSTS::NX - NG -
                   1; // Upper boundary (i == EULERCFD::CONSTS::NX - 1)
      for (std::size_t q = 0; q < N; ++q) {
        T scale = (q == 1) ? -1.0f : 1.0f;
        src[q][id(i, j, k)] = scale * src[q][id(jump_index, j, k)];
      }
    } else if constexpr (D == 4) {
      if (((j < EULERCFD::CONSTS::NY - NG))) {
        return;
      }
      jump_index = EULERCFD::CONSTS::NY - NG - 1;
      for (std::size_t q = 0; q < N; ++q) {
        T scale = (q == 2) ? -1.0f : 1.0f;
        src[q][id(i, j, k)] = scale * src[q][id(i, jump_index, k)];
      }
    } else if constexpr (D == 5) {
      if ((((k < EULERCFD::CONSTS::NZ - NG)))) {
        return;
      }
      jump_index = EULERCFD::CONSTS::NZ - NG - 1;
      for (std::size_t q = 0; q < N; ++q) {
        T scale = (q == 3) ? -1.0f : 1.0f;
        src[q][id(i, j, k)] = scale * src[q][id(i, j, jump_index)];
      }
    }
  }
}

template <typename T, T DS, std::size_t N, std::size_t N2>
__global__ void kernel_calc_xtr(std::array<T *, N> primitives,
                                std::array<T *, N> primitives_xtr,
                                std::array<T *, N2> gradients, std::size_t len,
                                T dt) {

  const std::size_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid >= len) {
    return;
  }
  std::size_t i, j, k;
  _1d23dindex_(tid, i, j, k);

  auto id = [](std::size_t i, std::size_t j, std::size_t k) -> std::size_t {
    return id_f(i, j, k);
  };

  if (!(i >= 1 && i < EULERCFD::CONSTS::NX - 1 && j >= 1 &&
        j < EULERCFD::CONSTS::NY - 1 && k >= 1 &&
        k < EULERCFD::CONSTS::NZ - 1)) {
    return;
  }

  T rho = primitives[0][id(i, j, k)];
  T vx = primitives[1][id(i, j, k)];
  T vy = primitives[2][id(i, j, k)];
  T vz = primitives[3][id(i, j, k)];
  T p = primitives[4][id(i, j, k)];

  T drho_dx = gradients[0][id(i, j, k)];
  T drho_dy = gradients[1][id(i, j, k)];
  T drho_dz = gradients[2][id(i, j, k)];

  T dvx_dx = gradients[3][id(i, j, k)];
  T dvx_dy = gradients[4][id(i, j, k)];
  T dvx_dz = gradients[5][id(i, j, k)];

  T dvy_dx = gradients[6][id(i, j, k)];
  T dvy_dy = gradients[7][id(i, j, k)];
  T dvy_dz = gradients[8][id(i, j, k)];

  T dvz_dx = gradients[9][id(i, j, k)];
  T dvz_dy = gradients[10][id(i, j, k)];
  T dvz_dz = gradients[11][id(i, j, k)];

  T dp_dx = gradients[12][id(i, j, k)];
  T dp_dy = gradients[13][id(i, j, k)];
  T dp_dz = gradients[14][id(i, j, k)];

  T *rho_xtr = primitives_xtr[0];
  T *vx_xtr = primitives_xtr[1];
  T *vy_xtr = primitives_xtr[2];
  T *vz_xtr = primitives_xtr[3];
  T *p_xtr = primitives_xtr[4];

  rho_xtr[id(i, j, k)] = rho - 0.5f * dt *
                                   (vx * drho_dx + rho * dvx_dx + vy * drho_dy +
                                    rho * dvy_dy + vz * drho_dz + rho * dvz_dz);

  vx_xtr[id(i, j, k)] =
      vx - 0.5f * dt *
               (vx * dvx_dx + vy * dvx_dy + vz * dvx_dz + (1.0f / rho) * dp_dx);

  vy_xtr[id(i, j, k)] =
      vy - 0.5f * dt *
               (vx * dvy_dx + vy * dvy_dy + vz * dvy_dz + (1.0f / rho) * dp_dy);

  vz_xtr[id(i, j, k)] =
      vz - 0.5f * dt *
               (vx * dvz_dx + vy * dvz_dy + vz * dvz_dz + (1.0f / rho) * dp_dz);

  p_xtr[id(i, j, k)] =
      p - 0.5f * dt *
              (EULERCFD::CONSTS::GAMMA * p * (dvx_dx + dvy_dy + dvz_dz) +
               vx * dp_dx + vy * dp_dy + vz * dp_dz);
}

template <typename T, T DS>
__global__ void
kernel_calc_fluxes(T *mass_flux_x, T *momentum_x_flux_x, T *momentum_y_flux_x,
                   T *momentum_z_flux_x, T *energy_flux_x, T *drho, T *dvx,
                   T *dvy, T *dvz, T *dp, T *rho, T *vx, T *vy, T *vz, T *p,
                   std::size_t len, std::array<std::size_t, 3> offsets,
                   std::size_t dim) {

  const std::size_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid >= len) {
    return;
  }

  std::size_t i, j, k;
  _1d23dindex_(tid, i, j, k);

  auto id = [](std::size_t i, std::size_t j, std::size_t k) -> std::size_t {
    return id_f(i, j, k);
  };

  if (!(i >= 1 && i < EULERCFD::CONSTS::NX - 1 && j >= 1 &&
        j < EULERCFD::CONSTS::NY - 1 && k >= 1 &&
        k < EULERCFD::CONSTS::NZ - 1)) {
    return;
  }

  auto van_leer = [](T a, T b) -> T {
    if (a * b <= 0.0) {
      return 0.0;
    } else {
      return (2.0 * a * b) /
             (a + b + 1e-8); // Add small epsilon to prevent division by zero
    }
  };

  T rl = rho[id(i + offsets[0], j + offsets[1], k + offsets[2])] -
         drho[id(i + offsets[0], j + offsets[1], k + offsets[2])] * (DS / 2.0f);

  T rr = rho[id(i, j, k)] + drho[id(i, j, k)] * (DS / 2.0f);

  T vxl = vx[id(i + offsets[0], j + offsets[1], k + offsets[2])] -
          dvx[id(i + offsets[0], j + offsets[1], k + offsets[2])] * (DS / 2.0f);

  T vxr = vx[id(i, j, k)] + dvx[id(i, j, k)] * (DS / 2.0f);

  T vyl = vy[id(i + offsets[0], j + offsets[1], k + offsets[2])] -
          dvy[id(i + offsets[0], j + offsets[1], k + offsets[2])] * (DS / 2.0f);

  T vyr = vy[id(i, j, k)] + dvy[id(i, j, k)] * (DS / 2.0f);

  T vzl = vz[id(i + offsets[0], j + offsets[1], k + offsets[2])] -
          dvz[id(i + offsets[0], j + offsets[1], k + offsets[2])] * (DS / 2.0f);

  T vzr = vz[id(i, j, k)] + dvz[id(i, j, k)] * (DS / 2.0f);

  T pl = p[id(i + offsets[0], j + offsets[1], k + offsets[2])] -
         dp[id(i + offsets[0], j + offsets[1], k + offsets[2])] * (DS / 2.0f);

  T pr = p[id(i, j, k)] + dp[id(i, j, k)] * (DS / 2.0f);

  T en_left = 0.5f * (rl * (vxl * vxl + vyl * vyl + vzl * vzl)) +
              pl / (EULERCFD::CONSTS::GAMMA - 1.0f);

  T en_right = 0.5f * (rr * (vxr * vxr + vyr * vyr + vzr * vzr)) +
               pr / (EULERCFD::CONSTS::GAMMA - 1.0f);

  T rho_star = (rl + rr) / 2.0f;

  T momentum_x_star = (rl * vxl + rr * vxr) / 2.0f;
  T momentum_y_star = (rl * vyl + rr * vyr) / 2.0f;
  T momentum_z_star = (rl * vzl + rr * vzr) / 2.0f;

  T en_star = 0.5f * (en_left + en_right);

  T p_star = (EULERCFD::CONSTS::GAMMA - 1.0f) *
             (en_star - 0.5f *
                            (momentum_x_star * momentum_x_star +
                             momentum_y_star * momentum_y_star +
                             momentum_z_star * momentum_z_star) /
                            rho_star);

  mass_flux_x[id(i, j, k)] = momentum_x_star;

  momentum_x_flux_x[id(i, j, k)] =
      (momentum_x_star * momentum_x_star / rho_star) + p_star;

  momentum_y_flux_x[id(i, j, k)] =
      (momentum_x_star * momentum_y_star) / rho_star;

  momentum_z_flux_x[id(i, j, k)] =
      (momentum_x_star * momentum_z_star) / rho_star;

  energy_flux_x[id(i, j, k)] =
      (en_star + p_star) * (momentum_x_star / rho_star);

  // Sources
  if constexpr (EULERCFD::CONSTS::bcs[0] == EULERCFD::BC::CONDUCTING_WALL) {
    if (i == 1) {
      const T t_in =
          p[id(i, j, k)] / (EULERCFD::CONSTS::RS * rho[id(i, j, k)]);
      const T hc_out = 10.0f - EULERCFD::CONSTS::VOUT_MS +
                       10.0f * std::sqrt(EULERCFD::CONSTS::VOUT_MS);
      constexpr T S = EULERCFD::CONSTS::DELTA * EULERCFD::CONSTS::DELTA;
      constexpr T V = EULERCFD::CONSTS::DELTA * EULERCFD::CONSTS::DELTA *
                      EULERCFD::CONSTS::DELTA;
      constexpr T d = 0.008f;
      constexpr T kappa = 0.2f;
      const T Q = (EULERCFD::CONSTS::T_OUT_KELVIN - t_in) *
                  (S / ((d / kappa) + (1.0f / hc_out)));
      energy_flux_x[id(i, j, k)] += Q / V;
    }
  }
  
  if  constexpr(EULERCFD::CONSTS::bcs[1] == EULERCFD::BC::CONDUCTING_WALL) {
    if (j == 1) {
      const T t_in =
          p[id(i, j, k)] / (EULERCFD::CONSTS::RS * rho[id(i, j, k)]);
      const T hc_out = 10.0f - EULERCFD::CONSTS::VOUT_MS +
                       10.0f * std::sqrt(EULERCFD::CONSTS::VOUT_MS);
      constexpr T S = EULERCFD::CONSTS::DELTA * EULERCFD::CONSTS::DELTA;
      constexpr T V = EULERCFD::CONSTS::DELTA * EULERCFD::CONSTS::DELTA *
                      EULERCFD::CONSTS::DELTA;
      constexpr T d = 0.008f;
      constexpr T kappa = 0.2f;
      const T Q = (EULERCFD::CONSTS::T_OUT_KELVIN - t_in) *
                  (S / ((d / kappa) + (1.0f / hc_out)));
      energy_flux_x[id(i, j, k)] += Q / V;
    }
  }
  
  if constexpr (EULERCFD::CONSTS::bcs[2] == EULERCFD::BC::CONDUCTING_WALL) {
    if (k == 1) {
      const T t_in =
          p[id(i, j, k)] / (EULERCFD::CONSTS::RS * rho[id(i, j, k)]);
      const T hc_out = 10.0f - EULERCFD::CONSTS::VOUT_MS +
                       10.0f * std::sqrt(EULERCFD::CONSTS::VOUT_MS);
      constexpr T S = EULERCFD::CONSTS::DELTA * EULERCFD::CONSTS::DELTA;
      constexpr T V = EULERCFD::CONSTS::DELTA * EULERCFD::CONSTS::DELTA *
                      EULERCFD::CONSTS::DELTA;
      constexpr T d = 0.008f;
      constexpr T kappa = 0.2f;
      const T Q = (EULERCFD::CONSTS::T_OUT_KELVIN - t_in) *
                  (S / ((d / kappa) + (1.0f / hc_out)));
      energy_flux_x[id(i, j, k)] += Q / V;
    }
  }
 
  if constexpr (EULERCFD::CONSTS::bcs[3] == EULERCFD::BC::CONDUCTING_WALL) {
    if (i == EULERCFD::CONSTS::NX-2) {
      const T t_in =
          p[id(i, j, k)] / (EULERCFD::CONSTS::RS * rho[id(i, j, k)]);
      const T hc_out = 10.0f - EULERCFD::CONSTS::VOUT_MS +
                       10.0f * std::sqrt(EULERCFD::CONSTS::VOUT_MS);
      constexpr T S = EULERCFD::CONSTS::DELTA * EULERCFD::CONSTS::DELTA;
      constexpr T V = EULERCFD::CONSTS::DELTA * EULERCFD::CONSTS::DELTA *
                      EULERCFD::CONSTS::DELTA;
      constexpr T d = 0.008f;
      constexpr T kappa = 0.2f;
      const T Q = (EULERCFD::CONSTS::T_OUT_KELVIN - t_in) *
                  (S / ((d / kappa) + (1.0f / hc_out)));
      energy_flux_x[id(i, j, k)] += Q / V;
    }
  }
  
  if  constexpr(EULERCFD::CONSTS::bcs[4] == EULERCFD::BC::CONDUCTING_WALL) {
    if (j == EULERCFD::CONSTS::NY-2) {
      const T t_in =
          p[id(i, j, k)] / (EULERCFD::CONSTS::RS * rho[id(i, j, k)]);
      const T hc_out = 10.0f - EULERCFD::CONSTS::VOUT_MS +
                       10.0f * std::sqrt(EULERCFD::CONSTS::VOUT_MS);
      constexpr T S = EULERCFD::CONSTS::DELTA * EULERCFD::CONSTS::DELTA;
      constexpr T V = EULERCFD::CONSTS::DELTA * EULERCFD::CONSTS::DELTA *
                      EULERCFD::CONSTS::DELTA;
      constexpr T d = 0.008f;
      constexpr T kappa = 0.2f;
      const T Q = (EULERCFD::CONSTS::T_OUT_KELVIN - t_in) *
                  (S / ((d / kappa) + (1.0f / hc_out)));
      energy_flux_x[id(i, j, k)] += Q / V;
    }
  }
  
  if constexpr (EULERCFD::CONSTS::bcs[5] == EULERCFD::BC::CONDUCTING_WALL) {
    if (k == EULERCFD::CONSTS::NZ-2) {
      const T t_in =
          p[id(i, j, k)] / (EULERCFD::CONSTS::RS * rho[id(i, j, k)]);
      const T hc_out = 10.0f - EULERCFD::CONSTS::VOUT_MS +
                       10.0f * std::sqrt(EULERCFD::CONSTS::VOUT_MS);
      constexpr T S = EULERCFD::CONSTS::DELTA * EULERCFD::CONSTS::DELTA;
      constexpr T V = EULERCFD::CONSTS::DELTA * EULERCFD::CONSTS::DELTA *
                      EULERCFD::CONSTS::DELTA;
      constexpr T d = 0.008f;
      constexpr T kappa = 0.2f;
      const T Q = (EULERCFD::CONSTS::T_OUT_KELVIN - t_in) *
                  (S / ((d / kappa) + (1.0f / hc_out)));
      energy_flux_x[id(i, j, k)] += Q / V;
    }
  }

  // Gravity term if enabled
  if (offsets[2] == 1 && EULERCFD::CONSTS::GRAVITY) {
    T h = (EULERCFD::CONSTS::NZ - 2.0f * EULERCFD::CONSTS::NGHOSTS) *
              EULERCFD::CONSTS::DELTA -
          (k - 2.0f) * EULERCFD::CONSTS::DELTA;
    momentum_x_flux_x[id(i, j, k)] += 1.0f * rho_star * EULERCFD::CONSTS::G * h;
    energy_flux_x[id(i, j, k)] +=
        0.5f * rho_star * (vxl + vxr) * EULERCFD::CONSTS::G * h;
  }

  T c_l = std::sqrt(EULERCFD::CONSTS::GAMMA * pl / rl) + std::abs(vxl);
  T c_r = std::sqrt(EULERCFD::CONSTS::GAMMA * pr / rr) + std::abs(vxr);
  T c_star = std::max(c_l, c_r);

  mass_flux_x[id(i, j, k)] = momentum_x_star - 0.5 * c_star * (rl - rr);
  momentum_x_flux_x[id(i, j, k)] -= c_star * (rl * vxl - rr * vxr) / 2.0f;
  momentum_y_flux_x[id(i, j, k)] -= c_star * (rl * vyl - rr * vyr) / 2.0f;
  momentum_z_flux_x[id(i, j, k)] -= c_star * (rl * vzl - rr * vzr) / 2.0f;
  energy_flux_x[id(i, j, k)] -= c_star * (en_left - en_right) / 2.0f;

  // Apply scaling factor `ds`
  mass_flux_x[id(i, j, k)] *= DS;
  momentum_x_flux_x[id(i, j, k)] *= DS;
  momentum_y_flux_x[id(i, j, k)] *= DS;
  momentum_z_flux_x[id(i, j, k)] *= DS;
  energy_flux_x[id(i, j, k)] *= DS;
}

template <typename T, T DS, std::size_t N, std::size_t N2>
__global__ void kernel_calc_addfluxes(std::array<T *, N2> fluxes,
                                      std::array<T *, N> conserved,
                                      std::size_t len, T dt) {

  const std::size_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid >= len) {
    return;
  }

  std::size_t i, j, k;
  _1d23dindex_(tid, i, j, k);

  auto id = [](std::size_t i, std::size_t j, std::size_t k) -> std::size_t {
    return id_f(i, j, k);
  };

  if (!(i >= 1 && i < EULERCFD::CONSTS::NX - 1 && j >= 1 &&
        j < EULERCFD::CONSTS::NY - 1 && k >= 1 &&
        k < EULERCFD::CONSTS::NZ - 1)) {
    return;
  }

  T *mass = conserved[0];
  T *momentum_x = conserved[1];
  T *momentum_y = conserved[2];
  T *momentum_z = conserved[3];
  T *energy = conserved[4];

  T *mass_flux_x = fluxes[0];
  T *mass_flux_y = fluxes[1];
  T *mass_flux_z = fluxes[2];
  T *momentum_x_flux_x = fluxes[3];
  T *momentum_x_flux_y = fluxes[4];
  T *momentum_x_flux_z = fluxes[5];
  T *momentum_y_flux_x = fluxes[6];
  T *momentum_y_flux_y = fluxes[7];
  T *momentum_y_flux_z = fluxes[8];
  T *momentum_z_flux_x = fluxes[9];
  T *momentum_z_flux_y = fluxes[10];
  T *momentum_z_flux_z = fluxes[11];
  T *energy_flux_x = fluxes[12];
  T *energy_flux_y = fluxes[13];
  T *energy_flux_z = fluxes[14];

  mass[id(i, j, k)] -=
      dt * DS *
      (mass_flux_x[id(i, j, k)] - mass_flux_x[id(i - 1, j, k)] +
       mass_flux_y[id(i, j, k)] - mass_flux_y[id(i, j - 1, k)] +
       mass_flux_z[id(i, j, k)] - mass_flux_z[id(i, j, k - 1)]);

  momentum_x[id(i, j, k)] -=
      dt * DS *
      (momentum_x_flux_x[id(i, j, k)] - momentum_x_flux_x[id(i - 1, j, k)] +
       momentum_x_flux_y[id(i, j, k)] - momentum_x_flux_y[id(i, j - 1, k)] +
       momentum_x_flux_z[id(i, j, k)] - momentum_x_flux_z[id(i, j, k - 1)]);

  momentum_y[id(i, j, k)] -=
      dt * DS *
      (momentum_y_flux_x[id(i, j, k)] - momentum_y_flux_x[id(i - 1, j, k)] +
       momentum_y_flux_y[id(i, j, k)] - momentum_y_flux_y[id(i, j - 1, k)] +
       momentum_y_flux_z[id(i, j, k)] - momentum_y_flux_z[id(i, j, k - 1)]);

  momentum_z[id(i, j, k)] -=
      dt * DS *
      (momentum_z_flux_x[id(i, j, k)] - momentum_z_flux_x[id(i - 1, j, k)] +
       momentum_z_flux_y[id(i, j, k)] - momentum_z_flux_y[id(i, j - 1, k)] +
       momentum_z_flux_z[id(i, j, k)] - momentum_z_flux_z[id(i, j, k - 1)]);

  energy[id(i, j, k)] -=
      dt * DS *
      (energy_flux_x[id(i, j, k)] - energy_flux_x[id(i - 1, j, k)] +
       energy_flux_y[id(i, j, k)] - energy_flux_y[id(i, j - 1, k)] +
       energy_flux_z[id(i, j, k)] - energy_flux_z[id(i, j, k - 1)]);
}
