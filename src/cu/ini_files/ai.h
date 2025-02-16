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
#include <array>

// device callable stuff
#define host_only __host__
#define dev_host __host__ __device__
// allocators
#define DEVICE_MATRIX_MALLOC cudaMalloc
#define DEVICE_MATRIX_FREE cudaFree
#define HOST_MATRIX_MALLOC malloc
#define HOST_MATRIX_FREE free
// profilers
#define PROFILE_RUN
// #define USE_MORTON

#ifdef PROFILE_RUN
#include <nvToolsExt.h>
#define PROFILE_START(msg) nvtxRangePushA((msg))
#define PROFILE_END() nvtxRangePop()
#else
#define PROFILE_START(msg)
#define PROFILE_END()
#endif

namespace EULERCFD {
enum SETUP { KHI, TRB,SOD,GREENHOUSE };
enum class BC { PERIODIC, WALL, CONDUCTING_WALL,GREENHOUSE, OUTFLOW ,KHI,SOD };

namespace CONSTS {

inline constexpr float GAMMA = 5.0 / 3.0;
inline constexpr float CFL = 0.2;
inline constexpr float DELTA = 1.0 / 128;
inline constexpr std::size_t NX = 128;
inline constexpr std::size_t NY = 6;
inline constexpr std::size_t NZ = 128;
inline constexpr float LX = NX * DELTA;
inline constexpr float LY = NY * DELTA;
inline constexpr float LZ = NZ * DELTA;
inline constexpr std::size_t NGHOSTS = 2;
inline constexpr float TMAX = 2;
inline constexpr float TOUT = 0.0000000001;
inline constexpr std::size_t MAXSTEPS = 100000000000;

inline constexpr float G = 9.81;
inline constexpr float RS = 287.0;
inline constexpr bool GRAVITY = false;
inline constexpr float T_OUT_KELVIN = 273.0f;
inline constexpr float VOUT_MS = 0.0f;
constexpr SETUP RUN_SETUP=SETUP::KHI;
inline constexpr std::array<BC, 6> bcs = {BC::PERIODIC, BC::PERIODIC,
                                          BC::PERIODIC, BC::PERIODIC,
                                          BC::PERIODIC, BC::PERIODIC};

} // namespace CONSTS

namespace DEVICE_PARAMETERS {
inline constexpr std::size_t MAX_BLOCKSIZE = 256;
inline constexpr std::size_t WARPSIZE = 32;
} // namespace DEVICE_PARAMETERS

} // namespace EULERCFD

__host__ __device__ static inline std::size_t id_f(std::size_t i, std::size_t j,
                                                   std::size_t k) {
  return i * (EULERCFD::CONSTS::NY * EULERCFD::CONSTS::NZ) +
         j * EULERCFD::CONSTS::NZ + k;
}

inline __device__ void _1d23dindex_(std::size_t tid, std::size_t &i,
                                    std::size_t &j, std::size_t &k) noexcept {
  i = tid / (EULERCFD::CONSTS::NY * EULERCFD::CONSTS::NZ);
  j = (tid % (EULERCFD::CONSTS::NY * EULERCFD::CONSTS::NZ)) /
      EULERCFD::CONSTS::NZ;
  k = tid % EULERCFD::CONSTS::NZ;
}
