/* File:   eulernv.cu
* Authors: Kostis Papadakis and Adam Kit (2024)
* Description: A program that acts as a heat source to heat up our green house

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
USA.
* */

#include "include/constants.h"
#include "include/grid.hpp"
#include "include/matrix3d.hpp"
#include "include/physics.hpp"
#include <array>

using namespace EULERCFD;
using namespace EULERCFD::CONSTS;
using namespace EULERCFD::DEVICE_PARAMETERS;

void set_log_level() {
  spdlog::set_level(spdlog::level::info);
  if (const char *env_p = std::getenv("DEBUG")) {
    if (strncmp(env_p, "1", 1) == 0) {
      spdlog::set_level(spdlog::level::debug);
      return;
    }
  }
  if (const char *env_p = std::getenv("INFO")) {
    if (strncmp(env_p, "0", 0) == 0) {
      spdlog::set_level(spdlog::level::off);
    }
  }
}

template <typename T> consteval auto get_init_function() {
  if constexpr (RUN_SETUP == EULERCFD::SETUP::KHI) {
    return &init_khi<T>;
  } else if constexpr (RUN_SETUP == EULERCFD::SETUP::TRB) {
    return &init_trb<T>;
  } else if constexpr (RUN_SETUP == EULERCFD::SETUP::GREENHOUSE) {
    return &init_greenhouse<T>;
  } else if constexpr (RUN_SETUP == EULERCFD::SETUP::SOD) {
    return &init_sod<T>;
  } else {
    static_assert(sizeof(T) < 0,
                  "Setup not supported"); // hack I know :) but compilers should
                                          // support my logic here
  }
}

int main() {
  using type_t = float;
  spdlog::stopwatch sw;
  set_log_level();
  spdlog::info("Starting Simulation!");
  EULERCFD::compute(
      EULERCFD::Grid<type_t, EULERCFD::GridInfo<type_t>{NX, NY, NZ, NGHOSTS, LX, LY, LZ},
           BACKEND::DEVICE>{},
      type_t(EULERCFD::CONSTS::TMAX), EULERCFD::CONSTS::MAXSTEPS,
      type_t(EULERCFD::CONSTS::TOUT), get_init_function<type_t>());
  spdlog::info("Simulation done in {0:f} seconds!", sw);
  return 0;
}
