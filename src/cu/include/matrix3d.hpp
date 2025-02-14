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
#include "spdlog/spdlog.h"
#include <array>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <cuda_runtime.h>
#include <driver_types.h>
#include <iostream>
#include <stdexcept>
#include <type_traits>

namespace EULERCFD {

inline consteval std::array<std::size_t, 2> launch_params(
    std::size_t arraySize,
    std::size_t blockSize = EULERCFD::DEVICE_PARAMETERS::MAX_BLOCKSIZE) {
  std::size_t gridSize = (arraySize + blockSize - 1) / blockSize;
  return {gridSize, blockSize};
}

template <typename T, std::array<BC, 6> bcs = EULERCFD::CONSTS::bcs>
struct GridInfo {
  std::size_t _nx, _ny, _nz, _nghosts;
  T _lx, _ly, _lz;
  dev_host constexpr std::size_t nx() const noexcept { return _nx; }
  dev_host constexpr std::size_t ny() const noexcept { return _ny; }
  dev_host constexpr std::size_t nz() const noexcept { return _nz; }
  dev_host constexpr std::size_t size() const noexcept { return nx()*ny()*nz(); }
  dev_host constexpr T lx() const noexcept { return _lx; }
  dev_host constexpr T ly() const noexcept { return _ly; }
  dev_host constexpr T lz() const noexcept { return _lz; }
  dev_host constexpr T dsx() const noexcept { return EULERCFD::CONSTS::DELTA; }
  dev_host constexpr T dsy() const noexcept { return EULERCFD::CONSTS::DELTA; }
  dev_host constexpr T dsz() const noexcept { return EULERCFD::CONSTS::DELTA; }
  dev_host constexpr T dv() const noexcept { return dsx() * dsy() * dsz(); }
};

enum BACKEND { HOST, DEVICE };

template <typename T, GridInfo<T> Info, BACKEND backend = BACKEND::HOST>
class Matrix3d {
public:
  Matrix3d() {
    _data = allocate(size());
    if (!_data) {
      spdlog::critical("ERROR: Failed to allocate");
      abort();
    }
  }

  Matrix3d(const Matrix3d &other) {
    if (nx() != other.nx() || ny() != other.ny() || nz() != other.nz()) {
      spdlog::critical("ERROR: Failed to construct due to size mismatch");
      throw std::runtime_error("Size Mismatch");
    }
    _data = allocate(size());
    if constexpr (backend == BACKEND::HOST) {
      std::memcpy(_data, other._data, size() * sizeof(T));
    }
    if constexpr (backend == BACKEND::DEVICE) {
      cudaMemcpy(_data, other._data, size() * sizeof(T),
                 cudaMemcpyDeviceToDevice);
    }
  }

  Matrix3d<T, Info, BACKEND::DEVICE> &
  operator=(const Matrix3d<T, Info, BACKEND::HOST> &other) {
    if (nx() != other.nx() || ny() != other.ny() || nz() != other.nz()) {
      spdlog::critical("ERROR: Failed to assign due to size mismatch");
      throw std::runtime_error("Size Mismatch");
    }
    cudaMemcpy(_data, other._data, size() * sizeof(T), cudaMemcpyHostToDevice);
    return *this;
  }

  ~Matrix3d() {
    if (_data) {
      deallocate(_data, size());
    }
  }

  host_only void export_to_host(T *buffer) const noexcept {
    if constexpr (backend == BACKEND::HOST) {
      std::memcpy(buffer, _data, size() * sizeof(T));
    }
    if constexpr (backend == BACKEND::DEVICE) {
      cudaMemcpy(buffer, _data, size() * sizeof(T), cudaMemcpyDeviceToHost);
    }
  }

  T &operator()(std::size_t i, std::size_t j, std::size_t k) noexcept {
    return _data[index(i, j, k)];
  }

  T &operator()(std::size_t i, std::size_t j, std::size_t k) const noexcept {
    return _data[index(i, j, k)];
  }

  const T *data() const noexcept { return _data; }

  T *data() noexcept { return _data; }

   constexpr std::size_t size() const noexcept {
    return Info._nx * Info._ny * Info._nz;
  }
  dev_host constexpr std::size_t nx() const noexcept { return Info._nx; }
  dev_host constexpr std::size_t ny() const noexcept { return Info._ny; }
  dev_host constexpr std::size_t nz() const noexcept { return Info._nz; }

  dev_host void fill(T &&val) noexcept {

    if constexpr (backend == BACKEND::HOST) {
      std::fill(_data, &_data[size() + 1], val);
    }
    if constexpr (backend == BACKEND::DEVICE) {
      abort();
    }
  }

  // private:
  dev_host constexpr std::size_t index(std::size_t i, std::size_t j,
                                       std::size_t k) const noexcept {
    return id_f(i,j,k);
  }

  host_only T *allocate(std::size_t elements) {
    T *ptr = nullptr;
    if constexpr (backend == BACKEND::HOST) {
      spdlog::debug("Host allocation of [{0:d} elements and  {1:d}] bytes ", elements,elements * sizeof(T));
      ptr = (T *)HOST_MATRIX_MALLOC(elements * sizeof(T));
      memset(ptr, 0, elements * sizeof(T));
    }
    if constexpr (backend == BACKEND::DEVICE) {
      spdlog::debug("Cuda allocation of [{0:d} elemenents and {1:d}] bytes ", elements,elements * sizeof(T));
      DEVICE_MATRIX_MALLOC(&ptr, elements * sizeof(T));
      cudaMemset(ptr, 0, elements * sizeof(T));
    }
    if (!ptr){
      spdlog::critical("ERROR: Failed to allocate [{0:d}] bytes",elements*sizeof(T));
      throw std::bad_alloc();
    }
    return ptr;
  }

  host_only void deallocate(T *ptr, std::size_t elements) {
    if constexpr (backend == BACKEND::HOST) {
      spdlog::debug("Host deallocation of [{0:d}] bytes ",
                    elements * sizeof(T));
      HOST_MATRIX_FREE(ptr);
    }
    if constexpr (backend == BACKEND::DEVICE) {
      spdlog::debug("Cuda deallocation of [{0:d}] bytes ",
                    elements * sizeof(T));
      DEVICE_MATRIX_FREE(ptr);
    }
  }

  T *_data;
};
} // namespace EULERCFD
