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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 * */
#pragma once
#include "matrix3d.hpp"
#include <array>
#include <cstddef>
#include <highfive/highfive.hpp>
#include <vector>

namespace EULERCFD {

template <typename T, GridInfo<T> Info, BACKEND Backend> class Grid {
public:

  Grid(){}
  Grid(const Grid& other)=delete;
  Grid(Grid&& other)=delete;
  Grid& operator=(const Grid& other)=delete;
  Grid& operator=(Grid&& other)=delete;
  ~Grid(){}
  
  Matrix3d<T, Info, Backend> rho, vx, vy, vz, p;
  Matrix3d<T, Info, Backend> rho_xtr, vx_xtr, vy_xtr, vz_xtr, p_xtr;
  Matrix3d<T, Info, Backend> drho_x, drho_y, drho_z, dvx_x, dvx_y, dvx_z, dvy_x,
      dvy_y, dvy_z, dvz_x, dvz_y, dvz_z, dp_x, dp_y, dp_z;
  Matrix3d<T, Info, Backend> mass, momx, momy, momz, e;
  Matrix3d<T, Info, Backend> fmass_x, fmass_y, fmass_z, fmomx_x, fmomx_y,
      fmomx_z, fmomy_x, fmomy_y, fmomy_z, fmomz_x, fmomz_y, fmomz_z, fe_x, fe_y,
      fe_z;

  dev_host constexpr std::size_t size() const noexcept {
    return Info._nx * Info._ny * Info._nz;
  }
  dev_host constexpr std::size_t nx() const noexcept { return Info._nx; }
  dev_host constexpr std::size_t ny() const noexcept { return Info._ny; }
  dev_host constexpr std::size_t nz() const noexcept { return Info._nz; }
  dev_host constexpr T lx() const noexcept { return Info._lx; }
  dev_host constexpr T ly() const noexcept { return Info._ly; }
  dev_host constexpr T lz() const noexcept { return Info._lz; }
  dev_host constexpr T dsx() const noexcept { return EULERCFD::CONSTS::DELTA; }
  dev_host constexpr T dsy() const noexcept { return EULERCFD::CONSTS::DELTA; }
  dev_host constexpr T dsz() const noexcept { return EULERCFD::CONSTS::DELTA; }
  dev_host constexpr T dv() const noexcept { return dsx() * dsy() * dsz(); }

  constexpr std::array<std::size_t, 3> get_dims() const noexcept {
    return {nx(), ny(), nz()};
  }

  constexpr std::array<T *, 5> get_primitive_pointers() noexcept {
    return {rho.data(), vx.data(), vy.data(), vz.data(), p.data()};
  }

  constexpr std::array<T *, 5> get_primitive_xtr_pointers() noexcept {
    return {rho_xtr.data(), vx_xtr.data(), vy_xtr.data(), vz_xtr.data(),
            p_xtr.data()};
  }

  constexpr std::array<T *, 15> get_gradients_pointers() noexcept {
    return {drho_x.data(), drho_y.data(), drho_z.data(), dvx_x.data(),
            dvx_y.data(),  dvx_z.data(),  dvy_x.data(),  dvy_y.data(),
            dvy_z.data(),  dvz_x.data(),  dvz_y.data(),  dvz_z.data(),
            dp_x.data(),   dp_y.data(),   dp_z.data()};
  }

  constexpr std::array<T *, 15> get_fluxes_pointers() noexcept {
    return {fmass_x.data(), fmass_y.data(), fmass_z.data(), fmomx_x.data(),
            fmomx_y.data(), fmomx_z.data(), fmomy_x.data(), fmomy_y.data(),
            fmomy_z.data(), fmomz_x.data(), fmomz_y.data(), fmomz_z.data(),
            fe_x.data(),    fe_y.data(),    fe_z.data()};
  }

  constexpr std::array<T *, 5> get_conserved_pointers() noexcept {
    return {mass.data(), momx.data(), momy.data(), momz.data(), e.data()};
  }

  void store(const char *filename) const {
    using namespace HighFive;
    File file(filename, File::Truncate);
    std::vector<T> buffer(size(), 0);

    rho.export_to_host(buffer.data());
    auto dset = file.createDataSet<T>("primitives/rho", DataSpace(get_dims()));
    dset.write_raw(buffer.data());

    vx.export_to_host(buffer.data());
    dset = file.createDataSet<T>("primitives/vx", DataSpace(get_dims()));
    dset.write_raw(buffer.data());

    vy.export_to_host(buffer.data());
    dset = file.createDataSet<T>("primitives/vy", DataSpace(get_dims()));
    dset.write_raw(buffer.data());

    vz.export_to_host(buffer.data());
    dset = file.createDataSet<T>("primitives/vz", DataSpace(get_dims()));
    dset.write_raw(buffer.data());

    p.export_to_host(buffer.data());
    dset = file.createDataSet<T>("primitives/p", DataSpace(get_dims()));
    dset.write_raw(buffer.data());

    // mass.export_to_host(buffer.data());
    // dset = file.createDataSet<T>("conserved/mass", DataSpace(get_dims()));
    // dset.write_raw(buffer.data());

    // momx.export_to_host(buffer.data());
    // dset = file.createDataSet<T>("conserved/momx", DataSpace(get_dims()));
    // dset.write_raw(buffer.data());

    // momy.export_to_host(buffer.data());
    // dset = file.createDataSet<T>("conserved/momy", DataSpace(get_dims()));
    // dset.write_raw(buffer.data());

    // momz.export_to_host(buffer.data());
    // dset = file.createDataSet<T>("conserved/momz", DataSpace(get_dims()));
    // dset.write_raw(buffer.data());

    // e.export_to_host(buffer.data());
    // dset = file.createDataSet<T>("conserved/e", DataSpace(get_dims()));
    // dset.write_raw(buffer.data());

    // drho_x.export_to_host(buffer.data());
    // dset = file.createDataSet<T>("gradients/rhox", DataSpace(get_dims()));
    // dset.write_raw(buffer.data());

    // drho_y.export_to_host(buffer.data());
    // dset = file.createDataSet<T>("gradients/rhoy", DataSpace(get_dims()));
    // dset.write_raw(buffer.data());

    // drho_z.export_to_host(buffer.data());
    // dset = file.createDataSet<T>("gradients/rhoz", DataSpace(get_dims()));
    // dset.write_raw(buffer.data());
    

    fmass_x.export_to_host(buffer.data());
    dset = file.createDataSet<T>("fluxes/fmassx", DataSpace(get_dims()));
    dset.write_raw(buffer.data());

    fmass_y.export_to_host(buffer.data());
    dset = file.createDataSet<T>("fluxes/fmassy", DataSpace(get_dims()));
    dset.write_raw(buffer.data());

    fmass_z.export_to_host(buffer.data());
    dset = file.createDataSet<T>("fluxes/fmassz", DataSpace(get_dims()));
    dset.write_raw(buffer.data());


    fmomx_x.export_to_host(buffer.data());
    dset = file.createDataSet<T>("fluxes/fmomx_x", DataSpace(get_dims()));
    dset.write_raw(buffer.data());

    fmomx_y.export_to_host(buffer.data());
    dset = file.createDataSet<T>("fluxes/fmomx_y", DataSpace(get_dims()));
    dset.write_raw(buffer.data());

    fmomx_z.export_to_host(buffer.data());
    dset = file.createDataSet<T>("fluxes/fmomx_z", DataSpace(get_dims()));
    dset.write_raw(buffer.data());
    
    fmomy_x.export_to_host(buffer.data());
    dset = file.createDataSet<T>("fluxes/fmomy_x", DataSpace(get_dims()));
    dset.write_raw(buffer.data());

    fmomy_y.export_to_host(buffer.data());
    dset = file.createDataSet<T>("fluxes/fmomy_y", DataSpace(get_dims()));
    dset.write_raw(buffer.data());

    fmomy_z.export_to_host(buffer.data());
    dset = file.createDataSet<T>("fluxes/fmomy_z", DataSpace(get_dims()));
    dset.write_raw(buffer.data());
    
    fmomz_x.export_to_host(buffer.data());
    dset = file.createDataSet<T>("fluxes/fmomz_x", DataSpace(get_dims()));
    dset.write_raw(buffer.data());

    fmomz_y.export_to_host(buffer.data());
    dset = file.createDataSet<T>("fluxes/fmomz_y", DataSpace(get_dims()));
    dset.write_raw(buffer.data());

    fmomz_z.export_to_host(buffer.data());
    dset = file.createDataSet<T>("fluxes/fmomz_z", DataSpace(get_dims()));
    dset.write_raw(buffer.data());
    
    fe_x.export_to_host(buffer.data());
    dset = file.createDataSet<T>("fluxes/fe_x", DataSpace(get_dims()));
    dset.write_raw(buffer.data());

    fe_y.export_to_host(buffer.data());
    dset = file.createDataSet<T>("fluxes/fe_y", DataSpace(get_dims()));
    dset.write_raw(buffer.data());

    fe_z.export_to_host(buffer.data());
    dset = file.createDataSet<T>("fluxes/fe_z", DataSpace(get_dims()));
    dset.write_raw(buffer.data());

  }
  
};
} // namespace EULERCFD
