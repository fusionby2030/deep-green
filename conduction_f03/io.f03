! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
module io
   use types_and_kinds
   use simgrid
   use hdf5
   IMPLICIT NONE 
   private :: add_field_h5, create_hdf5file, open_hdf5file, close_hdf5file
   public :: statewrite
contains
   !only for double 8 and same 3D shape for all arrays
   !Credits to https://gist.github.com/dmentipl/ed609e7278050bd6de0016c09155e039
   subroutine add_field_h5(id, name, x, error, compression_level)
      real(rk), intent(in)  :: x(:, :, :)
      character(*), intent(in)  :: name
      integer(hid_t), intent(in)  :: id
      integer(4), intent(in)::compression_level
      integer, intent(out) :: error

      integer, parameter :: ndims = 3
      integer(hsize_t)   :: xshape(ndims)
      integer(hsize_t)   :: chunk(ndims)
      integer(hid_t)     :: dspace_id
      integer(hid_t)     :: dset_id
      integer(hid_t)     :: prop_id
      integer(hid_t)     :: dtype_id

      xshape = shape(x)
      chunk = shape(x)
      dtype_id = h5t_native_double

      ! create dataspace
      call h5screate_simple_f(ndims, xshape, dspace_id, error)
      if (error /= 0) then
         write (*, '("cannot create hdf5 dataspace",/)')
         return
      end if

      call h5pcreate_f(h5p_dataset_create_f, prop_id, error)
      call h5pset_deflate_f(prop_id, compression_level, error)
      call h5pset_chunk_f(prop_id, ndims, chunk, error)
      if (error /= 0) then
         write (*, '("cannot create hdf5 property list",/)')
         return
      end if

      ! create dataset in file
      call h5dcreate_f(id, name, dtype_id, dspace_id, dset_id, error, prop_id)
      if (error /= 0) then
         write (*, '("cannot create hdf5 dataset",/)')
         return
      end if

      ! write to file
      call h5dwrite_f(dset_id, dtype_id, x, xshape, error)
      if (error /= 0) then
         write (*, '("cannot write to hdf5 file",/)')
         return
      end if

      ! close property list
      call h5pclose_f(prop_id, error)
      if (error /= 0) then
         write (*, '("cannot close hdf5 property list",/)')
         return
      end if

      ! close dataset
      call h5dclose_f(dset_id, error)
      if (error /= 0) then
         write (*, '("cannot close hdf5 dataset",/)')
         return
      end if

      ! close dataspace
      call h5sclose_f(dspace_id, error)
      if (error /= 0) then
         write (*, '("cannot close hdf5 dataspace",/)')
         return
      end if

   end subroutine add_field_h5

   subroutine create_hdf5file(filename, file_id, error)
      character(len=*), intent(in)  :: filename
      integer(hid_t), intent(out) :: file_id
      integer, intent(out) :: error
      integer :: filter_info
      integer :: filter_info_both
      logical :: avail

      ! initialise hdf5
      call h5open_f(error)
      if (error /= 0) then
         write (*, '("cannot initialise hdf5",/)')
         return
      end if

      ! check if gzip compression is available.
      call h5zfilter_avail_f(h5z_filter_deflate_f, avail, error)
      if (.not. avail) then
         write (*, '("gzip filter not available.",/)')
         return
      end if
      call h5zget_filter_info_f(h5z_filter_deflate_f, filter_info, error)
      filter_info_both = ior(h5z_filter_encode_enabled_f, h5z_filter_decode_enabled_f)
      if (filter_info /= filter_info_both) then
         write (*, '("gzip filter not available for encoding and decoding.",/)')
         return
      end if

      ! create file
      call h5fcreate_f(filename, h5f_acc_trunc_f, file_id, error)
      if (error /= 0) then
         write (*, '("cannot create hdf5 file",/)')
         return
      end if

   end subroutine create_hdf5file

   !Credits to https://gist.github.com/dmentipl/ed609e7278050bd6de0016c09155e039
   subroutine open_hdf5file(filename, file_id, error)
      character(len=*), intent(in)  :: filename
      integer(hid_t), intent(out) :: file_id
      integer, intent(out) :: error

      ! initialise hdf5
      call h5open_f(error)
      if (error /= 0) then
         write (*, '("cannot initialise hdf5",/)')
         return
      end if

      ! open file
      call h5fopen_f(filename, h5f_acc_rdwr_f, file_id, error)
      if (error /= 0) then
         write (*, '("cannot open hdf5 file",/)')
         return
      end if

   end subroutine open_hdf5file

   !Credits to https://gist.github.com/dmentipl/ed609e7278050bd6de0016c09155e039
   subroutine close_hdf5file(file_id, error)
      integer(hid_t), intent(in)  :: file_id
      integer, intent(out) :: error

      ! close file
      call h5fclose_f(file_id, error)
      if (error /= 0) then
         write (*, '("cannot close hdf5 file",/)')
         return
      end if

      ! close hdf5
      call h5close_f(error)
      if (error /= 0) then
         write (*, '("cannot close hdf5",/)')
         return
      end if

   end subroutine close_hdf5file

   subroutine statewrite(time, tstep, grid)
      integer(ik), intent(in) :: tstep
      real(rk), intent(in) :: time
      type(SimulationGrid), pointer :: grid
      integer(4):: compression_level = 0
      character(32) :: fname
      integer(hid_t) :: file_id
      integer:: error
      write (fname, '(a,"000",i7.7,a)') "state.", tstep, ".h5"
      call create_hdf5file(fname, file_id, error)
      call add_field_h5(file_id, "T", grid%primitives%T, error, compression_level)
      call add_field_h5(file_id, "TNEW", grid%primitives%TNEW, error, compression_level)
      call add_field_h5(file_id, "alpha", grid%coefficients%alpha, error, compression_level)
      call add_field_h5(file_id, 'source1', grid%sources%s1, error, compression_level)
      ! call add_field_h5(file_id, "time", time, error, compression_level)
      ! call add_field_h5(file_id, "rho", rho, error, compression_level)
      ! call add_field_h5(file_id, "vx", vx, error, compression_level)
      ! call add_field_h5(file_id, "vy", vy, error, compression_level)
      ! call add_field_h5(file_id, "vz", vz, error, compression_level)
      ! call add_field_h5(file_id, "p", p, error, compression_level)
      ! call add_field_h5(file_id, "mass", mass, error, compression_level)
      ! call add_field_h5(file_id, "momentum_x", momentum_x, error, compression_level)
      ! call add_field_h5(file_id, "momentum_y", momentum_y, error, compression_level)
      ! call add_field_h5(file_id, "momentum_z", momentum_z, error, compression_level)
      ! call add_field_h5(file_id, "energy", energy, error, compression_level)
      ! call add_field_h5(file_id, "temperature", temp, error, compression_level)
      call close_hdf5file(file_id, error)
   end subroutine statewrite
end module io
