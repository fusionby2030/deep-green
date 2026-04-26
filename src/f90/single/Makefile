FC = gfortran
FOPT = -fopenmp -ffree-line-length-none -Wno-unused
H5FC = /opt/homebrew/Cellar/hdf5/2.1.1/bin/h5fc

# BYO HDF5!!! 
HDF5 = /opt/homebrew/Cellar/hdf5/2.1.1
HDF5_LIBLINK = -Wl,-rpath,$(HDF5)/lib
HDF5_LIB = -L$(HDF5)/lib -lhdf5_hl_fortran -lhdf5_fortran -lhdf5
HDF5_INC = -I$(HDF5)/include 

all:
	$(FC) -c types_and_kinds.f90
#   $(H5FC) -c io.f90
	$(FC) -c io.f90 $(HDF5_INC)
	$(FC) $(FOPT) -c global.f90 initialization.f90 physics.f90 	 
#	$(H5FC) *.o euler_cfd.f90 $(FOPT) -o ./main 
	$(FC) $(HDF5_LIBLINK) *.o euler_cfd.f90 $(FOPT) -o ./main $(HDF5_INC) $(HDF5_LIB)	
	mkdir ./results # because we write state there  
clean:
	rm *.mod *.o ./main
	rm -rf ./results
