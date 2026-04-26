FC = gfortran
FOPT = -fopenmp -ffree-line-length-none -Wno-unused

H5FC = /opt/homebrew/Cellar/hdf5/2.1.1/bin/h5fc

all:
	$(FC) -c types_and_kinds.f90
	$(H5FC) -c io.f90
	$(FC) $(FOPT) -c global.f90 initialization.f90 physics.f90 	 
	$(H5FC) *.o euler_cfd.f90 $(FOPT) -o ./main 
	
clean:
	rm *.mod *.o ./main
