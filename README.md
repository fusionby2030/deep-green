# deep-green

To build and run the project, we use [meson](https://mesonbuild.com/SimpleStart.html). 

## Build, Run, and

### Conduction Standalone 

0. `cd conduction`
1. `meson setup build`
2. `meson compile -C build`
3. `./build/conduction ./run/conductive_input.nml`

The output files will be dumped to `/results/fname.dat`. 

The input file can be modified in `run/conductive_input.nml`

## Convection Standalone 

0. `cd convection`
1. `meson setup build --buildtype==release`
2. `meson compile -C build`
3. `./build/convection`

## Coupled 

0. `cd coupled && mkdir results`
1. `python3 -m numpy.f2py  -c drivers.f90 -m driver --fcompiler=gfortran --f90flags='-O3 -ffree-line-length-none'`
2. `python3 main.py run/*.yaml`

### Scalings 

0. `cd scalings`
1. `'./chmod +x generate_configs.sh`
2. `python3 -m numpy.f2py  -c drivers.f90 -m driver --fcompiler=gfortran --f90flags='-O3 -ffree-line-length-none'`
2. `mkdir results`
3. 