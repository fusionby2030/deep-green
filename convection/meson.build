project('ghouse', 'fortran',
  version : '0.1',
  default_options : ['warning_level=3'])

run_command('mkdir', 'results')
hdf5 = dependency('hdf5', language: 'fortran',required : true)
add_global_arguments('-fopenmp','-ffree-line-length-none','-Wno-unused', language : 'fortran')

executable('convection',
           'euler_cfd.f90',
           'types_and_kinds.f90',
           'physics.f90',
           'global.f90',
           'initialization.f90',
           'io.f90',
           link_args:['-fopenmp'],
           dependencies : hdf5,
           install : true)
