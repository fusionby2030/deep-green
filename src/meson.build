project('ghouse', 'fortran',
  version : '0.1',
  default_options : ['warning_level=3'])

hdf5 = dependency('hdf5', language: 'fortran',required : true)
add_global_arguments('-ffree-line-length-none', language : 'fortran')

executable('ghouse',
           'euler_cfd.f90',
           'types_and_kinds.f90',
           'physics.f90',
           'global.f90',
           'initialization.f90',
           'io.f90',
           dependencies : hdf5,
           install : true)
