project('conductive', 'fortran',
  version : '0.1',
  default_options : ['warning_level=3'])


hdf5 = dependency('hdf5', language: 'fortran',required : true)
add_global_arguments('-fopenmp','-ffree-line-length-none','-Wno-unused', language : 'fortran')

executable('conductive',
           'types_and_kinds.f03',
           'io.f03',
           'globals.f03',
           'grid.f03',
           'initialization.f03',
           'physics.f03',
           'main.f03',
           link_args:['-fopenmp'],
           dependencies : hdf5,
           install : true)
