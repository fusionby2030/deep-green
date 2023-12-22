project('ghouse', 'fortran',
  version : '0.1',
  default_options : ['warning_level=3'])

add_global_arguments('-ffree-line-length-none', language : 'fortran')
executable('ghouse',
           'euler_cfd.f90',
           install : true)
