project('conduction', 'fortran',
  version : '0.1',
  default_options : ['warning_level=3'])

run_command('mkdir', 'results')
add_global_arguments('-O3', language : 'fortran')
executable('conduction',
           'heat_equation.f90',
           install : true)
