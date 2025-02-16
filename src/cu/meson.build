project(
    'eulernv',
    ['cuda', 'cpp'],
    version: '3.0',
    default_options: ['warning_level=3', 'cpp_std=c++20'],
)
add_global_arguments(
    [
        '-Xcompiler',
        '-Wno-pedantic',
        '--std=c++20',
        '-Xcudafe',
        '--diag_suppress=20208',
    ],
    language: 'cuda',
)

CUSTOM_INSTALL_PATH = '/home/kstppd/software/'

spdlog_include = include_directories(
    CUSTOM_INSTALL_PATH + '/spdlog/include/',
    is_system: true,
)
h5_include = include_directories(
    CUSTOM_INSTALL_PATH + '/HighFive/include',
    is_system: true,
)
libspdlog = meson.get_compiler('cuda').find_library(
    'spdlog',
    dirs: [CUSTOM_INSTALL_PATH + 'spdlog/build/'],
)
cuda_dep = dependency('cuda', version: '>=11', modules: ['nvToolsExt'])
hdf5_dep = dependency('hdf5', language: 'cpp')


includes = [spdlog_include, h5_include]
deps = [cuda_dep, hdf5_dep]
add_project_arguments(
    [
        '--expt-relaxed-constexpr',
        '-std=c++20',
        '--expt-extended-lambda',
        '-Xcompiler',
        '-fopenmp',
    ],
    language: 'cuda',
)
executable(
    'eulernv',
    sources: ['eulernv.cu'],
    include_directories: includes,
    link_args: '-fopenmp',
    dependencies: deps,
    install: true,
)



