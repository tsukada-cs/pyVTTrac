project(
    'pyvttrac',
    'c', 'fortran',
    version: '2.0.0',
    license: 'BSD-2-Clause',
    meson_version: '>=1.1.0',
    default_options: [
        'buildtype=release',
        'fortran_std=legacy',
    ],
)

py_mod = import('python')
py = py_mod.find_installation(pure: false)

fc = meson.get_compiler('fortran')

omp_dep = dependency('openmp', language: 'fortran', required: get_option('use_openmp'))

fortran_args = ['-cpp']
if get_option('native_optimization')
    fortran_args += fc.get_supported_arguments('-march=native')
endif

incdir_numpy = run_command(py,
    ['-c', 'import numpy; print(numpy.get_include())'],
    check: true
).stdout().strip()

incdir_f2py = run_command(py,
    ['-c', 'import numpy.f2py; print(numpy.f2py.get_include())'],
    check: true
).stdout().strip()

inc_np = include_directories(incdir_numpy, incdir_f2py)

subdir('src/pyvttrac')
