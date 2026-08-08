project(
    'pyvttrac',
    'c', 'fortran',
    version: '2.2.0',
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
if get_option('debug_bounds_check')
    fortran_args += fc.get_supported_arguments('-fcheck=bounds', '-fbacktrace')
endif

# Never add -ffast-math/-Ofast (or any flag implying -fno-signed-zeros/
# -ffinite-math-only) to fortran_args: core.f90's find_score_peak() relies on
# IEEE-compliant NaN comparisons (NaN >= x is false) to skip missing scores,
# and store_score()'s division-by-zero-variance paths rely on IEEE-compliant
# 0/0 -> NaN rather than a trap or an optimizer-assumed-unreachable branch.

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
