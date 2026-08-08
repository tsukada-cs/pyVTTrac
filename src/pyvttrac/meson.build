# Build the `pyvttrac._core` extension module from the Fortran computational
# core (core.f90) via a hand-maintained f2py signature file (_core.pyf). 
# f2py itself only generates the C/Fortran wrapper glue here (via a custom_target); 
# actual compilation goes through meson so that OpenMP, optimization flags, 
# and cross-platform wheel builds are handled uniformly with the rest of the project.

core_pyf = files('_core.pyf')
core_f90 = files('core.f90')

_core_sources = custom_target('_core_sources',
    input: core_pyf,
    output: ['_coremodule.c', '_core-f2pywrappers2.f90'],
    command: [
        py, '-m', 'numpy.f2py', '@INPUT@',
        '--build-dir', '@OUTDIR@',
    ],
)

# `_core` mixes Fortran (core.f90, only source with `!$omp` pragmas) and C
# (f2py's wrapper glue) sources. OpenMP must NOT be attached to the target
# via `dependencies:` -- that applies its compile args to *every* source
# regardless of language, including `-fopenmp` on the C files, which Apple
# clang (the default `cc` on macOS) rejects outright ("unsupported option
# '-fopenmp'"); plain GCC on Linux silently tolerates the same mistake,
# which is why this only broke on macOS. Instead, add the flag only to
# `fortran_args` (compilation) and `link_args` (the final link step, which
# runs through gfortran because of `link_language: 'fortran'` below).
omp_fortran_args = []
omp_link_args = []
if omp_dep.found()
    omp_fortran_args = ['-fopenmp']
    omp_link_args = ['-fopenmp']
endif

py.extension_module(
    '_core',
    [core_f90, _core_sources, incdir_f2py / 'fortranobject.c'],
    fortran_args: fortran_args + omp_fortran_args,
    link_args: omp_link_args,
    include_directories: inc_np,
    link_language: 'fortran',
    install: true,
    subdir: 'pyvttrac',
)

python_sources = [
    '__init__.py',
    'api.py',
    'bidirectional.py',
    'config.py',
    'grid.py',
    'result.py',
    'seeds.py',
    'status.py',
]

py.install_sources(
    python_sources,
    subdir: 'pyvttrac',
)
