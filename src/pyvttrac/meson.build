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

py.extension_module(
    '_core',
    [core_f90, _core_sources, incdir_f2py / 'fortranobject.c'],
    fortran_args: fortran_args,
    include_directories: inc_np,
    dependencies: omp_dep.found() ? [omp_dep] : [],
    link_language: 'fortran',
    install: true,
    subdir: 'pyvttrac',
)

python_sources = [
    '__init__.py',
    'api.py',
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
