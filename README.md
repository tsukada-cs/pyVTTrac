# pyVTTrac

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://tsukada-cs.github.io/pyVTTrac/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://tsukada-cs.github.io/pyVTTrac/dev) -->
<!-- [![Build Status](https://github.com/tsukada-cs/pyVTTrac/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/tsukada-cs/pyVTTrac/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/tsukada-cs/pyVTTrac/branch/main/graph/badge.svg)](https://codecov.io/gh/tsukada-cs/pyVTTrac) -->

## pyVTTrac: Python implementation for Velocimetry by Template Tracking
This library provides the python implementation for `VTTrac.jl` (https://github.com/tsukada-cs/VTTrac.jl). It does not use module variables, so it should be good for parallel execution.

The algorithm used in this library is the simple template matching of PIV (particle image velocimetry) for monochromatic image-like data, but the matching is conducted multiple times in a Lagrangian manner as in PTV (particle tracking velocimetry) over a number of times specified by the parameter named `ntrac`. The default scoring method for template matching is the cross correlation coefficient, as in the basic PIV. Both forward and backward tracking is available. Use the parameter `itstep`; tracking is backward along time sequence, if it is negative.

### Available scoring methods
* XCOR: cross-correlation, `cov(x',y')/sig(x)/sig(y)`, where `x` represents the template sub-image and `y` represents the target sub-image that are slided.
* NCOV: normalized covariance, `cov(x',y')/sig(x)^2`: covariance normalized the variance of the template sub-image `x`.

### Screening
A check (result screening) based on velocity change along trajectory is available (by using the threshold parameter named `vxch` and `vych`), so it is recommended to always set `ntrac >= 2`. Further screening is available for initial templates (e.g., in terms of the complexity and contrast) and the quality of the results (e.g., score threshold); see the source code.

### Dimensions
Spatial coordinates are based on array indices, with the distance between adjacent grid points always being 1, so they are non-dimensional. The velocities are based on non-dimensional spacial displacement over time difference, where, time can either be dimensional or non-dimensional.

### Related packages
* `VTTrac.jl` by Taiga Tsukada: https://github.com/tsukada-cs/VTTrac.jl (submodule of this module)
* `VTTrac` by Takeshi Horinouchi: https://github.com/thorinouchi/VTTrac (`VTTrac.jl` were made as julia-language implementation for it)

### References
* Horinouchi, T., S. Tsujino, M. Hayashi, U. Shimada, W. Yanase, A. Wada, and H. Yamada, 2023: Stationary and Transient Asymmetric Features in Tropical Cyclone Eye with Wavenumber-1 Instability: Case Study for Typhoon Haishen (2020) with Atmospheric Motion Vectors from 30-Second Imaging. Monthly Weather Review, 151, 253–273, https://doi.org/10.1175/MWR-D-22-0179.1.
* Tsukada, T., T. Horinouchi, and S. Tsujino, 2024: Wind Distribution in the Eye of Tropical Cyclone Revealed by a Novel Atmospheric Motion Vector Derivation. JGR Atmospheres, 129, e2023JD040585, https://doi.org/10.1029/2023JD040585.


## Installation and Test
### How to install
```shell
$ git clone --recurse-submodules https://github.com/tsukada-cs/pyVTTrac.git
$ cd pyVTTrac
$ pip install .
```
### How to test
```shell
$ pip install pytest
$ pytest -s
```
If pytest raise an error like "ERROR: Unable to load dependent library", you may need to add your library directory to the `LD_LIBRARY_PATH` envrionmental variable.
For example, if you use anaconda3 or miniconda3, you can add the following line to your `.bashrc` for bash user or `.zshrc` for zsh user.
```shell
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/anaconda3/lib
```

### How to run
```shell
$ python example/sample.py
```