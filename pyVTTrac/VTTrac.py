from pathlib import Path

from julia.api import Julia
jl = Julia(compiled_modules=False)

from julia import Main, Pkg
path = Path(__file__).parent.parent / "VTTrac.jl/src/VTTrac.jl"
Pkg.activate(str(path.parent.parent))
Main.include(str(path))

import numpy as np

class VTT:
    """
    Say hello using Julia.
    """
    def __init__(self, z, tax=None, zmiss=-999.0, fmiss=-999.0, imiss=-999):
        z = np.array(z, np.float32)
        self.vtt = Main.VTTrac.VTT(z, tax, zmiss, fmiss, imiss)
    
