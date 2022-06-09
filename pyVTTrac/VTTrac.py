from pathlib import Path

import julia
julia.install()

from julia.api import Julia
jl = Julia(compiled_modules=False)

from julia import Main, Pkg
path = Path(__file__).parent.parent / "VTTrac.jl/src/VTTrac.jl"
Pkg.activate(str(path.parent.parent))
Main.include(str(path))

import numpy as np
import xarray as xr

class VTT:
    def __init__(self, z, t=None, zmiss=None, fmiss=-999.0, imiss=-999):
        """
        Object initialization (VTT)

        You need to call #setup subsequently before calling #trac

        Parameters
        ----------
        z: array-like
            The "images" used for tracking.
        t: array-like or None
            Time for each image (if None, [0,1,2,..]).
        zmiss: float
            Missing value in z (None if there is no missing).
        """
        z = np.array(z, np.float32)
        self.o = Main.VTTrac.VTT(z, t=t, zmiss=zmiss, fmiss=fmiss, imiss=imiss)
        
        self.attrs = {}
        self.attrs["nt"] = self.o.nt
        self.attrs["ny"] = self.o.ny
        self.attrs["nx"] = self.o.nx
        self.attrs["dtmean"] = self.o.dtmean
        self.attrs["zmiss"] = self.o.zmiss
        self.attrs["chk_zmiss"] = self.o.chk_zmiss
        self.attrs["fmiss"] = self.o.fmiss
        self.attrs["imiss"] = self.o.imiss

    def __repr__(self):
        return repr(self.attrs)    

    @property
    def nt(self):
        return self.attrs["nt"]
    @property
    def ny(self):
        return self.attrs["ny"]
    @property
    def nx(self):
        return self.attrs["nx"]
    @property
    def dtmean(self):
        return self.attrs["dtmean"]
    
    def __setitem__(self, key, value):
        self.attrs[key] = value
        exec(f"Main.VTTrac.set_{key}(self.o, value)")

    def setup(self, nsx, nsy, vxhw=None, vyhw=None, ixhw=None, iyhw=None, subgrid=True, subgrid_gaus=False,
        itstep=1, ntrac=2, score_method="xcor", score_th0=0.8, score_th1=0.7, vxch=None, vych=None,
        peak_inside_th=None, min_contrast=None, use_init_temp=False):
        """
        Setup for tracking.

        Parameters
        ----------
        o: VTT
            The object.
        nsx: int
            Submimage x sizes.
        nsy: int
            Submimage y sizes.
        vxch: float or None
            (either `v[xy]hw` or `i[xy]hw` are MANDATORY).
            the dimensions along which to perform the computation.
            search velocity range half sizes to set `i[xy]hw`.
            Seach at least to cover +-v?hw around the first guess or previous step.
            (the result can be outside the range.)
        vyhw: float or None
            See `vxhw` document.
        ixhw: int or None
            (either `v[xy]hw` or `i[xy]hw` are MANDATORY)
            Max displacement fro template match (can be set indirecly through `v[xy]hw`).
        iyhw: int or None
            See `ixhw` document.
        subgrid: bool, default True
            Whether to conduct subgrid tracking.
        subgrid_gaus: bool, default True
            Whether subgrid peak finding is by gaussian.
        itstep: int, default 1
            Step of `t`'s used (skip if >1).
        ntrack: int, default 2
            Max tracking times from initial loc.
        score_method: str, default "xcor"`
            `"xcor"` for cross-correlation, `"ncov"` for normalized covariance.
        score_th0: float, default 0.8
            The minimum score required for the 1st tracking.
        score_th1: float, default 0.7
            The minimum score required for subsequent tracking.
        vxch: float, optional
            If non-`nothing`, the max tolerant vx
            change between two consecutive tracking.
        vych: float, optional
            If non-`nothing`, the max tolerant vy
            change between two consecutive tracking.
        peak_inside_th: float, optional
            If non-`nothing`, an initial template is used only when it is peaked (max or min) inside,
            exceeding the max or min along the sides by the ratio specified by its value.
        min_contrast: float, optional
            If non-`nothing`, an initial template is used only when 
            it has a difference in max and min greater than its value.
        """

        Main.VTTrac.setup(
                self.o, nsx, nsy,
                vxhw=vxhw, vyhw=vyhw,
                ixhw=ixhw, iyhw=iyhw, subgrid=subgrid, subgrid_gaus=subgrid_gaus, itstep=itstep, ntrac=ntrac, score_method=score_method,
                score_th0=score_th0, score_th1=score_th1, vxch=vxch, vych=vych, peak_inside_th=peak_inside_th,
                min_contrast=min_contrast, use_init_temp=use_init_temp
            )

        self.attrs["nsx"] = self.o.nsx
        self.attrs["nsy"] = self.o.nsy
        self.attrs["vxhw"] = self.o.vxhw
        self.attrs["vyhw"] = self.o.vyhw
        self.attrs["ixhw"] = self.o.ixhw
        self.attrs["iyhw"] = self.o.iyhw
        self.attrs["subgrid"] = self.o.subgrid
        self.attrs["subgrid_gaus"] = self.o.subgrid_gaus
        self.attrs["itstep"] = self.o.itstep
        self.attrs["ntrac"] = self.o.ntrac
        self.attrs["score_method"] = self.o.score_method
        self.attrs["score_th0"] = self.o.score_th0
        self.attrs["score_th1"] = self.o.score_th1
        self.attrs["vxch"] = self.o.vxch
        self.attrs["vych"] = self.o.vych
        self.attrs["peak_inside_th"] = self.o.peak_inside_th
        self.attrs["min_contrast"] = self.o.min_contrast
        self.attrs["use_init_temp"] = self.o.use_init_temp
    
    @property
    def zmiss(self):
        return self.attrs["zmiss"]
    @property
    def fmiss(self):
        return self.attrs["fmiss"]
    @property
    def imiss(self):
        return self.attrs["imiss"]
    @property
    def nsx(self):
        return self.attrs["nsx"]
    @property
    def nsy(self):
        return self.attrs["nsy"]
    @property
    def vxhw(self):
        return self.attrs["vxhw"]
    @property
    def vyhw(self):
        return self.attrs["vyhw"]
    @property
    def ixhw(self):
        return self.attrs["ixhw"]
    @property
    def iyhw(self):
        return self.attrs["iyhw"]
    @property
    def subgrid(self):
        return self.attrs["subgrid"]
    @property
    def subgrid_gaus(self):
        return self.attrs["subgrid_gaus"]
    @property
    def itstep(self):
        return self.attrs["itstep"]
    @property
    def ntrac(self):
        return self.attrs["ntrac"]
    @property
    def score_method(self):
        return self.attrs["score_method"]
    @property
    def score_th0(self):
        return self.attrs["score_th0"]
    @property
    def score_th1(self):
        return self.attrs["score_th1"]
    @property
    def vxch(self):
        return self.attrs["vxch"]
    @property
    def vych(self):
        return self.attrs["vych"]
    @property
    def peak_inside_th(self):
        return self.attrs["peak_inside_th"]
    @property
    def min_contrast(self):
        return self.attrs["min_contrast"]
    @property
    def use_init_temp(self):
        return self.attrs["use_init_temp"]
    
    def trac(self, tid0, x0, y0, vxg0=None, vyg0=None, out_subimage=False, out_score_ary=False, asxarray=True):
        """
        Conduct tracking.

        Parameters
        ----------
        tid0: array-like
            Tracking initial time indices.
        x0: array-like
            Tracking initial template-center x location (index-based; non-integer for subgrid).
        y0: array-like
            Tracking initial template-center y location (index-based; non-integer for subgrid).
        vxg0: array-like
            First guess of vx (to search around it). Can be 0.
        vyg0: array-like
            First guess of vy (to search around it). Can be 0.
        out_subimage: bool, default False
            Whether output subimages.
        out_score_ary: bool, default False
            Whether output score arrays.
        asxarray: bool, default True
            Whether output as `xarray.Dataset`

        Returns
        -------
        count: np.ndarray
            The number of successful tracking for each initial template. Shape: [len]
        tid: np.ndarray
            time index of the trajectories (tid0 and subsequent ones). Shape: [ntrac+1, len]
        x: np.ndarray
            x locations of the trajectories (x0 and derived ones). Shape: [ntrac+1, len]
        y: np.ndarray
            y locations of trajectories (x0 and derived ones). Shape: [ntrac+1, len]
        vx: np.ndarray
            Derived x-velocity. Shape: [ntrac, len]
        vy: np.ndarray
            Derived y-velocity. Shape: [ntrac, len]
        score: np.ndarray
            Scores along the trajectory (max values, possibly at subgrid). Shape: [ntrac, len]
        zss: np.ndarray or None
            If `out_subimage` is `True`, the subimages along the track. Shape: [nsx, nsy, ,ntrac+1, len]
        score_arry: np.ndarray or None
            If `out_score_ary` is `True`, the entire scores. Shape: [(x-sliding size, y-sliding size, ntrac+1, len]

        Notes
        -----
        The shapes of `tid`, `x`, `y`, `vxg`, `vyg` must be the same (`ndim` can be > 1).
        This shape shall be expressed as sh in what follows.
        e.g. if the shape of the input arrays are [k,l], the shapes of
        the outputs shall be like, `count`: [k,l], `tid`: [ntrac+1,k,l]
        
        If `asxarray` is `True`, the return values are combined into a single `xr.Dataset`.
        """

        if np.isscalar(tid0):
            tid0 = np.full(x0.shape, tid0)
        
        # Julia is 1-based, so +1
        tid0 = np.asarray(tid0) + 1
        x0 = np.asarray(x0) + 1
        y0 = np.asarray(y0) + 1

        results = Main.VTTrac.trac(self.o, tid0, x0, y0, vxg=vxg0, vyg=vyg0, out_subimage=out_subimage, out_score_ary=out_score_ary, to_missing=False)
        count, status, tid, x, y, vx, vy, score, zss, score_ary = results

        # Python is 0-based, so -1
        tid -= 1
        x -= 1
        y -= 1

        # assign missing value
        count = np.array(count)
        status = np.array(status)
        tid = np.array(tid).astype(float)
        tid[tid==self.imiss] = np.nan
        x = np.array(x)
        x[x==self.fmiss] = np.nan
        y = np.array(y)
        y[y==self.fmiss] = np.nan
        vx = np.array(vx)
        vx[vx==self.fmiss] = np.nan
        vy = np.array(vy)
        vy[vy==self.fmiss] = np.nan
        score = np.array(score)
        score[score==self.fmiss] = np.nan
        if out_subimage:
            zss = np.array(zss)
            zss[zss==self.fmiss] = np.nan
        if out_score_ary:
            score_ary = np.array(score_ary)
            score_ary[score_ary==self.fmiss] = np.nan

        if asxarray:
            ds = self.to_xarray(count, status, tid, x, y, vx, vy, score, zss, score_ary, out_subimage=out_subimage, out_score_ary=out_score_ary)
            return ds
        else:
            return count, status, tid, x, y, vx, vy, score, zss, score_ary

    def to_xarray(self, count, status, tid, x, y, vx, vy, score, zss=None, score_ary=None, out_subimage=False, out_score_ary=False):
        sh = count.shape
        dims = [None]*len(sh)
        dimnames = [""]*len(sh)
        for i, length in enumerate(sh):
            dims[i] = np.arange(length)
            dimnames[i] = f"n{i}"
        ds = xr.Dataset(
            data_vars=dict(
                z = (["t", "y", "x"], self.o.z),
                count = (dimnames, count),
                status = (dimnames, status),
                tid = (["it_rel", *dimnames], tid),
                xloc = (["it_rel", *dimnames], x),
                yloc = (["it_rel", *dimnames], y),
                vx = (["it_rel_v", *dimnames], vx),
                vy = (["it_rel_v", *dimnames], vy),
                score = (["it_rel_v", *dimnames], score),
            ),
            coords={
                "t": (["t"], self.o.t),
                "it_rel": np.arange(self.ntrac+1)*self.itstep,
                "it_rel_v": np.arange(self.ntrac)*self.itstep + 0.5*np.sign(self.itstep)
            }
        )
        if out_subimage:
            ds = ds.assign_coords({"sx":np.arange(self.nsx), "sy": np.arange(self.nsy)})
            ds["zss"] = (["sx", "sy", "it_rel", *dimnames], zss)
        if out_score_ary:
            ds = ds.assign_coords({"scx":np.arange(self.ixhw*2 + 1), "scy": np.arange(self.iyhw*2 + 1)})
            ds["score_ary"] = (["scx", "scy", "it_rel_v", *dimnames], score_ary)
        
        ds.attrs = self.attrs
        return ds

    def set_grid_par(self, x0, y0, dx, dy, ucfact, ucufact):
        self.x0 = x0
        self.y0 = y0
        self.dx = abs(dx)
        self.dy = abs(dy)
        self.ucfact = ucfact
        self.ucufact = ucufact

    def setup_eq_grid(self, **kwargs):
        for key in ("vxhw", "vxch"):
            if kwargs[key]:
                if self.ucfact:
                    kwargs[key] = kwargs[key] / (self.dx*self.ucfact)
                else:
                    kwargs[key] = kwargs[key] / self.dx
        for key in ("vyhw", "vych"):
            if kwargs[key]:
                if self.ucfact:
                    kwargs[key] = kwargs[key] / (self.dy*self.ucfact)
                else:
                    kwargs[key] = kwargs[key] / self.dy
        self.setup(**kwargs)

    def trac_eq_grid(self, tid0, x, y, **opts):
        x = (x-self.x0)/self.dx
        y = (y-self.y0)/self.dy

        opts["asxarray"] = False
        count, status, tid, x, y, vx, vy, score, zss, score_ary = self.trac(tid0, x, y, **opts)

        x = x * self.dx + self.x0
        y = y * self.dy + self.y0
        if self.ucfact is None:
            vx = vx*self.dx
            vy = vy*self.dy
        else:
            vx = vx * (self.ucfact*self.dx)
            vy = vy * (self.ucfact*self.dy)
        ds = self.to_xarray(count, status, tid, x, y, vx, vy, score, zss, score_ary)
        return ds