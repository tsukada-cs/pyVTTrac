import math

import numpy as np
import xarray as xr

import juliacall
from pathlib import Path
path = Path(__file__).parent.parent / "VTTrac.jl/src/VTTrac.jl"
jl = juliacall.newmodule("pyVTTrac")
jl.include(str(path))
jl.eval("using VTTrac")


class VTT:
    def __init__(self, z, t=None, mask=None, zmiss=None, fmiss=-999.0, imiss=-999):
        """
        Object initialization (VTT)

        You need to call #setup subsequently before calling #trac

        Parameters
        ----------
        z: array-like
            The "images" used for tracking.
        t: array-like or None
            Time for each image (if None, [0,1,2,..]).
        mask: array-like or None
            The "masks" used for masking. True positions are ignored when calculate score.
        zmiss: float
            Missing value in z (None if there is no missing).
        fmiss: float
            Missing float value.
        imiss: int
            Missing int value.
        """
        z = np.array(z, np.float32)
        self.o = jl.VTTrac.VTT(z, t=t, mask=mask, zmiss=zmiss, fmiss=fmiss, imiss=imiss)
        
        self.attrs = {}
        self.attrs["nt"] = self.o.nt
        self.attrs["ny"] = self.o.ny
        self.attrs["nx"] = self.o.nx
        self.attrs["zmiss"] = self.o.zmiss
        self.attrs["fmiss"] = self.o.fmiss
        self.attrs["imiss"] = self.o.imiss
        self.attrs["chk_zmiss"] = self.o.chk_zmiss
        self.attrs["chk_mask"] = self.o.chk_mask

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
    def chk_zmiss(self):
        return self.attrs["chk_zmiss"]
    @property
    def chk_mask(self):
        return self.attrs["chk_mask"]
    
    def __setitem__(self, key, value):
        self.attrs[key] = value
        exec(f"jl.VTTrac.set_{key}(self.o, value)")
    def __getitem__(self, key):
        return self.attrs[key]

    def setup(self, nsx, nsy, vxhw=None, vyhw=None, ixhw=None, iyhw=None, subgrid=True, subgrid_gaus=False,
        itstep=1, ntrac=2, score_method="xcor", Sth0=0.8, Sth1=0.7, vxch=None, vych=None, 
        peak_inside_th=None, Cth=None, use_init_temp=False, min_samples=1):
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
        Sth0: float, default 0.8
            The minimum score required for the 1st tracking.
        Sth1: float, default 0.7
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
        Cth: float, optional
            If non-`nothing`, an initial template is used only when 
            it has a difference in max and min greater than its value.
        min_samples: int, default 1
            Minimum number of visible values to calculate score when `chk_mask` is True.
        """

        jl.VTTrac.setup(
                self.o, nsx, nsy,
                vxhw=vxhw, vyhw=vyhw,
                ixhw=ixhw, iyhw=iyhw, subgrid=subgrid, subgrid_gaus=subgrid_gaus, itstep=itstep, ntrac=ntrac, score_method=score_method,
                Sth0=Sth0, Sth1=Sth1, vxch=vxch, vych=vych, peak_inside_th=peak_inside_th,
                Cth=Cth, use_init_temp=use_init_temp, min_samples=min_samples
            )
        
        self.attrs["dtmean"] = self.o.dtmean
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
        self.attrs["Sth0"] = self.o.Sth0
        self.attrs["Sth1"] = self.o.Sth1
        self.attrs["vxch"] = self.o.vxch
        self.attrs["vych"] = self.o.vych
        self.attrs["peak_inside_th"] = self.o.peak_inside_th
        self.attrs["Cth"] = self.o.Cth
        self.attrs["use_init_temp"] = self.o.use_init_temp
        self.attrs["min_samples"] = self.o.min_samples
    
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
    def Sth0(self):
        return self.attrs["Sth0"]
    @property
    def Sth1(self):
        return self.attrs["Sth1"]
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
    def Cth(self):
        return self.attrs["Cth"]
    @property
    def use_init_temp(self):
        return self.attrs["use_init_temp"]
    @property
    def min_samples(self):
        return self.attrs["min_samples"]
    
    def calc_ixyhw_from_v(self, vxhw, vyhw, dt):
        ixhw = math.ceil(abs(vxhw * dt)) + 1 # max displacement
        iyhw = math.ceil(abs(vyhw * dt)) + 1 # +1 is margin to find peak
        return ixhw, iyhw

    def calc_ixyhw_from_v_eq_grid(self, vxhw, vyhw, dt):
        if self.ucfact:
            vxhw = vxhw / (self.dx*self.ucfact)
        else:
            vxhw = vxhw / self.dx
        if self.ucfact:
            vyhw = vyhw / (self.dy*self.ucfact)
        else:
            vyhw = vyhw / self.dy
        return self.calc_ixyhw_from_v(vxhw, vyhw, dt)
    
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
        status: np.ndarray
            Tracking status. Shape: [len]
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

        results = jl.VTTrac.trac(self.o, tid0, x0, y0, vxg=vxg0, vyg=vyg0, out_subimage=out_subimage, out_score_ary=out_score_ary, to_missing=False)
        count, status, tid, x, y, vx, vy, score, zss, score_ary = results

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

        # Python is 0-based, so -1
        tid -= 1
        x -= 1
        y -= 1

        if asxarray:
            ds = self.to_xarray(count, status, tid, x, y, vx, vy, score, zss, score_ary)
            return ds
        else:
            return count, status, tid, x, y, vx, vy, score, zss, score_ary

    def to_xarray(self, count, status, tid, x, y, vx, vy, score, zss=None, score_ary=None):
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
        if zss is not None:
            ds = ds.assign_coords({"sy": np.arange(self.nsy), "sx":np.arange(self.nsx)})
            ds["zss"] = (["it_rel", "sy", "sx", *dimnames], zss)
        if score_ary is not None:
            ds = ds.assign_coords({"scy": np.arange(self.iyhw*2 + 1), "scx":np.arange(self.ixhw*2 + 1)})
            ds["score_ary"] = (["it_rel_v", "scy", "scx", *dimnames], score_ary)
        
        ds.attrs = self.attrs
        return ds

    def set_grid_par(self, x0, y0, dx, dy, ucfact=1, ucufact=None):
        """
        Setup grid parameter

        Parameters
        ----------
        x0: int or float
            Minimum x
        y0: int or float
            Minimum y
        dx: int or float
            Delta x
        dy: int or float
            Delta y
        ucfact: int or float, optional
            unit change factor
        ucfact: int or float, optional
            unit change unit factor
        """
        self.x0 = x0
        self.y0 = y0
        self.dx = abs(dx)
        self.dy = abs(dy)
        self.ucfact = ucfact
        self.ucufact = ucufact

    def setup_eq_grid(self, **kwargs):
        """
        Setup for tracking on eq grid.

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
        Sth0: float, default 0.8
            The minimum score required for the 1st tracking.
        Sth1: float, default 0.7
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
        Cth: float, optional
            If non-`nothing`, an initial template is used only when 
            it has a difference in max and min greater than its value.
        """
        for key in ("vxhw", "vxch"):
            if key in kwargs:
                kwargs[key] = kwargs[key] / (self.dx*self.ucfact)
        for key in ("vyhw", "vych"):
            if key in kwargs:
                kwargs[key] = kwargs[key] / (self.dy*self.ucfact)
        self.setup(**kwargs)

    def trac_eq_grid(self, tid0, x, y, **opts):
        """
        Conduct tracking on eq grid.

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
        status: np.ndarray
            Tracking status. Shape: [len]
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
        x = (x-self.x0)/self.dx
        y = (y-self.y0)/self.dy

        opts["asxarray"] = False
        count, status, tid, x, y, vx, vy, score, zss, score_ary = self.trac(tid0, x, y, **opts)

        x = x * self.dx + self.x0
        y = y * self.dy + self.y0
        if self.ucfact is None:
            vx = vx * self.dx
            vy = vy * self.dy
        else:
            vx = vx * (self.ucfact*self.dx)
            vy = vy * (self.ucfact*self.dy)

        ds = self.to_xarray(count, status, tid, x, y, vx, vy, score, zss, score_ary)
        return ds