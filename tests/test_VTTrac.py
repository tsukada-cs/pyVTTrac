import pytest

import numpy as np
import xarray as xr

from pyVTTrac import VTTrac


class TestVTTrac:
    """Test the FastLineDetector class."""
    nt = 20
    ny = 80
    nx = 100
    tax = np.arange(nt, dtype=np.float64)
    tg = tax[:,None,None]
    yg = np.arange(ny)[None,:,None]
    xg = np.arange(nx)[None,None,:]
    k = 2 * np.pi / 10
    cx = 1.2
    cy = 1.2
    z = np.sin(k*(xg-cx*tg))*np.cos(k*(yg-cy*tg))
    vtt = VTTrac.VTT(z, tax)

    def test_init(self):
        assert self.vtt.attrs["nt"] == self.nt
        assert self.vtt.attrs["ny"] == self.ny
        assert self.vtt.attrs["nx"] == self.nx
    
    def test_setup(self):
        nsx = 5
        nsy = 5
        ntrac = self.nt-1
        self.vtt.setup(nsx, nsy, vxhw=1.8, vyhw=1.8, ntrac=ntrac)

        assert self.vtt.attrs["dtmean"] == self.vtt.attrs["itstep"]*(self.tax[-1]-self.tax[0])/(self.nt-1)
        assert self.vtt.attrs["nsx"] == nsx
        assert self.vtt.attrs["nsy"] == nsy
        assert self.vtt.attrs["vxhw"] == 1.8
        assert self.vtt.attrs["vyhw"] == 1.8
        assert self.vtt.attrs["ixhw"] == 3
        assert self.vtt.attrs["iyhw"] == 3
        assert self.vtt.attrs["ntrac"] == ntrac

    def test_trac(self):
        n = 6
        tid0 = np.zeros(n).astype(np.int64)
        x0 = 2.5*np.arange(n) + 7.5
        y0 = 1.0*np.arange(n) + 10.5

        ds = self.vtt.trac(tid0, x0, y0, out_subimage=True, out_score_ary=True, asxarray=True)
        assert type(ds) == xr.Dataset
        
        results = self.vtt.trac(tid0, x0, y0, out_subimage=True, out_score_ary=True, asxarray=False)
        assert type(results) == tuple