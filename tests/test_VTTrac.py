import unittest
from unittest import mock

import numpy as np
import xarray as xr

from pyVTTrac import VTTrac

class TestVTTrac(unittest.TestCase):
    """Test the FastLineDetector class."""
    def test_VTT(self):
        nt = 20
        ny = 80
        nx = 100
        t = np.arange(nt, dtype=np.float64)
        tg = t[:,None,None]
        yg = np.arange(ny)[None,:,None]
        xg = np.arange(nx)[None,None,:]
        k = 2 * np.pi / 10
        cx = 1.2
        cy = 1.2
        z = np.sin(k*(xg-cx*tg))*np.cos(k*(yg-cy*tg))
        
        self.vtt = VTTrac.VTT(z, t)

        self.assertEqual(self.vtt.attrs["nt"], nt)
        self.assertEqual(self.vtt.attrs["ny"], ny)
        self.assertEqual(self.vtt.attrs["nx"], nx)
        

        nsx = 5
        nsy = 5
        ntrac = nt-1
        self.vtt.setup(nsx, nsy, vxhw=1.8, vyhw=1.8, ntrac=ntrac)

        self.assertEqual(self.vtt.attrs["dtmean"], self.vtt.attrs["itstep"]*(t[-1]-t[0])/(nt-1))
        self.assertEqual(self.vtt.attrs["nsx"], nsx)
        self.assertEqual(self.vtt.attrs["nsy"], nsy)
        self.assertEqual(self.vtt.attrs["vxhw"], 1.8)
        self.assertEqual(self.vtt.attrs["vyhw"], 1.8)
        self.assertEqual(self.vtt.attrs["ixhw"], 3)
        self.assertEqual(self.vtt.attrs["iyhw"], 3)
        self.assertEqual(self.vtt.attrs["ntrac"], ntrac)


        n = 6
        tid0 = np.zeros(n).astype(np.int64)
        x0 = 2.5*np.arange(n) + 7.5
        y0 = 1.0*np.arange(n) + 10.5

        ds = self.vtt.trac(tid0, x0, y0, out_subimage=True, out_score_ary=True, asxarray=True)
        self.assertTrue(type(ds) == xr.Dataset)
        
        results = self.vtt.trac(tid0, x0, y0, out_subimage=True, out_score_ary=True, asxarray=False)
        self.assertTrue(type(results) == tuple)