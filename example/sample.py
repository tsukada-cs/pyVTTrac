#%%
import numpy as np

from pyVTTrac import VTTrac
# %%
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
#%%
vtt = VTTrac.VTT(z, t)
# %%
nsx = 5
nsy = 5
ntrac = nt-1
vtt.setup(nsx, nsy, vxhw=1.8, vyhw=1.8, ntrac=ntrac)
# %%
n = 6
tid0 = np.zeros(n).astype(np.int64)
x0 = 2.5*np.arange(n) + 7.5
y0 = 1.0*np.arange(n) + 10.5

ds = vtt.trac(tid0, x0, y0, out_subimage=True, out_score_ary=True)

print(ds)