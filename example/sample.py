#%%
import numpy as np
import matplotlib.pyplot as plt

from pyVTTrac import VTTrac

#%% Generate a synthetic dataset
nt = 20
ny = 80
nx = 100

t = np.arange(nt, dtype=np.float64)

tg = t[:,None,None]
yg = np.arange(ny)[None,:,None]
xg = np.arange(nx)[None,None,:]
k = 2 * np.pi / 12
cx = 1.2
cy = 1.2

z = np.sin(k*(xg-cx*tg))*np.cos(k*(yg-cy*tg))
#%% Initialize VTTrac object (first to use julia here)
vtt = VTTrac.VTT(z, t)
# %% Traking setup
nsx = 5
nsy = 5
ntrac = nt-1
vtt.setup(nsx, nsy, vxhw=1.8, vyhw=1.8, ntrac=ntrac)
# %% Make starting points and conduct tracking
n = 5
tid0 = np.zeros(n).astype(np.int64)
x0 = 2.5 * np.arange(n) + 7.5
y0 = 1.0 * np.arange(n) + 10.5

ds = vtt.trac(tid0, x0, y0, out_subimage=True, out_score_ary=True)
ds
# %% Plot the tracking result
it = 0
fig, ax = plt.subplots()
ax.imshow(ds.z[it], cmap="RdBu_r")
ax.scatter(ds.xloc[it], ds.yloc[it], c="k", ec="w", s=20)
ax.set(xlabel="x", ylabel="y", aspect="equal")
plt.show()
# %% Print the tracking result
print(f"cx = {cx}, averaged tracking result = {ds.vx.mean().item()}")
print(f"cy = {cy}, averaged tracking result = {ds.vy.mean().item()}")

# %%
