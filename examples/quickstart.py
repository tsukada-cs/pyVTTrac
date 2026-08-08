#%%
import matplotlib.pyplot as plt
import numpy as np

import pyvttrac as vt

#%% Generate a synthetic dataset: a traveling wave pattern
nt, ny, nx = 20, 80, 100
t = np.arange(nt, dtype=np.float64)

tg = t[:, None, None]
yg = np.arange(ny)[None, :, None]
xg = np.arange(nx)[None, None, :]
k = 2 * np.pi / 30
cx, cy = 1.2, 1.2

z = np.sin(k * (xg - cx * tg)) * np.cos(k * (yg - cy * tg))
z = z.astype(np.float32)

#%% Seed points and tracking
x0, y0 = vt.seed_grid(z.shape, spacing=(5.0, 2.5), margin=(10.5, 7.5))

res = vt.track(
    z, x0, y0, t0=0, time=t,
    template=(5, 5),
    search_velocity=(1.8, 1.8),
    nsteps=nt - 1,
    fixed_template=False,
    diagnostics=True,
)
res

#%% Plot the tracking result
it = 0
fig, ax = plt.subplots()
ax.imshow(z[it], cmap="RdBu_r", origin="lower")
ax.scatter(res.x[it], res.y[it], c="k", ec="w", s=20)
ax.set(xlabel="x", ylabel="y", aspect="equal", title=f"t index = {it}")
plt.show()

#%% Print the tracking result
print(f"cx = {cx}, mean tracked vx = {np.nanmean(res.vx):.4f}")
print(f"cy = {cy}, mean tracked vy = {np.nanmean(res.vy):.4f}")
print(f"tracked {res.ok.sum()} / {res.ok.size} points to completion")

# Prefer an xarray.Dataset or a tidy DataFrame? (optional dependencies)
# ds = res.to_xarray()
# df = res.to_dataframe()

# %%
