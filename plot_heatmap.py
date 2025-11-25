# plot_heatmap.py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# --- config ---
csv_path = "heatmap.csv"
output_png = "heatmap.png"
cmap = "inferno"     # choose "viridis", "plasma", "inferno", etc.
points_interp = 300  # resolution of output grid (300 -> higher quality)

# optional: polygon boundary you used (feet). Edit if you want overlay
poly = [
    (56.509, 6.20),
    (156.641, 95.380),
    (125.857, 170.050),
    (56.931, 162.723)
]

# --- load ---
df = pd.read_csv(csv_path)
xs = df['x_ft'].values
ys = df['y_ft'].values
zs = df['lux'].values

# make a regular grid for plotting
xi = np.linspace(xs.min(), xs.max(), points_interp)
yi = np.linspace(ys.min(), ys.max(), points_interp)
xi_grid, yi_grid = np.meshgrid(xi, yi)

# interpolate scattered xy->z onto grid; method='linear' is good; 'cubic' smoother but slower
zi = griddata((xs, ys), zs, (xi_grid, yi_grid), method='linear')

plt.figure(figsize=(9,8))
# show interpolated heatmap; use origin='lower' so y increases up
plt.imshow(zi, extent=(xi.min(), xi.max(), yi.min(), yi.max()), origin='lower', aspect='equal', cmap=cmap)
plt.colorbar(label='Illuminance (lux)')
# overlay polygon
poly_x = [p[0] for p in poly] + [poly[0][0]]
poly_y = [p[1] for p in poly] + [poly[0][1]]
plt.plot(poly_x, poly_y, color='white', linewidth=1.5)
plt.scatter(xs, ys, s=6, c='white', alpha=0.6)   # raw sample points (optional)
plt.title("Illuminance heatmap (lux) — interpolated")
plt.xlabel("x (ft)"); plt.ylabel("y (ft)")
plt.tight_layout()
plt.savefig(output_png, dpi=200)
plt.show()
print("Saved", output_png)

