import numpy as np
import matplotlib.pyplot as plt

# The grain ID values to use, as strings (must be strings! use quotes)
inner_grain_ID = "1"
outer_grain_ID = "2"

# Grid domain, in voxels
# Lui 2021 has domain of 120 um, so for 32 voxels, one voxel is 120/32 = 3.72 um
size_x = 1024
size_y = 1024
size_z = 4

# Circular grain diameter, in voxels
circle_radius = 384    # 90 um

# Name for the file to create
# Warning! The file will be overwritten if it already exists!
output_filename = "grainID_circle_{:d}_radius.txt".format(circle_radius)

# Values for the calculation (do not edit below this point)
circle_radius_squared = circle_radius**2
center_x = size_x / 2
center_y = size_y / 2
center_z = size_z / 2

grain_ids = np.full((size_x, size_y, size_z), int(outer_grain_ID), dtype=int)

for x in np.arange(size_x):
  for y in np.arange(size_y):
    dist_squared = (x - center_x)**2 + (y - center_y)**2
    if dist_squared < circle_radius_squared:
      grain_ids[x, y, :] = int(inner_grain_ID)

with open(output_filename, "w") as f:
  # Write header
  f.write("Grain ID file, circle of radius {:d} voxels\n".format(circle_radius))
  for x in np.arange(size_x):
    for y in np.arange(size_y):
      row = [str(grain_ids[x, y, z]) for z in np.arange(size_z)]
      f.write(" ".join(row) + "\n")

plot_filename_base = output_filename.rsplit(".", 1)[0]

for z in np.arange(size_z):
  fig, ax = plt.subplots()
  im = ax.imshow(grain_ids[:, :, z].T, origin="lower", cmap="viridis")
  ax.set_title("Grain IDs at z = {}".format(z))
  ax.set_xlabel("x")
  ax.set_ylabel("y")
  fig.colorbar(im, ax=ax, label="Grain ID")
  fig.tight_layout()
  fig.savefig("{}_z{:03d}.png".format(plot_filename_base, z), dpi=300)
  plt.close(fig)
