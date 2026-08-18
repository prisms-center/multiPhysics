# The grain ID values to use, as strings (must be string literals! use quotes!)
inner_grain_ID = "1"
outer_grain_ID = "2"

# Grid domain size, in voxels
size_x = 240
size_y = 120
size_z = 2

# Name for the file to create
# Warning! The file will be overwritten if it already exists!!
output_filename = "grainID.txt"



# Values for the calculation (do not edit below this point)
center_x = size_x / 2

with open(output_filename, "w") as f:
  # Write header
  f.write("Grain ID file, planar bicrystal for twin transmission\n")
  for x in range(size_x):
    for y in range(size_y):
      row = []
      for z in range(size_z):
        if x < center_x:
          row.append(inner_grain_ID)
        else:
          row.append(outer_grain_ID)
      f.write(" ".join(row) + "\n")

print("Grain ID file written to {}".format(output_filename))

