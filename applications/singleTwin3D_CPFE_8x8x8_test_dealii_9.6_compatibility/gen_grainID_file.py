def generate_grainID_file(filename, M, N, P):
    with open(filename, 'w') as f:
        # Write the first row
        f.write("xxxxxxxxxxxxxxx\n")
        
        # Write the MxN rows with P ones separated by spaces
        row = " ".join(["1"] * P) + "\n"
        for _ in range(M * N):
            f.write(row)

# Set the dimensions
M, N, P = 8, 8, 8  # Change as needed

# Generate the file
generate_grainID_file("grainID_single_8x8x8.txt", M, N, P)
