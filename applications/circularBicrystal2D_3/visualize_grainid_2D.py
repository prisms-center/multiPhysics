#!/usr/bin/env python3
"""
Script to visualize 2D grain ID data as a pseudocolor plot.

Usage:
    python visualize_grainid_2D.py <nx> <ny> <filename>

Arguments:
    nx: Number of columns in the output array
    ny: Number of rows in the output array
    filename: Path to the text file containing numerical data
"""

import sys
import numpy as np
import matplotlib.pyplot as plt


def main():
    # Check command line arguments
    if len(sys.argv) != 4:
        print("Usage: python visualize_grainid_2D.py <nx> <ny> <filename>")
        sys.exit(1)
    
    try:
        nx = int(sys.argv[1])
        ny = int(sys.argv[2])
        filename = sys.argv[3]
    except ValueError:
        print("Error: nx and ny must be integers")
        sys.exit(1)
    
    # Read the file, skipping the first line
    try:
        data = np.loadtxt(filename, skiprows=1)
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found")
        sys.exit(1)
    except ValueError as e:
        print(f"Error reading file: {e}")
        sys.exit(1)
    
    # Reshape the data into nx x ny array
    expected_size = nx * ny
    if data.size != expected_size:
        print(f"Error: Expected {expected_size} values ({nx}x{ny}), but got {data.size}")
        sys.exit(1)
    
    data_array = data.reshape((ny, nx))
    
    # Create pseudocolor plot
    plt.figure(figsize=(10, 8))
    plt.pcolormesh(data_array, cmap='viridis', shading='auto')
    plt.colorbar(label='Grain ID')
    plt.xlabel('X (nx)')
    plt.ylabel('Y (ny)')
    plt.title('2D Grain ID Visualization')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()
