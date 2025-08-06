import numpy as np
import argparse
from scipy.spatial.transform import Rotation

def rodrigues_to_euler(rodrigues_vector, degrees=True):
    """
    Convert Rodrigues vector to Bunge Euler angles (ZXZ convention).

    Parameters:
    rodrigues_vector : np.ndarray - Rodrigues vector (3 elements)
    degrees : bool - If True, output angles in degrees (default: True)

    Returns:
    euler_angles : np.ndarray - Bunge Euler angles (phi1, Phi, phi2)
    """
    # Compute rotation angle theta
    norm_r = np.linalg.norm(rodrigues_vector)
    if np.isclose(norm_r, 0):
        return np.array([0.0, 0.0, 0.0])  # No rotation case

    theta = 2 * np.arctan(norm_r)

    # Compute rotation axis
    r = rodrigues_vector / norm_r

    # Convert rotation vector to rotation matrix
    r_matrix = Rotation.from_rotvec(r * theta).as_matrix()

    # Convert rotation matrix to Euler angles (ZXZ convention)
    euler_angles = Rotation.from_matrix(r_matrix).as_euler('ZXZ', degrees=degrees)

    return euler_angles

def main():
    parser = argparse.ArgumentParser(description="Convert Rodrigues vector to Bunge Euler angles.")
    parser.add_argument("rx", type=float, help="Rodrigues vector component x")
    parser.add_argument("ry", type=float, help="Rodrigues vector component y")
    parser.add_argument("rz", type=float, help="Rodrigues vector component z")
    parser.add_argument("--radians", action="store_true", help="Set this flag if you want Euler angles in radians")

    args = parser.parse_args()

    rod_vector = np.array([args.rx, args.ry, args.rz])
    euler_angles = rodrigues_to_euler(rod_vector, degrees=not args.radians)
    print("Euler angles (Bunge ZXZ):", euler_angles)

if __name__ == "__main__":
    main()
