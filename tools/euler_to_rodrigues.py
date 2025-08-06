import numpy as np
import argparse
from scipy.spatial.transform import Rotation

def euler_to_rodrigues(phi1, Phi, phi2, degrees=True):
    """
    Convert Bunge Euler angles (ZXZ convention) to Rodrigues vector.
    
    Parameters:
    phi1 : float - First rotation (about Z-axis)
    Phi  : float - Second rotation (about X-axis)
    phi2 : float - Third rotation (about Z-axis)
    degrees : bool - If True, angles are given in degrees (default: True)
    
    Returns:
    rodrigues_vector : np.ndarray - Rodrigues vector (3 elements)
    """
    
    # Convert Euler angles to rotation matrix using SciPy
    r = Rotation.from_euler('ZXZ', [phi1, Phi, phi2], degrees=degrees)
    R = r.as_matrix()

    # Compute rotation angle theta
    theta = np.arccos((np.trace(R) - 1) / 2)

    if np.isclose(theta, 0):
        return np.array([0.0, 0.0, 0.0])  # No rotation case

    # Compute rotation axis
    rx = (R[2,1] - R[1,2]) / (2 * np.sin(theta))
    ry = (R[0,2] - R[2,0]) / (2 * np.sin(theta))
    rz = (R[1,0] - R[0,1]) / (2 * np.sin(theta))
    r = np.array([rx, ry, rz])

    # Compute Rodrigues vector
    rodrigues_vector = r * np.tan(theta / 2)

    return rodrigues_vector

def main():
    parser = argparse.ArgumentParser(description="Convert Bunge Euler angles to Rodrigues vector.")
    parser.add_argument("phi1", type=float, help="First Euler angle (phi1)")
    parser.add_argument("Phi", type=float, help="Second Euler angle (Phi)")
    parser.add_argument("phi2", type=float, help="Third Euler angle (phi2)")
    parser.add_argument("--radians", action="store_true", help="Set this flag if input angles are in radians")

    args = parser.parse_args()

    rod_vector = euler_to_rodrigues(args.phi1, args.Phi, args.phi2, degrees=not args.radians)
    print("Rodrigues vector:", rod_vector)

if __name__ == "__main__":
    main()
