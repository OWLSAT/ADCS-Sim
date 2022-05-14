"""
Various helper functions for the OwlSat ADCS Simulation and a few constants too
"""

from math import cos, sin
import numpy as np
import math
import random
from scipy.spatial.transform import Rotation as Rot

import matplotlib.pyplot as plt

import pickle

pi = np.pi

YEAR_MS = 31556952000
SIDEREAL_DAY_MS = 86164091


def rotate_rotv(v, rot):
    """
    Rotates a vector v in 3D given a rotation vector
    """
    R = Rot.from_rotvec(rot)
    v_rot = R.apply(v)
    return v_rot

def rotate_angles(v, yaw, pitch, roll):
    """
    Rotates a vector v in 3D space given yaw pitch and roll angles
    """
    R = Rot.from_euler('xyz', [yaw, pitch, roll], degrees=False)
    v_rot = R.apply(v)
    return v_rot

def attitude_from_vectors(v1, v2):
    """
    Computes the attitude angles (yaw, pitch, roll) based on the pair of orthonormal vectors
    """
    R, error = Rot.align_vectors([v1,v2], [np.array([1,0,0]),np.array([0,1,0])])
    angles = R.as_euler('xyz', degrees=False)
    return angles

def xyz_to_rtp(x, y, z):
    """
    Converts cartesian coordinates into our custom spherical coordinates.
    x, y, and z are cartesian position in meters
    """
    h = math.sqrt(x ** 2 + y ** 2)
    r = math.sqrt(h ** 2 + z ** 2)
    theta = math.atan2(float(y), float(x))
    phi = math.atan2(float(z), float(h))
    return r, theta, phi


def rtp_to_xyz(r, theta, phi):
    """
    Converts our custom spherical coordinatesinto cartesian coordinates.
    """
    x = r * cos(theta) * cos(phi)
    y = r * sin(theta) * cos(phi)
    z = r * sin(phi)
    return x, y, z


def neu_to_xyz(x, y, z, n, e, u):
    """
    Converts magnetic field vector from north, east, down to x, y, z
    (units of magnetic flux density are nanoteslas)
    """
    r, theta, phi = xyz_to_rtp(x, y, z)

    conversion_matrix = np.array([[-sin(phi) * cos(theta), -sin(theta), x / r],
                                  [-sin(phi) * sin(theta), cos(theta), y / r],
                                  [cos(phi), 0, z / r]])
    matrix_thing = np.array([[n, e, u]]).transpose()
    xyz = np.matmul(conversion_matrix, matrix_thing)
    return xyz[0, 0], xyz[1, 0], xyz[2, 0]


def vec_projection(u, v):
    """
    u, v are numpy vectors
    Projects u ONTO v
    """
    proj = (np.dot(u, v) / np.linalg.norm(v) ** 2) * v
    return proj


def vec_unit(v):
    return v / np.linalg.norm(v)

def get_noise_V(noise):
    """
    Returns a random vector of noise +-noise
    """
    x_dev = ((2 * random.random()) - 1) * noise
    y_dev = ((2 * random.random()) - 1) * noise
    z_dev = ((2 * random.random()) - 1) * noise
    dev_V = np.array([x_dev, y_dev, z_dev])
    return dev_V

def angle(v1, v2, acute=True):
    """
    Gives angle between two vectors
    """
    # v1 is your firsr vector
    # v2 is your second vector
    norm_value = (np.linalg.norm(v1) * np.linalg.norm(v2))
    if norm_value == 0:
        return 0
    theta = np.arccos(np.dot(v1, v2) / norm_value)
    if (acute == True):
        return theta
    else:
        return 2 * np.pi - theta

### DATA PLOTTING ###
def linegraph(data, xlbl, ylbl, xmin=0, xmax=0, ymin=0, ymax=0, size=(5,3)):
    f = plt.figure(figsize=size, dpi=100)
    plt.plot(data)
    plt.xlabel(xlbl)
    plt.ylabel(ylbl)
    if xmin != xmax:
        if ymin != ymax:
            plt.axis([xmin,xmax,ymin,ymax])
    plt.show()

def multiLinegraph(dataSet, xlbl, ylbl, xmin=0, xmax=0, ymin=0, ymax=0, size=(5,3)):
    f = plt.figure(figsize=size, dpi=100)
    for data in dataSet:
        plt.plot(data)
    plt.xlabel(xlbl)
    plt.ylabel(ylbl)
    if xmin != xmax:
        if ymin != ymax:
            plt.axis([xmin,xmax,ymin,ymax])
    plt.show()


### Provide functions to save and load the computed arrays to a file on disk

def save_arrays(pos_array, mag_array, dt_val, p_val):
    # Group the data into a single object:
    DATA_PACKAGE = tuple([pos_array, mag_array, dt_val, p_val])
    # Pickle it into a file
    with open("orbitArrays.txt", 'wb') as f:
        pickle.dump(DATA_PACKAGE, f)
    f.close()

def load_arrays():
    # Un-Pickle data into object
    with open("orbitArrays.txt", 'rb') as f:
        DATA_PACKAGE = pickle.load(f)
    f.close()
    pos_array = DATA_PACKAGE[0]
    mag_array = DATA_PACKAGE[1]
    dt_val = DATA_PACKAGE[2]
    p_val = DATA_PACKAGE[3]
    return pos_array, mag_array, dt_val, p_val