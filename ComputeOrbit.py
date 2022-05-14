"""
Code to precompute the orbital position values and the corresponding magnetic field values at those positions
and store everything in a file to load into the simulation
"""

from Helpers import *
import numpy as np
from datetime import datetime
import ppigrf

class PseudoEnvironment():
    def __init__(self):
        self.time = 0  # in ms # summer solstice when sun is over prime meridian
        self.G = 6.67430E-11  # m^3 kg^-1 s^-2
        self.earth_radius = 6.3781E6  # m
        self.earth_rotation = ((((2 * pi) / SIDEREAL_DAY_MS) * self.time) % 2 * pi)  # an angle, in radians
        self.earth_tilt = 0.4084  # rad
        self.earth_mass = 5.97237E24  # kg
        self.earth_center = np.array([0, 0, 0])
        self.earth_dipole_moment = 8.22E22  # A * m^2

    def longitude(self, t, theta):
        """
        Returns longitude at time t and angle theta
        """
        return (theta + pi - (((2 * pi) / SIDEREAL_DAY_MS) * t)) % 2 * pi - pi

    def mag_field(self, pos_xyz):
        """
        Arguments:
            pos_xyz: a numpy array representing the cartesian coordinates of the satellite's position

        Returns:
            Cartesian magnetic flux density vector B representing the magnetic field at the point pos_xyz
        """
        x = pos_xyz[0]
        y = pos_xyz[1]
        z = pos_xyz[2]

        r, t, p = xyz_to_rtp(x, y, z)
        h = (r - self.earth_radius) / 1000.0  # Its supposed to be in km avobe sea-level

        date = datetime(2021, 3, 28)

        Be, Bn, Bu = ppigrf.igrf(self.longitude(self.time, t), p, h, date)  # returns east, north, up

        # Convert nT to T
        Be = Be[0, 0] / (10 ** 9)
        Bn = Bn[0, 0] / (10 ** 9)
        Bu = Bu[0, 0] / (10 ** 9)

        # Convert to xyz frame of reference
        Bx, By, Bz = neu_to_xyz(x, y, z, Bn, Be, Bu)

        return np.array([Bx, By, Bz])


env = PseudoEnvironment()


def precompute_orbit_params(p, h, i):
    """
    p is the number of simulation points in a single orbit
    h is the altitude (above sea level)
    i is the inclination (degrees from prograde-equatorial)
    """
    positions = []
    mag_fields = []
    p = p
    r = h + env.earth_radius
    M = env.earth_mass
    G = env.G
    dpos = 2 * np.pi / p
    v = np.sqrt(G * M / r)
    T = 2 * np.pi * r / v
    dt = T / p
    t = 0
    last_percent = -1
    for k in range(p):
        x = r * cos(t)
        y = r * sin(t) * cos(i * np.pi / 180)
        z = r * sin(t) * sin(i * np.pi / 180)
        t += dpos
        positions.append(np.array([x, y, z]))
        mag_fields.append(env.mag_field(np.array([x, y, z])))
        percent = round(k/p*100)
        if percent % 1 == 0:
            if percent != last_percent:
                print("Computing... [" + str(percent) + "%]")
                last_percent = percent
    print("Timestep length: " + str(dt))
    return positions, mag_fields, dt, p


### COMPUTE ORBIT ARRAYS
print("=====================================================================================")
print("Computing orbital array values...")
print("-------------------------------------------------------------------------------------")
POSITION_ARRAY, MAGFIELD_ARRAY, Gdt, Gp = precompute_orbit_params(50000, 565000, 97.5)
print("-------------------------------------------------------------------------------------")
print("Computation finished!")
print("=====================================================================================")
print("Saving to file...")
save_arrays(POSITION_ARRAY, MAGFIELD_ARRAY, Gdt, Gp)
print("Saved!")
print("=====================================================================================")

