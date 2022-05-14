"""
Contains the class for the environment that the satellite is simulated in
"""

from Helpers import *
import numpy as np
import ppigrf
from datetime import datetime
import math


class Environment:
    def __init__(self):
        self.time = 0  # in ms # summer solstice when sun is over prime meridian
        self.G = 6.67430E-11  # m^3 kg^-1 s^-2
        self.earth_radius = 6.3781E6  # m
        self.earth_rotation = ((((2 * pi) / SIDEREAL_DAY_MS) * self.time) % 2 * pi)  # an angle, in radians
        self.earth_tilt = 0.4084  # rad
        self.earth_mass = 5.97237E24  # kg
        self.earth_center = np.array([0, 0, 0])
        self.earth_dipole_moment = 8.22E22  # A * m^2

    def update_time(self, new_time):
        self.time = new_time

    def longitude(self, t, theta):
        """
        Returns longitude at time t and angle theta
        """
        return (theta + pi - (((2 * pi) / SIDEREAL_DAY_MS) * t)) % 2 * pi - pi

    def compute_sat_torque(self, m, B):
        T = np.cross(m, B)
        # print("Torque = ", T)
        return np.array(T)

    def compute_net_force(self, sat_pos_xyz, sat_mass):
        h = np.linalg.norm(sat_pos_xyz)
        F = self.gravity(h) * sat_mass
        return F * -1 * vec_unit(sat_pos_xyz)

    def gravity(self, h):
        # Returns acceleration due to gravity G * m1 * m2 / r**2
        return (self.G * self.earth_mass) / h ** 2

    def atmo_density(self, pos_xyz):
        return 0

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

    def surface_velocity(self, pos_xyz, vel_xyz):
        """
        pos_xyz: a numpy array representing the position of the satellite
        vel_xyz: a numpy array representing the velocity of the satellite

        Returns: surface_vel: a numpy vector representing the velocity of the satellite on the surface of the earth.
        """
        equatorial_rotational_vel = self.earth_radius * self.earth_rotation * np.array(
            [-pos_xyz[1], pos_xyz[0]]) * (1 / math.sqrt((pos_xyz[0]) ** 2 + (pos_xyz[1]) ** 2))
        surface_vel = (np.linalg.norm(pos_xyz) / self.earth_radius) * vec_projection(vel_xyz,
                                                                                     pos_xyz) - equatorial_rotational_vel * math.sqrt(
            pos_xyz[0] ** 2 + pos_xyz[1] ** 2) / math.sqrt(pos_xyz[0] ** 2 + pos_xyz[1] ** 2 + pos_xyz[2] ** 2)

        return surface_vel

    def is_sun_visible(self, t, pos_xyz):
        """
        t: time
        pos_xyz: a numpy vector representing the cartesian coordinates of the satellite at time t
        """
        # Returns true if sun is visible from position pos at time t
        sun_theta = sin((2 * pi * t) / YEAR_MS)
        sun_phi = self.earth_tilt * cos((2 * pi * t) / YEAR_MS)
        sun_r = 1
        sun_x, sun_y, sun_z = rtp_to_xyz(sun_r, sun_theta, sun_phi)

        n = np.array([-sun_x, -sun_y, -sun_z])

        cos_angle = (np.dot(n, pos_xyz) / np.linalg.norm(n)) / np.linalg.norm(pos_xyz)

        if cos_angle > 0:
            return False

        a = self.earth_center  # Center of Earth
        p = pos
        d = np.linalg.norm((p - a) - n * np.dot((p - a), n))

        if d <= self.earth_radius:
            return False
        else:
            return True