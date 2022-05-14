"""
Contains the construction for the OwlSat Cube Satellite
"""

from SatBuilder import *

### OwlSat Satellite Instance Declaration ###

OWLSAT_EDGE_DIMENSION = 0.1 # m
OWLSAT_MASS = 1.387

owlsat_sensors = {"Accelerometer": Accelerometer(),
                  "Gyroscope": Gyroscope(),
                  "Magnetometer": Magnetometer(),
                  "EUV_x+": EUV(),
                  "EUV_x-": EUV(),
                  "EUV_y+": EUV(),
                  "EUV_y-": EUV(),
                  "EUV_z+": EUV(),
                  "EUV_z-": EUV(),
                  "GPS_pos": GNSS_POS(),
                  "GPS_vel": GNSS_VEL()}
owlsat_actuators = {"Magnetorquer": Magnetorquer()}

owlsat_physics = {"dimensions": OWLSAT_EDGE_DIMENSION,
                  "mass": OWLSAT_MASS,
                  "moi": 1.0/6.0 * OWLSAT_MASS}

owlsat = Satellite(owlsat_sensors, owlsat_actuators, owlsat_physics)
