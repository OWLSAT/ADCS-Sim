"""
Contains the classes used to define and construct the OwlSat instance as well as its associated sensors and actuators
"""

import numpy as np
from Helpers import get_noise_V, rotate_angles, vec_unit
import random

class Satellite:
    def __init__(self, sensors, actuators, physics):
        """
        A class representing the satellite.
        sensors: a dict of instances of class Sensor() to add to the satellite
        actuators: a dict of instances of class Actuator() to add to the satellite
        physics: a dictionary of values specifying the physical properties of the satellite
                 (Satellite treated as a point mass inside a cube volume)
        """
        self.sensors = sensors
        self.sensor_data = {}
        for i in self.sensors:
            self.sensor_data[i] = 0
        self.actuators = actuators
        self.state = 0
        # State 0: Deploying
        # State 1: Initializing
        # State 2: Detumbling
        # State 3: Connecting
        # State 4: Mission/Stabilize
        # State 5: Error

        self.dimensions = physics["dimensions"]  # Length of cube side in m
        self.volume = self.dimensions ** 3  # Satellite volume in m^3
        self.mass = physics["mass"]  # Satellite mass in kg
        self.moi = physics["moi"]  # Satellite moment of inertia in kg*m^2

        self.power = {"solar": 2.7,  # Solar power generation in W
                      "battery": 0,  # Battery charge in A/hr
                      "usage": 0}  # Power usage in W

        self.time = 0  # Time in ms since poweron

        self.sim_state = {"pos":np.array([0,0,0]), # XYZ position coordinates
                          "vel":np.array([0,0,0]), # Velocity vector
                          "att":np.array([0,0,0]), # Attitude angles (wrt. earth-sun vector)
                          "rot":np.array([0,0,0])} # Angular velocity vector

        self.BDOT_K = 2
        self.BDOT_vals = []
        self.dt_val = 0

    def update_state(self, timestep, sim):
        self.time += timestep
        self.dt_val = sim.dt
        for i in self.sensors:
            self.sensor_data[i] = self.sensors[i].measure(sim)
        self.BDOT_vals.append(self.sensor_data["Magnetometer"])
        if len(self.BDOT_vals) > self.BDOT_K:
            self.BDOT_vals.pop(0)

    def compute_actuator_output(self):
        if self.state == 0:  # Detumbling
            ### Uses control scheme based on sensor_data to produce final net magnetic moment
            if len(self.BDOT_vals) == self.BDOT_K:
                voltages = self.actuators["Magnetorquer"].Bdot(self.BDOT_vals, self.dt_val)
            else:
                voltages = np.array([0,0,0])
            #voltages = self.actuators["Magnetorquer"].gyro_detumble(self.sensor_data["Gyroscope"], self.sensor_data["Magnetometer"])
        else:
            voltages = np.array([0, 0, 0])
        # Translate voltages into moments for simulation
        sat_mag_moment = self.actuators["Magnetorquer"].moment_from_voltage(voltages)
        return sat_mag_moment

class Sensor:
    """
    A class representing a generic sensor.
    """
    def __init__(self):
        self.noise = 0

    def measure(self, sim):
        return 0


# Subclasses (types) of Sensors below

class Clock(Sensor): # Subclass of sensor
    def measure(self, sim):
        return sim.time


class Accelerometer(Sensor): # Subclass of Sensor
    def __init__(self):
        # Noise: 180ug
        self.noise = 180*(10**-6)

    def measure(self, sim):
        position = sim.get_current_state()[0]
        mass = sim.sat.mass
        F = sim.env.compute_net_force(position, mass)
        a = F / mass
        return a + get_noise_V(self.noise)


class Gyroscope(Sensor): # Subclass of Sensor
    def __init__(self):
        # Noise: 0.008º/s
        self.noise = 0.008*np.pi/180

    def measure(self, sim):
        real_val = sim.get_current_state()[3]
        return real_val + get_noise_V(self.noise)


class Magnetometer(Sensor): # Subclass of Sensor
    def __init__(self):
        # Noise: 0.3uT
        #self.noise = 0.3*(10**-6)
        self.noise = 100*(10**-9)

    def measure(self, sim):
        rot = sim.get_current_state()[2]
        sat_yaw = rot[0]
        sat_pit = rot[1]
        sat_rol = rot[2]

        magfield = sim.read_mag_field()

        # Get magnetic field vector as seen by the magnetometer
        referenced_field = rotate_angles(magfield, -sat_yaw, -sat_pit, -sat_rol)

        return referenced_field + get_noise_V(self.noise)


class EUV(Sensor): # Subclass of Sensor
    def measure(self, sim):
        return 0


class GNSS_POS(Sensor): # Subclass of Sensor
    def __init__(self):
        # Noise: 8m
        self.noise = 8
        # Accuracy of position measurement is <8m
        #   - Max of 8 meters of deviation in total

    def measure(self, sim):
        real_val = sim.get_current_state()[0]
        length = random.random() * self.noise
        total_deviation = vec_unit(get_noise_V(1)) * length
        return real_val + total_deviation


class GNSS_VEL(Sensor): # Subclass of Sensor
    def __init__(self):
        # Noise: 0.1m/s
        self.noise = 0.1

    def measure(self, sim):
        real_val = sim.get_current_state()[1]
        return real_val + get_noise_V(self.noise)


class Magnetorquer:
    def moment_from_voltage(self, voltage):
        """
        Computes and returns the magnetic moment vector produced by the input voltages to the coils,
        in satellite frame of reference
        """
        moments = [0, 0, 0]
        # Equations for X - Y torque coils (0 at 0V - 0.3 at 3.3V)
        moments[0] = np.interp(voltage[0], [-3.3, 3.3], [-0.3, 0.3])
        moments[1] = np.interp(voltage[1], [-3.3, 3.3], [-0.3, 0.3])
        # Equation for Z torque coil (0 at 0V - 0.34 at 3.3V)
        moments[2] = np.interp(voltage[2], [-3.3, 3.3], [-0.34, 0.34])
        return moments

    def gyro_detumble(self, gyro, magfield):
        k = 0.0001
        moment = -k * np.cross(gyro, magfield) / np.dot(magfield, magfield)
        voltages = [0,0,0]
        max_moment = np.max(np.abs(moment))
        for i in range(3):
            voltages[i] = 3.3*(moment[i]/float(max_moment))
        return np.array(voltages)

    def Bdot(self, previousB, dt):
        # Separate into xyz values
        Bx = []
        By = []
        Bz = []
        for i in range(len(previousB)):
            Bx.append(previousB[i][0])
            By.append(previousB[i][1])
            Bz.append(previousB[i][2])

        # Compute avg derivative
        dBx = np.diff(Bx) / dt
        dBy = np.diff(By) / dt
        dBz = np.diff(Bz) / dt
        avg_dBx = np.sum(dBx) / len(dBx)
        avg_dBy = np.sum(dBy) / len(dBy)
        avg_dBz = np.sum(dBz) / len(dBz)
        avg_dBx = dBx[0]
        avg_dBy = dBy[0]
        avg_dBz = dBz[0]

        # Actuate coils depending on the sign of the derivatives
        V = 3.3
        T = 0.01 * (10 ** -6)
        voltages = [0, 0, 0]
        moment = [-avg_dBx, -avg_dBy, -avg_dBz]
        max_moment = np.max(np.abs(moment))
        for i in range(3):
            voltages[i] = 3.3 * (moment[i] / float(max_moment))

        for i in range(3):
            if moment[i] > 0:
                voltages[i] = 3.3
            if moment[i] < 0:
                voltages[i] = -3.3

        return np.array(voltages)
