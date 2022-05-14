"""
Contains the code that runs the simulation for OwlSat ADCS systems
"""

from Environment import Environment
from Helpers import *


#### SIMULATION ###

class Simulation:
    def __init__(self, sat):
        self.time = 0
        self.timestep = 0.01  ## time step in s

        self.inclination = 0

        self.env = Environment()

        self.sat = sat
        # Simulation values
        self.pos = np.array([0, 0, 0])
        self.vel = np.array([0, 0, 0])
        self.att = [np.array([1, 0, 0]),
                    np.array([0, 1, 0])] # Instead of angles, use a pair of orthonormal vectors
        self.rot = np.array([0, 0, 0])

        self.sat_angles = np.array([0, 0, 0]) # Used to store the conversion of both vectors onto attitude angles

        self.last_torque = np.array([0, 0, 0])

        self.positions = []
        self.mag_fields = []
        self.current_pos = -1
        self.dt = 0
        self.points = 0

        # Big dictionary to store simulation data
        self.state_data = {"Xpos":[],
                           "Ypos":[],
                           "Zpos":[],
                           "Xvel":[],
                           "Yvel":[],
                           "Zvel":[],
                           "Yaw":[],
                           "Pitch":[],
                           "Roll":[],
                           "XrVel":[],
                           "YrVel":[],
                           "ZrVel":[],
                           "Xmag":[],
                           "Ymag":[],
                           "Zmag":[],
                           "Xtorque":[],
                           "Ytorque":[],
                           "Ztorque":[],
                           "Xsatmag":[],
                           "Ysatmag":[],
                           "Zsatmag":[]}

    def read_mag_field(self):
        return self.mag_fields[self.current_pos]

    def set_initial_params(self, pos, vel, rot, rvl):
        self.pos = pos
        self.vel = vel
        self.att = rot
        self.rot = rvl
        self.sat_angles = attitude_from_vectors(self.att[0],self.att[1])

    def get_current_state(self):
        return tuple([self.pos, self.vel, self.sat_angles, self.rot, self.att])

    def run_simulation(self, timespan):
        self.state_data = {"Xpos": [], "Ypos": [], "Zpos": [],
                           "Xvel": [], "Yvel": [], "Zvel": [],
                           "Yaw": [], "Pitch": [], "Roll": [],
                           "XrVel": [], "YrVel": [], "ZrVel": [],
                           "Xmag": [], "Ymag": [], "Zmag": [],
                           "Xtorque": [], "Ytorque": [], "Ztorque": [],
                           "Xsatmag": [], "Ysatmag": [], "Zsatmag": []}
        self.sat.update_state(self.dt, self)
        timestep_count = int(timespan // self.dt)
        tps = int(timestep_count // timespan)
        print("Total timesteps to simulate: " + str(timestep_count))
        for i in range(timestep_count):
            if i % tps == 0:
                self.iterate(collect=True)
            else:
                self.iterate(collect=False)
            if i % 1000 == 0:
                print("[" + str(round(i / float(timestep_count) * 100)) + "%]")
        return self.state_data

    def iterate(self, collect=True):
        # Update satellite position: ------------------------------------------------------------
        self.current_pos = (self.current_pos + 1) % self.points # Update rails position index
        self.pos = self.positions[self.current_pos] # Set position to next rails position
        if collect == True:
            self.state_data["Xpos"].append(self.pos[0])
            self.state_data["Ypos"].append(self.pos[1])
            self.state_data["Zpos"].append(self.pos[2])

        # Update satellite velocity vector: -----------------------------------------------------
        self.vel = np.array([0, 0, 0]) # Currently no need to do this since satellite is on rails
        if collect == True:
            self.state_data["Xvel"].append(self.vel[0])
            self.state_data["Yvel"].append(self.vel[1])
            self.state_data["Zvel"].append(self.vel[2])

        # Update satellite attitude: ------------------------------------------------------------
        self.att[0] = rotate_rotv(self.att[0], self.rot*self.dt)
        self.att[1] = rotate_rotv(self.att[1], self.rot*self.dt)
        self.sat_angles = attitude_from_vectors(self.att[0],self.att[1])
        if collect == True:
            self.state_data["Yaw"].append(self.sat_angles[0])
            self.state_data["Pitch"].append(self.sat_angles[1])
            self.state_data["Roll"].append(self.sat_angles[2])

        # Update satellite angular velocity: ----------------------------------------------------
        magfield = self.read_mag_field()
        sat_magfield = rotate_angles(magfield, -self.sat_angles[0], -self.sat_angles[1], -self.sat_angles[2])
        if collect == True:
            self.state_data["Xsatmag"].append(sat_magfield)
            self.state_data["Ysatmag"].append(sat_magfield)
            self.state_data["Zsatmag"].append(sat_magfield)
        sat_mag_moment = self.sat.compute_actuator_output()
        actual_moment = rotate_angles(sat_mag_moment, self.sat_angles[0], self.sat_angles[1], self.sat_angles[2])
        actual_torque = self.env.compute_sat_torque(actual_moment, self.mag_fields[self.current_pos])
        if collect == True:
            self.state_data["Xtorque"].append(actual_torque[0])
            self.state_data["Ytorque"].append(actual_torque[1])
            self.state_data["Ztorque"].append(actual_torque[2])
        self.rot = self.rot + (actual_torque / self.sat.moi)*self.dt
        if collect == True:
            self.state_data["XrVel"].append(self.rot[0])
            self.state_data["YrVel"].append(self.rot[1])
            self.state_data["ZrVel"].append(self.rot[2])

        # Update remaining state items
        self.sat.update_state(self.dt, self)
        self.time += self.dt
        self.env.update_time(self.time)

    def load_parameters(self, POSITION_ARRAY, MAGFIELD_ARRAY, dt_params, num_points):
        self.positions = POSITION_ARRAY
        self.mag_fields = MAGFIELD_ARRAY
        self.dt = dt_params
        self.points = num_points
