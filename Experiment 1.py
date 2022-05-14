### EXPERIMENT 1 ###

from Simulation import Simulation
from OwlSat import owlsat
import numpy as np

sim = Simulation(owlsat)

sim_length_s = 5540
data_density = 100

earth_radius =  6.3781E6 #m

pos = np.array([420000+earth_radius,0,0])
vel = np.array([0,7657.419503231067,0])
rot = np.array([0,0,0])
rvl = np.array([0,0,0])

sim.set_initial_params(pos, vel, rot, rvl)

sim.run_simulation(sim_length_s, data_density, mode=0)
