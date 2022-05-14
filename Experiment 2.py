### EXPERIMENT 2 ###

from Simulation import Simulation
from OwlSat import owlsat
import numpy as np
from Helpers import linegraph, multiLinegraph, load_arrays

sim = Simulation(owlsat)

POSITION_ARRAY, MAGFIELD_ARRAY, Gdt, Gp = load_arrays()
sim.load_parameters(POSITION_ARRAY, MAGFIELD_ARRAY, Gdt, Gp)

sim_length_s = 20000

earth_radius = 6.3781E6 #m

pos = np.array([565000+earth_radius,0,0])
vel = np.array([0,7657.419503231067,0])
rot = [np.array([1, 0, 0]),
       np.array([0, 1, 0])]

rvl = np.array([-0.1,0,0.1])

sim.set_initial_params(np.array([1,0,0]), vel, rot, rvl)

SD = sim.run_simulation(sim_length_s)

"""
self.state_data = {"Xpos": [], "Ypos": [], "Zpos": [],
                   "Xvel": [], "Yvel": [], "Zvel": [],
                   "Yaw": [], "Pitch": [], "Roll": [],
                   "XrVel": [], "YrVel": [], "ZrVel": [],
                   "Xmag": [], "Ymag": [], "Zmag": [],
                   "Xtorque": [], "Ytorque": [], "Ztorque": [],
                   "Xsatmag": [], "Ysatmag": [], "Zsatmag": []}
"""

zero = np.zeros(len(SD["Xpos"]))

"""
linegraph(Xrvl, "Time", "X rotation velocity")
linegraph(Yrvl, "Time", "Y rotation velocity")
linegraph(Zrvl, "Time", "Z rotation velocity")
linegraph(XT, "Time", "X Torque", size=(15,3))
linegraph(YT, "Time", "Y Torque", size=(15,3))
linegraph(ZT, "Time", "Z Torque", size=(15,3))
linegraph(XB, "Time", "X MagField", size=(15,3))
linegraph(YB, "Time", "Y MagField", size=(15,3))
linegraph(ZB, "Time", "Z MagField", size=(15,3))
linegraph(satXB, "Time", "Sat X MagField", size=(15,3))
linegraph(satYB, "Time", "Sat Y MagField", size=(15,3))
linegraph(satZB, "Time", "Sat Z MagField", size=(15,3))
"""

#multiLinegraph([np.array(Xpos)*(180/np.pi),np.array(Ypos)*180/np.pi,np.array(Zpos)*180/np.pi], "Time", "Angle", size=(15,3))
multiLinegraph([SD["XrVel"],SD["YrVel"],SD["ZrVel"],zero], "Time", "Angular velocity", size=(15,3))
#multiLinegraph([VX,VY,VZ,zero], "Time", "Coil Voltages", size=(15,3))
#multiLinegraph([SD["Xtorque"],SD["Ytorque"],SD["Ztorque"],zero], "Time", "Torque", size=(15,3))
#multiLinegraph([SD["Xmag"],SD["Ymag"],SD["Zmag"],zero], "Time", "MagField", size=(15,3))
#multiLinegraph([SD["Xsatmag"],SD["Ysatmag"],SD["Zsatmag"],zero], "Time", "Sat MagField", size=(15,3))
