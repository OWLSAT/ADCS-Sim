# ADCS-Sim
ADCS satellite detumbling, orientation determination and orbit simulation code

Should be able to run detumbling simulations from the Experiment2.py file.

Orbit position and magnetic field arrays are precomputed. If they need to be recalculated use the code in the ComputeOrbit.py file

Simulation files are:
 - Environment.py (Provides the Low-Earth-Orbit environment)
 - SatBuilder.py (Provides the classes to construct satellite objects + define sensors and actuators - CONTAINS THE DETUMBLING ALGORITHMS)
 - OwlSat.py (Uses the SatBuilder to create a simulated version of OwlSat - Physical properties as well as sensor/actuator set)
 - Simulation.py (Processes and executes the simulation of the satellite and collects the data)
 - Helpers.py (Contains various coordinate transformation functions + linear algebra helpers)

Additional files:
 - ComputeOrbit.py (Used to precompute position and magnetic field values for a given orbit - Saves incredible amounts of computation time when simulating)
 - orbitArrays.txt (File where the precomputed orbit position and magnetic field values are stored - Used by the Simulation.py code)
 - Experiment1.py & Experiment2.py (Use cases of the simulation to run various experiments - Experiment1 might be obsolete)
 - ppigrf.py (International Geomagnetic Reference Field library - Provides magnetif field model of the Earth)
 - IGRF13.sch (Database used by the ppigrf library to compute magnetic field values)
