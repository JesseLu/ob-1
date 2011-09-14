import sim_2dte

dims = (40, 30)
sim = sim_2dte.simulate(dims, 0.4)
sim.insert_structure()
b_val = sim.da.getVecArray(sim.b) # Obtain access to elements of b.
b_val[dims[0]/2, dims[1]/2] = 1; # Define location of point source.

sim.solve('bcgs')
# print sim.A.getValues(range(9), range(9))
