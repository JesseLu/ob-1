from em_2dte import sim

dims = (60, 40)
sim = sim.FieldSolve(dims, 0.6)
sim.insert_structure()
b_val = sim.da.getVecArray(sim.b) # Obtain access to elements of b.
b_val[dims[0]/2, dims[1]/2] = 1; # Define location of point source.

sim.solve()
# print sim.A.getValues(range(9), range(9))

# pylab.contourf(sim.da.getVecArray(sim.x)[:])
# pylab.show()
