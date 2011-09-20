from em_2dte import sim
import numpy as np
from plot import gnuplotter

gplot = gnuplotter.create(debug_option=1)

epsilon = np.concatenate((np.ones((40, 25)), 4.25*np.ones((40, 10)), 
                            np.ones((40, 25))), axis=1)
sim = sim.FieldSolve(epsilon, 0.6)
sim.insert_structure()
b_val = sim.da.getVecArray(sim.b) # Obtain access to elements of b.
b_val[epsilon.shape[0]/2, epsilon.shape[1]/2] = 1; # Define location of point source.

sim.solve()
sim.plot(gplot)
raw_input('Please press return to continue...\n')
# print sim.A.getValues(range(9), range(9))

# pylab.contourf(sim.da.getVecArray(sim.x)[:])
# pylab.show()
