from em_2dte import sim
import numpy as np
from plot import gnuplotter

gplot = gnuplotter.create(debug_option=1)

epsilon = np.concatenate((np.ones((60, 15)), 1*np.ones((60, 11)), 
                            np.ones((60, 15))), axis=1)
sim = sim.solve(0.2, epsilon, d_pml=10, mode_pos=30)
# sim.insert_structure(0.1)

# sim.solve(0.1)
sim.plot(gplot)
# gplot.save('test.ps')
raw_input('Please press return to continue...\n')
# print sim.A.getValues(range(9), range(9))

# pylab.contourf(sim.da.getVecArray(sim.x)[:])
# pylab.show()
