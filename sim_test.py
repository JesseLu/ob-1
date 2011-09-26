from em_2dte import sim
import numpy as np
from plot import gnuplotter

gplot = gnuplotter.create(debug_option=0)

epsilon = np.concatenate((np.ones((60, 15)), 12.25*np.ones((60, 9)), 
                            np.ones((60, 15))), axis=1)
sim = sim.solve(0.2, epsilon)
sim.plot(gplot)
raw_input('Please press return to continue...\n')
