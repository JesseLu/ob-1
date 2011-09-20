import numpy as np
from em_2dte import wgmode
from plot import gnuplotter

gplot = gnuplotter.create(debug_option=1)
epsilon = np.concatenate((np.ones(40), 12.25*np.ones(10), np.ones(40)))
wg_mode = wgmode.solve(0.1, epsilon, 3)
# print wg_mode.field
wg_mode.plot(gplot)
raw_input('Please press return to continue...\n')
