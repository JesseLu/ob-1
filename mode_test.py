import numpy as np
from em_2dte import wgmode
from plot import gnuplotter

gplot = gnuplotter.create(debug_option=1)
epsilon = np.concatenate((np.ones(20), 12.25*np.ones(11), np.ones(20)))
wg_mode = wgmode.solve(0.1, epsilon, 1)
# print wg_mode.field
wg_mode.plot(gplot)
print wg_mode.beta
raw_input('Please press return to continue...\n')
