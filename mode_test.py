import numpy as np
from em_2dte import wgmode

epsilon = np.concatenate((np.ones(40), 12.25*np.ones(10), np.ones(40)))
wg_mode = wgmode.solve(0.1, epsilon, 2)
print wg_mode.field
