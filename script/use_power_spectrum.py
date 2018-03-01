"""
Example using lsslognormal.PowerSpectrum class
"""

import lsslognormal

ps = lsslognormal.PowerSpectrum('../data/planck_matterpower.dat')

print('Pk = %e' % ps(0.1)) # P(k = 0.1)

