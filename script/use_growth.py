"""
Example using lsslognormal.growth module

growth.D(omega_m, z=z, a=a)
growth.f(omega_m, z=z, a=a)

Provide either redshift z or scale factor a
"""

import lsslognormal

omega_m = 0.308
z = 1.0


# growth factor
d = lsslognormal.growth.D(omega_m, z=1.0)

print('D(z=%.1f) = %.4f' % (z, d))

# growth rate
f = lsslognormal.growth.f(omega_m, z=1.0)

print('f(z=%.1f) = %.4f' % (z, f))


