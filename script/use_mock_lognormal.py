"""
Example using lsslognormal
"""

import lsslognormal

omega_m = 0.308
z = 1.0


# growth factor
#d = lsslognormal.growth.D(omega_m, z=1.0)

# growth rate
#f = lsslognormal.growth.f(omega_m, z=1.0)

boxsize = 1000.0    # box size on a side
nc = 64
n = 10              # number of output particles (mean)
z = 1.0

mock = lsslognormal.generate_lognormal_mock("../data/planck_matterpower.dat",
                                            nc, boxsize, n)
#, *,
#                            redshift=0.0,
#                            seed=None,
#                            rsd=False, omega_m=None, bias=1.0,
#                            fix_amplitude=False):


print(mock.asarray())


