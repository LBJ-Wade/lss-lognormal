import random
import numpy as np
import lsslognormal
import lsslognormal._lsslognormal as c

def generate_lognormal_mock(ps, nc, boxsize, n, *,
                            redshift=0.0,
                            seed=None,
                            rsd=False, omega_m=None, bias=1.0,
                            fix_amplitude=False):
    """
    Args:
      ps (str or lsslognormal.PowerSpectrum): file name or PowerSpectrum
      nc (int): number of grids per dimension
      boxsize (float): periodic box length on a side
      n (int): number of output particles
      
      seed (int): random seed (default is random.random)
      rsd (bool): apply redshift-space distortion (False -> real space)
      omega_m (float): Omega matter
      bias (float): linear bias
      fix_amplitude: fix |delta(k)| instead of Gaussian random noise

    Note:
      Particle distribution is uniform within the grid of size boxsize/nc.
    """
    # power spectrum 
    if isinstance(ps, str):
        ps = lsslognormal.PowerSpectrum(ps)
        
    if not isinstance(ps, lsslognormal.PowerSpectrum):
        raise TypeError('ps is neither file name or PowerSpectrum')
    
    particles = lsslognormal.Particles()

    # growth factor D(z)
    if redshift == 0.0:
        growth = 1.0
    else:
        growth = lsslognormal.growth.D(omega_m, z=redshift)

    # growth rate f = dln D/dln a
    if rsd:
        f = lsslognormal.growth.f(omega_m, z=redshift)
    else:
        f = 0.0

    # random seed
    if seed is None:
        seed = random.randint(1, 2**31)
        
    mock = c._mock_lognormal_generate(ps._ps,
                                      nc, boxsize, n, seed,
                                      growth, bias, f, int(fix_amplitude),
                                      particles._particles)
    
    return particles

    
