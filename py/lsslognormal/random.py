import lsslognormal
import lsslognormal._lsslognormal as c

def generate_octant(nc, boxsize, nbar, *,
                     seed=None,
                     r_min=None, r_max=None):
    """
    Args:
      nc (int): number of grids per dimension
      boxsize (float): periodic box length on a side
      nbar (float): number density
      
      seed (int): random seed (default is random.random)

      r_min (float): r_min for kind='octant-shell'
      r_max (float): r_max for kind='octant-shell'

    Returns:
      Random particles in an octant shell
    """
    
    particles = lsslognormal.Particles()

    # random seed
    if seed is None:
        seed = random.randint(1, 2**31)

    c._random_generate_octant(nc, boxsize, nbar, seed,
                              r_min, r_max,
                              particles._particles)
    
    return particles
