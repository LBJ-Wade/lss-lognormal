import lsspert._lsspert as c

def D(omega_m, *, a=None, z=None):
    """
    Compute linear growth factor D(z)

    Args:
       omega_m (float): Omega_m0
       a (float): scale factor
       z (float): redshift
    """

    if a is None and z is None:
        raise ValueError('Need to give either a or z')

    if z is not None:
        if a is not None:
            print('value error?', a, z)
            raise ValueError('Cannot set both a and z')
        
        a = 1.0/(1.0 + z)

    #print('a=', a)
    #print('omega_m', omega_m)

    return c._growth_D(a, omega_m)


def f(omega_m, *, a=None, z=None):
    """
    Compute linear growth rate f(z)

    Args:
       omega_m (float): Omega_m0
       a (float): scale factor
       z (float): redshift    
    """

    if a is None and z is None:
        raise ValueError('Need to give either a or z')

    if z is not None:
        if a is not None:
            raise ValueError('Cannot set both a and z')
        a = 1.0/(1.0 + z)

    return c._growth_f(a, omega_m)
