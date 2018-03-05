import lsslognormal._lsslognormal as c

class PowerSpectrum:
    def __init__(self, filename):
        try:
            self._ps = c._power_spectrum_alloc(filename)
        except FileNotFoundError:
            raise FileNotFoundError('Unable to open file: %s' % filename)
        except OSError:
            raise OSError('Error reading file: %s' % filename)


    def __call__(self, k):
        """P(k)"""
        
        return c._power_spectrum_call(self._ps, k)
