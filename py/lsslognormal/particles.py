import lsslognormal._lsslognormal as c

class Particles:
    def __init__(self):
        self._particles = c._particles_alloc()

    def asarray(self):
        return c._particles_asarray(self._particles)

    def clear(self):
        pass # not implemented
    
