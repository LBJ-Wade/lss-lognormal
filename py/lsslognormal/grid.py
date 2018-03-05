import lssps
import lssps._lssps as c
from numbers import Number
import warnings

class Grid:
    """Grid is a 3-dimensional cubic grid with nc points per dimension

    Methods:
        grid[:]: return grid data as an array
        clear(): reset all data to 0
        fft_forward():   Fast Fourier Transform to Fourier space

    Properties:
        grid.boxsize:
        grid.x0:

    """

    def __init__(self, nc, boxsize, x0=None, offset=None):
        # _grid is the pointer to a C++ Grid object
        self._grid = c._grid_alloc(nc)
        self.boxsize = boxsize

        if x0 is not None:
            self.x0 = x0

        if offset is not None:
            self.offset = offset

    def __getitem__(self, index):
        """grid[ix, iy, iz]: value at grid point (ix, iy, iz)
           grid[:]         : whole grid as a numpy.array"""
        if self.mode == 'real-space':
            return c._grid_fx_asarray(self._grid)[index]
        elif self.mode == 'fourier-space':
            return c._grid_fk_asarray(self._grid)[index]

    def clear(self):
        """Reset the grid with zeros"""
        c._grid_clear(self._grid)
        self.mode = 'real-space'
        self.interlacing = None
        
        if self.shifted is not None:
            self.shifted.clear()

        return self

    def fft_forward(self):
        """Fast Fourier Transform from real space to Fourier space"""
        c._grid_fft_forward(self._grid)

        return self

    @property
    def mode(self):
        """Current mode of the grid:
        'real-space' or 'fourier-space'
        """
        return c._grid_get_mode(self._grid)

    @mode.setter
    def mode(self, m):
        if m == 'real-space':
            c._grid_set_mode(self._grid, 1)
        elif m == 'fourier-space':
            c._grid_set_mode(self._grid, 2)
        else:
            raise TypeError('grid.mode must be real-space or fourier-space')

    @property
    def nc(self):
        """Number of grid points per dimension"""
        return c._grid_nc(self._grid)

    @property
    def boxsize(self):
        """Length of the cubic box on a side"""
        return c._grid_get_boxsize(self._grid)

    @boxsize.setter
    def boxsize(self, value):
        c._grid_set_boxsize(self._grid, value)

    @property
    def x0(self):
        """Coordinate of the box corner"""
        return c._grid_get_x0(self._grid)

    @x0.setter
    def x0(self, x0):
        """Coordinate of the box corner
           grid.x0 = (0.0, 0.0, 0.0)"""
        c._grid_set_x0(self._grid, x0[0], x0[1], x0[2])


def zeros(nc, boxsize, x0=None, offset=0.0, *, interlacing=False):
    """Return a new empty grid filled with zeros
    Args:
        nc (int): number of grids per dimension
        boxsize: length of the box on a side
        x0: corner of the box
        offset: offset of the grid points from the box corner
                give a tuple of 2 floats to set both offsets for
                the main and shifted grid
        interlacing: attach shifted grid which is shifted by half gridspacing
    """

    grid = Grid(nc, boxsize, x0, offset)
    grid.clear()

    return grid

