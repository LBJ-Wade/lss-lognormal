from distutils.core import setup, Extension
import numpy as np

setup(name='lsslognormal',
      version='0.0.1',
      description='Lognormal mock for Large-scale structure cosmology',
      author='Jun Koda',
      py_modules=['lsslognormal.mock_lognormal', 'lsslognormal.growth',
      ],
      ext_modules=[
          Extension('lsslognormal._lsslognormal',
                    ['py_package.cpp',
                     'py_growth.cpp', 'py_power_spectrum.cpp',
                     'py_particles.cpp', 'py_grid.cpp',
                     'py_mock_lognormal.cpp', 
                     'power_spectrum.cpp', 'growth.cpp',
                     'grid.cpp',
                     'mock_lognormal.cpp',
                    ],
                    include_dirs = [np.get_include(),],
                    libraries = ['gsl', 'gslcblas', 'fftw3'],
                    undef_macros = ['NDEBUG'],
          )
      ],
      packages=['lsslognormal'],
)


