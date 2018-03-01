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
                     'py_mock_lognormal.cpp', 'py_growth.cpp',
                     'power_spectrum.cpp', 'growth.cpp',
                    ],
                    libraries = ['gsl', 'gslcblas'],
                    undef_macros = ['NDEBUG'],
          )
      ],
      packages=['lsslognormal'],
)


