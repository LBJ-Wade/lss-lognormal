from distutils.core import setup, Extension
import numpy as np

import os


# Remove the "-Wstrict-prototypes" compiler option, which isn't valid for C++.
# https://stackoverflow.com/questions/8106258
import distutils.sysconfig
cfg_vars = distutils.sysconfig.get_config_vars()
for key, value in cfg_vars.items():
    if type(value) == str:
        cfg_vars[key] = value.replace("-Wstrict-prototypes", "")

# directories for include -I(idir)
idirs = os.environ["IDIRS"]
if idirs:
    idirs = idirs.split()
else:
    idirs = []

idirs = [np.get_include(),] + idirs # path for numpy

# directories for libraries -L(dir)
ldirs = os.environ["LDIRS"]
if ldirs:
    ldirs = ldirs.split()
else:
    ldirs = []

#ldirs = ['/storage/local/home/astrorm3/koda/miniconda3/lib', ] + ldirs
libs = os.environ['LIBS'].split()

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
                    include_dirs = idirs,
                    #extra_compile_args = [os.environ['OPT']],
                    library_dirs = ldirs,
                    libraries = libs,
                    #libraries = ['gsl', 'gslcblas', 'fftw3'],
                    undef_macros = ['NDEBUG'],
          )
      ],
      packages=['lsslognormal'],
)


