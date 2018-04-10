all: install

# compiler
OPENMP    := -fopenmp
CXXOPENMP ?= g++
CXX       ?= $(CXXOPENMP) $(OPENMP)
CC        := $(CXX)

# library paths
FFTW3_DIR ?= #e.g. /Users/jkoda/Research/opt/gcc/fftw3
GSL_DIR   ?= 
HDF5_DIR  ?=

DIR_PATH   = $(FFTW3_DIR) $(GSL_DIR) $(HDF5_DIR)

IDIRS    += $(foreach dir, $(DIR_PATH), $(dir)/include)
LDIRS    += $(foreach dir, $(DIR_PATH), $(dir)/lib)

# linking libraries
LIBS += fftw3 gsl gslcblas

ifdef OPENMP
  LIBS += fftw3$(FFTWSUF)_omp
  #LIBS += fftw3$(FFTWSUF)_threads # for thread parallelization instead of omp
endif

export CC CXX IDIRS LDIRS LIBS


SRC := growth.cpp power.cpp

build: setup.py $(SRC)
	python3 setup.py build_ext --inplace

install:
	python3 setup.py install

.PHONY: clean

clean:
	rm -rf build fs/*.so
