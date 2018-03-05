
all: install

SRC := growth.cpp power.cpp

build: setup.py $(SRC)
	python3 setup.py build_ext --inplace

install:
	python3 setup.py install

.PHONY: clean

clean:
	rm -rf build fs/*.so
