.PHONY: all build-dep build clean

build:
	python3 ./src/setup.py build_ext --build-lib .

build-dep:
	python3 -m pip install -r ./requirements.txt

test:
	python3 -m unittest discover ./unittests

all: clean build-dep build test

help:
	python3 ./make.py --help

clean:
	rm -rf build/*
	rm -rf src/cython/build/*
	rm -f src/py/*.pyd
	rm -f src/py/*/*.pyd
	rm -f src/py/*/*/*.pyd
	rm -f src/py/*.so
	rm -f src/py/*/*.so
	rm -f src/py/*/*/*.so
