.PHONY: all build-dep build clean

build:
	python -m pip install .

build-dep:
	python -m pip install -r ./requirements.txt

test:
	python -m unittest discover ./tests

all: clean build-dep build test

help:
	python ./make.py --help

clean:
	rm -rf build/*
	rm -rf _cython/build/*
	rm -f src/*/*.pyd
	rm -f src/*/*/*.pyd
	rm -f src/*/*/*/*.pyd
	rm -f src/*/*.so
	rm -f src/*/*/*.so
	rm -f src/*/*/*/*.so
