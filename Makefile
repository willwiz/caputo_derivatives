
all: clean build-dep build unittest

build-dep:
	python3 -m pip install -r ./requirements.txt

build:
	python3 ./src/setup.py build_ext --build-lib .

unittest:
	python3 -m unittest discover ./unittests

buildtest: build test

clean:
	rm -rf build/*
	rm -rf src/cython/build/*
	rm -f src/py/*.pyd
	rm -f src/py/*/*.pyd
	rm -f src/py/*/*/*.pyd
	rm -f src/py/*.so
	rm -f src/py/*/*.so
	rm -f src/py/*/*/*.so
