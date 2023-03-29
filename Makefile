.PHONY: all build-dep build clean
all: clean build-dep build

build-dep:
	python3 -m pip install -r ./requirements.txt

build:
	python3 ./src/setup.py build_ext --build-lib .

clean:
	rm -rf build/*
	rm -rf src/cython/build/*
	rm -f src/py/*.pyd
	rm -f src/py/*/*.pyd
	rm -f src/py/*/*/*.pyd
	rm -f src/py/*.so
	rm -f src/py/*/*.so
	rm -f src/py/*/*/*.so
