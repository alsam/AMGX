# several sample targets

debug:
	mkdir -p build_debug && cd build_debug && cmake .. -DCMAKE_BUILD_TYPE=Debug -DCUDA_ARCH="35 52 60" -G Ninja && ninja -j4
