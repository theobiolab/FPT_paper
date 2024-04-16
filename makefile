CC = c++
mac_CFLAGS = -O2 -shared -std=c++14 -fPIC -undefined dynamic_lookup
linux_CFLAGS = -O2 -shared -std=c++14 -fPIC

# eigen and boost paths in lib folder (assumes make dependencies has been executed)
eigen = -I $(shell pwd)/lib/eigen-3.4.0
boost = -I $(shell pwd)/lib/boost_1_81_0

# gmp, mpfr, mpc and pybind11 paths in lib folder (assumes make dependencies has been executed)
lgmp = -L$(shell pwd)/lib/gmp-6.3.0/.libs -lgmp
igmp = -I$(shell pwd)/lib/gmp-6.3.0
lmpfr = -L$(shell pwd)/lib/mpfr-4.2.1/src/.libs -lmpfr
impfr = -I$(shell pwd)/lib/mpfr-4.2.1/src
lmpc = -L$(shell pwd)/lib/mpc-1.3.1/src/.libs -lmpc
impc = -I$(shell pwd)/lib/mpc-1.3.1/src
pybind11 = `python -m pybind11 --includes`

# checks for system specific compile flags
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	CFLAGS = $(linux_CFLAGS)
else ifeq ($(UNAME_S),Darwin)
	CFLAGS = $(mac_CFLAGS)
else
	$(error Unsupported operating system: $(UNAME_S))
endif

dependencies_docker: boost_and_eigen
ladders_docker: Ladder_3_dock Ladder_4_dock Ladder_5_dock Ladder_6_dock

dependencies: boost_and_eigen libgmp libmpfr libmpc
tests: test_pybind test_boost test_mpfr test_mpc test_gmp
ladders: Ladder_3 Ladder_4 Ladder_5 Ladder_6

boost_and_eigen:
	@mkdir -p ./lib 
	@cd lib

	wget https://boostorg.jfrog.io/artifactory/main/release/1.81.0/source/boost_1_81_0.tar.gz
	@tar -xf boost_1_81_0.tar.gz
	@rm boost_1_81_0.tar.gz

	wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
	@tar -xf eigen-3.4.0.tar.gz
	@rm eigen-3.4.0.tar.gz

	@mv -v boost_1_81_0 lib
	@mv -v eigen-3.4.0 lib

libgmp:
	wget https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz
	@tar -xf gmp-6.3.0.tar.xz
	@rm gmp-6.3.0.tar.xz
	@mv -v gmp-6.3.0 lib
	@cd lib/gmp-6.3.0 && ./configure && make

libmpfr:
	wget https://www.mpfr.org/mpfr-current/mpfr-4.2.1.tar.xz
	@tar -xf mpfr-4.2.1.tar.xz
	@rm mpfr-4.2.1.tar.xz
	@mv -v mpfr-4.2.1 lib
	@cd lib/mpfr-4.2.1 && \
		LDFLAGS="-L$(shell pwd)/lib/gmp-6.3.0/.libs" \
		CPPFLAGS="-I$(shell pwd)/lib/gmp-6.3.0" \
		./configure --with-gmp=$(shell pwd)/lib/gmp-6.3.0 && \
		make

libmpc:
	wget https://ftp.gnu.org/gnu/mpc/mpc-1.3.1.tar.gz
	@tar -xf mpc-1.3.1.tar.gz
	@rm mpc-1.3.1.tar.gz
	@mv -v mpc-1.3.1 lib
	@cd lib/mpc-1.3.1 && \
		LDFLAGS="-L$(shell pwd)/lib/mpfr-4.2.1/src/.libs" \
		CPPFLAGS="-I$(shell pwd)/lib/mpfr-4.2.1/src" \
		./configure --with-mpfr=$(shell pwd)/lib/mpfr-4.2.1 && \
		make

test_pybind: tests/pybind_test.cpp
	$(CC) $(CFLAGS) $(pybind11) $(eigen) $(boost) $(igmp) $(impfr) $(impc) tests/pybind_test.cpp -o tests/pybind_test.so ${lgmp} ${lmpfr} ${lmpc} 

test_boost: tests/boost_test.cpp
	$(CC) $(CFLAGS) $(pybind11) $(eigen) $(boost) $(igmp) $(impfr) $(impc) tests/boost_test.cpp -o tests/boost_test.so ${lgmp} ${lmpfr} ${lmpc} 

test_mpfr: tests/mpfr_test.cpp
	$(CC) $(CFLAGS) $(pybind11) $(eigen) $(boost) $(igmp) $(impfr) $(impc) tests/mpfr_test.cpp -o tests/mpfr_test.so ${lgmp} ${lmpfr} ${lmpc} 

test_mpc: tests/mpc_test.cpp
	$(CC) $(CFLAGS) $(pybind11) $(eigen) $(boost) $(igmp) $(impfr) $(impc) tests/mpc_test.cpp -o tests/mpc_test.so ${lgmp} ${lmpfr} ${lmpc} 

test_gmp: tests/gmp_test.cpp
	$(CC) $(CFLAGS) $(pybind11) $(eigen) $(boost) $(igmp) $(impfr) $(impc) tests/gmp_test.cpp -o tests/gmp_test.so ${lgmp} ${lmpfr} ${lmpc} 

Ladder_3: src/Ladder_3_v3.cpp
	$(CC) $(CFLAGS) $(pybind11) $(eigen) $(boost) $(igmp) $(impfr) $(impc) src/Ladder_3_v3.cpp -o bin/Ladder_3_v3.so ${lgmp} ${lmpfr} ${lmpc} 

Ladder_4: src/Ladder_4.cpp
	$(CC) $(CFLAGS) $(pybind11) $(eigen) $(boost) $(igmp) $(impfr) $(impc) src/Ladder_4.cpp -o bin/Ladder_4.so ${lgmp} ${lmpfr} ${lmpc} 

Ladder_5: src/Ladder_5.cpp
	$(CC) $(CFLAGS) $(pybind11) $(eigen) $(boost) $(igmp) $(impfr) $(impc) src/Ladder_5.cpp -o bin/Ladder_5.so ${lgmp} ${lmpfr} ${lmpc} 

Ladder_6: src/Ladder_6_prec_100_v2.cpp
	$(CC) $(CFLAGS) $(pybind11) $(eigen) $(boost) $(igmp) $(impfr) $(impc) src/Ladder_6_prec_100_v2.cpp -o bin/Ladder_6_prec_100_v2.so ${lgmp} ${lmpfr} ${lmpc} 

Ladder_3_dock: src/Ladder_3_v3.cpp
	$(CC) $(CFLAGS) $(pybind11) $(eigen) $(boost) src/Ladder_3_v3.cpp -o bin/Ladder_3_v3.so -lmpfr -lmpc

Ladder_4_dock: src/Ladder_4.cpp
	$(CC) $(CFLAGS) $(pybind11) $(eigen) $(boost) src/Ladder_4.cpp -o bin/Ladder_4.so -lmpfr -lmpc

Ladder_5_dock: src/Ladder_5.cpp
	$(CC) $(CFLAGS) $(pybind11) $(eigen) $(boost) src/Ladder_5.cpp -o bin/Ladder_5.so -lmpfr -lmpc

Ladder_6_dock: src/Ladder_6_prec_100_v2.cpp
	$(CC) $(CFLAGS) $(pybind11) $(eigen) $(boost) src/Ladder_6_prec_100_v2.cpp -o bin/Ladder_6_prec_100_v2.so -lmpfr -lmpc


test_ladders: tests/test_calculations.py
	source config.sh && cd tests && python test_calculations.py
