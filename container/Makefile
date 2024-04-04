CC = c++
CFLAGS = -O2 -Wall -shared -std=c++11 -fPIC 
eigen = -I $(shell pwd)/lib/eigen-3.4.0
boost = -I $(shell pwd)/lib/boost_1_81_0
pybind11 = `python3 -m pybind11 --includes`
libs = -lmpfr -lmpc

all_container: test_pybind test_boost test_mpfr Ladder_3 Ladder_6 Ladder_6_prec_100 Ladder_6_prec_50 triangle_graph
ubuntu_setup: dependencies_ubuntu boost eigen test_pybind test_boost test_mpfr
all_ubuntu: Ladder_3 Ladder_6 Ladder_6_prec_100 Ladder_6_prec_50 triangle_graph

dependencies_ubuntu: 
	sudo apt-get install -y libmpfr-dev
	sudo apt update && sudo apt upgrade -y
	sudo apt-get install -y build-essential
	sudo apt-get install -y libmpfr-dev
	sudo apt-get install -y libmpc-dev
	sudo apt-get install -y autotools-dev libicu-dev libbz2-dev libboost-all-dev
	sudo apt-get install -y wget tar unzip git 
	sudo add-apt-repository -y ppa:deadsnakes/ppa && sudo apt install -y python3.8 python3-pip
	pip3 install pybind11
	pip3 install scipy
	pip3 install matplotlib
	pip3 install pandas

boost: 
	mkdir -p ./lib 
	cd lib
	wget https://boostorg.jfrog.io/artifactory/main/release/1.81.0/source/boost_1_81_0.tar.gz
	tar -xf boost_1_81_0.tar.gz
	rm boost_1_81_0.tar.gz
	mv -v boost_1_81_0 lib

eigen: 
	mkdir -p ./lib 
	cd lib
	wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
	tar -xf eigen-3.4.0.tar.gz
	rm eigen-3.4.0.tar.gz
	mv -v eigen-3.4.0 lib

test_pybind: test/pybind_test.cpp
	$(CC) $(CFLAGS) $(pybind11) test/pybind_test.cpp -o test/pybind_test.so

test_boost: test/boost_test.cpp
	$(CC) $(boost) test/boost_test.cpp -o test/boost_test.so $(libs)

test_mpfr: test/mpfr_test.cpp
	$(CC) test/mpfr_test.cpp -o test/mpfr_test.so $(libs)

Ladder_3: src/Ladder_3.cpp
	$(CC) $(CFLAGS) $(pybind11) $(eigen) $(boost) src/Ladder_3.cpp -o Ladder_3.so $(libs)

Ladder_6 : src/Ladder_6.cpp
	$(CC) $(CFLAGS) $(pybind11) $(eigen) $(boost) src/Ladder_6.cpp -o Ladder_6.so $(libs)

Ladder_6_prec_100 : src/Ladder_6_prec_100.cpp
	$(CC) $(CFLAGS) $(pybind11) $(eigen) $(boost) src/Ladder_6_prec_100.cpp -o Ladder_6_prec_100.so $(libs)

Ladder_6_prec_50 : src/Ladder_6_prec_50.cpp
	$(CC) $(CFLAGS) $(pybind11) $(eigen) $(boost) src/Ladder_6_prec_50.cpp -o Ladder_6_prec_50.so $(libs)

triangle_graph : src/triangle_graph.cpp
	$(CC) $(CFLAGS) $(pybind11) $(eigen) $(boost) src/triangle_graph.cpp -o triangle_graph.so $(libs)

install: 
	cp Ladder_3.so $(shell pwd)/bin
	cp Ladder_6.so $(shell pwd)/bin
	cp Ladder_6_prec_100.so $(shell pwd)/bin
	cp Ladder_6_prec_50.so $(shell pwd)/bin
	cp triangle_graph.so $(shell pwd)/bin

install_container: 
	cp Ladder_3.so /home/fpt/bin
	cp Ladder_6.so /home/fpt/bin
	cp Ladder_6_prec_100.so /home/fpt/bin
	cp Ladder_6_prec_50.so /home/fpt/bin
	cp triangle_graph.so /home/fpt/bin

clean: 
	rm Ladder_3.so
	rm Ladder_6.so
	rm Ladder_6_prec_100.so
	rm Ladder_6_prec_50.so
	rm triangle_graph.so
	