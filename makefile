CC = c++
CFLAGS = -O2 -Wall -shared -std=c++11 -fPIC 
eigen = -I /usr/lib/eigen-3.4.0
boost = -I /usr/lib/boost_1_71_0
pybind11 = `python3 -m pybind11 --includes`
libs = -lmpfr -lmpc

all : test_pybind test_boost test_mpfr Ladder_3 Ladder_6 Ladder_6_prec_100 Ladder_6_prec_50 triangle_graph

test_pybind: test/pybind_test.cpp
	$(CC) $(CFLAGS) $(pybind11) test/pybind_test.cpp -o test/pybind_test.so

test_boost: test/boost_test.cpp
	$(CC) $(CFLAGS) test/boost_test.cpp -o test/boost_test.so $(libs)

test_mpfr: test/mpfr_test.cpp
	$(CC) $(CFLAGS) test/mpfr_test.cpp -o test/mpfr_test.so $(libs)

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
	cp Ladder_3.so /home/fpt/bin
	cp Ladder_6.so /home/fpt/bin
	cp Ladder_6_prec_100.so /home/fpt/bin
	cp Ladder_6_prec_50.so /home/fpt/bin
	cp triangle_graph.so /home/fpt/bin
	
	cp Ladder_3.so /usr/lib/python3/dist-packages
	cp Ladder_6.so /usr/lib/python3/dist-packages
	cp Ladder_6_prec_100.so /usr/lib/python3/dist-packages
	cp Ladder_6_prec_50.so /usr/lib/python3/dist-packages
	cp triangle_graph.so /usr/lib/python3/dist-packages

clean: 
	rm Ladder_3.so
	rm Ladder_6.so
	rm Ladder_6_prec_100.so
	rm Ladder_6_prec_50.so
	rm triangle_graph.so
