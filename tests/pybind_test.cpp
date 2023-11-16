#include <pybind11/pybind11.h>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/mpc.hpp>
#include <Eigen/Dense>

int add(int i, int j) {
    return i + j;
}

PYBIND11_MODULE(pybind_test, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("add", &add, "A function that adds two numbers");
}