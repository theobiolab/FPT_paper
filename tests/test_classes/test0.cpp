#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>   //return matrices in python
#include <pybind11/embed.h>
#include <pybind11/stl.h> 
#include <cstdlib>

#include <iostream>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/mpc.hpp>
#include <Eigen/Dense>

using namespace Eigen;
namespace py=pybind11;

template <typename T, typename IOType>
void modify_class_attribute(T &attribute){ //important to pass by reference
        std::cout << "before" << attribute << "\n";
        std::cout.flush();
        attribute = 3;
        std::cout << "after" << attribute << "\n";
        std::cout.flush();
    } 


template<typename T, typename IOType>
class Base {
    protected:
        T attribute;

    public: 
        void set_attribute_class(){
            std::cout << "parent function\n";
            std::cout.flush();
        }

};

template <typename T, typename IOType> 
class child: public Base<T, IOType>{

    public:        
    
        
        void set_attribute_class(){
            std::cout << "Child function\n";
            std::cout.flush();
            modify_class_attribute<T, IOType>(this->attribute);
            std::cout << "Printing attribute\n";
            std::cout << this->attribute << "\n";
            std::cout.flush();
        }

      
};

const int MANTISSA_PRECISION = 20;
typedef boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<MANTISSA_PRECISION>> PreciseType;


PYBIND11_MODULE(test0, m){
    py::class_<child<long double, long double>>(m, "test_class0", py::module_local())
        .def(py::init())
        .def("set_attribute_class", &child<long double, long double>::set_attribute_class);
}

