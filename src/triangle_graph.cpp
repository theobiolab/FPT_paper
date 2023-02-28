#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>   //return matrices in python
#include <pybind11/embed.h>
#include <pybind11/stl.h>     // py::cast

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
Matrix<T,Dynamic,Dynamic> set_laplacian_TF(IOType TF_conc, Matrix<T,Dynamic,Dynamic> L){
    
        int n_states = L.cols();
        int n_states_line = (n_states-1)/2;
         
        Matrix<T, Dynamic, Dynamic> laplacian_mult = L; 
        
        for (int i = n_states_line; i<n_states-1; i++){
            T TF_conc_casted = (T) TF_conc;
            laplacian_mult(i,i-n_states_line) = laplacian_mult(i,i-n_states_line)*TF_conc; 
        }
        
        
        for (int i=0; i<n_states; i++){
            laplacian_mult(i,i) = 0.0; 
            laplacian_mult(i,i) = -(laplacian_mult.col(i).sum()); 
        }
            
        return laplacian_mult.template cast<T>(); 
    }
    

template <typename T, typename IOType>
IOType compute_response(Matrix<IOType,Dynamic,1> ss_vector){  
    IOType response; 
    //ladder 3 response = Pr^ss(state 3) + Pr^ss(state 3T)
    response = ss_vector[2] + ss_vector[5];
    return response;

}

#include "../include/scoring2D.hpp"

const int MANTISSA_PRECISION = 100;
typedef boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<MANTISSA_PRECISION>> PreciseType;

template <typename T, typename IOType> 
class childScoring: public scoring2D<T, IOType>{
    
    public:
    
    void set_laplacian_from_child_class(Matrix<IOType,Dynamic,1> pset){
    
        ///error for the lenght of partset
    
        Matrix<T, 4, 4> lap; 
        
        Matrix<T,Dynamic,1> parset = pset.template cast<T>();  
        
        
        T l12 = parset(0);
        T l13 = parset(1);
        T l21 = parset(2);
        T l23 = parset(3);
        T l31 = parset(4);
        T l32 = parset(5);
        T l34 = parset(6);
        
        lap <<    0,     l21,     l31,      0, 
                l12,       0,     l32,      0, 
                l13,     l23,       0,      0, 
                  0,       0,     l34,      0; 
               
        for (int i=0; i<4; i++){
            lap(i, i) = -(lap.col(i).sum());
        }

        this->laplacian = lap; 
    }       
};



PYBIND11_MODULE(triangle_graph, m){
    py::class_<childScoring<PreciseType, long double>>(m, "triangle_graph")
        .def(py::init())
        .def("getLaplacian",          &childScoring<PreciseType, long double>::get_laplacian)
        .def("setTriangleLaplacian",  &childScoring<PreciseType, long double>::set_laplacian_from_child_class)
        .def("get_FPT_stat",          &childScoring<PreciseType, long double>::get_Mean_CV_FPT)
        .def("getSteadyState",        &childScoring<PreciseType, long double>::get_steady_state)
        .def("compute_response",      &childScoring<PreciseType, long double>::Response_class)  
        .def("mult_laplacian_TF",     &childScoring<PreciseType, long double>::set_laplacian_TF_class) 
        .def("getLaplacianTF",        &childScoring<PreciseType, long double>::get_laplacian_TF)
        .def("simple_score",          &childScoring<PreciseType, long double>::SimpleScore, 
                                        py::arg("TF_range_for_score"), py::arg("n_TF_conc_points"))
        .def("monotonicity",          &childScoring<PreciseType, long double>::monotonicity_test,
                                        py::arg("range_to_check"), py::arg("n_points_for_ss"), py::arg("delta_TH"))        
        .def("score_no_response",     &childScoring<PreciseType, long double>::scoreNoResponse, 
                                        py::arg("TF_range_for_score"), py::arg("n_points_for_ss"), py::arg("delta_TH"))      
        .def("score_response",        &childScoring<PreciseType, long double>::scoreWithResponse,
                                        py::arg("n_points_for_ss"), py::arg("delta_TH"))
        .def("score_2Dgaus_delta",    &childScoring<PreciseType, long double>::score2DgaussianDelta, 
                                        py::arg("starting_TF_range"), py::arg("TF_range_no_response"), 
                                        py::arg("n_points_score"), py::arg("n_points_ss"), py::arg("delta_TH")); 
}


