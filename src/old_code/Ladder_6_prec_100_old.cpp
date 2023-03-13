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
Matrix<T,Dynamic,Dynamic> set_laplacian_TF(IOType TF_conc, const Matrix<T,Dynamic,Dynamic> L){
    
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
    response = ss_vector[5] + ss_vector[11];
    return response;

}

#include "../include/scoring2D_old.hpp"

const int MANTISSA_PRECISION = 100;
typedef boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<MANTISSA_PRECISION>> PreciseType;

template <typename T, typename IOType> 
class childScoring: public scoring2D<T, IOType>{
    
    public:
    
    void set_ladder_3_laplacian(Matrix<IOType,Dynamic,1> pset){
    
        ///error for the lenght of partset
    
        Matrix<T, 13, 13> lap; 
        
        Matrix<T,Dynamic,1> parset = pset.template cast<T>();  
        
        T a12  = parset(0);
        T a23  = parset(1);
        T a34  = parset(2);
        T a45  = parset(3);
        T a56  = parset(4);
        
        T b12  = parset(5); 
        T b23  = parset(6);
        T b34  = parset(7);  
        T b45  = parset(8);
        T b56  = parset(9);
        
        T a12T  = parset(10);
        T a23T  = parset(11);
        T a34T  = parset(12);
        T a45T  = parset(13);
        T a56T  = parset(14);
        
        T b12T  = parset(15); 
        T b23T  = parset(16);
        T b34T  = parset(17);  
        T b45T  = parset(18);
        T b56T  = parset(19);
        
        T kon  = parset(20);
        T koff = parset(21);
        T P    = parset(22);
        T PT   = parset(23);
        
        //        [0]     [1]    [2]      [3]      [4]      [5]    [6]    [7]    [8]     [9]   [10]    [11]  [12]
        //         1       2      3        4        5        6     1T      2T     3T     4T     5T      6T    M
        
        lap <<     0,    b12,     0,       0,       0,       0,  koff,     0,      0,     0,     0,     0,    0,
                 a12,      0,   b23,       0,       0,       0,     0,  koff,      0,     0,     0,     0,    0,
                   0,    a23,     0,     b34,       0,       0,     0,     0,   koff,     0,     0,     0,    0,
                   0,      0,   a34,       0,     b45,       0,     0,     0,      0,  koff,     0,     0,    0,
                   0,      0,     0,     a45,       0,     b56,     0,     0,      0,     0,  koff,     0,    0,
                   0,      0,     0,       0,     a56,       0,     0,     0,      0,     0,     0,  koff,    0,
                 kon,      0,     0,       0,       0,       0,     0,  b12T,      0,     0,     0,     0,    0,
                   0,    kon,     0,       0,       0,       0,  a12T,     0,   b23T,     0,     0,     0,    0,
                   0,      0,   kon,       0,       0,       0,     0,  a23T,      0,  b34T,     0,     0,    0,
                   0,      0,     0,     kon,       0,       0,     0,     0,   a34T,     0,  b45T,     0,    0,
                   0,      0,     0,       0,     kon,       0,     0,     0,      0,  a45T,     0,  b56T,    0,
                   0,      0,     0,       0,       0,     kon,     0,     0,      0,     0,  a56T,     0,    0,
                   0,      0,     0,       0,       0,       P,     0,     0,      0,     0,     0,    PT,    0; 
                       
        for (int i=0; i<13; i++){
            lap(i, i) = -(lap.col(i).sum());
        }

        this->laplacian = lap; 
    }       
};



PYBIND11_MODULE(Ladder_6_prec_100, m){
    py::class_<childScoring<PreciseType, long double>>(m, "Ladder_6_prec_100")
        .def(py::init())
        .def("getLaplacian",          &childScoring<PreciseType, long double>::get_laplacian)
        .def("setLaplacian_ladder_6", &childScoring<PreciseType, long double>::set_ladder_3_laplacian)
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


