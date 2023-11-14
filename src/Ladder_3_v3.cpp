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

#include "../include/scoring2D_v3.hpp"

using namespace Eigen;
namespace py=pybind11;

template <typename T, typename IOType>
void set_laplacian(Matrix<IOType,Dynamic,1> parset, Matrix<T, Dynamic, Dynamic> &lap){ //modifying laplacian attribute directly
    
        
        //Matrix<T,Dynamic,1> parset = pset.template cast<T>();  
        
        T a12  = parset(0);
        T a23  = parset(1);
        T b12  = parset(2); // check indexing 
        T b23  = parset(3);
        T a12T = parset(4);
        T a23T = parset(5);
        T b12T = parset(6);
        T b23T = parset(7);
        T kon  = parset(8);
        T koff = parset(9);
        T P    = parset(10);
        T PT   = parset(11);
        
        
        lap <<     0,    b12,     0,      koff,        0,       0,     0, 
                 a12,      0,   b23,          0,    koff,       0,     0, 
                   0,    a23,     0,          0,       0,    koff,     0,
                 kon,      0,     0,          0,    b12T,       0,     0, 
                   0,    kon,     0,       a12T,       0,    b23T,     0,
                   0,      0,   kon,          0,    a23T,       0,     0,
                   0,      0,     P,          0,       0,      PT,     0;
               
        for (int i=0; i<7; i++){
            lap(i, i) = -(lap.col(i).sum());
        }

    }  

template <typename T, typename IOType>
void modify_TF_conc(IOType TF_conc, const Matrix<T,Dynamic,Dynamic> &L, Matrix<T,Dynamic,Dynamic> &laplacian_mult ){
    
        int n_states = L.cols();
        int n_states_line = (n_states-1)/2; 
        
        for (int i = n_states_line; i<n_states-1; i++){
            T TF_conc_casted = (T) TF_conc;
            laplacian_mult(i,i-n_states_line) = L(i,i-n_states_line)*TF_conc; 
        }
        
        
        for (int i=0; i<n_states; i++){
            laplacian_mult(i,i) = 0.0; 
            laplacian_mult(i,i) = -(laplacian_mult.col(i).sum()); 
        }
            
        //return laplacian_mult.template cast<T>(); 
    }  
    

template <typename T, typename IOType>
IOType response_function(Matrix<IOType,Dynamic,1> &ss_vector){  //maybe ss_vector should be T,Dynamic,1 and response should also be T ?
    IOType response; 
    //ladder 3 response = Pr^ss(state 3) + Pr^ss(state 3T)
    response = ss_vector[2] + ss_vector[5];
    return response;
}


const int MANTISSA_PRECISION = 100;
typedef boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<MANTISSA_PRECISION>> PreciseType;



template <typename T, typename IOType> 
class childscoring_ladder3: public scoring2D_v3<T, IOType>{
    
    public:

        childscoring_ladder3(){
            this->laplacian.resize(7,7);
            this->laplacian_TF.resize(7,7);
        }

        void set_laplacian_class(Matrix<IOType,Dynamic,1> pset){
            set_laplacian<T, IOType>( pset, this->laplacian);
            this->laplacian_TF=this->laplacian; //the first time, assign values, then when set_laplacian_TF is called, only the cells with TF conc and the diagonals need to be modified
        }

        void set_laplacian_TF(IOType TF_conc){
            modify_TF_conc<T, IOType>(TF_conc, this->laplacian, this->laplacian_TF);
        }

        IOType compute_response(Matrix<IOType,Dynamic,1> ss_vector){ //should be T
            return response_function<T, IOType>(ss_vector);
        }
    
};



PYBIND11_MODULE(Ladder_3_v3, m){
    py::class_<childscoring_ladder3<PreciseType, long double>>(m, "Ladder_3_v3", py::module_local())
        .def(py::init())
        .def("setLaplacian", &childscoring_ladder3<PreciseType, long double>::set_laplacian_class)
        .def("getLaplacian",          &childscoring_ladder3<PreciseType, long double>::get_laplacian)
        .def("getSteadyState",        &childscoring_ladder3<PreciseType, long double>::get_steady_state)
        .def("get_FPT_stat",          &childscoring_ladder3<PreciseType, long double>::get_Mean_CV_FPT)
        .def("set_TF_and_compute_response",      &childscoring_ladder3<PreciseType, long double>::set_conc_and_compute_response) 
        .def("mult_laplacian_TF",     &childscoring_ladder3<PreciseType, long double>::set_laplacian_TF) 
        .def("getLaplacianTF",        &childscoring_ladder3<PreciseType, long double>::get_laplacian_TF)
        .def("simple_score",          &childscoring_ladder3<PreciseType, long double>::SimpleScore, 
                                        py::arg("TF_range_for_score"), py::arg("n_TF_conc_points"), py::arg("low_acc"), py::arg("up_acc"), py::arg("computemeanscore")=true, py::arg("computecvscore")=true, py::arg("sqpar")=0.01)
        .def("triming",               &childscoring_ladder3<PreciseType, long double>::trim_range_response, 
                                        py::arg("TF_range"), py::arg("n_points"), py::arg("deltaTH"))
        .def("redef_score_RESP",      &childscoring_ladder3<PreciseType, long double>::redef_score_RESP, 
                                        py::arg("n_points"), py::arg("deltaTH"), py::arg("low_acc"), py::arg("up_acc"), py::arg("computemeanscore")=true, py::arg("computecvscore")=true,py::arg("sqpar")=0.01)
        .def("redef_score_noRESP",    &childscoring_ladder3<PreciseType, long double>::redef_score_noRESP, 
                                        py::arg("TF_range_noresp"), py::arg("n_points"), py::arg("deltaTH"), py::arg("low_acc"), py::arg("up_acc"), py::arg("computemeanscore")=true, py::arg("computecvscore")=true,py::arg("sqpar")=0.01)
        .def("monotonicity",          &childscoring_ladder3<PreciseType, long double>::monotonicity_test,
                                        py::arg("range_to_check"), py::arg("n_points_for_ss"), py::arg("delta_TH"));        
        //.def("score_no_response",     &childscoring_ladder3<PreciseType, long double>::scoreNoResponse, 
        //                                py::arg("TF_range_for_score"), py::arg("n_points_for_ss"), py::arg("delta_TH"), py::arg("low_acc"), py::arg("up_acc"))      
        //.def("score_response",        &childscoring_ladder3<PreciseType, long double>::scoreWithResponse,
        //                                py::arg("n_points_for_ss"), py::arg("delta_TH"), py::arg("low_acc"), py::arg("up_acc"));
}

//set_TF_and_compute_response before was compute_response in python, and Response_class in c

