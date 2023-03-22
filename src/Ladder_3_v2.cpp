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
    response = ss_vector[2] + ss_vector[5];
    return response;

}

#include "../include/scoring2D_v2.hpp"

const int MANTISSA_PRECISION = 100;
typedef boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<MANTISSA_PRECISION>> PreciseType;

template <typename T, typename IOType> 
class childscoringv2: public scoring2D_v2<T, IOType>{
    
    public:
    
    void set_ladder_3_laplacian(Matrix<IOType,Dynamic,1> pset){
    
        ///error for the lenght of partset
    
        Matrix<T, 7, 7> lap; 
        
        Matrix<T,Dynamic,1> parset = pset.template cast<T>();  
        
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

        this->laplacian = lap; 
    }  

    void set_ladder_3_ss_laplacian_TF(Matrix<IOType,Dynamic,1> pset, long TF_conc){
    
        ///laplacian with TF concentration, for ss calculation. Keeping the last row and column so the rest of the code works well
    
        Matrix<T, 7, 7> lap; 
        
        Matrix<T,Dynamic,1> parset = pset.template cast<T>();  
        
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
                 kon*TF_conc,      0,     0,          0,    b12T,       0,     0, 
                   0,    kon*TF_conc,     0,       a12T,       0,    b23T,     0,
                   0,      0,   kon*TF_conc,          0,    a23T,       0,     0,
                   0,      0,     0,          0,       0,      0,     0;
               
        for (int i=0; i<7; i++){
            lap(i, i) = -(lap.col(i).sum());
        }

        this->laplacian_TF = lap; 
    }       
};



PYBIND11_MODULE(Ladder_3_v2, m){
    py::class_<childscoringv2<PreciseType, long double>>(m, "Ladder_3_v2", py::module_local())
        .def(py::init())
        .def("getLaplacian",          &childscoringv2<PreciseType, long double>::get_laplacian)
        .def("getLaplacianred",          &childscoringv2<PreciseType, long double>::get_laplacian_red)
        .def("setLaplacian_ladder_3", &childscoringv2<PreciseType, long double>::set_ladder_3_laplacian)
        .def("setLaplacian_ladder_3_TF", &childscoringv2<PreciseType, long double>::set_ladder_3_ss_laplacian_TF)
        .def("get_FPT_stat",          &childscoringv2<PreciseType, long double>::get_Mean_CV_FPT)
        .def("getSteadyState",        &childscoringv2<PreciseType, long double>::get_steady_state)
        .def("compute_response",      &childscoringv2<PreciseType, long double>::Response_class) 
        .def("compute_response_noTF",      &childscoringv2<PreciseType, long double>::Response_class_noTF)  
        .def("mult_laplacian_TF",     &childscoringv2<PreciseType, long double>::set_laplacian_TF_class) 
        .def("getLaplacianTF",        &childscoringv2<PreciseType, long double>::get_laplacian_TF)
        .def("simple_score",          &childscoringv2<PreciseType, long double>::SimpleScore, 
                                        py::arg("TF_range_for_score"), py::arg("n_TF_conc_points"), py::arg("low_acc"), py::arg("up_acc"), py::arg("computemeanscore")=true, py::arg("computecvscore")=true)
        .def("triming",               &childscoringv2<PreciseType, long double>::trim_range_response, 
                                        py::arg("TF_range"), py::arg("n_points"), py::arg("deltaTH"))
        .def("redef_score_RESP",      &childscoringv2<PreciseType, long double>::redef_score_RESP, 
                                        py::arg("n_points"), py::arg("deltaTH"), py::arg("low_acc"), py::arg("up_acc"), py::arg("computemeanscore")=true, py::arg("computecvscore")=true)
        .def("redef_score_noRESP",    &childscoringv2<PreciseType, long double>::redef_score_noRESP, 
                                        py::arg("TF_range_noresp"), py::arg("n_points"), py::arg("deltaTH"), py::arg("low_acc"), py::arg("up_acc"), py::arg("computemeanscore")=true, py::arg("computecvscore")=true)
        .def("monotonicity",          &childscoringv2<PreciseType, long double>::monotonicity_test,
                                        py::arg("range_to_check"), py::arg("n_points_for_ss"), py::arg("delta_TH"));        
        //.def("score_no_response",     &childscoringv2<PreciseType, long double>::scoreNoResponse, 
        //                                py::arg("TF_range_for_score"), py::arg("n_points_for_ss"), py::arg("delta_TH"), py::arg("low_acc"), py::arg("up_acc"))      
        //.def("score_response",        &childscoringv2<PreciseType, long double>::scoreWithResponse,
        //                                py::arg("n_points_for_ss"), py::arg("delta_TH"), py::arg("low_acc"), py::arg("up_acc"));
}



