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
void set_ladder3_laplacian(Matrix<IOType,12,1> parset, IOType TF_conc, Matrix<T,Dynamic, Dynamic> &laplacian){ 
      //adding TF_conc and modifying laplacian attribute directly, passed by reference. Since this function is made for each graph, it is ok to specify the matrix.
    
        ///error for the lenght of partset
    
        
        
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
        T kon_TF = kon*TF_conc;
        
        laplacian <<     0,    b12,     0,      koff,        0,       0,     0, 
                 a12,      0,   b23,          0,    koff,       0,     0, 
                   0,    a23,     0,          0,       0,    koff,     0,
                 kon_TF,      0,     0,          0,    b12T,       0,     0, 
                   0,    kon_TF,     0,       a12T,       0,    b23T,     0,
                   0,      0,   kon_TF,          0,    a23T,       0,     0,
                   0,      0,     P,          0,       0,      PT,     0;
               
        for (int i=0; i<7; i++){
            laplacian(i, i) = -(laplacian.col(i).sum());
        }

        //laplacian = lap; 
    } 

template <typename T, typename IOType>
void modify_TF_conc_ladder3(T TF_conc, const Matrix<T,Dynamic,Dynamic> &L, Matrix<T,Dynamic,Dynamic> &laplacian_mult ){
        //laplacian_mult is laplacian_TF attribute of the class
        std::cout << "inside" << "\n";
        std::cout.flush();

        int n_states = 7; //L.cols();
        int n_states_line = 3; //(n_states-1)/2;
         
        laplacian_mult = L; 
        T TF_conc_casted = (T) TF_conc;

        for (int i = n_states_line; i<n_states-1; i++){
            laplacian_mult(i,i-n_states_line) = laplacian_mult(i,i-n_states_line)*TF_conc_casted; 
        }
        
        
        for (int i=0; i<n_states; i++){
            laplacian_mult(i,i) = 0.0; 
            laplacian_mult(i,i) = -(laplacian_mult.col(i).sum()); 
        }
            
        //return laplacian_mult.template cast<T>(); 
    }  


template<typename T, typename IOType>
class Base {
    public: 
        
        Matrix<T, Dynamic, Dynamic> laplacian; //laplacian matrix   
        Matrix<T, Dynamic, Dynamic>  laplacian_TF; 
        Matrix<T, Dynamic, Dynamic> laplacian_red; //the reduced version, for ss calculation

        virtual void modify_TF_conc_class(T TF_conc)=0; //virtual function so that it must be defined in derived classes
        void test_modify_TF_conc(){
            //test function that uses the modify_TF_conc_class function defined in child class
            //https://stackoverflow.com/questions/4869216/calling-derived-class-function-from-base-class
            std::cout << "testing" << "\n";
            std::cout.flush();
            std::cout << "before" << "\n";
            std::cout << this->laplacian_TF(3,0) << "\n";
            modify_TF_conc_class(10);
            std::cout << "after" << "\n";
            std::cout << this->laplacian_TF(3,0)<< "\n";
            std::cout.flush();
        }

        //// Getter for the laplacian: 
        Matrix<IOType, Dynamic, Dynamic> get_laplacian(){
        ////void get_laplacian(){
            return this->laplacian.template cast<IOType>();
        }
        //
        //// Getter laplacian TF: 
        Matrix<IOType, Dynamic, Dynamic> get_laplacian_TF(){
            return this->laplacian_TF.template cast<IOType>();
        }
        
       
};

template <typename T, typename IOType> 
class Child: public Base<T, IOType>{
    public:
    Child(){ //initialise the laplacian and laplaican_TF matrices so they can be passed by reference and modified when desired
        this->laplacian.resize(7,7);
        this->laplacian_TF.resize(7,7);
    }

    void set_laplacian_class(Matrix<IOType,Dynamic,1> pset, IOType TF_conc = 1){
           set_ladder3_laplacian<T, IOType>( pset, TF_conc, this->laplacian);
           std::cout << "Printing laplacian element after setting laplacian\n";
           std::cout << (IOType) this->laplacian(0,0) << "\n";
           std::cout.flush();
       }
    void modify_TF_conc_class(T TF_conc){
            modify_TF_conc_ladder3<T, IOType>(TF_conc, this->laplacian, this->laplacian_TF);
    }
};
//    this->laplacian(0,1)=2;
    //this->laplacian_TF=Matrix<T, 7, 7>::Zero();
//    public:
       
//      
//};

const int MANTISSA_PRECISION = 20;
typedef boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<MANTISSA_PRECISION>> PreciseType;


PYBIND11_MODULE(test1, m){
    py::class_<Child<PreciseType, long double>>(m, "test_class1", py::module_local())
        .def(py::init())
        .def("setLaplacian", &Child<PreciseType, long double>::set_laplacian_class, py::arg("pset"), py::arg("TF_conc")=1)
        .def("call_test", &Child<PreciseType, long double>::test_modify_TF_conc)
        .def("getLaplacian", &Child<PreciseType, long double>::get_laplacian)
        //.def("print_laplacian_element", &childscoringv2<PreciseType, long double>::print_laplacian_element)
        .def("getLaplacian_TF", &Child<PreciseType, long double>::get_laplacian_TF);
}

