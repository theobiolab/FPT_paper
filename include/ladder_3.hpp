#ifndef LADDER_3_H
#define LADDER_3_H


/* This header file must be included in the superclass scoring2D file
*  In this header are defined this two functions: 
*       1) IOType compute_response(double TF_conc);
*       2) void set_laplacian_TF(IOType TF_conc);
*/


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

#include "../scoring2D.hpp"

# endif // LADDER_3_H