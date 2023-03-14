#ifndef MAIN_CLASS_HPP
#define MAIN_CLASS_HPP

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>   //return matrices in python
#include <pybind11/embed.h>
#include <pybind11/stl.h> 

#include <iostream>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>

#include "KBNSum.hpp"
#include "kahan.hpp"
#include "math.hpp"

using namespace Eigen;
namespace py=pybind11;

/*
* This header contains a base class for performing the computations of the steady state
* input-output response and the first passage time scorings 
*
*
*
*/

template <typename T, typename IOType>
class scoring2D {
    
    // protected members can be accessed from inheritaded classes
    protected:        
    
        Matrix<T, Dynamic, Dynamic>  laplacian;      // laplacian matrix   
        Matrix<T, Dynamic, Dynamic>  laplacian_TF; 
        
        py::module scipy = py::module::import("scipy.interpolate");
        py::function interp1d = scipy.attr("interp1d"); 
        
        IOType gaussian_score(IOType new_point, IOType initial_point, IOType std)
        {
            IOType rel_new_point = (new_point-initial_point)/initial_point;
            IOType exp_fact = -(1.0/2.0)*((rel_new_point*rel_new_point)/(std*std));            
            return exp(exp_fact);          
        }
        
        IOType gaussian_2D_score(IOType new_point_mean, IOType initial_point_mean, 
                                 IOType new_point_CV, IOType initial_point_CV,
                                 IOType std)
        {
            IOType rel_MEAN = (new_point_mean-initial_point_mean)/initial_point_mean;
            IOType rel_CV   = (new_point_CV-initial_point_CV)/initial_point_CV;
            
            IOType exp_fact_MEAN = (1.0/2.0)*((rel_MEAN*rel_MEAN)/(std*std));
            IOType exp_fact_CV   = (1.0/2.0)*((rel_CV*rel_CV)/(std*std));
           
            return exp(-(exp_fact_MEAN+exp_fact_CV));          
        }
        
        std::vector<IOType> derivative(std::vector<IOType> x, std::vector<IOType> y){
            
            IOType dx = x[1]-x[0];  
            std::vector<IOType> dy_dx(x.size()-1); 
            
            for (int i=1; i<x.size(); i++){
                IOType dy = y[i]-y[i-1];
                dy_dx[i-1] = dy/dx; 
            }   
            
            return dy_dx; 
        }
           
    public: 
        
        // declaration of the set_laplacian function to be defined in the child class
        void set_laplacian_from_child_class(Matrix<IOType, Dynamic, Dynamic> lap); 
        
        // declaration of the compute_response function to be defined in the child class
        IOType compute_response_class(Matrix<IOType,Dynamic,1> ss_vector){
            return compute_response<T,IOType>(ss_vector); 
        }
        
        // declaration of the laplacian_TF to be defined in the child class
        void set_laplacian_TF_class(IOType TF_conc){
            this->laplacian_TF = set_laplacian_TF<T,IOType>(TF_conc, this->laplacian); 
        }
    
        // Getter for the laplacian: 
        Matrix<IOType, Dynamic, Dynamic> get_laplacian(){
            return laplacian.template cast<IOType>();
        }
        
        // Getter laplacian TF: 
        Matrix<IOType, Dynamic, Dynamic> get_laplacian_TF(){
            return laplacian_TF.template cast<IOType>();
        }
        
        IOType Response_class(IOType TF_conc){
            set_laplacian_TF_class(TF_conc); 
            Matrix<IOType,Dynamic,1> ss_vector = get_steady_state(true); 
            return compute_response<T,IOType>(ss_vector); 
        }
  
        
        // Method for the calculation of mean first passage time and CV first passage time
        Matrix<IOType, Dynamic, Dynamic> get_Mean_CV_FPT(bool use_laplacian_TF = false) {
            
            Matrix<T,Dynamic,Dynamic> L; 
            
            //define L with an if 
            if (use_laplacian_TF){
                L = this->laplacian_TF; 
            } else {
                L = this->laplacian; 
            }
             
            int t = L.rows()-1; 
            
            Matrix<T, Dynamic, Dynamic> sublaplacian = L.block(0, 0, t, t);
            Matrix<T, Dynamic, 1> b = Matrix<T, Dynamic, 1>::Zero(t);

            // Get the left-hand matrix in the first-passage time linear system
            Matrix<T, Dynamic, Dynamic> A_Mu  = sublaplacian.transpose()*sublaplacian.transpose();
            Matrix<T, Dynamic, Dynamic> A_Var = -(sublaplacian.transpose()*sublaplacian.transpose()*sublaplacian.transpose());       

            // Get the right-hand side vector in the first-passage time linear system
            for (int i=0; i<L.cols()-1; i++){
                b(i) = L(t, i); 
            }

            // Solve the linear system with QR decomposition 
            Matrix<T, Dynamic, 1> solution_Mu  = solveByQRD<T>(A_Mu, b);  
            Matrix<T, Dynamic, 1> solution_sec_moment = solveByQRD<T>(A_Var, b);    
            
            T mean_FPT = solution_Mu(0);
            T var_FPT  = 2*solution_sec_moment(0) - mean_FPT*mean_FPT;
            T CV_FPT   = sqrt(var_FPT)/mean_FPT;
            
            Matrix<T, 2, 1> output;
            output(0) = mean_FPT;
            output(1) = CV_FPT;
            
            // return a vector with the mean first passage time and the CV   
            return output.template cast<IOType>();  
        }
        
        
        // Get steady state
        
        Matrix<IOType, Dynamic, 1> get_steady_state(bool use_laplacian_TF = false){
            
            Matrix<T,Dynamic,Dynamic> L;
            
            if (use_laplacian_TF){
                L = this->laplacian_TF; 
            } else {
                L = this->laplacian; 
            }

            int n_states = L.rows(); 
            Matrix<T, Dynamic, 1> steady_state;
            Matrix<T, Dynamic, Dynamic> L_red = L.block(0, 0, n_states-1, n_states-1);  

            // add the last two elements of L to the corresponding diagonal elements of L_red

            for (int i; i<n_states-1; i++)
                L_red(i,i) = L_red(i,i) + L(n_states-1, i);


            try{
                steady_state = getOneDimNullspaceFromSVD<T>(L_red);
            }
            catch (const std::runtime_error& e){
                throw;
            }

            T norm = steady_state.sum();
            for (int i=0; i<steady_state.size();i++){
                steady_state[i] = steady_state[i] / norm;
            }

            return steady_state.template cast<IOType>();
        }       
        
        
        // Simple mean-CV scoring
        Matrix<IOType,2,1> SimpleScore(Matrix<IOType,2,1> TF_range, IOType n_points){
        
            IOType    low = TF_range(0); 
            IOType    up  = TF_range(1);  
            IOType    delta = abs(low-up)/(n_points-1);  
            
            // Compute the initial point mean-CV -> first TF concentration
            set_laplacian_TF_class(pow(10, low)); 
            Matrix<IOType,2,1> InitialPoint = get_Mean_CV_FPT(true);
            
            // Filter for mean values in our range of interest  
            // MOST OF THE POINTS ARE DISCARDED IN THIS STEP
            // if (InitialPoint(0)<400 || InitialPoint(0)>5000){
            
            if ((InitialPoint(0)<200)||(InitialPoint(0)>5000)){   // just filter to avoid artifacts due to low mean values 
                 Matrix<IOType,2,1> neg_results; 
                 neg_results << -1.0, -1.0; 
                 return neg_results; 
            } 
            
            // compute the scores
            IOType      mean_score;             
            IOType TF_conc = low+delta; 
            int    accepted = 0;
            
            for (int i = 1; i<n_points; i++){
            
                // initialize the variable mean_score
                if (i == 1 || accepted == 0){                
                    set_laplacian_TF_class(pow(10, TF_conc));
                    
                    try{ 
                
                        Matrix<IOType,2,1> new_point = get_Mean_CV_FPT(true);       
                        IOType tmp_mean_scores = this->gaussian_score(new_point(0), InitialPoint(0), 0.01);
                        
                        mean_score = tmp_mean_scores;
        
                        TF_conc  += delta;
                        accepted += 1;
                        continue; 
                        
                    } catch (...){
                
                        TF_conc += delta;
                        continue; 
                    }
                }
                
                // compute the scores for each other TF concentration and keep the lowerst mean_score
                
                set_laplacian_TF_class(pow(10, TF_conc)); 
                
                try{ 
                
                    Matrix<IOType,2,1> new_point = get_Mean_CV_FPT(true);       
                    IOType tmp_mean_scores = this->gaussian_score(new_point(0), InitialPoint(0), 0.01); 
                    
                    // keep only the lowerst mean scores
                    if (tmp_mean_scores<mean_score){
                        mean_score = tmp_mean_scores; 
                        TF_conc += delta;
                    } else {
                        TF_conc += delta;
                        continue;
                    }
                    
                } catch (...){
                
                    TF_conc += delta;
                    continue; 

                }
            }  //end of the for loop          
            
            
            // compute CV score            
            set_laplacian_TF_class(pow(10, up));
            Matrix<IOType,2,1> new_point = get_Mean_CV_FPT(true);       
            IOType CV_score = this->gaussian_score(new_point(1), InitialPoint(1), 0.05);

            // check that at least one mean_score as been accepted
            if (accepted==0){ 
                Matrix<IOType,2,1> neg_results; 
                neg_results << -1.0, -1.0; 
                return neg_results;  
            }
            
            // return the score            
            Matrix<IOType,2,1> results; 
            results << mean_score, CV_score;            
            return results; 
        }
        
        
        // Monotonicity check and trimmed range_calculation
        std::vector<IOType> monotonicity_test(Matrix<IOType,Dynamic,1>TF_range,
                                              const int n_points, const IOType deltaTH){
        
            IOType low = TF_range(0);
            IOType up  = TF_range(1);
            IOType delta = abs(low-up)/(n_points-1.0); 
            
            // initialize a response vector and a concentration vector
            std::vector<IOType> resp(n_points); 
            std::vector<IOType> conc_array(n_points); 
            IOType conc = low;
            for (int i=0; i<n_points; i++){
                conc_array[i] = conc; 
                resp[i] = Response_class(pow(10, conc)); 
                conc += delta; 
            }
            
            // Evaluate if the response vector is monotonic, return 0.0 if it is non-monotonic
            IOType previous_diff = resp[1]-resp[0]; 
            for (int i=2; i<n_points; i++){
                IOType diff = resp[i]-resp[i-1]; 
                
                if (diff!=0){
                    if (not((diff>0)and(previous_diff>=0))){
                        if (not((diff<0)and(previous_diff<=0))){
                            std::vector<IOType> res{0.0}; 
                            return res;
                        }
                    }
                }                
                previous_diff = diff; 
            }

            // calculate the delta_response
            IOType deltaResp = resp[n_points-1] - resp[0];

            // check response delta is greater(smaller) than deltaTH(-deltaTH)
            if ((deltaResp >= deltaTH)||(deltaResp <= -deltaTH)){
                std::vector<IOType> dresp_dTF = this->derivative(conc_array, resp); //logspace derivative
                
                IOType new_low;  
                IOType new_up; 
                bool   min_acc = false; 
                
                for (int i=1; i<dresp_dTF.size(); i++){
                    if ((abs(dresp_dTF[i])>1e-5)and(abs(dresp_dTF[i-1])<=1e-5)){
                        if (!min_acc){
                            new_low = conc_array[i+1];
                            min_acc = true; 
                        }
                    } else if ((abs(dresp_dTF[i-1])>1e-5)and(abs(dresp_dTF[i])<=1e-5)){
                        new_up  = conc_array[i+1];
                    } else {
                        continue;                     
                    }
                }
                
                std::vector<IOType> res{new_low, new_up, deltaResp}; 
                return res;    
            } else {
                std::vector<IOType> res{low, up, deltaResp};
                return res;    
            }   
        }
        
        
        // compute scoring with no response -deltaTF<delta<deltaTH
        Matrix<IOType,2,1> scoreNoResponse(Matrix<IOType,Dynamic,1>TF_range_score,
                                              const int n_points, const IOType deltaTH){
                                                                      
            IOType low = -20.0;
            IOType up  = 20.0;
            IOType delta = abs(low-up)/(n_points-1.0); 

            //evaluate delta response, return [-1,-1] if "there is a response" in the range considered; 
            IOType resp_low = Response_class(pow(10, low));
            IOType resp_up  = Response_class(pow(10, up));
            IOType deltaResp = resp_up - resp_low;
                       
            if ((deltaResp >= deltaTH)||(deltaResp <= -deltaTH)){
                Matrix<IOType,2,1> ret;
                ret << -1.0, -1.0; 
                return ret;   
            }
            
            // initialize a response vector and a concentration vector
            std::vector<IOType> resp(n_points); 
            std::vector<IOType> conc_array(n_points); 
            IOType conc = low;
            for (int i=0; i<n_points; i++){
                conc_array[i] = conc; 
                resp[i] = Response_class(pow(10, conc)); 
                conc += delta; 
            }

            // Evaluate if the response vector is monotonic, return [-1,-1] if non-monotonic
            IOType previous_diff = resp[1]-resp[0]; 
            for (int i=2; i<n_points; i++){
                IOType diff = resp[i]-resp[i-1]; 

                if (diff!=0){
                    if (not((diff>0)and(previous_diff>=0))){
                        if (not((diff<0)and(previous_diff<=0))){
                            Matrix<IOType,2,1> ret; 
                            ret << -1.0,-1.0; 
                            return ret;
                        }
                    }
                }                
                previous_diff = diff; 
            }
            
            // costume TF_range 
            return SimpleScore(TF_range_score, 15);         
        }
        
        
        // compute scoring with response -deltaTF<delta<deltaTH
        Matrix<IOType,2,1> scoreWithResponse(const int n_points, const IOType deltaTH){
                                              
            IOType low = -20.0;
            IOType up  = 20.0;
            IOType delta = abs(low-up)/(n_points-1.0); 
            Matrix<IOType,2,1> TF_range; 
            TF_range << low, up; 
            
            //evaluate delta response, return [-1,-1] if "there is no response" in the range considered; 
            IOType resp_low = Response_class(pow(10, low));
            IOType resp_up  = Response_class(pow(10, up));
            IOType deltaResp = resp_up - resp_low;
            if ((deltaResp <= deltaTH)and(deltaResp >= -deltaTH)){ 
                Matrix<IOType,2,1> ret;
                ret << -1.0, -1.0; 
                return ret;
            }
            
            // trim the concentration range
            std::vector<IOType> trimmed_range = monotonicity_test(TF_range, n_points, deltaTH);
            
            // check for monotonicity, if non-monotonic return [-1,-1]
            if (trimmed_range.size()==1){
                Matrix<IOType,2,1> ret;
                ret << -1.0, -1.0; 
                return ret;
            }          
            
            // compute the score with the trimmed_range 
            Matrix<IOType,2,1> new_TF_range; 
            new_TF_range << trimmed_range[0], trimmed_range[1]; 
            
            // double check that there is a delta response after the trimming
            
            return SimpleScore(new_TF_range,15);                                                                         
        }
        

        // compute 2D gaussian scoring
        IOType helper_2Dgauss_score(Matrix<IOType,Dynamic,1> TF_range, IOType n_points){
            
            IOType    low = TF_range(0); 
            IOType    up  = TF_range(1);  
            IOType    delta = abs(low-up)/(n_points-1);  
            
            // Compute the initial point mean-CV -> first TF concentration
            set_laplacian_TF_class(pow(10, low)); 
            Matrix<IOType,2,1> InitialPoint = get_Mean_CV_FPT(true);
            
            // Filter for mean values in our range of interest  
            // MOST OF THE POINTS ARE DISCARDED IN THIS STEP
            //if (InitialPoint(0)<400 || InitialPoint(0)>5000){
            if (InitialPoint(0)<200){ 
                 return -1.0;  
            }
            
            // compute the scores
            IOType scores;             
            IOType TF_conc=low+delta; 
            int    accepted=0;
            
            for (int i = 1; i<n_points; i++){
            
                // initialize the variable mean_score
                if (i == 1 || accepted == 0){                
                    set_laplacian_TF_class(pow(10, TF_conc));
                    
                    try{ 
                
                        Matrix<IOType,2,1> new_point = get_Mean_CV_FPT(true);       
                        IOType tmp_scores = this->gaussian_2D_score(new_point(0), InitialPoint(0), 
                                                                    new_point(1), InitialPoint(1), 0.01);
                        
                        scores = tmp_scores;
        
                        TF_conc  += delta;
                        accepted += 1;
                        continue; 
                        
                    } catch (...){
                
                        TF_conc += delta;
                        continue; 
                    }
                }
                
                // compute the scores for each other TF concentration and keep the lowerst mean_score
                
                set_laplacian_TF_class(pow(10, TF_conc)); 
                
                try{ 
                
                    Matrix<IOType,2,1> new_point = get_Mean_CV_FPT(true);       
                    IOType tmp_scores = this->gaussian_2D_score(new_point(0), InitialPoint(0), 
                                                                new_point(1), InitialPoint(1), 0.01);
                    
                    // keep only the lowerst mean scores
                    if (tmp_scores<scores){
                        scores = tmp_scores; 
                        TF_conc += delta;
                    } else {
                        TF_conc += delta;
                        continue;
                    }
                    
                } catch (...){
                
                    TF_conc += delta;
                    continue; 

                }
            }  //end of the for loop 
            
            // check that at least one mean_score as been accepted
            if (accepted==0){ 
                return -1.0; 
            }
            
            // return the gaussian score
            return scores;      
        }
        
        
        // compute 2D gaussian - delta response scoring
        Matrix<IOType,2,1> score2DgaussianDelta(Matrix<IOType,Dynamic,1> starting_TF_range,
                                                      Matrix<IOType,Dynamic,1> TF_range_no_response, 
                                                      IOType n_points_score, IOType n_points_ss, IOType delta_TH){
                                                      
            // compute the delta response, trim TF concentration range, check for monotonicity
            std::vector<IOType> new_range_mono = monotonicity_test(starting_TF_range, n_points_ss, delta_TH);
            
            // return [-1,-1] if non monotonic
            if (new_range_mono.size()==1){
                Matrix<IOType,2,1> ret; 
                ret << -1.0, -1.0; 
                return ret; 
            }
        
            // compute scoring depending on the response delta 
            if ((new_range_mono[2] < delta_TH)and(new_range_mono[2] > -delta_TH)){
            
                // compute gaussian scoring with TF_range_no_response
                Matrix<IOType,2,1> res; 
                
                IOType gauss = helper_2Dgauss_score(TF_range_no_response, n_points_score);  
                
                if (gauss<0){
                    res << -1.0, -1.0; 
                    return res; 
                } else {
                    res << gauss, new_range_mono[2]; 
                    return res;
                }
                
            } else if ((new_range_mono[2] >= delta_TH)||(new_range_mono[2] <= -delta_TH)){
            
                // compute gaussian scoring with new_range_mono
                Matrix<IOType,2,1> new_range_TF; 
                Matrix<IOType,2,1> res; 
                
                new_range_TF << new_range_mono[0], new_range_mono[1]; 
                IOType gauss = helper_2Dgauss_score(new_range_TF, n_points_score);
                
                if (gauss<0){
                    res << -1.0, -1.0; 
                    return res; 
                } else {
                    res << gauss, new_range_mono[2]; 
                    return res;
                }
            }
        }      
}; 

#endif 




















