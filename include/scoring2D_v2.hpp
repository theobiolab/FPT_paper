#ifndef MAIN_CLASS_HPP
#define MAIN_CLASS_HPP

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
class scoring2D_v2 {
    
    // protected members can be accessed from inheritaded classes
    protected:        
    
        Matrix<T, Dynamic, Dynamic>  laplacian;      // laplacian matrix   
        Matrix<T, Dynamic, Dynamic>  laplacian_TF; 
        Matrix<T, Dynamic, Dynamic> laplacian_red; //the reduced version, for ss calculation
        IOType gaussian_score(IOType new_point, IOType initial_point, IOType std_sq)
        {
            IOType rel_new_point = (new_point-initial_point)/initial_point;
            IOType exp_fact = -(0.5)*((rel_new_point*rel_new_point)/(std_sq));            
            return exp(exp_fact);          
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
        
        std::vector<IOType> derivative_abs(std::vector<IOType> x, std::vector<IOType> y){
            
            IOType dx = x[1]-x[0];  
            std::vector<IOType> dy_dx(x.size()-1); 
            
            for (int i=1; i<x.size(); i++){
                IOType dy = y[i]-y[i-1];
                dy_dx[i-1] = std::abs(dy/dx); 
            }   
            
            return dy_dx; 
        }
           
    public: 
        
        // declaration of the set_laplacian function to be defined in the child class
        //void set_laplacian_from_child_class(Matrix<IOType, Dynamic, Dynamic> lap); 

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

        // Getter laplacian TF red: 
        Matrix<IOType, Dynamic, Dynamic> get_laplacian_red(){
            return laplacian_red.template cast<IOType>();
        }
        
        IOType Response_class(IOType TF_conc){
            set_laplacian_TF_class(TF_conc); 
            Matrix<IOType,Dynamic,1> ss_vector = get_steady_state(true); 
            return compute_response<T,IOType>(ss_vector); 
        }

        IOType Response_class_noTF(){
            //set_laplacian_TF_class(TF_conc); 
            Matrix<IOType,Dynamic,1> ss_vector = get_steady_state(true); 
            return compute_response<T,IOType>(ss_vector); 
        }
  
        
        // Method for the calculation of mean first passage time and/or CV first passage time
        Matrix<IOType, Dynamic, Dynamic> get_Mean_CV_FPT(bool use_laplacian_TF = false, bool cv = true) { //adding option to not compute cv score to save cost 
            
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
            Matrix<T, Dynamic, Dynamic> sublaplacian_t=sublaplacian.transpose();
            
            // Get the left-hand matrix in the first-passage time linear system
            Matrix<T, Dynamic, Dynamic> A_Mu  = sublaplacian_t*sublaplacian_t;

            // Get the right-hand side vector in the first-passage time linear system
            for (int i=0; i<L.cols()-1; i++){
                b(i) = L(t, i); 
            }

            // Solve the linear system with QR decomposition 
            Matrix<T, Dynamic, 1> solution_Mu  = solveByQRD<T>(A_Mu, b);
            T mean_FPT = solution_Mu(0);
            
            T CV_FPT;
            if (cv){
                Matrix<T, Dynamic, Dynamic> A_Var = -(A_Mu*sublaplacian_t);       
                Matrix<T, Dynamic, 1> solution_sec_moment = solveByQRD<T>(A_Var, b);    
                T var_FPT  = 2*solution_sec_moment(0) - mean_FPT*mean_FPT;
                CV_FPT   = sqrt(var_FPT)/mean_FPT;
            }else{
                CV_FPT=-1;
            }
            
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

            for (int i=0; i<n_states-1; i++){
                L_red(i,i) = L_red(i,i) + L(n_states-1, i);
            }
            this -> laplacian_red = L_red;

            try{
                steady_state = getOneDimNullspaceFromSVD<T>(L_red);
            }
            catch (const std::runtime_error& e){
                throw;
            }

            for (int i=0; i<steady_state.size();i++){
                //py::print(i, "before:", steady_state[i]);
                steady_state[i]=abs(steady_state[i]);
                //py::print("after:", steady_state[i]);
            }

            T norm = steady_state.sum();
            for (int i=0; i<steady_state.size();i++){
                //std::cout << i << "," << steady_state[i];
                steady_state[i] = steady_state[i] / norm;
            }

            return steady_state.template cast<IOType>();
        }       
        
        
        // Simple mean-CV scoring
        Matrix<IOType,2,1> SimpleScore(Matrix<IOType,2,1> TF_range, IOType n_points,
                                        IOType low_acc, IOType up_acc, bool computemeanscore = true, bool computecvscore = true){
        
            IOType    low = TF_range(0); 
            IOType    up  = TF_range(1);  
            IOType    delta = abs(low-up)/(n_points-1);  
            
            // Compute the initial point mean-CV -> first TF concentration
            set_laplacian_TF_class(pow(10, low)); 
            Matrix<IOType,2,1> InitialPoint = get_Mean_CV_FPT(true, true); 
            
            // MOST OF THE POINTS ARE DISCARDED IN THIS STEP
            // Mean filter to avoid numerical artifacts due to too small/big mean values
            if ((InitialPoint(0)<low_acc)||(InitialPoint(0)>up_acc)){   
                 Matrix<IOType,2,1> neg_results; 
                 neg_results << -1.0, -1.0; 
                 return neg_results; 
            } 
            
            // compute the scores. For each TF concentration the mean and CV scores are computed.
            // only the lowerst (worst) scores are kept. 
            IOType      mean_score;     
            IOType      cv_score;
            IOType      TF_conc = low+delta; 
            int         accepted = 0;
            IOType std_sq = 0.01 * 0.01;
            IOType tmp_mean_scores =-1;
            IOType tmp_cv_scores = -1;
            Matrix<IOType,2,1> new_point;
            for (int i = 1; i<n_points; i++){
            
                // initialize the variables mean_score and cv_score
                if (i == 1 || accepted == 0){                
                    set_laplacian_TF_class(pow(10, TF_conc));
                    
                    try{ 

                        if ((computemeanscore && computecvscore)||(computecvscore)){
                    
                            new_point = get_Mean_CV_FPT(true, true);
                        
                        } else { //meanonly
                            new_point = get_Mean_CV_FPT(true, false);

                        } 

                        if (computemeanscore){ 
                            tmp_mean_scores = this->gaussian_score(new_point(0), InitialPoint(0), std_sq);
                        } else {
                            tmp_mean_scores=-1;
                        }
                        if (computecvscore) {  
                            tmp_cv_scores  = this->gaussian_score(new_point(1), InitialPoint(1), std_sq);
                        }else{
                            tmp_cv_scores=-1;

                        }
                        
                        cv_score   = tmp_cv_scores;
                        mean_score = tmp_mean_scores;
        
                        TF_conc  += delta;
                        accepted += 1;
                        continue; 
                        
                    } catch (...){
                    
                        TF_conc += delta;
                        continue; 
                        
                    }
                }
                
                // compute the scores for each other TF concentration and keep the lowerst scores
                
                set_laplacian_TF_class(pow(10, TF_conc)); 
                
                try{ 

                
                    if ((computemeanscore && computecvscore)||(computecvscore)){
                    
                        new_point = get_Mean_CV_FPT(true, true);
                        
                    } else { //meanonly
                        
                        new_point = get_Mean_CV_FPT(true, false);

                    } 


                    

                    if (computemeanscore){ 
                            tmp_mean_scores = this->gaussian_score(new_point(0), InitialPoint(0), std_sq);
                        } else {
                            tmp_mean_scores=-1;
                        }
                    if (computecvscore) {  
                            tmp_cv_scores  = this->gaussian_score(new_point(1), InitialPoint(1), std_sq);
                        }else{
                            tmp_cv_scores=-1;

                    }
                        
                   // keep only the lowerst mean scores
                    if (tmp_mean_scores<mean_score){
                        mean_score = tmp_mean_scores; 
                    } 
                    // keep only the lowest cv score
                    if (tmp_cv_scores<cv_score){
                        cv_score = tmp_cv_scores; 
                    }
                    
                    TF_conc += delta;
                    
                } catch (...){
                
                    TF_conc += delta;
                    continue; 

                }
            }  //end of the for loop          
                       
            // return the score            
            Matrix<IOType,2,1> results; 
            results << mean_score, cv_score;            
            return results; 
        }
        
        
        // compute Delta Response and trim range
        std::vector<IOType> trim_range_response(Matrix<IOType,Dynamic,1>TF_range,
                                                 const int n_points, const IOType deltaTH){


                        
        
            IOType low = TF_range(0);
            IOType up  = TF_range(1);
            IOType delta = abs(low-up)/(n_points-1.0); 
            
            // initialize a response vector and a concentration vector and compute the delta
            std::vector<IOType> resp(n_points); 
            std::vector<IOType> conc_array(n_points);             
            IOType conc = low;
            IOType min_resp = Response_class(pow(10, conc)); 
            IOType max_resp = Response_class(pow(10, conc));  
            
            for (int i=0; i<n_points; i++){
                conc_array[i] = conc; 
                resp[i] = Response_class(pow(10, conc)); 
                conc += delta; 
                
                if (resp[i] < min_resp){
                    min_resp = resp[i];
                }
                if (resp[i] > max_resp){
                    max_resp = resp[i]; 
                }
            } 
            
            IOType deltaResp = max_resp - min_resp;
            
            // if the Delta filter is passsed trim the range and return the trimmed range otherwise return the whole range
            if (deltaResp >= deltaTH){
                std::vector<IOType> dresp_dTF = this->derivative_abs(conc_array, resp); //logspace derivative
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
                
                // compute the new delta for the trimmed range
                std::vector<IOType> new_resp_range; 
                for (int i=0; i<conc_array.size(); i++){
                    if ((new_low<=conc_array[i])and(conc_array[i]<=new_up)){
                        IOType resp_tmp = resp[i];                    
                        new_resp_range.push_back(resp_tmp);   
                    }
                }
                
                IOType new_max = *max_element(new_resp_range.begin(), new_resp_range.end());
                IOType new_min = *min_element(new_resp_range.begin(), new_resp_range.end());
                IOType new_deltaResp = new_max-new_min; 
                
                std::vector<IOType> res{new_low, new_up, new_deltaResp}; 
                return res;      
            } else {
                std::vector<IOType> res{low, up, deltaResp};
                return res;  
            }
        }    

        // compute scoring with redefined scoring and trim_range_response 
        Matrix<IOType,2,1> redef_score_RESP(const int n_points, const IOType deltaTH, 
                                            const IOType low_acc, const IOType up_acc, bool computemeanscore = true, bool computecvscore = true){
        
            // initial range
            IOType low = -30.0;
            IOType up  = 30.0;
            Matrix<IOType,2,1> TF_range; 
            TF_range << low, up; 
            
            // range trimming
            std::vector<IOType> trimmed_range = trim_range_response(TF_range, n_points, deltaTH);
            
            if (trimmed_range[2] >= deltaTH){
                Matrix<IOType,2,1> new_TF_range; 
                new_TF_range << trimmed_range[0], trimmed_range[1];
                return SimpleScore(new_TF_range,15,low_acc,up_acc, computemeanscore, computecvscore);
            } else {
                Matrix<IOType,2,1> ret;
                ret << -1.0, -1.0; 
                return ret;
            }                                                                            
        }
        
        Matrix<IOType,2,1> redef_score_noRESP(Matrix<IOType,Dynamic,1>TF_range_noresp, 
                                                const int n_points, const IOType deltaTH, 
                                                const IOType low_acc, const IOType up_acc, bool computemeanscore = true, bool computecvscore = true){
            IOType low = -30.0;
            IOType up  = 30.0;
            Matrix<IOType,2,1> TF_range; 
            TF_range << low, up; 
            
            std::vector<IOType> trimmed_range = trim_range_response(TF_range, n_points, deltaTH);
            
            if (trimmed_range[2] < deltaTH){
                return SimpleScore(TF_range_noresp,15,low_acc,up_acc, computemeanscore, computecvscore);
            } else {
                Matrix<IOType,2,1> ret;
                ret << -1.0, -1.0; 
                return ret;
            }                                                                            
        }
              
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////// OLD FUNCTIONS ///////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
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
        
        // // compute scoring with no response -deltaTF<delta<deltaTH
        // Matrix<IOType,2,1> scoreNoResponse(Matrix<IOType,Dynamic,1>TF_range_score,
        //                                       const int n_points, const IOType deltaTH, 
        //                                       const IOType low_acc, const IOType up_acc){
                                                                      
        //     IOType low = -20.0;
        //     IOType up  = 20.0;
        //     IOType delta = abs(low-up)/(n_points-1.0); 

        //     //evaluate delta response, return [-1,-1] if "there is a response" in the range considered; 
        //     IOType resp_low = Response_class(pow(10, low));
        //     IOType resp_up  = Response_class(pow(10, up));
        //     IOType deltaResp = resp_up - resp_low;
                       
        //     if ((deltaResp >= deltaTH)||(deltaResp <= -deltaTH)){
        //         Matrix<IOType,2,1> ret;
        //         ret << -1.0, -1.0; 
        //         return ret;   
        //     }
            
        //     // initialize a response vector and a concentration vector
        //     std::vector<IOType> resp(n_points); 
        //     std::vector<IOType> conc_array(n_points); 
        //     IOType conc = low;
        //     for (int i=0; i<n_points; i++){
        //         conc_array[i] = conc; 
        //         resp[i] = Response_class(pow(10, conc)); 
        //         conc += delta; 
        //     }

        //     // Evaluate if the response vector is monotonic, return [-1,-1] if non-monotonic
        //     IOType previous_diff = resp[1]-resp[0]; 
        //     for (int i=2; i<n_points; i++){
        //         IOType diff = resp[i]-resp[i-1]; 

        //         if (diff!=0){
        //             if (not((diff>0)and(previous_diff>=0))){
        //                 if (not((diff<0)and(previous_diff<=0))){
        //                     Matrix<IOType,2,1> ret; 
        //                     ret << -1.0,-1.0; 
        //                     return ret;
        //                 }
        //             }
        //         }                
        //         previous_diff = diff; 
        //     }
            
        //     // costume TF_range 
        //     return SimpleScore(TF_range_score, 15, low_acc,up_acc);         
        // }
        
        
        // // compute scoring with response -deltaTF<delta<deltaTH
        // Matrix<IOType,2,1> scoreWithResponse(const int n_points, const IOType deltaTH, 
        //                                      const IOType low_acc, const IOType up_acc){
                                              
        //     IOType low = -20.0;
        //     IOType up  = 20.0;
        //     IOType delta = abs(low-up)/(n_points-1.0); 
        //     Matrix<IOType,2,1> TF_range; 
        //     TF_range << low, up; 
            
        //     //evaluate delta response, return [-1,-1] if "there is no response" in the range considered; 
        //     IOType resp_low = Response_class(pow(10, low));
        //     IOType resp_up  = Response_class(pow(10, up));
        //     IOType deltaResp = resp_up - resp_low;
        //     if ((deltaResp <= deltaTH)and(deltaResp >= -deltaTH)){ 
        //         Matrix<IOType,2,1> ret;
        //         ret << -1.0, -1.0; 
        //         return ret;
        //     }
            
        //     // trim the concentration range
        //     std::vector<IOType> trimmed_range = monotonicity_test(TF_range, n_points, deltaTH);
            
        //     // check for monotonicity, if non-monotonic return [-1,-1]
        //     if (trimmed_range.size()==1){
        //         Matrix<IOType,2,1> ret;
        //         ret << -1.0, -1.0; 
        //         return ret;
        //     }          
            
        //     // compute the score with the trimmed_range 
        //     Matrix<IOType,2,1> new_TF_range; 
        //     new_TF_range << trimmed_range[0], trimmed_range[1]; 
            
        //     // double check that there is a delta response after the trimming
            
        //     return SimpleScore(new_TF_range,15,low_acc,up_acc);                                                                         
        // }
    
}; 
#endif 





