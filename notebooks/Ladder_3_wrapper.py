import numpy as np
import sys, os
sys.path.insert(0, "../bin")
from Ladder_3 import Ladder_3


class L3(Ladder_3):
    """ 
    Child class with decorators to make the code easier to understand
    """

    def setLaplacian_ladder_3(self, converted_pset: np.ndarray): 
        """
        Given a an array of parameters set the laplacian matrix
            
            Parameters: 
            
                converted_pset: np.ndarray parameters converted following the
                                Ladder_3 scheme.  
        
        """ 
        super().setLaplacian_ladder_3(converted_pset)
    
    
    def mult_laplacian_TF(self, tf_conc: float): 
        """
        Setter to multiply laplacian matrix by TF concentration (the laplacian is not overwritten)
            
            Parameters: 
            
                tf_conc: the concentration of transcription factor that will 
                         be the multiplicative factor for the k_on (see implementation)
                         
        """
        super().mult_laplacian_TF(tf_conc)
        
        
    def getLaplacian(self) -> np.ndarray: 
        "Getter for the laplacian matrix"
        return super().getLaplacian()
    
    
    def getLaplacianTF(self) -> np.ndarray: 
        "Getter for the laplacian matrix multiplied by the TF concentration"
        return super().getLaplacianTF()
    
    
    def compute_response(self, tf_conc: float) -> float:
        """
        Given a TF concentration compute the steady state response
            
            Parameters: 
            
                tf_conc: TF concentration (the mult_laplacian_TF step is done authomatically)
            
            Returns: 
            
                A float representing the steady state response value at the provided TF concentration
                see implementation for the definition of steady state response
            
        """
        return super().compute_response(tf_conc)
    
    
    def get_FPT_stat(self, use_TF_mult_laplacian: bool) -> np.ndarray:
        """
        Compute FPT stat (mean FPT and CV FPT)
        
            Parameters: 
            
                use_TF_mult_laplacian: if true the laplacian multiplied by the TF concentration
                                       is used for the computation of the FPT statistics
                                       
            returns: 
            
                np.ndarray with the mean FPT and the CV of the FPT distribution between state 1 and 
                state N in the model
        
        """
        return super().get_FPT_stat(use_TF_mult_laplacian)
    
    
    def getSteadyState(self, use_TF_mult_laplacian: bool) -> np.ndarray: 
        """
        Compute steady state vector
        
            Parameters: 
            
                use_TF_mult_laplacian: if true the laplacian multiplied by the TF concentration
                                       is used for the computation of the steady state vector
                                       
            returns: 
            
                np.ndarray with the steady state probability distribution across 
                all the states of the system
        
        """
        
        return super().getSteadyState(use_TF_mult_laplacian)
    
    
    def triming(self, initial_tf_range: list, n_points: int, delta_response_TH: float) -> np.ndarray: 
        """
        Trim TF concentration range
        
            Parameters: 
            
                initial_tf_range: initial TF concentration range where to compute 
                                  the steady state response vector
                
                n_points: length of the concentration array
                
                delta_response_TH: dynamic range treshold for the selection,
                                   if |Delta Response| > delta_response_TH the trimming is preformed, 
                                   otherwise initial_tf_range is returned
                                   
            Returns: 
            
                np.ndarray with the lower and upper bounds of the new range and the absolute
                value of the dynamic range.
                
        """
        return super().triming(initial_tf_range, n_points, delta_response_TH)
    
    
    def simple_score(self, tf_range: list, n_points: int, low_mean: float, up_mean: float):
        """
        Compute the simple gaussian scoring on mean and CV of FPT within the TF conc range defined. 
        
            Parameters:
            
                tf_range: on which the score is computed
                
                n_points: length of the TF concentration array 
                
                low_mean: points with initial mean FPT (mean FPT at the lowest concentration) 
                          below low_mean are not used for the score computation
                
                up_mean: points with initial mean FPT (mean FPT at the lowest concentration) 
                         above up_mean are not used for the score computation
        
            Returns: 
            
                [-1, -1]: if the initial mean FPT is lower than low_mean or greather than up_mean, 
                          or if an error occurs during the computation. 
                          
                [mean_score, CV_score]: if the selection filters are passed
        """
        return super().simple_score(tf_range, n_points, low_mean, up_mean)
    
   