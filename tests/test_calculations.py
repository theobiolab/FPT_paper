import sympy as sym
import numpy as np
import logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
from sympy import *
import pickle
import tqdm

import sys, os
sys.path.insert(0, "../bin")
from Ladder_3_v3 import Ladder_3_v3
from Ladder_6_prec_100_v2 import Ladder_6_prec_100_v2
from Ladder_4 import Ladder_4
from Ladder_5 import Ladder_5

from sympy import N as prec
# obj = Ladder_3()

a12 = sym.Symbol("a_12")
a23 = sym.Symbol("a_23")
a34 = sym.Symbol("a_34")
a45 = sym.Symbol("a_45")
a56 = sym.Symbol("a_56")

b12 = sym.Symbol("b_12")
b23 = sym.Symbol("b_23")
b34 = sym.Symbol("b_34")
b45 = sym.Symbol("b_45")
b56 = sym.Symbol("b_56")

a12T = sym.Symbol("a_12T")
a23T = sym.Symbol("a_23T")
a34T = sym.Symbol("a_34T")
a45T = sym.Symbol("a_45T")
a56T = sym.Symbol("a_56T")

b12T = sym.Symbol("b_12T")
b23T = sym.Symbol("b_23T")
b34T = sym.Symbol("b_34T")
b45T = sym.Symbol("b_45T")
b56T = sym.Symbol("b_56T")

kon  = sym.Symbol("kon")
koff = sym.Symbol("koff")
P = sym.Symbol("P")
PT = sym.Symbol("PT")


def Q_matrices(LG: np.ndarray, N: int, cQ:np.ndarray,
               cSigma:int = 1, k:int = 0):
    """ Evaluation of "matrices of in-forests"
            LG -> is minus the transposed laplacian
                  matrix of the graph G
            N  -> is the number of vertex in the graph G
            cQ -> is the NxN identity
    """
    if k == N-1:
        return []
    else: 
        sigma_p1 = trace(LG@cQ)/(k+1)
        Q_p1     =  -LG@cQ + sigma_p1*sym.eye(N)
        return [Q_p1] + Q_matrices(LG, N, cQ = Q_p1, cSigma = sigma_p1, k = k+1)
    
def MeanFPT(Qs, N): 
    """ Evaluate the Mean of the First passage time distribution given an array of N-1 "matrices of in-forests".
            Qs -> array of N-1 matrices in-forests
            N  -> number of states in the graph
    """
    
    Q_Nmin1 = Qs[-1]    # Q_matrices computes the matrix for k in [0,N-1], so the last one is N-1
    Q_Nmin2 = Qs[-2]    # While this one is N-2
    
    # Q_Nmin2[0,:-1] from the first to N-1 of the first column Q_Nmin1[N-1, N-1] the last entry in the matrix.
    num = 0
    for i in Q_Nmin2[0, :-1]:
        num = i + num
    den = Q_Nmin1[N-1, N-1]
    return num/den

def ss(Qs, N): 
    """ Evaluate the Mean of the First passage time distribution given an array of N-1 "matrices of in-forests".
            Qs -> array of N-1 matrices in-forests
            N  -> number of states in the graph
    """
    ss_vect = [Qs[-1][i,i]/trace(Qs[-1]) for i in range(N)]
    return ss_vect

def L3_matrices(pset): 

    Lap_3_ss = sym.Matrix([[0,    b12,     0,      koff,        0,       0], 
                            [a12,      0,   b23,          0,    koff,       0], 
                            [0,    a23,     0,          0,       0,    koff],
                            [kon,      0,     0,          0,    b12T,       0], 
                            [0,    kon,     0,       a12T,       0,    b23T],
                            [0,      0,   kon,          0,    a23T,       0]])
        
    Lap_3 = sym.Matrix([[0,    b12,     0,      koff,        0,       0,     0], 
                        [a12,      0,   b23,          0,    koff,       0,     0], 
                        [0,    a23,     0,          0,       0,    koff,     0],
                        [kon,      0,     0,          0,    b12T,       0,     0], 
                        [0,    kon,     0,       a12T,       0,    b23T,     0],
                        [0,      0,   kon,          0,    a23T,       0,     0],
                        [0,      0,     P,          0,       0,      PT,     0]])


    for i in range(Lap_3.shape[0]):
        Lap_3[i,i] = -np.sum(Lap_3[:,i])
    for i in range(Lap_3_ss.shape[0]):
        Lap_3_ss[i,i] = -np.sum(Lap_3_ss[:,i])

    return Lap_3.evalf(100, subs = pset), Lap_3_ss.evalf(100, subs = pset)

def L4_matrices(pset): 

    Lap_4 = sym.Matrix([[0,    b12,     0,       0,    koff,      0,     0,     0,     0],
                        [a12,      0,   b23,       0,       0,   koff,     0,     0,     0],
                        [0,    a23,     0,     b34,       0,      0,  koff,     0,     0],
                        [0,      0,   a34,       0,       0,      0,     0,  koff,     0],
                        [kon,      0,     0,       0,       0,   b12T,     0,     0,     0],
                        [0,    kon,     0,       0,    a12T,      0,  b23T,     0,     0],
                        [0,      0,   kon,       0,       0,   a23T,     0,  b34T,     0],
                        [0,      0,     0,     kon,       0,      0,  a34T,     0,     0],
                        [0,      0,     0,       P,       0,      0,     0,     P,     0]])

    Lap_4_ss = sym.Matrix([[0,    b12,     0,       0,    koff,      0,     0,     0],
                           [a12,    0,    b23,      0,       0,   koff,     0,     0],
                           [0,    a23,     0,     b34,       0,      0,  koff,     0],
                           [0,      0,   a34,       0,       0,      0,     0,  koff],
                           [kon,    0,     0,       0,       0,   b12T,     0,     0],
                           [0,    kon,     0,       0,    a12T,      0,  b23T,     0],
                           [0,      0,   kon,       0,       0,   a23T,     0,  b34T],
                           [0,      0,     0,     kon,       0,      0,  a34T,     0]])


    for i in range(Lap_4.shape[0]):
        Lap_4[i,i] = -np.sum(Lap_4[:,i])
    for i in range(Lap_4_ss.shape[0]):
        Lap_4_ss[i,i] = -np.sum(Lap_4_ss[:,i])

    return Lap_4.evalf(100, subs = pset), Lap_4_ss.evalf(100, subs = pset)

def L5_matrices(pset): 

    Lap_5 = sym.Matrix([[   0,    b12,     0,       0,      0,    koff,       0,     0,     0,     0,     0], 
                           [ a12,      0,   b23,       0,      0,       0,    koff,     0,     0,     0,     0],
                           [   0,    a23,     0,     b34,      0,       0,       0,  koff,     0,     0,     0],
                           [   0,      0,   a34,       0,    b45,       0,       0,     0,  koff,     0,     0],
                           [   0,      0,     0,     a45,      0,       0,       0,     0,     0,  koff,     0],
                           [ kon,      0,     0,       0,      0,       0,    b12T,     0,     0,     0,     0],
                           [   0,    kon,     0,       0,      0,    a12T,       0,  b23T,     0,     0,     0],
                           [   0,      0,   kon,       0,      0,       0,    a23T,     0,  b34T,     0,     0],
                           [   0,      0,     0,     kon,      0,       0,       0,  a34T,     0,  b45T,     0],
                           [   0,      0,     0,       0,    kon,       0,       0,     0,  a45T,     0,     0],
                           [   0,      0,     0,       0,      P,       0,       0,     0,     0,     P,     0]])
        
    Lap_5_ss = sym.Matrix([[   0,    b12,     0,       0,      0,    koff,       0,     0,     0,     0], 
                        [ a12,      0,   b23,       0,      0,       0,    koff,     0,     0,     0],
                        [   0,    a23,     0,     b34,      0,       0,       0,  koff,     0,     0],
                        [   0,      0,   a34,       0,    b45,       0,       0,     0,  koff,     0],
                        [   0,      0,     0,     a45,      0,       0,       0,     0,     0,  koff],
                        [ kon,      0,     0,       0,      0,       0,    b12T,     0,     0,     0],
                        [   0,    kon,     0,       0,      0,    a12T,       0,  b23T,     0,     0],
                        [   0,      0,   kon,       0,      0,       0,    a23T,     0,  b34T,     0],
                        [   0,      0,     0,     kon,      0,       0,       0,  a34T,     0,  b45T],
                        [   0,      0,     0,       0,    kon,       0,       0,     0,  a45T,     0]])

    for i in range(Lap_5.shape[0]):
        Lap_5[i,i] = -np.sum(Lap_5[:,i])
    for i in range(Lap_5_ss.shape[0]):
        Lap_5_ss[i,i] = -np.sum(Lap_5_ss[:,i])

    return Lap_5.evalf(100, subs = pset), Lap_5_ss.evalf(100, subs = pset)

def L6_matrices(pset): 

    Lap_6 = sym.Matrix([ [0,    b12,     0,       0,       0,       0,  koff,     0,      0,     0,     0,     0,    0],
                            [a12,    0,   b23,       0,       0,       0,     0,  koff,      0,     0,     0,     0,    0],
                            [0,    a23,     0,     b34,       0,       0,     0,     0,   koff,     0,     0,     0,    0],
                            [0,      0,   a34,       0,     b45,       0,     0,     0,      0,  koff,     0,     0,    0],
                            [0,      0,     0,     a45,       0,     b56,     0,     0,      0,     0,  koff,     0,    0],
                            [0,      0,     0,       0,     a56,       0,     0,     0,      0,     0,     0,  koff,    0],
                            [kon,    0,     0,       0,       0,       0,     0,  b12T,      0,     0,     0,     0,    0],
                            [0,    kon,     0,       0,       0,       0,  a12T,     0,   b23T,     0,     0,     0,    0],
                            [0,      0,   kon,       0,       0,       0,     0,  a23T,      0,  b34T,     0,     0,    0],
                            [0,      0,     0,     kon,       0,       0,     0,     0,   a34T,     0,  b45T,     0,    0],
                            [0,      0,     0,       0,     kon,       0,     0,     0,      0,  a45T,     0,  b56T,    0],
                            [0,      0,     0,       0,       0,     kon,     0,     0,      0,     0,  a56T,     0,    0],
                            [0,      0,     0,       0,       0,       P,     0,     0,      0,     0,     0,    PT,    0]])
        
    Lap_6_ss = sym.Matrix([[0,    b12,     0,       0,       0,       0,  koff,     0,      0,     0,     0,     0],
                        [a12,    0,   b23,       0,       0,       0,     0,  koff,      0,     0,     0,     0],
                        [0,    a23,     0,     b34,       0,       0,     0,     0,   koff,     0,     0,     0],
                        [0,      0,   a34,       0,     b45,       0,     0,     0,      0,  koff,     0,     0],
                        [0,      0,     0,     a45,       0,     b56,     0,     0,      0,     0,  koff,     0],
                        [0,      0,     0,       0,     a56,       0,     0,     0,      0,     0,     0,  koff],
                        [kon,    0,     0,       0,       0,       0,     0,  b12T,      0,     0,     0,     0],
                        [0,    kon,     0,       0,       0,       0,  a12T,     0,   b23T,     0,     0,     0],
                        [0,      0,   kon,       0,       0,       0,     0,  a23T,      0,  b34T,     0,     0],
                        [0,      0,     0,     kon,       0,       0,     0,     0,   a34T,     0,  b45T,     0],
                        [0,      0,     0,       0,     kon,       0,     0,     0,      0,  a45T,     0,  b56T],
                        [0,      0,     0,       0,       0,     kon,     0,     0,      0,     0,  a56T,     0]])


    for i in range(Lap_6.shape[0]):
        Lap_6[i,i] = -np.sum(Lap_6[:,i])
    for i in range(Lap_6_ss.shape[0]):
        Lap_6_ss[i,i] = -np.sum(Lap_6_ss[:,i])

    return Lap_6.evalf(100, subs = pset), Lap_6_ss.evalf(100, subs = pset)


if __name__=="__main__": 

        logging.info('Running test for Ladder 3')

        obj = Ladder_3_v3()
        sample = 10**np.random.uniform(-1,1,12)
        list_symbols = [a12, a23, b12, b23, a12T, a23T, b12T, b23T, kon, koff, P, PT]
        dict_eval = {list_symbols[i]:p for i,p in enumerate(sample)}

        # recurrence calculations
        lap_fpt, lap_ss = L3_matrices(dict_eval)
        Q_mat_ss  = Q_matrices(-lap_ss.T, N = lap_ss.shape[0], cQ = sym.eye(lap_ss.shape[0]))
        Q_mat_fpt = Q_matrices(-lap_fpt.T, N = lap_fpt.shape[0], cQ = sym.eye(lap_fpt.shape[0]))
        ss_vector = ss(Q_mat_ss, N = lap_ss.shape[0])
        ss_recurrence = ss_vector[2]+ss_vector[5]
        mean_recurrence = MeanFPT(Q_mat_fpt, N = lap_fpt.shape[0])

        # SVD calculations
        obj.setLaplacian(sample)
        ss_svd_vect = obj.getSteadyState(False) 
        ss_svd = ss_svd_vect[2]+ss_svd_vect[5]
        mean_fpt = obj.get_FPT_stat(False, False)[0][0]

        logging.info(f'Steady State recurrence = {ss_recurrence}')
        logging.info(f'Steady State svd = {ss_svd}')
        logging.info(f'Mean FPT recurrence = {mean_recurrence}')
        logging.info(f'Mean FPT svd = {round(mean_fpt,5)}')

        tollerance_ss = np.isclose(float(ss_recurrence), float(ss_svd))
        tollearnce_mean_fpt = np.isclose(float(mean_recurrence), float(mean_fpt))

        if tollerance_ss and tollearnce_mean_fpt:
            logging.info(f'Ladder 3 test PASSED \n')
        else: 
            logging.info(f'Ladder 3 test FAILED \n') 

        

        logging.info('Running test for Ladder 4')
        obj = Ladder_4()
        sample = 10**np.random.uniform(-1,1,16)
        list_symbols = [a12, a23, a34, b12, b23, b34, a12T, a23T, a34T, b12T, b23T, b34T, kon, koff, P, PT]
        dict_eval = {list_symbols[i]:p for i,p in enumerate(sample)}

        # recurrence calculations
        lap_fpt, lap_ss = L4_matrices(dict_eval)
        Q_mat_ss  = Q_matrices(-lap_ss.T, N = lap_ss.shape[0], cQ = sym.eye(lap_ss.shape[0]))
        Q_mat_fpt = Q_matrices(-lap_fpt.T, N = lap_fpt.shape[0], cQ = sym.eye(lap_fpt.shape[0]))
        ss_vector = ss(Q_mat_ss, N = lap_ss.shape[0])
        ss_recurrence = ss_vector[3]+ss_vector[7]
        mean_recurrence = MeanFPT(Q_mat_fpt, N = lap_fpt.shape[0])

        # SVD calculations
        obj.setLaplacian(sample)
        ss_svd_vect = obj.getSteadyState(False) 
        ss_svd = ss_svd_vect[3]+ss_svd_vect[7]
        mean_fpt = obj.get_FPT_stat(False, False)[0][0]

        logging.info(f'Steady State recurrence = {ss_recurrence}')
        logging.info(f'Steady State svd = {ss_svd}')
        logging.info(f'Mean FPT recurrence = {mean_recurrence}')
        logging.info(f'Mean FPT svd = {mean_fpt}')

        tollerance_ss = np.isclose(float(ss_recurrence), float(ss_svd))
        tollearnce_mean_fpt = np.isclose(float(mean_recurrence), float(mean_fpt))

        if tollerance_ss and tollearnce_mean_fpt:
            logging.info(f'Ladder 4 test PASSED \n')
        else: 
            logging.info(f'Ladder 4 test FAILED \n') 

        logging.info('Running test for Ladder 5')
        obj = Ladder_5()
        sample = 10**np.random.uniform(-1,1,20)
        list_symbols = [a12, a23, a34, a45, b12, b23, b34, b45, a12T, a23T, a34T, a45T, b12T, b23T, b34T, b45T, kon, koff, P, PT]
        dict_eval = {list_symbols[i]:p for i,p in enumerate(sample)}

        # recurrence calculations
        lap_fpt, lap_ss = L5_matrices(dict_eval)
        Q_mat_ss  = Q_matrices(-lap_ss.T, N = lap_ss.shape[0], cQ = sym.eye(lap_ss.shape[0]))
        Q_mat_fpt = Q_matrices(-lap_fpt.T, N = lap_fpt.shape[0], cQ = sym.eye(lap_fpt.shape[0]))
        ss_vector = ss(Q_mat_ss, N = lap_ss.shape[0])
        ss_recurrence = ss_vector[4]+ss_vector[9]
        mean_recurrence = MeanFPT(Q_mat_fpt, N = lap_fpt.shape[0])

        # SVD calculations
        obj.setLaplacian(sample)
        ss_svd_vect = obj.getSteadyState(False) 
        ss_svd = ss_svd_vect[4]+ss_svd_vect[9]
        mean_fpt = obj.get_FPT_stat(False, False)[0][0]

        logging.info(f'Steady State recurrence = {ss_recurrence}')
        logging.info(f'Steady State svd = {ss_svd}')
        logging.info(f'Mean FPT recurrence = {mean_recurrence}')
        logging.info(f'Mean FPT svd = {mean_fpt}')

        tollerance_ss = np.isclose(float(ss_recurrence), float(ss_svd))
        tollearnce_mean_fpt = np.isclose(float(mean_recurrence), float(mean_fpt))

        if tollerance_ss and tollearnce_mean_fpt:
            logging.info(f'Ladder 5 test PASSED \n')
        else: 
            logging.info(f'Ladder 5 test FAILED \n') 


        logging.info('Running test for Ladder 6')
        obj = Ladder_6_prec_100_v2()
        sample = 10**np.random.uniform(-1,1,24)
        list_symbols = [a12, a23, a34, a45, a56, b12, b23, b34, b45, b56, a12T, a23T, a34T, a45T, a56T, b12T, b23T, b34T, b45T, b56T, kon, koff, P, PT]
        dict_eval = {list_symbols[i]:p for i,p in enumerate(sample)}

        # recurrence calculations
        lap_fpt, lap_ss = L6_matrices(dict_eval)
        Q_mat_ss  = Q_matrices(-lap_ss.T, N = lap_ss.shape[0], cQ = sym.eye(lap_ss.shape[0]))
        Q_mat_fpt = Q_matrices(-lap_fpt.T, N = lap_fpt.shape[0], cQ = sym.eye(lap_fpt.shape[0]))
        ss_vector = ss(Q_mat_ss, N = lap_ss.shape[0])
        ss_recurrence = ss_vector[5]+ss_vector[11]
        mean_recurrence = MeanFPT(Q_mat_fpt, N = lap_fpt.shape[0])

        # SVD calculations
        obj.setLaplacian(sample)
        ss_svd_vect = obj.getSteadyState(False) 
        ss_svd = ss_svd_vect[5]+ss_svd_vect[11]
        mean_fpt = obj.get_FPT_stat(False, False)[0][0]

        logging.info(f'Steady State recurrence = {ss_recurrence}')
        logging.info(f'Steady State svd = {ss_svd}')
        logging.info(f'Mean FPT recurrence = {mean_recurrence}')
        logging.info(f'Mean FPT svd = {mean_fpt}')

        tollerance_ss = np.isclose(float(ss_recurrence), float(ss_svd))
        tollearnce_mean_fpt = np.isclose(float(mean_recurrence), float(mean_fpt))

        if tollerance_ss and tollearnce_mean_fpt:
            logging.info(f'Ladder 5 test PASSED \n')
        else: 
            logging.info(f'Ladder 5 test FAILED \n') 

