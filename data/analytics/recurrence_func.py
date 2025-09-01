import sympy as sym
import numpy as np
import matplotlib.pyplot as plt
from sympy import degree
from sympy import *
from sympy.printing import mathematica
import itertools 
from collections import OrderedDict

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
        sigma_p1 = sym.simplify(np.trace(LG@cQ)/(k+1))
        Q_p1     = sym.simplify(-LG@cQ + sigma_p1*sym.eye(N))
        return [Q_p1] + Q_matrices(LG, N, cQ = Q_p1, cSigma = sigma_p1, k = k+1) 
    
def MeanFPT(Qs, N): 
    """ Evaluate the Mean of the First passage time distribution given an array of N-1 "matrices of in-forests".
            Qs -> array of N-1 matrices in-forests
            N  -> number of states in the graph
    """
    
    Q_Nmin1 = Qs[-1]    # Q_matrices computes the matrix for k in [0,N-1], so the last one is N-1
    Q_Nmin2 = Qs[-2]    # While this one is N-2
    
    # Q_Nmin2[0,:-1] from the first to N-1 of the first column Q_Nmin1[N-1, N-1] the last entry in the matrix.
    ratio = Q_Nmin2[0, :-1]/Q_Nmin1[N-1, N-1]
    return np.sum(ratio)


def CV_FPT(Qs, N, mean): 
    """ Evaluate the Std of the First passage time distribution given an array of N-1 "matrices of in-forests".
            Qs -> array of N-1 matrices in-forests
            N  -> number of states in the graph
    """
    
    Q_Nmin1 = Qs[-1]
    Q_Nmin2 = Qs[-2]
    
    mean_squared = mean**2
    
    first_term = []
    for j in range(N-1):
        for k in range(N-1):
            num = Q_Nmin2[0, j]*Q_Nmin2[j, k]
            den = (Q_Nmin1[N-1, N-1])**2
            
            first_term.append(num/den)
    
    std = sym.sqrt(2*np.sum(first_term) - mean_squared)
    
    return std/mean
    


def collect_terms(der_num_x, alpha, beta, x):

    der_num_x = sym.Poly(der_num_x, x) 

    all_terms=[]
    for cn,coeff in enumerate(der_num_x.coeffs()[::-1]):
        termsc=[] 
        for monomial in sym.Add.make_args(coeff): 
            noneps=[]
            eps=[]
            number=1
            eps_base=[] 
            noneps_base=[] 
            args= sym.Mul.make_args(monomial)
            for arg_ in args:
                if type(arg_)==sym.core.power.Pow:
                    arg=sym.core.power.Pow.as_base_exp(arg_)[0] 
                else:
                    arg=arg_
                if  arg in [alpha, beta]:
                    eps.append(arg_)
                    eps_base.append(arg.name)
                else:
                    if type(arg)!=sym.core.numbers.Integer and type(arg)!=sym.core.numbers.NegativeOne:
                        noneps_base.append(arg.name)
                        noneps.append(arg_)
                    else:
                        number=arg

            noneps_o=np.array(noneps)[np.argsort(np.array(noneps_base))]
            eps_o=np.array(eps)[np.argsort(np.array(eps_base))]
            termsc.append([number,noneps_o,eps_o])
        all_terms.append(termsc)

    return all_terms

def write_in(all_terms, filename): 

    math_output=open(f"{filename}.in","w")

    # identify the unique terms
    for i, term_list in enumerate(all_terms):

        factor_list = [k[1] for k in term_list]
        baits = set(map(tuple, factor_list))

        
        while len(baits)>0: 
            bait_fact = np.array(baits.pop()) 

            eq = 0
            for j in term_list: 
                if set(j[1]) == set(bait_fact): 
                    if len(j[-1])>0:
                        eq += j[0]*np.prod(j[-1])
                    else:
                        eq += j[0]
            
            math_output.write(str(i)+";"+mathematica.mathematica_code(eq).replace("_","")+"\n")

    math_output.close()
