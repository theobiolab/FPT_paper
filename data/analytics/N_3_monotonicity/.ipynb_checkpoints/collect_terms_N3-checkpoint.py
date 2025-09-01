import sympy as sym
import numpy as np 
from sympy import *
from sympy.printing import mathematica
from collections import OrderedDict

import sys
sys.path.insert(0, "..")
import recurrence_func

import logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)


l12  = sym.symbols('l_{12}')
l23  = sym.symbols('l_{23}')
l21  = sym.symbols('l_{21}')
l32  = sym.symbols('l_{32}')
l12T = sym.symbols('l_{12}^T')
l23T = sym.symbols('l_{23}^T')
l21T = sym.symbols('l_{21}^T')
l32T = sym.symbols('l_{32}^T')
kon  = sym.symbols('k_{on}')
koff  = sym.symbols('k_{off}')
r    = sym.symbols('r')
x    = sym.symbols("x")

        
lap = sym.Matrix([[0,    l21,     0,      koff,        0,       0,     0], 
                  [l12,    0,   l32,          0,    koff,       0,     0], 
                  [0,    l23,     0,          0,       0,    koff,     0],
                  [kon*x,    0,     0,          0,    l21T,       0,     0], 
                  [0,    kon*x,     0,       l12T,       0,    l32T,     0],
                  [0,      0,   kon*x,          0,    l23T,       0,     0],
                  [0,      0,     r,          0,       0,       r,     0]])

lap_ss = sym.Matrix([[0,    l21,     0,      koff,        0,      0], 
                    [l12,    0,   l32,          0,    koff,       0], 
                    [0,    l23,     0,          0,       0,    koff],
                    [kon*x,    0,     0,          0,    l21T,       0], 
                    [0,    kon*x,     0,       l12T,       0,    l32T],
                    [0,      0,   kon*x,          0,    l23T,       0]])


for i in range(lap.shape[0]): 
    lap[i,i] = -np.sum(lap[:,i])
LG_fpt = -lap.T

for i in range(lap_ss.shape[0]): 
    lap_ss[i,i] = -np.sum(lap_ss[:,i])
LG_ss = -lap_ss.T


####################### compute Q matrices general case ###########################

logging.info('Computing Q matrices for the FPT laplacian')
qc = recurrence_func.Q_matrices(LG = LG_fpt, N = LG_fpt.shape[0], cQ = sym.eye(LG_fpt.shape[0]))
logging.info('Computing Q matrices for the SS laplacian')
q_ss = recurrence_func.Q_matrices(LG = LG_ss, N = LG_ss.shape[0], cQ = sym.eye(LG_ss.shape[0]))

####################### Computing the steady state and MFPT and the derivatives numerator #################

logging.info('Computing steady state and mfpt derivatives')

p_ss    = sym.simplify((q_ss[-1][2,2]+q_ss[-1][5,5])/(np.trace(q_ss[-1])))
mean_fpt = sym.simplify(recurrence_func.MeanFPT(qc, LG_fpt.shape[0]))

num_p_ss, den_p_ss = sym.fraction(p_ss)
num_mean_fpt, den_mean_fpt = sym.fraction(mean_fpt)

num_d_ss_dx   = sym.simplify(sym.diff(num_p_ss, x)*den_p_ss - num_p_ss*sym.diff(den_p_ss, x))
num_d_mfpt_dx = sym.simplify(sym.expand(sym.diff(num_mean_fpt, x)*den_mean_fpt - num_mean_fpt*sym.diff(den_mean_fpt, x)))

# ##################### Collect terms for all forward ###############################

# logging.info('Factoring out the regulatory terms for g12, g23')

# g12    = sym.symbols('g12')
# g23    = sym.symbols('g23')

# forward_dss   = sym.simplify(num_d_ss_dx.subs([(l12T, l12*g12), (l23T, l23*g23), (l21T, l21), (l32T, l32)]).expand())    # subs regulatory parameters
# forward_dmfpt = sym.simplify(num_d_mfpt_dx.subs([(l12T, l12*g12), (l23T, l23*g23), (l21T, l21), (l32T, l32)]).expand())  # subs regulatory parameters

# collected_forward_ss = recurrence_func.collect_terms(forward_dss, alpha = g12, beta = g23, x = x)
# collected_forward_mfpt = recurrence_func.collect_terms(forward_dmfpt, alpha = g12, beta = g23, x = x)

# recurrence_func.write_in(all_terms = collected_forward_mfpt, filename = "mfpt_g12_g23")
# recurrence_func.write_in(all_terms = collected_forward_ss, filename = "ss_g12_g23")

# ##################### Collect terms for all backward ###############################

# logging.info('Factoring out the regulatory terms for g21, g32')

# g21    = sym.symbols('g21')
# g32    = sym.symbols('g32')

# bakward_dss   = sym.simplify(num_d_ss_dx.subs([(l12T, l12), (l23T, l23), (l21T, l21*g21), (l32T, l32*g32)]).expand())    # subs regulatory parameters
# bakward_dmfpt = sym.simplify(num_d_mfpt_dx.subs([(l12T, l12), (l23T, l23), (l21T, l21*g21), (l32T, l32*g32)]).expand())  # subs regulatory parameters

# collected_bakward_ss = recurrence_func.collect_terms(bakward_dss, alpha = g21, beta = g32, x = x)
# collected_bakward_mfpt = recurrence_func.collect_terms(bakward_dmfpt, alpha = g21, beta = g32, x = x)

# recurrence_func.write_in(all_terms = collected_bakward_mfpt, filename = "mfpt_g21_g32")
# recurrence_func.write_in(all_terms = collected_bakward_ss, filename = "ss_g21_g32")

# ##################### Collect terms for b1 and f2 ###############################

# logging.info('Factoring out the regulatory terms for g21, g23')

# g21    = sym.symbols('g21')
# g23    = sym.symbols('g23')

# b1f2_dss   = sym.simplify(num_d_ss_dx.subs([(l12T, l12), (l23T, l23*g23), (l21T, l21*g21), (l32T, l32)]).expand())    # subs regulatory parameters
# b1f2_dmfpt = sym.simplify(num_d_mfpt_dx.subs([(l12T, l12), (l23T, l23*g23), (l21T, l21*g21), (l32T, l32)]).expand())  # subs regulatory parameters

# collected_b1f2_ss = recurrence_func.collect_terms(b1f2_dss, alpha = g21, beta = g23, x = x)
# collected_b1f2_mfpt = recurrence_func.collect_terms(b1f2_dmfpt, alpha = g21, beta = g23, x = x)

# recurrence_func.write_in(all_terms = collected_b1f2_mfpt, filename = "mfpt_g21_g23")
# recurrence_func.write_in(all_terms = collected_b1f2_ss, filename = "ss_g21_g23")


# ##################### Collect terms for f1 and b2 ###############################

# logging.info('Factoring out the regulatory terms for g12, g32')

# g12    = sym.symbols('g12')
# g32    = sym.symbols('g32')

# f1b2_dss   = sym.simplify(num_d_ss_dx.subs([(l12T, l12*g12), (l23T, l23), (l21T, l21), (l32T, l32*g32)]).expand())    # subs regulatory parameters
# f1b2_dmfpt = sym.simplify(num_d_mfpt_dx.subs([(l12T, l12*g12), (l23T, l23), (l21T, l21), (l32T, l32*g32)]).expand())  # subs regulatory parameters

# collected_f1b2_ss = recurrence_func.collect_terms(f1b2_dss, alpha = g12, beta = g32, x = x)
# collected_f1b2_mfpt = recurrence_func.collect_terms(f1b2_dmfpt, alpha = g12, beta = g32, x = x)

# recurrence_func.write_in(all_terms = collected_f1b2_mfpt, filename = "mfpt_g12_g32")
# recurrence_func.write_in(all_terms = collected_f1b2_ss, filename = "ss_g12_g32")


##################### Collect terms for f1 ###############################

logging.info('Factoring out the regulatory terms for g12')

g12    = sym.symbols('g12')

f1b2_dss   = sym.simplify(num_d_ss_dx.subs([(l12T, l12*g12), (l23T, l23), (l21T, l21), (l32T, l32)]).expand())    # subs regulatory parameters
f1b2_dmfpt = sym.simplify(num_d_mfpt_dx.subs([(l12T, l12*g12), (l23T, l23), (l21T, l21), (l32T, l32)]).expand())  # subs regulatory parameters

collected_f1b2_ss = recurrence_func.collect_terms(f1b2_dss, alpha = g12, beta = None, x = x)
collected_f1b2_mfpt = recurrence_func.collect_terms(f1b2_dmfpt, alpha = g12, beta = None, x = x)

recurrence_func.write_in(all_terms = collected_f1b2_mfpt, filename = "mfpt_g12")
recurrence_func.write_in(all_terms = collected_f1b2_ss, filename = "ss_g12")


##################### Collect terms for f2 ###############################

logging.info('Factoring out the regulatory terms for g23')

g23    = sym.symbols('g23')

f1b2_dss   = sym.simplify(num_d_ss_dx.subs([(l12T, l12), (l23T, l23*g23), (l21T, l21), (l32T, l32)]).expand())    # subs regulatory parameters
f1b2_dmfpt = sym.simplify(num_d_mfpt_dx.subs([(l12T, l12), (l23T, l23*g23), (l21T, l21), (l32T, l32)]).expand())  # subs regulatory parameters

collected_f1b2_ss = recurrence_func.collect_terms(f1b2_dss, alpha = g23, beta = None, x = x)
collected_f1b2_mfpt = recurrence_func.collect_terms(f1b2_dmfpt, alpha = g23, beta = None, x = x)

recurrence_func.write_in(all_terms = collected_f1b2_mfpt, filename = "mfpt_g23")
recurrence_func.write_in(all_terms = collected_f1b2_ss, filename = "ss_g23")


##################### Collect terms for b1 ###############################

logging.info('Factoring out the regulatory terms for g21')

g21    = sym.symbols('g21')

f1b2_dss   = sym.simplify(num_d_ss_dx.subs([(l12T, l12), (l23T, l23), (l21T, l21*g21), (l32T, l32)]).expand())    # subs regulatory parameters
f1b2_dmfpt = sym.simplify(num_d_mfpt_dx.subs([(l12T, l12), (l23T, l23), (l21T, l21*g21), (l32T, l32)]).expand())  # subs regulatory parameters

collected_f1b2_ss = recurrence_func.collect_terms(f1b2_dss, alpha = g21, beta = None, x = x)
collected_f1b2_mfpt = recurrence_func.collect_terms(f1b2_dmfpt, alpha = g21, beta = None, x = x)

recurrence_func.write_in(all_terms = collected_f1b2_mfpt, filename = "mfpt_g21")
recurrence_func.write_in(all_terms = collected_f1b2_ss, filename = "ss_g21")

##################### Collect terms for b2 ###############################

logging.info('Factoring out the regulatory terms for g32')

g32    = sym.symbols('g32')

f1b2_dss   = sym.simplify(num_d_ss_dx.subs([(l12T, l12), (l23T, l23), (l21T, l21), (l32T, l32*g32)]).expand())    # subs regulatory parameters
f1b2_dmfpt = sym.simplify(num_d_mfpt_dx.subs([(l12T, l12), (l23T, l23), (l21T, l21), (l32T, l32*g32)]).expand())  # subs regulatory parameters

collected_f1b2_ss = recurrence_func.collect_terms(f1b2_dss, alpha = g32, beta = None, x = x)
collected_f1b2_mfpt = recurrence_func.collect_terms(f1b2_dmfpt, alpha = g32, beta = None, x = x)

recurrence_func.write_in(all_terms = collected_f1b2_mfpt, filename = "mfpt_g32")
recurrence_func.write_in(all_terms = collected_f1b2_ss, filename = "ss_g32")






