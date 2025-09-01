#!/bin/bash


# forward activation g12 > 1 activation case
math -script test_exp_gt_0.wl mfpt_g12.in "g12>1 && g12>0"  
math -script test_exp_lt_0.wl mfpt_g12.in "g12>1 && g12>0"
math -script test_exp_gt_0.wl ss_g12.in "g12>1 && g12>0"  
math -script test_exp_lt_0.wl ss_g12.in "g12>1 && g12>0"

# forward repression g12 > 1 activation case
math -script test_exp_gt_0.wl mfpt_g12.in "g12<1 && g12>0"  
math -script test_exp_lt_0.wl mfpt_g12.in "g12<1 && g12>0"
math -script test_exp_gt_0.wl ss_g12.in "g12<1 && g12>0"  
math -script test_exp_lt_0.wl ss_g12.in "g12<1 && g12>0"


# backward activation g12 > 1 activation case
math -script test_exp_gt_0.wl mfpt_g21.in "g21>1 && g21>0"  
math -script test_exp_lt_0.wl mfpt_g21.in "g21>1 && g21>0"
math -script test_exp_gt_0.wl ss_g21.in "g21>1 && g21>0"  
math -script test_exp_lt_0.wl ss_g21.in "g21>1 && g21>0"

# backward repression g12 > 1 activation case
math -script test_exp_gt_0.wl mfpt_g21.in "g21<1 && g21>0"  
math -script test_exp_lt_0.wl mfpt_g21.in "g21<1 && g21>0"
math -script test_exp_gt_0.wl ss_g21.in "g21<1 && g21>0"  
math -script test_exp_lt_0.wl ss_g21.in "g21<1 && g21>0"
