#!/bin/bash

mkdir -p results

# g12 g23 > 1
math -script test_exp_gt_0.wl mfpt_g12_g23.in "g12>1 && g12>0 && g23>1 && g23>0"  
math -script test_exp_lt_0.wl mfpt_g12_g23.in "g12>1 && g12>0 && g23>1 && g23>0" 
math -script test_exp_gt_0.wl ss_g12_g23.in "g12>1 && g12>0 && g23>1 && g23>0" 
math -script test_exp_lt_0.wl ss_g12_g23.in "g12>1 && g12>0 && g23>1 && g23>0" 

# g12 g23 < 1
math -script test_exp_gt_0.wl mfpt_g12_g23.in "g12<1 && g12>0 && g23<1 && g23>0"  
math -script test_exp_lt_0.wl mfpt_g12_g23.in "g12<1 && g12>0 && g23<1 && g23>0" 
math -script test_exp_gt_0.wl ss_g12_g23.in "g12<1 && g12>0 && g23<1 && g23>0" 
math -script test_exp_lt_0.wl ss_g12_g23.in "g12<1 && g12>0 && g23<1 && g23>0" 

# g21 g32 > 1
math -script test_exp_gt_0.wl mfpt_g21_g32.in "g21>1 && g21>0 && g32>1 && g32>0"  
math -script test_exp_lt_0.wl mfpt_g21_g32.in "g21>1 && g21>0 && g32>1 && g32>0"
math -script test_exp_gt_0.wl ss_g21_g32.in "g21>1 && g21>0 && g32>1 && g32>0"  
math -script test_exp_lt_0.wl ss_g21_g32.in "g21>1 && g21>0 && g32>1 && g32>0"

# g21 g32 < 1
math -script test_exp_gt_0.wl mfpt_g21_g32.in "g21<1 && g21>0 && g32<1 && g32>0"
math -script test_exp_lt_0.wl mfpt_g21_g32.in "g21<1 && g21>0 && g32<1 && g32>0"
math -script test_exp_gt_0.wl ss_g21_g32.in "g21<1 && g21>0 && g32<1 && g32>0"
math -script test_exp_lt_0.wl ss_g21_g32.in "g21<1 && g21>0 && g32<1 && g32>0"
