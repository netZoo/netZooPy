import pytest
import numpy as np
import os
from netZooPy import dragon

def test_dragon():
    #1. test1
    print('Start Dragon run ...')
    n = 1000
    p1 = 500
    p2 = 100
    X1, X2, Theta, _ = dragon.simulate_dragon_data(eta11=0.005, eta12=0.005, eta22=0.05,
                                            p1=100, p2=500, epsilon=[0.1,0.1],
                                            n=n, seed=123)
    os.system('curl -O https://netzoo.s3.us-east-2.amazonaws.com/netZooPy/unittest_datasets/dragonx1.npy')
    os.system('curl -O https://netzoo.s3.us-east-2.amazonaws.com/netZooPy/unittest_datasets/dragonx2.npy')
    X1=np.load('dragonx1.npy')
    X2=np.load('dragonx2.npy')
    lambdas, lambdas_landscape = dragon.estimate_penalty_parameters_dragon(X1, X2)
    lambdasSingle=tuple([int(10*x)/10 for x in lambdas]) # 3 digit precision
    alamb=lambdas_landscape[1,1]
    assert(lambdasSingle == (0.9, 0.9))
    assert((alamb < 398.7*1.002) & (alamb > 398.7*0.998)) #0.2% of error
    assert(int(X1[1,1]*1000)/1000 ==0.880)
    assert(int(X2[1,1]*1000)/1000 ==0.664)
    
    #2. test2
    r = dragon.get_partial_correlation_dragon(X1, X2, lambdas)
    adj_p_vals, p_vals = dragon.estimate_p_values_dragon(r, n, p1, p2, lambdas)
    p_valstest=int(p_vals[2,1]*100)/100
    adj_p_valstest=int(adj_p_vals[2,1]*10)/10 # 3 digit precision
    assert(int(np.max(r)*100)/100 == 0.03)
    assert(int(r[1,2]*100000)/100000 == 0.00012)
    assert(p_valstest == 0.96)
    assert(adj_p_valstest == 0.9)

    #3. test monte carlo p-values
    p1 = 3
    p2 = 4
    n = 10
    lam = [0,0] # no shrinkage 
    test11 = np.array([[1,1/2.,-1/4.],
                   [1/2.,1,1/8.],
                   [1/4.,1/8.,1]])
    test12 = np.array([[-3/4.,1/2.,1/4.,0],
                  [-3/4.,1/2.,1/4.,0],
                  [-3/4.,1/2.,1/4.,0]])
    
    test21 = np.transpose(test12)

    test22 = np.array([[1,-3/4.,1/2.,-1/4.],
                  [-3/4.,1,1/8.,1/16.],
                  [1/2.,1/8.,1,1/32.],
                  [-1/4.,1/16.,1/32.,1]])
    test_mc_mat = np.identity(p1+p2)
    test_mc_mat[0:3,0:3] = test11
    test_mc_mat[0:3,3:7] = test12
    test_mc_mat[3:7,0:3] = test21
    test_mc_mat[3:7,3:7] = test22

    dragon_p_mc = dragon.dragon.estimate_p_values_mc(test_mc_mat,n,p1,p2,lam,seed=412)
    print(dragon_p_mc)
    ground_truth_mc_p = np.array([[0,0,1/3.,0,1/12.,7/12.,1],
                                 [0,0,2/3.,0,1/12.,7/12.,1],
                                 [1/3.,2/3.,0,0,1/12,7/12.,1],
                                 [0,0,0,0,0,0,5/6.],
                                 [1/12.,1/12.,1/12.,0,0,5/6.,5/6.],
                                 [7/12.,7/12.,7/12.,0,5/6.,0,5/6.],
                                 [1,1,1,5/6.,5/6.,5/6.,0]])
    print(ground_truth_mc_p)
    print(dragon_p_mc - ground_truth_mc_p)
    assert(np.array_equal(dragon_p_mc,ground_truth_mc_p))

test_dragon()