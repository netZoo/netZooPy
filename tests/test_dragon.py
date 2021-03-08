import pytest
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
    lambdas, lambdas_landscape = dragon.estimate_penalty_parameters_dragon(X1, X2)
    lambdas=tuple([int(1000*x)/1000 for x in lambdas]) # 3 digit precision
    alamb=int(lambdas_landscape[1,1]*10)/10
    assert(lambdas == (0.907, 0.913))
    assert(alamb == 398.7)

    #2. test2
    r = dragon.get_partial_correlation_dragon(X1, X2, lambdas)
    adj_p_vals, p_vals = dragon.estimate_p_values_dragon(r, n, p1, p2, lambdas)
    adj_p_vals=int(adj_p_vals[2,1]*1000)/1000 # 3 digit precision
    p_vals=int(p_vals[2,1]*100)/100
    assert(p_vals == 0.96)
    assert(adj_p_vals == 0.999)

     
