import pytest
import os
from netZooPy.lioness.lioness import Lioness
from netZooPy.panda.panda import Panda
import pandas as pd
import numpy as np
import glob

def test_lioness():
    print('Start lioness test!')
    #1. First generate temporary PANDA files as inputs for Lioness
    ppi            ='tests/puma/ToyData/ToyPPIData.txt'
    motif          ='tests/puma/ToyData/ToyMotifData.txt'
    expression_data='tests/puma/ToyData/ToyExpressionData.txt'
    lioness_file   ='Travis'
    rm_missing     =False
    output_file    ='panda.npy'
    gt_file        ='tests/panda/test_panda.txt'
    panda_obj      =Panda(expression_data, motif, ppi, save_tmp=True, remove_missing=rm_missing,
                      keep_expression_matrix=bool(lioness_file), modeProcess='legacy', save_memory=False)
    # Set parameters
    lioness_obj = Lioness(panda_obj, start=1, end=1)
    lioness_obj.save_lioness_results(lioness_file)
    # Read first lioness network
    res  = np.load('lioness_output/lioness.1.npy')
    gt = np.load('tests/lioness/lioness.1.npy')
    # Compare to ground truth
    assert(np.allclose(gt,res))

    #2. Testing Lioness with motif set to None to compute Lioness on coexpression networks
    motif          = None
    # Make sure to keep epxression matrix for next step
    panda_obj      = Panda(expression_data, motif, ppi, save_tmp=True, remove_missing=rm_missing,
                      keep_expression_matrix=True, modeProcess='legacy')
    lioness_obj    = Lioness(panda_obj, start=1, end=1)
    lioness_obj.save_lioness_results(lioness_file)
    # Read first lioness network
    res  = np.load('lioness_output/lioness.1.npy')
    gt   = np.load('tests/lioness/lionessCoexpression.1.npy')
    # Compare to ground truth
    assert(np.allclose(gt,res))

    #2. Testing Lioness in parallel
    c=np.randn(0,46)
    lioness_obj = Lioness(panda_obj, ncores=2,start=c,end=c+3)
    lioness_obj.save_lioness_results(lioness_file)
    
    # traces=glob.glob('lioness_output/*.npy')
    # res = pd.DataFrame()
    # for i,trace in enumerate(traces):
    #     data=np.load(trace)
    #     res=pd.concat([res,pd.DataFrame(data.flatten())],axis=1)
    # res=res.dropna().to_numpy()

    # np.save('lioness_output/lioness.all.npy',res)
    # Read first lioness network
    res  = np.load('lioness_output/lioness.'+str(c)+'.npy')
    gt = np.load('tests/lioness/lioness.'+str(c)+'.npy')
    # Compare to ground truth
    assert(np.allclose(gt,res))

