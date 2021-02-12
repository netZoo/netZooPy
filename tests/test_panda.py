import pytest
import os
from netZooPy.panda.panda import Panda
import pandas as pd
import numpy as np

def test_panda():
    #print(os.getcwd())
    print('Start Panda run ...')
    ppi            ='tests/puma/ToyData/ToyPPIData.txt'
    motif          ='tests/puma/ToyData/ToyMotifData.txt'
    expression_data='tests/puma/ToyData/ToyExpressionData.txt'
    lioness_file   =''
    rm_missing     = False
    output_file    ='travis_test_panda.txt'
    gt_file        ='tests/panda/union_test_panda.txt'
    gt_file_inter  ='tests/panda/inter_test_panda.txt'

    #0. Intersection
    panda_obj      = Panda(expression_data, motif, ppi, save_tmp=False, remove_missing=rm_missing,
                      keep_expression_matrix=bool(lioness_file), modeProcess='intersection')
    panda_obj.save_panda_results(output_file)
    res=pd.read_csv(output_file, sep=' ', header=None)
    gt =pd.read_csv(gt_file_inter, sep=' ', header=None)
    pd.testing.assert_frame_equal(res,gt,check_less_precise=False,check_exact=False)

    #0.1 Intersection via DataFrame
    expression = pd.read_csv(expression_data, sep='\t', header=None, index_col=0)
    motif_data = pd.read_csv(motif, sep='\t', names=['source','target','w'])
    ppi_data = pd.read_csv(ppi, sep='\t', header=None)
	    
    panda_obj      = Panda(expression, motif_data, ppi_data, save_tmp=False, remove_missing=rm_missing,
                      keep_expression_matrix=bool(lioness_file), modeProcess='intersection')
    panda_obj.save_panda_results(output_file)
    res=pd.read_csv(output_file, sep=' ', header=None)
    gt =pd.read_csv(gt_file_inter, sep=' ', header=None)
    pd.testing.assert_frame_equal(res,gt,check_less_precise=False,check_exact=False)

    #0.2 Intersection with symmetric PPI
    ppi_data = pd.read_csv(ppi, sep='\t', header=None)
    new_df   = pd.DataFrame(data={0:ppi_data[0],1:ppi_data[1],2:ppi_data[2]})
    ppi_data_symm  = pd.concat([ppi_data,new_df])
    panda_obj      = Panda(expression, motif_data, ppi_data_symm, save_tmp=False, remove_missing=rm_missing,
                      keep_expression_matrix=bool(lioness_file), modeProcess='intersection')
    panda_obj.save_panda_results(output_file)
    res=pd.read_csv(output_file, sep=' ', header=None)
    gt =pd.read_csv(gt_file_inter, sep=' ', header=None)
    pd.testing.assert_frame_equal(res,gt,check_less_precise=False,check_exact=False)

    #1. Union
    panda_obj      = Panda(expression_data, motif, ppi, save_tmp=False, remove_missing=rm_missing,
                      keep_expression_matrix=bool(lioness_file), modeProcess='union')
    panda_obj.save_panda_results(output_file)
    res=pd.read_csv(output_file, sep=' ', header=None)
    gt =pd.read_csv(gt_file, sep=' ', header=None)
    pd.testing.assert_frame_equal(res,gt,check_less_precise=False,check_exact=False)

    #2. In-degree and out-degree
    panda_obj      = Panda(expression_data, motif, ppi, save_tmp=False, remove_missing=rm_missing,
                       keep_expression_matrix=bool(lioness_file), modeProcess='union', save_memory=False)
    panda_obj.return_panda_indegree()
    panda_obj.return_panda_outdegree()
    # Lazy test
    assert (np.round(panda_obj.panda_indegree.iloc[0].loc['force'], 5) == 1.13971)
    assert (np.round(panda_obj.panda_outdegree.iloc[0].loc['force'], 5) == 1030.06840)

    #3. Legacy
    panda_obj = Panda(expression_data, motif, ppi, save_tmp=True, remove_missing=rm_missing,
                      keep_expression_matrix=True, save_memory=True, modeProcess='legacy')
    panda_obj.save_panda_results(output_file)
    gt_file         ='tests/panda/legacy_test_panda.txt'
    res = pd.read_csv(output_file, sep=' ', header=None)
    gt = pd.read_csv(gt_file, sep=' ', header=None)
    pd.testing.assert_frame_equal(res,gt,check_less_precise=False,check_exact=False)
    print('Test panda passed was successful!')

    #4. None Types
    i=0
    gt_test_panda = 'gt_panda'
    test_panda = 'test_panda'
    for modeProcess in ['legacy','union','intersection']:
        print(modeProcess)
        #Motif
        i=i+1
        panda_obj = Panda(expression_data, None, ppi, save_tmp=True, remove_missing=rm_missing,
                      keep_expression_matrix=True, save_memory=True, modeProcess=modeProcess)
        panda_obj.save_panda_results(test_panda + str(i) + '.txt')
        res = pd.read_csv(test_panda + str(i) + '.txt', sep=' ', header=None)
        os.system('curl -O https://netzoo.s3.us-east-2.amazonaws.com/netZooPy/tutorial_datasets/'+gt_test_panda + str(i) + '.txt')
        gt = pd.read_csv(gt_test_panda + str(i) + '.txt', sep=' ', header=None)
        pd.testing.assert_frame_equal(res,gt,check_less_precise=False,check_exact=False)
        #PPI
        i=i+1
        panda_obj = Panda(expression_data, motif, None, save_tmp=True, remove_missing=rm_missing,
                      keep_expression_matrix=True, save_memory=True, modeProcess=modeProcess)
        panda_obj.save_panda_results(test_panda + str(i) + '.txt')
        res = pd.read_csv(test_panda + str(i) + '.txt', sep=' ', header=None)
        os.system('curl -O https://netzoo.s3.us-east-2.amazonaws.com/netZooPy/tutorial_datasets/'+gt_test_panda + str(i) + '.txt')
        gt = pd.read_csv(gt_test_panda + str(i) + '.txt', sep=' ', header=None)
        pd.testing.assert_frame_equal(res,gt,check_less_precise=False,check_exact=False)
        #Expression
        i=i+1
        panda_obj = Panda(None, motif, ppi, save_tmp=True, remove_missing=rm_missing,
                      keep_expression_matrix=True, save_memory=True, modeProcess=modeProcess)
        panda_obj.save_panda_results(test_panda + str(i) + '.txt')
        res = pd.read_csv(test_panda + str(i) + '.txt', sep=' ', header=None)
        os.system('curl -O https://netzoo.s3.us-east-2.amazonaws.com/netZooPy/tutorial_datasets/'+gt_test_panda + str(i) + '.txt')
        gt = pd.read_csv(gt_test_panda + str(i) + '.txt', sep=' ', header=None)
        pd.testing.assert_frame_equal(res,gt,check_less_precise=False,check_exact=False)
        #Expression and PPI
        i=i+1
        panda_obj = Panda(None, motif, None, save_tmp=True, remove_missing=rm_missing,
                      keep_expression_matrix=True, save_memory=True, modeProcess=modeProcess)
        panda_obj.save_panda_results(test_panda + str(i) + '.txt')
        res = pd.read_csv(test_panda + str(i) + '.txt', sep=' ', header=None)
        os.system('curl -O https://netzoo.s3.us-east-2.amazonaws.com/netZooPy/tutorial_datasets/'+gt_test_panda + str(i) + '.txt')
        gt = pd.read_csv(gt_test_panda + str(i) + '.txt', sep=' ', header=None)
        pd.testing.assert_frame_equal(res,gt,check_less_precise=False,check_exact=False)
