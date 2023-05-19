import pytest
import os
from netZooPy import sambar
import pandas as pd
import numpy as np

def test_sambar():
    #print(os.getcwd())
    print('Start Sambar run ...')
    mut_file       ='tests/sambar/ToyData/mut.ucec.csv'
    cangenes       ='tests/sambar/ToyData/genes.txt'
    sign_file      ='tests/sambar/ToyData/h.all.v6.1.symbols.gmt'
    esize_file     ='tests/sambar/ToyData/esizef.csv'
    gt_file        ='tests/sambar/ToyData/sambar_gt.csv'

    #1. Sambar
    pathway_scores, cluster_groups = sambar.sambar(mut_file,esize_file,cangenes,sign_file) 
    #panda_obj.save_panda_results(output_file)
    #res=pd.read_csv(output_file, sep=' ', header=None)
    gt =pd.read_csv(gt_file, sep=',', header=0, index_col=0)
    pd.testing.assert_frame_equal(pathway_scores,gt,check_exact=False)

    #2. Sambar call without files
    a, b = sambar.sambar()
