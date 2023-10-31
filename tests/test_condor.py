import pytest
import os
from netZooPy import condor
import pandas as pd
import numpy as np
import random
import subprocess
import netZooPy.command_line as cmd

np.random.seed(0)
import random
random.seed(0)

def test_condor():
    print('Start Condor run ...')
    random.seed(10)
    network        ='tutorials/condor/toynetwork.csv'
    output_file1   ='tar_memb.txt'
    output_file2   ='reg_memb.txt'
    gt_file1       ='tests/condor/gh_tar_memb.txt'
    gt_file2       ='tests/condor/gh_reg_memb.txt'

    print('finished extra comparison')
    condor.run_condor(network)
    res1=pd.read_csv(output_file1, sep=',', header=None)
    res2=pd.read_csv(output_file2, sep=',', header=None)
    gt1 =pd.read_csv(gt_file1, sep=',', header=None)
    gt2 =pd.read_csv(gt_file2, sep=',', header=None)
    pd.testing.assert_frame_equal(res1,gt1,check_exact=False)
    pd.testing.assert_frame_equal(res2,gt2,check_exact=False)
    print('finished condor test')

    np.random.seed(0)
    random.seed(0)

    ## Test command line call
    result = subprocess.run(["netzoopy", "condor", "--help"], capture_output=False)
    assert result.returncode == 0

    # 1. Test command line
    cmd.condor.callback(network, tar_output="tar2.txt",reg_output="reg2.txt")
    res1=pd.read_csv('tar2.txt', sep=',', header=None)
    res2=pd.read_csv('reg2.txt', sep=',', header=None)
    gt1 =pd.read_csv(gt_file1, sep=',', header=None)
    gt2 =pd.read_csv(gt_file2, sep=',', header=None)

    print(res1.columns)
    print(gt1.columns)
    pd.testing.assert_frame_equal(res1,gt1,check_exact=False)
    pd.testing.assert_frame_equal(res2,gt2,check_exact=False)

    print('finished comparison with command line')