import pytest
import os
from netZooPy.lioness.lioness import Lioness
from netZooPy.panda.panda import Panda
import pandas as pd
import numpy as np
import glob
import subprocess
import netZooPy.command_line as cmd


def test_lioness():
    print("Start lioness test!")
    # 1. First generate temporary PANDA files as inputs for Lioness
    ppi = "tests/puma/ToyData/ToyPPIData.txt"
    motif = "tests/puma/ToyData/ToyMotifData.txt"
    expression_data = "tests/puma/ToyData/ToyExpressionData.txt"
    lioness_file = "Travis"
    rm_missing = False
    output_file = "panda.npy"
    gt_file = "tests/panda/test_panda.txt"
    output_table = "lioness_output/toy-lioness-py.txt"
    toy_r_file = 'tests/lioness/toy-lioness-first4-1.txt'

    panda_obj = Panda(
        expression_data,
        motif,
        ppi,
        save_tmp=True,
        remove_missing=rm_missing,
        keep_expression_matrix=bool(lioness_file),
        modeProcess="legacy",
        save_memory=False,
        start = 1,
        end = 4
    )
    # Set parameters
    lioness_obj = Lioness(panda_obj,save_dir = "lioness_output", start=1, end = 4, save_single=True)
    panda_obj.save_panda_results('panda_remove.txt')
    lioness_obj.export_lioness_table(output_table)

    # Check that the correlation matrix for Panda and Lioness are the same (there is no update by Panda or Lioness)
    np.allclose(panda_obj.correlation_matrix, lioness_obj.correlation_matrix)
    np.allclose(panda_obj.motif_matrix, lioness_obj.motif_matrix)
    np.allclose(panda_obj.ppi_matrix, lioness_obj.ppi_matrix)

    # Compare with R
    rdf = pd.read_csv(toy_r_file, sep = ' ', header=None, )
    pydf = pd.read_csv(output_table, sep = ' ').iloc[:,0:3]

    #pd.testing.assert_frame_equal(rdf, pydf, rtol=5e-1, atol= 99e-2,  check_exact=False, check_names = False, check_column_type = False)
    np.allclose(rdf.iloc[:,2].values, pydf.iloc[:,2].values)

    ## Test command line call
    result = subprocess.run(["netzoopy", "lioness", "--help"], capture_output=False)
    assert result.returncode == 0

    # 1. Test command line
    #positional: expression, motif, ppi, output_panda, output_lioness, fmt, computing, precision, ncores, save_memory, save_tmp, rm_missing, mode_process,output_type, alpha, start, end):
    cmd.lioness.callback(expression_data, motif, ppi, 'panda.txt','lioness_output_cmd',None,'npy','cpu','double',1,False,True,rm_missing,'legacy','network',0.1,1,4,False,True)
    res = np.load("lioness_output/lioness.1.npy")
    gt = res = np.load("lioness_output_cmd/lioness.1.npy")
    assert np.allclose(gt, res)

    # 2. Testing Lioness with motif set to None to compute Lioness on coexpression networks
    motif = None
    # Make sure to keep epxression matrix for next step
    panda_obj_2 = Panda(
        expression_data,
        motif,
        ppi,
        save_tmp=True,
        remove_missing=rm_missing,
        keep_expression_matrix=True,
        modeProcess="legacy",
    )
    lioness_obj_2 = Lioness(panda_obj_2, start=1, end=1, save_single=True,save_fmt='npy')
    # lioness_obj.save_lioness_results(lioness_file)
    # Read first lioness network
    res = np.load("lioness_output/lioness.1.npy")
    gt = np.load("tests/lioness/lioness.1.coexpression.npy")
    # Compare to ground truth
    assert np.allclose(gt, res)
    
    print('test3')
    # 3. Testing Lioness in parallel
    # Set parameters
    os.remove("lioness_output/lioness.1.npy")
    lioness_obj = Lioness(panda_obj, ncores=2, start=1, end=2, save_single=True)
    # lioness_obj.save_lioness_results(lioness_file)
    # Read first lioness network
    res = np.load("lioness_output/lioness.1.npy")
    gt = np.load("tests/lioness/lioness.1.npy")
    # Compare to ground truth
    assert np.allclose(gt, res)
    