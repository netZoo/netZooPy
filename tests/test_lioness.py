import pytest
import os
from netZooPy.lioness.lioness import Lioness
from netZooPy.panda.panda import Panda
import pandas as pd
import numpy as np
import glob


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
    lioness_obj = Lioness(panda_obj, start=1, end = 4)
    lioness_obj.export_lioness_table(output_table)
    # Check that the correlation matrix for Panda and Lioness are the same (there is no update by Panda or Lioness)
    np.allclose(panda_obj.correlation_matrix, lioness_obj.correlation_matrix)
    np.allclose(panda_obj.motif_matrix, lioness_obj.motif_matrix)
    np.allclose(panda_obj.ppi_matrix, lioness_obj.ppi_matrix)

    # Compare with R
    rdf = pd.read_csv(toy_r_file, sep = ' ', header=None)
    pydf = pd.read_csv(output_table, sep = ' ', header=None).iloc[:,0:3]

    pd.testing.assert_frame_equal(rdf, pydf, rtol=5e-1, atol= 99e-2,  check_exact=False)


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
    lioness_obj_2 = Lioness(panda_obj_2, start=1, end=1)
    # lioness_obj.save_lioness_results(lioness_file)
    # Read first lioness network
    res = np.load("lioness_output/lioness.1.npy")
    gt = np.load("tests/lioness/lioness.1.coexpression.npy")
    # Compare to ground truth
    assert np.allclose(gt, res)
    
    # 3. Testing Lioness in parallel
    # Set parameters
    os.remove("lioness_output/lioness.1.npy")
    lioness_obj = Lioness(panda_obj, ncores=2, start=1, end=2)
    # lioness_obj.save_lioness_results(lioness_file)
    # Read first lioness network
    res = np.load("lioness_output/lioness.1.npy")
    gt = np.load("tests/lioness/lioness.1.npy")
    # Compare to ground truth
    assert np.allclose(gt, res)
    