import pytest
import os
from netZooPy.ligress.ligress import Ligress
import pandas as pd
import numpy as np
import glob
import subprocess
import netZooPy.command_line as cmd


def test_lioness():
    print("Start ligress test!")
    # 1. First generate temporary PANDA files as inputs for Lioness
    ppi = "tests/puma/ToyData/ToyPPIData.txt"
    priors_table = "tests/ligress/prior_tab.csv"
    expression_data = "tests/ligress/expression_with_header_fm.txt"
    output_folder = "ligress/"
    c1_file = "tests/ligress/coexpression_c1.txt"
    c3_file = "tests/ligress/coexpression_c5.txt"
    p1_file = "tests/ligress/c1.txt"
    p3_file = "tests/ligress/c5.txt"
    output_c1 = "ligress/coexpression/coexpression_c1.txt"
    output_c3 = "ligress/coexpression/coexpression_c5.txt"
    output_p1 = "ligress/single_panda/c1.txt"
    output_p3 = "ligress/single_panda/c5.txt"


    # read expression data, prepare ppi+motif+expression universes
    ligress_obj = Ligress(expression_data, priors_table, ppi_file = ppi, output_folder=output_folder, mode_process='intersection')
    ligress_obj.run_ligress(keep_coexpression=True, save_memory=False, cores =1 ,computing_panda='cpu',alpha = 0.1, precision = 'double', th_motifs = 3)
    # Compare correlations
    c1df = pd.read_csv(c1_file, sep = ' ', index_col = 0)
    c3df = pd.read_csv(c3_file, sep = ' ', index_col = 0)
    c1odf = pd.read_csv(output_c1, sep = ' ', index_col = 0)
    c3odf = pd.read_csv(output_c3, sep = ' ', index_col = 0)

    pd.testing.assert_frame_equal(c1df, c1odf, rtol=5e-1, check_exact=False)
    pd.testing.assert_frame_equal(c3df, c3odf, rtol=5e-1, check_exact=False)

    ## Test command line call
    result = subprocess.run(["netzoopy", "ligress", "--help"], capture_output=False)
    assert result.returncode == 0

    # 1. Test command line
    #positional: expression, motif, ppi, output_panda, output_lioness, fmt, computing, precision, ncores, save_memory, save_tmp, rm_missing, mode_process,output_type, alpha, start, end):
    #cmd.lioness.callback(expression_data, motif, ppi, 'panda.txt','lioness_output_cmd','npy','cpu','double',1,False,True,rm_missing,'legacy','network',0.1,1,4)
    #res = np.load("lioness_output/lioness.1.npy")
    #gt = res = np.load("lioness_output_cmd/lioness.1.npy")
    #assert np.allclose(gt, res)

    # 2. Testing if panda are similar

    # Compare correlations
    p1df = pd.read_csv(p1_file, sep = '\t', index_col = 0)
    p3df = pd.read_csv(p3_file, sep = '\t', index_col = 0)
    p1odf = pd.read_csv(output_p1, sep = '\t', index_col = 0)
    p3odf = pd.read_csv(output_p3, sep = '\t', index_col = 0)

    pd.testing.assert_frame_equal(p1df, p1odf, rtol=5e-1,  atol= 99e-3, check_exact=False)
    pd.testing.assert_frame_equal(p3df, p3odf, rtol=5e-1,  atol= 99e-2, check_exact=False)

