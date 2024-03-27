import pytest
import os
from netZooPy.panda.panda import Panda
import pandas as pd
import numpy as np
import subprocess
import netZooPy.command_line as cmd


import logging

LOGGER = logging.getLogger(__name__)



def runPandatest(modeProcess,ppi,motif,expression_data,lioness_file,rm_missing,output_file,gt_file):
    panda_obj = Panda(
        expression_data,
        motif,
        ppi,
        save_tmp=False,
        remove_missing=rm_missing,
        keep_expression_matrix=bool(lioness_file),
        modeProcess=modeProcess,
    )
    res = panda_obj.panda_network
    gt = pd.read_csv(gt_file, sep=",", header=0, index_col=0)
    pd.testing.assert_frame_equal(res, gt, rtol=1e-7, atol=1e-7, check_exact=False, check_names=False)
    return

# PUMA toy data set has 1000 expression genes and 913 motif gene. The 913 motif genes are a subset of the 1000 expression
# genes
def test_panda():
    # print(os.getcwd())
    print("Start Panda run ...")
    ppi = "tests/puma/ToyData/ToyPPIData.txt"
    motif = "tests/puma/ToyData/ToyMotifData.txt"
    expression_data = "tests/puma/ToyData/ToyExpressionData.txt"
    lioness_file = ""
    rm_missing = False
    output_file = "travis_test_panda.txt"
    gt_file = "tests/panda/union_test_panda.txt"
    gt_file_inter = "tests/panda/inter_test_panda.txt"
    gt_file_rm = "tests/panda/rm_test_panda.txt"

    ## Test command line call
    result = subprocess.run(["netzoopy", "panda", "--help"], capture_output=False)
    assert result.returncode == 0
    
    LOGGER.warning('Test1')

    # 1. Intersection
    panda_obj = Panda(
        expression_data,
        motif,
        ppi,
        save_tmp=False,
        remove_missing=rm_missing,
        keep_expression_matrix=bool(lioness_file),
        modeProcess="intersection",
    )
   
    panda_obj.save_panda_results(output_file)
    res = pd.read_csv(output_file, sep=" ", header=None)
    gt = pd.read_csv(gt_file_inter, sep=" ", header=None)
    pd.testing.assert_frame_equal(res, gt, rtol=1e-12, check_exact=False)

    LOGGER.warning('Test1b')
    # 1.b Intersection from command line
    cmd.panda.callback(expression_data, motif, ppi, output_file, rm_missing = rm_missing, keep_expr=bool(lioness_file),mode_process='intersection', save_memory = True, old_compatible = True)
    res = pd.read_csv(output_file, sep=" ", header=None)
    gt = pd.read_csv(gt_file_inter, sep=" ", header=None)
    pd.testing.assert_frame_equal(res, gt, rtol=1e-12, check_exact=False)

    LOGGER.warning('Test1.1')
    # 1.1 Intersection via DataFrame
    expression = pd.read_csv(expression_data, sep="\t", header=None, index_col=0)
    motif_data = pd.read_csv(motif, sep="\t", names=["source", "target", "w"])
    ppi_data = pd.read_csv(ppi, sep="\t", header=None)

    panda_obj = Panda(
        expression,
        motif_data,
        ppi_data,
        save_tmp=False,
        remove_missing=rm_missing,
        keep_expression_matrix=bool(lioness_file),
        modeProcess="intersection",
    )
    panda_obj.save_panda_results(output_file)
    res = pd.read_csv(output_file, sep=" ", header=None)
    gt = pd.read_csv(gt_file_inter, sep=" ", header=None)
    pd.testing.assert_frame_equal(res, gt, rtol=1e-12, check_exact=False)
    LOGGER.warning('Test1.2')
    # 1.2 Intersection with symmetric PPI
    ppi_data = pd.read_csv(ppi, sep="\t", header=None)
    new_df = pd.DataFrame(data={0: ppi_data[0], 1: ppi_data[1], 2: ppi_data[2]})
    ppi_data_symm = pd.concat([ppi_data, new_df])
    panda_obj = Panda(
        expression,
        motif_data,
        ppi_data_symm,
        save_tmp=False,
        remove_missing=rm_missing,
        keep_expression_matrix=bool(lioness_file),
        modeProcess="intersection",
    )
    panda_obj.save_panda_results(output_file)
    res = pd.read_csv(output_file, sep=" ", header=None)
    gt = pd.read_csv(gt_file_inter, sep=" ", header=None)
    pd.testing.assert_frame_equal(res, gt, rtol=1e-12, check_exact=False)

    # 2. Union
    panda_obj = Panda(
        expression_data,
        motif,
        ppi,
        save_tmp=False,
        remove_missing=rm_missing,
        keep_expression_matrix=bool(lioness_file),
        modeProcess="union",
    )
    panda_obj.save_panda_results(output_file)
    res = pd.read_csv(output_file, sep=" ", header=None)
    gt = pd.read_csv(gt_file, sep=" ", header=None)
    pd.testing.assert_frame_equal(res, gt, rtol=1e-12, check_exact=False)


    # 2.b Union from command line

    cmd.panda.callback(expression_data, motif, ppi, output_file, rm_missing = rm_missing, keep_expr=bool(lioness_file),mode_process='union', save_memory = True, old_compatible = True)

    res = pd.read_csv(output_file, sep=" ", header=None)
    gt = pd.read_csv(gt_file, sep=" ", header=None)
    pd.testing.assert_frame_equal(res, gt, rtol=1e-12, check_exact=False)


    # 3. In-degree and out-degree
    panda_obj = Panda(
        expression_data,
        motif,
        ppi,
        save_tmp=False,
        remove_missing=rm_missing,
        keep_expression_matrix=bool(lioness_file),
        modeProcess="union",
        save_memory=False,
    )
    panda_obj.return_panda_indegree()
    panda_obj.return_panda_outdegree()
    # Lazy test
    assert (panda_obj.panda_indegree.iloc[0].loc["force"]*1e5).astype(int)/1e5 == 1.13969
    assert (panda_obj.panda_outdegree.iloc[0].loc["force"]*1e5).astype(int)/1e5 == 1030.06841

    # 4. Legacy
    panda_obj = Panda(
        expression_data,
        motif,
        ppi,
        save_tmp=True,
        remove_missing=rm_missing,
        keep_expression_matrix=True,
        save_memory=True,
        modeProcess="legacy",
    )
    #panda_obj.save_panda_results(output_file)
    #gt_file = "tests/panda/legacy_test_panda.txt"
    gt_file = "tests/panda/panda_gt_matlab.csv"
    res = panda_obj.panda_network
    gt = pd.read_csv(gt_file, sep=",", index_col=0, header=0)
    pd.testing.assert_frame_equal(res, gt, rtol=1e-12, atol=1e-12, check_exact=False, check_names=False)
    print("Test panda passed was successful!")

    # 3' Legacy with rm_missing=True
    panda_obj = Panda(
        expression_data,
        motif,
        ppi,
        save_tmp=True,
        remove_missing=True,
        keep_expression_matrix=True,
        save_memory=True,
        modeProcess="legacy",
    )
    panda_obj.save_panda_results(output_file)
    res = pd.read_csv(output_file, sep=" ", header=None)
    gt = pd.read_csv(gt_file_rm, sep=" ", header=None)
    pd.testing.assert_frame_equal(res, gt, rtol=1e-12, check_exact=False)
    print("Test panda passed was successful!")

    # 4. None Types
    i = 0
    gt_test_panda = "gt_panda"
    test_panda = "test_panda"
    for modeProcess in ["legacy", "union", "intersection"]:
        print(modeProcess)
        # Motif
        i = i + 1
        panda_obj = Panda(
            expression_data,
            None,
            ppi,
            save_tmp=True,
            remove_missing=rm_missing,
            keep_expression_matrix=True,
            save_memory=True,
            modeProcess=modeProcess,
        )
        panda_obj.save_panda_results(test_panda + str(i) + ".txt")
        res = pd.read_csv(test_panda + str(i) + ".txt", sep=" ", header=None)
        os.system(
            "curl -O https://netzoo.s3.us-east-2.amazonaws.com/netZooPy/tutorial_datasets/"
            + gt_test_panda
            + str(i)
            + ".txt"
        )
        gt = pd.read_csv(gt_test_panda + str(i) + ".txt", sep=" ", header=None)
        pd.testing.assert_frame_equal(res, gt, rtol=1e-12, check_exact=False)
        # PPI
        i = i + 1
        panda_obj = Panda(
            expression_data,
            motif,
            None,
            save_tmp=True,
            remove_missing=rm_missing,
            keep_expression_matrix=True,
            save_memory=True,
            modeProcess=modeProcess,
        )
        panda_obj.save_panda_results(test_panda + str(i) + ".txt")
        res = pd.read_csv(test_panda + str(i) + ".txt", sep=" ", header=None)
        os.system(
            "curl -O https://netzoo.s3.us-east-2.amazonaws.com/netZooPy/tutorial_datasets/"
            + gt_test_panda
            + str(i)
            + ".txt"
        )
        gt = pd.read_csv(gt_test_panda + str(i) + ".txt", sep=" ", header=None)
        pd.testing.assert_frame_equal(res, gt, rtol=1e-12, check_exact=False)
        # Expression
        i = i + 1
        panda_obj = Panda(
            None,
            motif,
            ppi,
            save_tmp=True,
            remove_missing=rm_missing,
            keep_expression_matrix=True,
            save_memory=True,
            modeProcess=modeProcess,
        )
        panda_obj.save_panda_results(test_panda + str(i) + ".txt")
        res = pd.read_csv(test_panda + str(i) + ".txt", sep=" ", header=None)
        os.system(
            "curl -O https://netzoo.s3.us-east-2.amazonaws.com/netZooPy/tutorial_datasets/"
            + gt_test_panda
            + str(i)
            + ".txt"
        )
        gt = pd.read_csv(gt_test_panda + str(i) + ".txt", sep=" ", header=None)
        pd.testing.assert_frame_equal(res, gt, rtol=1e-12, check_exact=False)
        # Expression and PPI
        i = i + 1
        panda_obj = Panda(
            None,
            motif,
            None,
            save_tmp=True,
            remove_missing=rm_missing,
            keep_expression_matrix=True,
            save_memory=True,
            modeProcess=modeProcess,
        )
        panda_obj.save_panda_results(test_panda + str(i) + ".txt")
        res = pd.read_csv(test_panda + str(i) + ".txt", sep=" ", header=None)
        os.system(
            "curl -O https://netzoo.s3.us-east-2.amazonaws.com/netZooPy/tutorial_datasets/"
            + gt_test_panda
            + str(i)
            + ".txt"
        )
        gt = pd.read_csv(gt_test_panda + str(i) + ".txt", sep=" ", header=None)
        pd.testing.assert_frame_equal(res, gt, rtol=1e-12, check_exact=False)

    # 5. pantests
    modeProcess = 'union'
    ## AgNet 0
    ppi = "tests/panda/ppi.txt"
    motif = "tests/panda/motif.txt"
    expression_data = "tests/panda/expression.txt"
    lioness_file = ""
    rm_missing = False
    output_file = "pan_test_panda.txt"
    gt_file = "tests/panda/uAgNet0.csv"
    runPandatest(modeProcess,ppi,motif,expression_data,lioness_file,rm_missing,output_file,gt_file)
    ## AgNet 1
    expression_data = "tests/panda/expression_down.txt"
    gt_file = "tests/panda/uAgNet1.csv"
    runPandatest(modeProcess,ppi,motif,expression_data,lioness_file,rm_missing,output_file,gt_file)
    ## AgNet 2
    expression_data = "tests/panda/expression.txt"
    motif = 'tests/panda/motif_down_gene.txt';
    gt_file = "tests/panda/uAgNet2.csv"
    runPandatest(modeProcess,ppi,motif,expression_data,lioness_file,rm_missing,output_file,gt_file)
    ## AgNet 3
    motif = 'tests/panda/motif_down_tf.txt';
    gt_file = "tests/panda/uAgNet3.csv"
    runPandatest(modeProcess,ppi,motif,expression_data,lioness_file,rm_missing,output_file,gt_file)
    ## AgNet 4
    motif = "tests/panda/motif.txt"
    ppi = "tests/panda/ppi_down.txt"
    gt_file = "tests/panda/uAgNet4.csv"
    runPandatest(modeProcess,ppi,motif,expression_data,lioness_file,rm_missing,output_file,gt_file)
    ## AgNet 5
    ppi = "tests/panda/ppi.txt"
    expression_data = "tests/panda/expression_up.txt"
    gt_file = "tests/panda/uAgNet5.csv"
    runPandatest(modeProcess,ppi,motif,expression_data,lioness_file,rm_missing,output_file,gt_file)
    ## AgNet 6
    expression_data = "tests/panda/expression.txt"
    motif = 'tests/panda/motif_up_gene.txt';
    gt_file = "tests/panda/uAgNet6.csv"
    runPandatest(modeProcess,ppi,motif,expression_data,lioness_file,rm_missing,output_file,gt_file)
    ## AgNet 7
    motif = 'tests/panda/motif_up_tf.txt';
    gt_file = "tests/panda/uAgNet7.csv"
    runPandatest(modeProcess,ppi,motif,expression_data,lioness_file,rm_missing,output_file,gt_file)
    ## AgNet 8
    motif = "tests/panda/motif.txt"
    ppi = "tests/panda/ppi_up.txt"
    gt_file = "tests/panda/uAgNet8.csv"
    runPandatest(modeProcess,ppi,motif,expression_data,lioness_file,rm_missing,output_file,gt_file)

    # 6. test square nonsymmetric matrices
    W = np.array([[1., 1., 0.],
                  [0., 0., 1.],
                  [1., 2., 3.]])
    W_gt = np.array([[ 1.        ,  0.5       , -1.75592895],
                     [-1.5       , -1.3660254 ,  0.81101776],
                     [-0.3660254 ,  0.8660254 ,  1.81093659]])
    W_res = Panda._normalize_network(self=[],x=W)
    assert(np.allclose(W_gt, W_res, rtol=1e-05, atol=1e-08))


def test_incorrect_gene_ids_raises_exception():
    '''
    This test checks that an appropriately descriptive exception
    is raised if one attempts to run PANDA on an expression matrix
    that does not have genes intersecting with those in the motif
    prior matrix
    '''

    # these numbers are arbitrary but don't need to be large
    # to trigger the exception since it would be raised prior
    # to any Panda iterations.
    num_genes = 20 # number of genes in the expression mtx
    num_samples = 30 # number of samples in the expression mtx
    num_tfs = 10 # number of transcription factors in the motif prior
    num_motif_genes = 10 # number of genes in the motif prior

    # create a dataframe of mock expression data
    exp_mtx_genes = [f'g{_}' for _ in range(num_genes)]
    exp_mtx_samples = [f's{_}' for _ in range(num_samples)]
    expression_data = pd.DataFrame(
        np.random.randint(0, 100, size=(num_genes, num_samples)),
        index=exp_mtx_genes,
        columns=exp_mtx_samples
    )

    # mock motif prior. Note that the genes (2nd col) start with 'G'
    # instead of 'g' in the expression matrix
    tf_list = [f'tf{_}' for _ in range(num_tfs)]
    motif_gene_list = [f'G{_}' for _ in range(num_motif_genes)]
    motif_data = pd.DataFrame(
        {
            0: np.repeat(tf_list, num_motif_genes),
            1: np.tile(motif_gene_list, num_tfs),
            2: np.random.random(num_tfs*num_motif_genes)
        }
    )

    # mock ppi prior
    ppi_data = pd.DataFrame(
        {
            0: np.repeat(tf_list, num_tfs),
            1: np.tile(tf_list, num_tfs),
            2: np.random.random(num_tfs**2)
        }
    )

    with pytest.raises(Exception,
                       match='Error when creating the motif network!.*'):
        panda_obj1 = Panda(
            expression_data,
            motif_data,
            ppi_data,
            modeProcess='legacy'
        )

    with pytest.raises(Exception,
                       match='Error when creating the motif network!.*'):
            panda_obj2 = Panda(
            expression_data,
            motif_data,
            ppi_data,
            modeProcess='intersection'
        )
