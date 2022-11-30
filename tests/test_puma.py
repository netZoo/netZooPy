import pytest
from netZooPy.puma.puma import Puma
import pandas as pd
from netZooPy.lioness.lioness_for_puma import LionessPuma


def test_puma():
    # print(os.getcwd())
    print("Start Puma run ...")
    ppi = "tests/puma/ToyData/ToyPPIData.txt"
    motif = "tests/puma/ToyData/ToyMotifData.txt"
    expression_data = "tests/puma/ToyData/ToyExpressionData.txt"
    mir_file = "tests/puma/ToyData/ToyMiRList.txt"
    lioness_file = ""
    rm_missing = False
    output_file = "travis_test_puma.txt"
    gt_file = "tests/puma/matlablike_test_puma.txt"

    # 1. Vanilla puma
    puma_obj = Puma(
        expression_data,
        motif,
        ppi,
        mir_file,
        save_tmp=False,
        remove_missing=rm_missing,
        keep_expression_matrix=bool(lioness_file),
        modeProcess="legacy"
    )
    puma_obj.save_puma_results(output_file)
    res = pd.read_csv(output_file, sep=" ", header=None)
    gt = pd.read_csv(gt_file, sep=" ", header=None)
    pd.testing.assert_frame_equal(res, gt, rtol=1e-5, check_exact=False)

    # 1.b Testing if passing a mir list works
    with open(mir_file, "r") as f:
        miR = f.read().splitlines()
    puma_obj = Puma(
        expression_data,
        motif,
        ppi,
        miR,
        save_tmp=False,
        remove_missing=rm_missing,
        keep_expression_matrix=bool(lioness_file),
        modeProcess="legacy"
    )
    puma_obj.save_puma_results(output_file)
    res = pd.read_csv(output_file, sep=" ", header=None)
    gt = pd.read_csv(gt_file, sep=" ", header=None)
    pd.testing.assert_frame_equal(res, gt, rtol=1e-5, check_exact=False)

    # 2. with argument values
    rm_missing = False
    print("Test puma was successful!")


    # Let's at least check that Lioness for puma runs
    puma_obj = Puma(
        expression_data,
        motif,
        ppi,
        miR,
        save_tmp=False,
        remove_missing=rm_missing,
        keep_expression_matrix=True,
        modeProcess="legacy"
    )
    lioness_obj = LionessPuma(puma_obj, computing='cpu', end = 3)
