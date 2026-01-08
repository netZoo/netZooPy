import pytest
import os
from netZooPy.panda.panda import Panda
from netZooPy.giraffe import giraffe
import pandas as pd
import numpy as np
import subprocess
import netZooPy.command_line as cmd

# Use the logger
import logging
LOGGER = logging.getLogger(__name__)

def test_giraffe():
    """First giraffe test with ToyData
    """
    LOGGER.warning('Test1')

    ppi_file = "tests/puma/ToyData/ToyPPIData.txt"
    motif_file = "tests/puma/ToyData/ToyMotifData.txt"
    expression_file = "tests/puma/ToyData/ToyExpressionData.txt"

    tfa_file = "tests/giraffe/Toygiraffe_TFA_hat.txt"
    r_file = "tests/giraffe/Toygiraffe_R_hat.txt"


    panda_obj = Panda(
            expression_file,
            motif_file,
            ppi_file,
            modeProcess="intersection",
            with_header=False,
            process_data_only = True)

    # get the preprocessed data from panda object
    expression = panda_obj.expression
    motif = panda_obj.motif_matrix_unnormalized
    ppi = panda_obj.ppi_matrix

    # Run GIRAFFE
    giraffe_model = giraffe.Giraffe(expression, motif.T, ppi)

    R_hat = giraffe_model.get_regulation() # Size (G, TF)
    TFA_hat = giraffe_model.get_tfa() # Size (TF, n)

    LOGGER.warning('Testing Rhat')
    Rres = pd.read_csv(r_file, sep="\t", index_col=0)
    # check Rres is same size
    assert Rres.shape == R_hat.shape
    # check values are close
    np.testing.assert_allclose(Rres.values, R_hat, atol=1e-5)

    LOGGER.warning('Testing TFAhat')
    TFAres = pd.read_csv(tfa_file, sep="\t", index_col=0, header = None)
    TFAres = TFAres
    # check TFAres is same size
    assert TFAres.shape == TFA_hat.shape
    # check values are close
    np.testing.assert_allclose(TFAres.values, TFA_hat, atol=1e-5)
