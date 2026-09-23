import os

import numpy as np
import pandas as pd
import pytest

from netZooPy.spider.spider import Spider


DATA = "tests/spider/ToyData"
EXPR = f"{DATA}/ToyExpressionData.txt"
MOTIF = f"{DATA}/ToyMotifData.txt"
EPIFILTER = f"{DATA}/ToyEpiFilterData.txt"
PPI = f"{DATA}/ToyPPIData.txt"


def test_spider_runs():
    """SPIDER produces an (nTFs x nGenes) network with no NaNs."""
    obj = Spider(EXPR, MOTIF, EPIFILTER, PPI,
                 save_memory=False, modeProcess="union")
    net = np.asarray(obj.spider_network)
    assert net.ndim == 2
    assert net.shape == (obj.num_tfs, obj.num_genes)
    assert not np.isnan(net).any()


def test_no_epifilter_equals_all_ones():
    """epifilter=None must equal an all-ones filter (SPIDER's null mask)."""
    motif = pd.read_csv(MOTIF, sep="\t", header=None)
    ones = motif.copy()
    ones.iloc[:, 2] = 1.0

    none_obj = Spider(EXPR, MOTIF, None, PPI,
                      save_memory=False, modeProcess="union")
    ones_obj = Spider(EXPR, MOTIF, ones, PPI,
                      save_memory=False, modeProcess="union")
    assert np.allclose(none_obj.spider_network, ones_obj.spider_network)


def test_epifilter_changes_network():
    """A non-trivial epifilter must change the inferred network."""
    with_f = Spider(EXPR, MOTIF, EPIFILTER, PPI,
                    save_memory=False, modeProcess="union")
    no_f = Spider(EXPR, MOTIF, None, PPI,
                  save_memory=False, modeProcess="union")
    assert not np.allclose(with_f.spider_network, no_f.spider_network)


def test_epifilter_size_mismatch_raises():
    """A filter that does not match the motif rows must raise."""
    motif = pd.read_csv(MOTIF, sep="\t", header=None)
    bad = motif.iloc[:-2].copy()
    with pytest.raises(ValueError):
        Spider(EXPR, MOTIF, bad, PPI, modeProcess="union")


def test_export_and_save(tmp_path):
    """Edge-list export exists (save_memory=False) and save works."""
    obj = Spider(EXPR, MOTIF, EPIFILTER, PPI,
                 save_memory=False, modeProcess="union")
    assert hasattr(obj, "export_spider_results")
    out = str(tmp_path / "spider_out.txt")
    obj.save_spider_results(out)
    assert os.path.exists(out)
