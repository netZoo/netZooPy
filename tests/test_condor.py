"""
Unit tests for CONDOR.

CONDOR seeds its community structure with a *randomised* igraph algorithm
(Leiden by default), so the integer community labels it produces are
arbitrary: the same partition can come back with its labels permuted from one
igraph release to the next (see issues #312 and #320). These tests therefore
never compare raw labels. They check

* the structure of the output (node sets, one community per node),
* that the partition matches the reference partition up to a relabelling
  (adjusted Rand index), and
* that the final bipartite modularity matches the reference,

and that results are stable across seeds and reproducible for a fixed seed.
"""
import os
import random
import subprocess

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

from netZooPy import condor
import netZooPy.command_line as cmd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
NETWORK = os.path.join(REPO, "tutorials", "condor", "toynetwork.csv")
GT_TAR = os.path.join(HERE, "condor", "gh_tar_memb.txt")
GT_REG = os.path.join(HERE, "condor", "gh_reg_memb.txt")

# Bipartite modularity of the reference partition stored in the ground-truth files.
REFERENCE_MODULARITY = 0.526667
MODULARITY_TOL = 1e-4
# Minimum adjusted Rand index for a partition to count as matching the reference.
MIN_ARI = 0.9


def _seed(seed):
    random.seed(seed)
    np.random.seed(seed)


def _run(seed, tar_out, reg_out):
    """Run CONDOR on the toy network with a fixed seed and return the condor object."""
    _seed(seed)
    return condor.run_condor(
        NETWORK, return_output=True, silent=True, tar_output=tar_out, reg_output=reg_out
    )


def _labels(memb, node_col):
    """Community labels ordered by node name, so two partitions line up node-for-node."""
    return memb.sort_values(node_col)["community"].to_numpy()


def _read_memb(path, node_col):
    df = pd.read_csv(path, index_col=0)
    assert list(df.columns) == [node_col, "community"]
    return df


def _assert_partition_matches(memb, gt_path, node_col):
    gt = _read_memb(gt_path, node_col)
    # Same nodes, each assigned to exactly one community.
    assert sorted(memb[node_col]) == sorted(gt[node_col])
    assert memb[node_col].is_unique
    # Same partition up to a relabelling of the communities.
    ari = adjusted_rand_score(_labels(gt, node_col), _labels(memb, node_col))
    assert ari >= MIN_ARI, f"partition differs from the reference (ARI={ari:.3f})"


def _assert_structure(co):
    net = pd.read_csv(NETWORK, index_col=0)
    reg_nodes = {"reg_" + str(n) for n in net.iloc[:, 0]}
    tar_nodes = {"tar_" + str(n) for n in net.iloc[:, 1]}
    assert set(co.reg_memb["reg"]) == reg_nodes
    assert set(co.tar_memb["tar"]) == tar_nodes
    for memb in (co.tar_memb, co.reg_memb):
        assert (memb["community"] >= 0).all()
        assert memb["community"].nunique() >= 2


def test_condor(tmp_path):
    tar_out, reg_out = str(tmp_path / "tar_memb.txt"), str(tmp_path / "reg_memb.txt")
    co = _run(10, tar_out, reg_out)

    _assert_structure(co)
    assert abs(co.modularity - REFERENCE_MODULARITY) <= MODULARITY_TOL
    _assert_partition_matches(co.tar_memb, GT_TAR, "tar")
    _assert_partition_matches(co.reg_memb, GT_REG, "reg")

    # The files written to disk carry the same membership as the returned object.
    _assert_partition_matches(_read_memb(tar_out, "tar"), GT_TAR, "tar")
    _assert_partition_matches(_read_memb(reg_out, "reg"), GT_REG, "reg")


def test_condor_reproducible_for_fixed_seed(tmp_path):
    runs = [_run(10, str(tmp_path / f"t{i}.txt"), str(tmp_path / f"r{i}.txt")) for i in range(2)]
    np.testing.assert_array_equal(_labels(runs[0].tar_memb, "tar"), _labels(runs[1].tar_memb, "tar"))
    np.testing.assert_array_equal(_labels(runs[0].reg_memb, "reg"), _labels(runs[1].reg_memb, "reg"))


def test_condor_stable_across_seeds(tmp_path):
    """The stochastic initialisation must not change the quality of the result."""
    runs = [_run(s, str(tmp_path / f"t{s}.txt"), str(tmp_path / f"r{s}.txt")) for s in (0, 42, 123)]
    for co in runs:
        assert abs(co.modularity - REFERENCE_MODULARITY) <= MODULARITY_TOL
    ref_tar, ref_reg = _labels(runs[0].tar_memb, "tar"), _labels(runs[0].reg_memb, "reg")
    for co in runs[1:]:
        assert adjusted_rand_score(ref_tar, _labels(co.tar_memb, "tar")) >= MIN_ARI
        assert adjusted_rand_score(ref_reg, _labels(co.reg_memb, "reg")) >= MIN_ARI


def test_condor_command_line(tmp_path):
    result = subprocess.run(["netzoopy", "condor", "--help"], capture_output=True)
    assert result.returncode == 0

    tar_out, reg_out = str(tmp_path / "tar2.txt"), str(tmp_path / "reg2.txt")
    _seed(10)
    # The CLI defaults to Louvain ("LCS"), whose result varies with the seed and the
    # igraph version; Leiden reaches the reference partition for every seed.
    cmd.condor.callback(NETWORK, initial_method="LDN", tar_output=tar_out, reg_output=reg_out)
    _assert_partition_matches(_read_memb(tar_out, "tar"), GT_TAR, "tar")
    _assert_partition_matches(_read_memb(reg_out, "reg"), GT_REG, "reg")
