from netZooPy import cobra
import numpy as np
import pandas as pd


def test_cobra():
    print("Start COBRA run ...")

    # Read COBRA input
    X_path = "tests/cobra/X.csv"
    X = pd.read_csv(X_path, index_col=0)
    expression_path = "tests/cobra/expression.csv"
    expression = pd.read_csv(expression_path, index_col=0)

    # Read ground-truth
    psi_gt_path = "tests/cobra/psi.csv"
    psi_gt = pd.read_csv(psi_gt_path, index_col=0)
    Q_gt_path = "tests/cobra/Q.csv"
    Q_gt = pd.read_csv(Q_gt_path, index_col=0)
    D_gt_path = "tests/cobra/D.csv"
    D_gt = pd.read_csv(D_gt_path, index_col=0)
    G_gt_path = "tests/cobra/G.csv"
    G_gt = pd.read_csv(G_gt_path, index_col=0)

    # Call COBRA
    psi, Q, D, G = cobra.cobra(X, expression)

    # Cast output to pandas
    psi = pd.DataFrame(psi, index=psi_gt.index, columns=psi_gt.columns)
    Q = pd.DataFrame(Q, index=Q_gt.index, columns=Q_gt.columns)
    D = pd.DataFrame(D, index=D_gt.index, columns=D_gt.columns)
    G = pd.DataFrame(G, index=G_gt.index, columns=G_gt.columns)

    # Test for equality
    pd.testing.assert_frame_equal(psi, psi_gt, rtol=1e-10, check_exact=False)
    pd.testing.assert_frame_equal(D, D_gt, rtol=1e-10, check_exact=False)
    pd.testing.assert_frame_equal(G, G_gt, rtol=1e-10, check_exact=False)

    q = psi.shape[0]
    for i in range(q):
        C = Q.to_numpy().dot(np.mean(X, axis=0)[i] * np.diag(psi.to_numpy()[i, :])).dot(Q.to_numpy().T)
        C_gt = Q_gt.to_numpy().dot(np.mean(X, axis=0)[i] * np.diag(psi_gt.to_numpy()[i, :])).dot(Q_gt.to_numpy().T)
        pd.testing.assert_frame_equal(pd.DataFrame(C), pd.DataFrame(C_gt), rtol=1e-10, check_exact=False)
