import numpy as np
from scipy.linalg import eigh,pinv


def cobra(X, expression, standardize=True):
    """
     COBRA decomposes a (partial) gene co-expression matrix as a
     linear combination of covariate-specific components.
     It can be applied for batch correction, differential co-expression
     analysis controlling for variables, and to understand the impact of
     variables of interest to the observed co-expression.
    Parameters
    -----------
        X               : array
            design matrix of size (n, q), n = number of samples, q = number of covariates
        expressionData  : array
            gene expression as a matrix of size (g, n), g = number of genes
        standardize     : bool
            flag to standardize the gene expression as a pre-processing step
    Returns
    ---------
        psi : array
            impact of each covariate on the eigenvalues as a matrix of size (q, n)
        Q   : array
            eigenvectors corresponding to non-zero eigenvalues as a matrix of size (g, n)
        D   : array
            list of length n containing the non-zero eigenvalues
        G   : array
            (standardized) gene expression as a matrix of size (g, n)
    """
    # Extract Shapes
    p, n = expression.shape
    assert p > n, "'expression is supposed to have higher number of genes (rows) than samples (columns)."
    assert X.shape[0] == n, "'expression' is of shape p*n, so design matrix 'X' should be of shape n*q."
    _, q = X.shape

    # Standardize Gene Expressions
    g = expression - expression.mean(axis=1)[:, None] if standardize else expression.copy()
    g = g / np.linalg.norm(g, axis=1)[:, None]

    # Co-expression Matrix
    c = np.dot(g, g.T)
    c_eigenvalues, c_eigenvectors = eigh(c, subset_by_index=(p - n, p - 1))  # gives eigenvalues in ascending order
    # Select Non-zero eigenvalues
    indices_nonzero = c_eigenvalues != 0.0
    assert np.any(indices_nonzero), "Co-expression matrix is zero (e.g. all eigen values are zero)."
    Q = c_eigenvectors[:, indices_nonzero].T[::-1].T

    #
    gtq = np.matmul(g.T, Q)
    d = c_eigenvalues[indices_nonzero][::-1]

    #
    xtx_inv = np.linalg.pinv(
        np.dot(X.T, X)
    )
    xtx_inv_xt = np.dot(
        xtx_inv, X.T
    )

    #
    psi = np.zeros((q, n))

    for i in range(q):
        for h in range(n):
            psi[i, h] = n * np.sum([
                (
                        xtx_inv_xt[i, k] * gtq[k, h] ** 2
                ) for k in range(n)
            ])

    return psi, Q, d, g