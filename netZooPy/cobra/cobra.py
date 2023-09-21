import numpy as np

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
    n = expression.shape[1]
    N = min(n, expression.shape[0])

    if standardize:
        expression = (expression.T - np.mean(expression, axis=1)).T
        g = (expression.T / np.sqrt(np.sum(expression ** 2, axis=1))).T
        Sigma = g.dot(g.T)

    Q, d, _ = np.linalg.svd(Sigma, full_matrices=True)
    Q = Q[:, [i for i in range(N)]]
    D = d[[i for i in range(N)]]

    hatmat = np.linalg.pinv(X.T.dot(X)).dot(X.T)

    psi = np.zeros((X.shape[1], N))

    for i in range(psi.shape[0]):
        psi[i, :] = (np.diag(Q.T.dot(g).dot(n * np.diag(hatmat[i,])).dot((Q.T.dot(g)).T)))

    return psi, Q, D, g