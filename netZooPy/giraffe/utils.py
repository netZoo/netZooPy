import numpy as np
import patsy


def check_symmetric(ppi):
    """
    Raises an exception if ppi does not have the expected structure.

    :param ppi: matrix to be checked
    """
    if ppi.shape[0] != ppi.shape[1]:
        raise Exception(
            "PPI must be a squared matrix. "
        )
    if np.diag(ppi).any() != 1:
        raise Exception(
            "PPI matrix must have ones on the diagonal. "
        )
    if np.any(np.abs(ppi - ppi.T) > 1e-8):
        raise Exception(
            "PPI matrix must be symmetrical. "
        )

def limma(pheno, exprs, covariate_formula, design_formula='1', rcond=1e-8):
    design_matrix = patsy.dmatrix(design_formula, pheno)
    
    design_matrix = design_matrix[:,1:]
    rowsum = design_matrix.sum(axis=1) -1
    design_matrix=(design_matrix.T+rowsum).T
    
    covariate_matrix = patsy.dmatrix(covariate_formula, pheno)
    design_batch = np.hstack((covariate_matrix,design_matrix))
    coefficients, res, rank, s = np.linalg.lstsq(design_batch, exprs.T, rcond=rcond)
    beta = coefficients[-design_matrix.shape[1]:]
    return exprs - design_matrix.dot(beta).T
