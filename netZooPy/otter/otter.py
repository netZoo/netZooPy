import numpy as np
from math import sqrt


def otter(W, P, C, lam=0.035, gamma=0.335, Iter=60, eta=0.00001, bexp=1):
    """
	Description:
	    OTTER infers gene regulatory networks using TF DNA binding
	    motif (W), TF PPI (P), and gene coexpression (C) through
	    minimzing the following objective:
	        min f(W)
	    with f(W) = (1-lam)/4*||WW' - P||^2 + (lam/4)*||W'W - C||^2 + (gamma/2)*||W||^2
	Inputs:
	    W     : TF-gene regulatory network based on TF motifs as a
	            matrix of size (t,g), g=number of genes, t=number of TFs
	    P     : TF-TF protein interaction network as a matrix of size (t,t)
	    C     : gene coexpression as a matrix of size (g,g)
	    lam   : it should be in [0,1].
	    gamma : penalization term
	    Iter  : number of iterations of the algorithm
	    eta   : the learning rate
	Outputs:
	    W     : Predicted TF-gene complete regulatory network as an adjacency matrix of size (t,g).
    Reference:
        Weighill, Deborah, et al. "Gene Regulatory Network Inference as Relaxed Graph Matching." BioRxiv (2020).
    """
    b1 = 0.9
    b2 = 0.999
    eps = 0.00000001
    b1t = b1 ** bexp
    b2t = b2 ** bexp

    t, g = W.shape

    C = C  / np.trace(C)
    W = W / np.sqrt(np.trace(W @ W.T))
    P = P / np.trace(P)

    P = P * (-(1 - lam)) + gamma * np.identity(t)
    C = C * (-lam)

    m = np.zeros((t, g))
    v = np.zeros((t, g))

    for i in range(Iter):
        grad = W @ W.T @ W + P @ W + W @ C
        m = b1 * m + (4 * (1 - b1)) * grad
        v = b2 * v + (16 * (1 - b2)) * grad ** 2
        b1t = b1t * b1
        b2t = b2t * b2
        alpha = eta * sqrt(1 - b2t) / (1 - b1t)
        epst = eps * sqrt((1 - b2t))
        W = W - alpha * (m / (epst + np.sqrt(v)))
    return W