import numpy as np
from math import sqrt

def otter(W, P, C, lam = 0.0035, gamma = 0.335, Iter = 300, eta = 0.00001, bexp = 1):
    b1 = 0.9
    b2 = 0.999
    eps = 0.00000001
    b1 = 0.9
    b2 = 0.999
    eps = 0.00000001
    b1t = b1**bexp
    b2t = b2**bexp

    t, g = W.shape
    P = P * (-(1-lam)/np.trace(P)) - (1-lam) * 0.0013
    C = C * (-lam /np.trace(C))
    W = P@W
    W = W/(-sqrt(np.trace(W @ W.T)))
    P = P + gamma*np.identity(t)
    m = np.zeros((t, g))
    v = np.zeros((t, g))

    for i in range(Iter):
        grad = W@W.T@W + P@W + W@C
        m = b1 * m + (4 * (1-b1)) * grad
        v = b2 * v + (16 * (1-b2)) * grad**2
        b1t = b1t * b1
        b2t = b2t * b2
        alpha = eta*sqrt(1-b2t)/(1-b1t)
        epst = eps * sqrt((1-b2t))
        W = W - alpha*(m/(epst + np.sqrt(v)))
    return W


