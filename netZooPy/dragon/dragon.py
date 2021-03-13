import numpy as np
import math
from pylab import meshgrid
from scipy.optimize import minimize
from scipy import optimize
import scipy.special as sc
import scipy.integrate as integrate
import statsmodels.stats.multitest as multi

def Scale(X):
    X_temp = X
    X_std = np.std(X_temp, axis=0)
    X_mean = np.mean(X_temp, axis=0)
    return (X_temp - X_mean) / X_std

def VarS(X):
    xbar = np.mean(X, 0)
    n = X.shape[0]
    x_minus_xbar = X - xbar
    a = x_minus_xbar*x_minus_xbar
    wbar = np.cov(X.T, bias=True)#x_minus_xbar.T@x_minus_xbar/n
    varS = a.T@a
    varS += - n*wbar**2
    varS *= n/(n-1)**3
    return(varS)

def EsqS(X):
    #xbar = np.mean(X, 0)
    n = X.shape[0]
    #x_minus_xbar = X - xbar
    wbar = np.cov(X.T, bias=True)#x_minus_xbar.T@x_minus_xbar/n
    ES2 = wbar**2*n**2/(n-1)**2
    return(ES2)

def estimate_penalty_GeneNet(X):
    vaS = VarS(X)
    vaS = np.sum(vaS) - np.trace(vaS)
    sqS = EsqS(X)
    sqS = np.sum(sqS) - np.trace(sqS)
    lam = vaS/sqS
    if lam > 1.:
        lam = 1.
    return(lam)

def estimate_penalty_parameters_dragon(X1, X2):    #X1 is n x p1 and X2 is n x p2
    n = X1.shape[0]
    p1 = X1.shape[1]
    p2 = X2.shape[1]
    X = np.append(X1, X2, axis=1)
    varS = VarS(X)
    eSqs = EsqS(X)
    IDs = np.cumsum([p1,p2])
    varS1 = varS[0:IDs[0],0:IDs[0]]
    varS12 = varS[0:IDs[0],IDs[0]:IDs[1]]
    varS2 = varS[IDs[0]:IDs[1],IDs[0]:IDs[1]]

    eSqs1 = eSqs[0:IDs[0],0:IDs[0]]
    eSqs12 = eSqs[0:IDs[0],IDs[0]:IDs[1]]
    eSqs2 = eSqs[IDs[0]:IDs[1],IDs[0]:IDs[1]]

    const = (np.sum(varS1) + np.sum(varS2) - 2.*np.sum(varS12)
            + 4.*np.sum(eSqs12))
    T1_1 = -2.*(np.sum(varS1) - np.trace(varS1) + np.sum(eSqs12))
    T1_2 = -2.*(np.sum(varS2) - np.trace(varS2) + np.sum(eSqs12))
    T2_1 = np.sum(eSqs1) - np.trace(eSqs1)
    T2_2 = np.sum(eSqs2) - np.trace(eSqs2)
    T3 = 2.*np.sum(eSqs12)
    T4 = 4.*(np.sum(varS12)-np.sum(eSqs12))

    def risk(lam):
        R = const + ((1.-lam[0]**2)*T1_1 + (1.-lam[1]**2)*T1_2
             + (1.-lam[0]**2)**2*T2_1 + (1.-lam[1]**2)**2*T2_2
             + (1.-lam[0]**2)*(1.-lam[1]**2)*T3 + lam[0]*lam[1]*T4)     #reparametrize lamx = 1-lamx**2 for better optimization
        return(R)

    x = np.arange(0., 1.01, 0.01)
    lamgrid = meshgrid(x, x)
    risk_grid = risk(lamgrid)
    indices = np.unravel_index(np.argmin(risk_grid.T, axis=None), risk_grid.shape)
    lams = [x[indices[0]],x[indices[1]]]

    res = minimize(risk, lams, method='L-BFGS-B',#'TNC',#'SLSQP',
                    tol=1e-12,
                    bounds = [[0.,1.],[0.,1.]])

    penalty_parameters = (1.-res.x[0]**2), (1.-res.x[1]**2)

    def risk_orig(lam):
        R = const + (lam[0]*T1_1 + lam[1]*T1_2
             + lam[0]**2*T2_1 + lam[1]**2*T2_2
             + lam[0]*lam[1]*T3 + np.sqrt(1-lam[0])*np.sqrt(1-lam[1])*T4)     #reparametrize lamx = 1-lamx**2 for better optimization
        return(R)
    risk_grid_orig = risk_orig(lamgrid)
    return(penalty_parameters, risk_grid_orig)

def get_shrunken_covariance_dragon(X1, X2, lambdas):
    n = X1.shape[0]
    p1 = X1.shape[1]
    p2 = X2.shape[1]
    p = p1 + p2
    X = np.append(X1, X2, axis=1)
    S = np.cov(X.T)
    Sigma = np.zeros((p,p))
    IDs = np.cumsum([0,p1,p2])
    for i in range(0,len(IDs)-1):
        T = np.diag(np.diag(S))
        Sigma[IDs[i]:IDs[i+1],IDs[i]:IDs[i+1]] = (
                    (1.-lambdas[i]) * S[IDs[i]:IDs[i+1],IDs[i]:IDs[i+1]]
                    + lambdas[i]*T[IDs[i]:IDs[i+1],IDs[i]:IDs[i+1]])
        for j in range(i+1,len(IDs)-1):
                Sigma[IDs[i]:IDs[i+1],IDs[j]:IDs[j+1]] = (
                                    np.sqrt((1.-lambdas[i])*(1.-lambdas[j]))
                                    * S[IDs[i]:IDs[i+1],IDs[j]:IDs[j+1]])
                Sigma[IDs[j]:IDs[j+1],IDs[i]:IDs[i+1]] = (
                                    Sigma[IDs[i]:IDs[i+1],IDs[j]:IDs[j+1]].T)
    return(Sigma)

def get_shrunken_covariance_GGM(X, lam):
    n = X.shape[0]
    p = X.shape[1]
    S = np.cov(X.T)
    T = np.diag(np.diag(S))
    Sigma = (1.-lam)*S + lam*T
    return(Sigma)

def get_precision_matrix_dragon(X1, X2, lambdas):
    Sigma = get_shrunken_covariance_dragon(X1, X2, lambdas)
    Theta = np.linalg.inv(Sigma)
    X = np.hstack((X1, X2))
    mu = np.mean(X, axis=0)
    return(Theta, mu)

def get_precision_matrix_GGM(X, lam):
    Sigma = get_shrunken_covariance_GGM(X, lam)
    Theta = np.linalg.inv(Sigma)
    mu = np.mean(X, axis=0)
    return(Theta, mu)

def get_partial_correlation_dragon(X1, X2, lambdas):
    Theta,_ = get_precision_matrix_dragon(X1, X2, lambdas)
    p = Theta.shape[0]
    A = np.sqrt(np.zeros((p,p))+np.diag(Theta))
    r = -Theta/A/A.T
    r = r - np.diag(np.diag(r))
    return(r)

def get_partial_correlation_from_precision(Theta):
    p = Theta.shape[0]
    A = np.sqrt(np.zeros((p,p))+np.diag(Theta))
    r = -Theta/A/A.T
    r = r - np.diag(np.diag(r))
    return(r)

def get_partial_correlation(X, lambda0):
    Theta, _ = get_precision_matrix_GGM(X, lambda0)
    p = Theta.shape[0]
    A = np.sqrt(np.zeros((p,p))+np.diag(Theta))
    r = -Theta/A/A.T
    r = r - np.diag(np.diag(r))
    return(r)

def estimate_kappa(n, p, lambda0, seed):
    Sigma = np.identity(p)
    np.random.seed(seed)
    X = np.random.multivariate_normal(np.zeros(p), Sigma, n)
    r_sim = get_partial_correlation(X, lambda0)
    r = r_sim[np.triu_indices(p,1)]
    logliterm = lambda x: 0.5*np.log(1. - x**2/(1.-lambda0)**2)
    term_Dlogli = np.sum(logliterm(r))
    Dlogli = lambda x: (0.5*len(r)*(sc.digamma(x/2)-sc.digamma((x-1)/2))
                        +term_Dlogli)
    DDlogli = lambda x: (1./4*len(r)*(sc.polygamma(1,x/2)-sc.polygamma(1,(x-1)/2)))
    #res = optimize.bisect(Dlogli, 1.001, 1000)#bracket=[1.001, 100.*n], x0=100,
                               #method='bisect')
    res = optimize.bisect(Dlogli, 1.001, 1000*n)#optimize.newton(Dlogli, n, DDlogli)
    return(res)

def estimate_p_values(r, n, lambda0, kappa='estimate', seed=1):    #here, r is the partial correlation matrix returned by get_partial_correlation
    p = r.shape[0]
    rvec = r[np.triu_indices(p,1)]
    if kappa=='analytic' and n - 1 - (p - 2) > 0:
        kappa = n - 1 - (p - 2)
        print('kappa n>>p: '+str(kappa))
    if kappa=='estimate' and lambda0<1.:
        kap_full = [estimate_kappa(n, p, lambda0, seedi) for seedi in range(10)]
        kappa = np.mean(kap_full)#estimate_kappa(n, p, lambda0, seed)
        print('kappa estimate: '+str(kappa))
    if not lambda0<1.:
        print('Warning: lambda >= 1')
        p_mat = np.ones((p,p))
        return(p_mat, p_mat)

    denom = sc.beta(1/2,(kappa-1)/2)*(1.-lambda0)
    integr = lambda x: x*sc.hyp2f1(1/2, (3-kappa)/2, 3/2, x**2/(lambda0-1)**2)

    def pval(r):
        result = 1. - 2.*integr(np.abs(r))/denom
        if result < 1e-08:
            f = lambda r0: 2.*(1. - np.abs(r0)**2/(1 - lambda0)**2)**(kappa/2 - 3/2)/denom
            result, _ = integrate.quad(lambda x: f(x), r, 1 - lambda0)
        return(result)

    def array_pval(r):
        return np.array([pval(ri) for ri in r])

    pvalues = array_pval(np.abs(rvec))
    #pvalues = 1. - 2.*integr(np.abs(rvec))/denom

    _, pvalues_adj ,_ ,_ = multi.multipletests(pvalues, alpha=0.05, method='fdr_bh')
    adj_pvalues_mat = np.zeros((p,p))
    adj_pvalues_mat[np.triu_indices(p,1)] = pvalues_adj
    adj_pvalues_mat += adj_pvalues_mat.T
    p_mat = np.zeros((p,p))
    p_mat[np.triu_indices(p,1)] = pvalues
    p_mat += p_mat.T
    return(adj_pvalues_mat, p_mat)

def estimate_kappa_dragon(n, p1, p2, lambdas, seed, simultaneous = False):
    Sigma1 = np.identity(p1)
    Sigma2 = np.identity(p2)
    np.random.seed(seed)
    X1 = np.random.multivariate_normal(np.zeros(p1), Sigma1, n)
    X2 = np.random.multivariate_normal(np.zeros(p2), Sigma2, n)
    r_sim = get_partial_correlation_dragon(X1, X2, lambdas)
    IDs = np.cumsum([0,p1,p2])
    r11 = r_sim[IDs[0]:IDs[1],IDs[0]:IDs[1]]
    r11 = r11[np.triu_indices(p1,1)]
    r22 = r_sim[IDs[1]:IDs[2],IDs[1]:IDs[2]]
    r22 = r22[np.triu_indices(p2,1)]
    r12 = r_sim[IDs[0]:IDs[1],IDs[1]:IDs[2]].flatten()

    logliterm = lambda x, lam: 0.5*np.log(1. - x**2/lam**2)
    if lambdas[0] == 1. or lambdas[1] == 1.:
        print("Warning: simultaneous fit not possible! Simultaneous set to False.")
        simultaneous = False
    if simultaneous is False:
        if lambdas[0] < 1.:
            term_Dlogli11 = np.sum(logliterm(r11, lam=1.-lambdas[0]))
            Dlogli11 = lambda x: (1./4*p1*(p1-1)
                            *(sc.digamma(x/2)-sc.digamma((x-1)/2))
                            +term_Dlogli11)
            kappa11 = optimize.bisect(Dlogli11, 1.001, 1000*n)
        else:
            kappa11 = np.nan

        if lambdas[1] < 1.:
            term_Dlogli22 = np.sum(logliterm(r22, lam=1.-lambdas[1]))
            Dlogli22 = lambda x: (1./4*p2*(p2-1)
                            *(sc.digamma(x/2)-sc.digamma((x-1)/2))
                            +term_Dlogli22)
            kappa22 = optimize.bisect(Dlogli22, 1.001, 1000*n)
        else:
            kappa22 = np.nan

        if lambdas[0] < 1. and lambdas[1] < 1.:
            term_Dlogli12 = np.sum(logliterm(r12, lam=np.sqrt(1.-lambdas[0])*np.sqrt(1.-lambdas[1])))
            Dlogli12 = lambda x: (1./2*p1*p2
                            *(sc.digamma(x/2)-sc.digamma((x-1)/2))
                            +term_Dlogli12)
            kappa12 = optimize.bisect(Dlogli12, 1.001, 1000*n)
        else:
            kappa12 = np.nan

    if simultaneous and lambdas[0] < 1. and lambdas[1] < 1.:
        logliterm = lambda x, lam: 0.5*np.log(1. - x**2/lam**2)
        term_Dlogli = (np.sum(logliterm(r11, lam=1.-lambdas[0]))
                   + np.sum(logliterm(r22, lam=1.-lambdas[1]))
                   + np.sum(logliterm(r12, lam=np.sqrt(1.-lambdas[0])
                     *np.sqrt(1.-lambdas[1])))
                   )

        Dlogli = lambda x: (1./4*(p1+p2)*(p1+p2-1)
                        *(sc.digamma(x/2)-sc.digamma((x-1)/2))
                        +term_Dlogli)
        DDlogli = lambda x: (1./8*(p1+p2)*(p1+p2-1)
                         *(sc.polygamma(1,x/2)-sc.polygamma(1,(x-1)/2)))
        res = optimize.bisect(Dlogli, 1.001, 1000*n)
        kappa11 = res
        kappa22 = res
        kappa12 = res
        return(kappa11, kappa22, kappa12)


    return kappa11, kappa22, kappa12


def estimate_p_values_dragon(r, n, p1, p2, lambdas, kappa='estimate', seed=1, simultaneous = False):    #here, r is the partial correlation matrix returned by get_partial_correlation
    lam = lambdas
    IDs = np.cumsum([0,p1,p2])
    p = p1+p2
    r11 = r[IDs[0]:IDs[1],IDs[0]:IDs[1]]
    r11 = r11[np.triu_indices(p1,1)]
    r22 = r[IDs[1]:IDs[2],IDs[1]:IDs[2]]
    r22 = r22[np.triu_indices(p2,1)]
    r12 = r[IDs[0]:IDs[1],IDs[1]:IDs[2]].flatten()

    if kappa=='estimate':
        if lam[0]==1. and lam[1]==1.:
            p_mat = np.ones((p,p))
            return(p_mat, p_mat)
        else:
            kap_full = [estimate_kappa_dragon(n, p1, p2, lam, seedi,
                                simultaneous=simultaneous) for seedi in range(10)]
            #print(kap_full)
            kappa = np.mean(kap_full, axis=0)
            kappa = kappa.tolist()#estimate_kappa_dragon(n, p1, p2, lam, seed,
                                #simultaneous=simultaneous)
        print('kappa estimate (k11, k22, k12): '+str(kappa))
    if kappa=='analytic':
        kappa = n - 1 - (p - 2)
        kappa = kappa, kappa, kappa
        print('kappa n>>p (k11, k22, k12): '+str(kappa))
    def denom(l,k):
        if k!=None and k==k:
            denom = sc.beta(1/2,(k-1)/2)*l
        else:
            denom = None
        return(denom)
    d11 = denom(1.-lam[0], kappa[0])
    d22 = denom(1.-lam[1], kappa[1])
    d12 = denom(np.sqrt(1.-lam[0])*np.sqrt(1.-lam[1]), kappa[2])

    pfunc = lambda x, l, k: (x*sc.hyp2f1(1/2, (3-k)/2, 3/2, (x/l)**2))

    n1 = int(p1*(p1-1)/2)
    n2 = int(p2*(p2-1)/2)
    n12 = int(p1*p2)

    def pval(r, l, k, d):
        result = 1. - 2.*pfunc(np.abs(r), l, k)/d
        if np.abs(result) < 1e-08:
            f = lambda r0: 2.*(1. - np.abs(r0)**2/l**2)**(k/2 - 3/2)/d
            result, _ = integrate.quad(lambda x: f(x), r, l)
        return(result)

    def array_pval(r, l, k, d):
        return np.array([pval(ri, l, k, d) for ri in r])

    if lam[0] < 1.:
        # p11 = 1. -2.*pfunc(np.abs(r11), 1.-lam[0], kappa[0])/d11
        p11 = array_pval(np.abs(r11), 1.-lam[0], kappa[0], d11)
    else:
        p11 = np.ones(n1)

    if lam[1] < 1.:
        # p22 = 1. -2.*pfunc(np.abs(r22), 1.-lam[1], kappa[1])/d22
        p22 = array_pval(np.abs(r22), 1.-lam[1], kappa[1], d22)
    else:
        p22 = np.ones(n2)

    if lam[0] < 1. and lam[1] < 1.:
        # p12 = 1. -2.*pfunc(np.abs(r12), np.sqrt(1.-lam[0])*np.sqrt(1.-lam[1]), kappa[2])/d12
        p12 = array_pval(np.abs(r12), np.sqrt(1.-lam[0])*np.sqrt(1.-lam[1]),
                         kappa[2], d12)
    else:
        p12 = np.ones(n12)

    IDs = np.cumsum([0, p1, p2])
    p_mat = np.zeros((p, p))

    p_mat[IDs[0]:IDs[1], IDs[0]:IDs[1]][np.triu_indices(p1, 1)] = p11
    p_mat[IDs[1]:IDs[2], IDs[1]:IDs[2]][np.triu_indices(p2, 1)] = p22
    p_mat[IDs[0]:IDs[1], IDs[1]:IDs[2]] = p12.reshape((p1, p2))
    p_mat += p_mat.T

    _, padj11, _, _ = multi.multipletests(p11, alpha=0.05, method='fdr_bh')
    _, padj22, _, _ = multi.multipletests(p22, alpha=0.05, method='fdr_bh')
    _, padj12, _, _ = multi.multipletests(p12, alpha=0.05, method='fdr_bh')
    adj_pvalues_mat = np.zeros((p, p))
    adj_pvalues_mat[IDs[0]:IDs[1], IDs[0]:IDs[1]][np.triu_indices(p1,1)] = padj11
    adj_pvalues_mat[IDs[1]:IDs[2], IDs[1]:IDs[2]][np.triu_indices(p2,1)] = padj22
    adj_pvalues_mat[IDs[0]:IDs[1], IDs[1]:IDs[2]] = padj12.reshape((p1,p2))

    adj_pvalues_mat += adj_pvalues_mat.T

    return(adj_pvalues_mat, p_mat)

def calculate_true_R(X1, X2, Sigma):
    x = np.arange(0., 1.01, 0.01)
    n_lams = len(x)
    lamgrid = meshgrid(x, x)
    R = np.zeros((n_lams,n_lams))
    for i in range(n_lams):
        print(i)
        for j in range(n_lams):
            R[i,j] = np.sum((Sigma
                                - get_shrunken_covariance_dragon(X1, X2,
                           [lamgrid[0][i,j],lamgrid[1][i,j]]))**2)
    return(R)


def simulate_dragon_data(eta11, eta12, eta22, p1, p2, epsilon, n, seed):
    np.random.seed(seed)
    Theta = np.identity(p1+p2)
    n11 = int(np.around(p1*(p1-1)/2*eta11))
    n12 = int(np.around(p1*p2*eta12))
    n22 = int(np.around(p2*(p2-1)/2*eta22))
    print("n11="+str(n11)+", n12="+str(n12)+", n22="+str(n22))
    IDs = np.cumsum([0,p1,p2])
    n11_IDs = np.random.choice(range(int(p1*(p1-1)/2)), size=n11, replace=False)
    n12_IDs = np.random.choice(range(int(p1*p2)), size=n12, replace=False)
    n22_IDs = np.random.choice(range(int(p2*(p2-1)/2)), size=n22, replace=False)
    Theta11 = Theta[IDs[0]:IDs[1],IDs[0]:IDs[1]]
    Theta12 = Theta[IDs[0]:IDs[1],IDs[1]:IDs[2]]
    Theta22 = Theta[IDs[1]:IDs[2],IDs[1]:IDs[2]]
    Theta11_vec = Theta11[np.triu_indices(p1,1)]
    Theta12_vec = Theta12.flatten()
    Theta22_vec = Theta22[np.triu_indices(p2,1)]

    Theta11_vec[n11_IDs] = np.random.uniform(-1.,1.,size=len(n11_IDs))
    Theta12_vec[n12_IDs] = np.random.uniform(-1.,1.,size=len(n12_IDs))
    Theta22_vec[n22_IDs] = np.random.uniform(-1.,1.,size=len(n22_IDs))

    Theta11[np.triu_indices(p1,1)] = Theta11_vec
    Theta22[np.triu_indices(p2,1)] = Theta22_vec
    Theta[IDs[0]:IDs[1],IDs[1]:IDs[2]] = Theta12_vec.reshape((p1,p2))
    Theta[IDs[0]:IDs[1],IDs[0]:IDs[1]] = Theta11
    Theta[IDs[1]:IDs[2],IDs[1]:IDs[2]] = Theta22
    Theta = Theta + Theta.T - np.identity(p1+p2)
    Theta = Theta - np.identity(p1+p2) + np.diag(np.sum(abs(Theta), axis=0)+ 0.0001)
    A = np.zeros((p1+p2,p1+p2)) +  np.sqrt(np.diag(Theta))
    Theta = Theta/A/A.T
    Sigma = np.linalg.inv(Theta)
    mu = np.zeros(p1+p2)
    X = np.random.multivariate_normal(mean=mu, cov=Sigma, size=n)

    noise1 = np.random.normal(0, epsilon[0], (n,p1))
    noise2 = np.random.normal(0, epsilon[1], (n,p2))

    X1 = X[:,IDs[0]:IDs[1]] + noise1
    X2 = X[:,IDs[1]:IDs[2]] + noise2
    return(X1, X2, Theta, Sigma)

def logli(X, Theta, mu):
    X_mu = X - mu
    n = X.shape[0]
    p = X.shape[1]
    S = X_mu.T@X_mu
    term1 = n/2.*np.linalg.slogdet(Theta)[1]     #det(A**(-1)) = 1/det(A)
    term2 = - 0.5*np.sum(S*Theta)     #np.trace(S@Theta)
    term3 = - 0.5*p*n*np.log(2.*math.pi)
    return((term1 + term2 + term3)/n)
