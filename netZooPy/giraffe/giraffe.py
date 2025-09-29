import random
#from sklearn.linear_model import LinearRegression
import torch
from torch import nn
from torch.functional import F
from typing import Callable
from netZooPy.cobra import *
from netZooPy.giraffe.utils import *
import logging
import pandas as pd

from torch.optim import Adam
from torch.optim.optimizer import required


class Giraffe(object):
    """
            GIRAFFE infers regulation and transcription factor activity using

            gene expression
            TF DNA binding motif
            TF interactions (PPI)

            while (optionally) controlling for covariates of interest (e.g. age).
            There are G genes, TF transcription factors, and n samples (cells, patients, ...)

            Parameters
            -----------
            expression : numpy array or pandas dataframe
                    Dimension G x n
            prior : numpy array or pandas dataframe
                    TF-gene regulatory network based on TF motifs as a
                    matrix of size (tf,g), g=number of genes, tf=number of TFs
            ppi : numpy array or pandas dataframe
                    TF-TF protein interaction network as a matrix of size (tf, tf)
                    Must be symmetrical and have ones on the diagonal
            design_matrix : np.ndarray, pd.DataFrame
                    COBRA design matrix of size (n, q), n = number of samples, q = number of covariates
            cobra_covariate_to_keep : int
                    Zero-indedex base of COBRA co-expression component to use
            iterations : int
                    number of iterations of the algorithm
            lr : float
                    the learning rate
                lam:

            Public Functions
            ---------
            get_regulation : numpy array
                    Returns the predicted TF-gene complete regulatory network as an adjacency matrix of size (tf,g)
            get_tfa : numpy array
                    Returns the predicted transcription factor activity  as an adjacency matrix of size (tf,n)

            References
            -------------
                Author: Soel Micheletti
        """

    def __init__(
            self,
            expression,
            prior,
            ppi,
            adjusting=None,
            design_matrix=None,
            cobra_covariate_to_keep=0,
            l1=0,
            max_iter=2000,
            min_iter=200,
            lr=1e-5,
            lam=None,
            balance_fn = None,
            save_computation = False,
            seed=42,  # For reproducibility
    ) -> None:

        torch.manual_seed(seed)
        random.seed(seed)
        np.random.seed(seed)
        self._save_computation = save_computation
        self.process_data(
            expression,
            prior,
            ppi,
            adjusting,
            design_matrix,
            cobra_covariate_to_keep
        )
        self._max_iter = max_iter
        self._min_iter = min_iter
        self._lr = lr
        self._l = l1
        self._lam = lam
        self._balance_fn = balance_fn
        if self._lam is not None and self._balance_fn is not None:
            raise ValueError("Can't set both custom weights and custom weight learning function! Please provide only one of both options. ")
        self._R, self._TFA = self._compute_giraffe()

    def process_data(
            self,
            expression,
            motif,
            ppi,
            adjusting,
            design_matrix,
            cobra_covariate_to_keep
    ):
        """
        Processes data files into normzalized data matrices.

        :param expression: numpy matrix or pandas dataframe containing gene expression data.
        :param motif: numpy matrix or pandas dataframe containing motif data.
        :param ppi: numpy or pandas dataframe containing PPI data. Must be symmetrical.
        :param adjusting: vector of covariates that needs to be adjusted.
        """
        self._expression = expression
        self._motif = motif
        self._ppi = ppi

        if isinstance(expression, pd.DataFrame):
            self._expression = expression.to_numpy()
        if isinstance(motif, pd.DataFrame):
            self._motif = motif.to_numpy()
        if isinstance(ppi, pd.DataFrame):
            self._ppi = ppi.to_numpy()

        check_symmetric(self._ppi) # Check that ppi has legal structure.

        if not isinstance(self._expression, np.ndarray):
            raise ValueError(
                "Error in processing expression data. Please provide a numpy array or a pandas dataframe. "
            )
        if not isinstance(self._motif, np.ndarray):
            raise ValueError(
                "Error in processing motif data. Please provide a numpy array or a pandas dataframe. "
            )
        if not isinstance(self._ppi, np.ndarray):
            raise ValueError(
                "Error in processing PPI data. Please provide a numpy array or a pandas dataframe. "
            )

        if adjusting is None:
            self._adjusting = None
        else:
            if len(adjusting.shape) == 1:
                adjusting = adjusting.reshape(len(adjusting), 1)
            self._adjusting = torch.Tensor(np.zeros(adjusting.shape))
            if isinstance(adjusting, pd.DataFrame):
                adjusting = adjusting.to_numpy()
            if adjusting is not None and not isinstance(adjusting, np.ndarray):
                raise ValueError(
                    "Error in processing adjusting covariates. Please provide a numpy array or a pandas dataframe. "
                )
            if adjusting is not None:
                for i in range(adjusting.shape[1]):
                    tmp = torch.Tensor(adjusting[:, i])
                    tmp /= torch.norm(tmp)
                    self._adjusting[:, i] = tmp

        # Normalize motif and PPI
        self._motif = motif / np.sqrt(np.trace(motif.T @ motif))
        self._ppi = self._ppi / np.trace(self._ppi)

        self._C = 0
        if not self._save_computation:
            # Compute co-expression matrix
            if (np.cov(self._expression).diagonal() == 0).any():
                self._expression += 1e-8

            if design_matrix is not None:
                psi, Q, d, g = cobra(design_matrix, self._expression)
                if cobra_covariate_to_keep < 0 or cobra_covariate_to_keep >= psi.shape[0]:
                    raise AttributeError(
                        "Invalid COBRA component! Valid COBRA components are in range " + str(0) + " - " + str(psi.shape[0] - 1)
                    )
                self._C = Q.dot(np.diag(psi[cobra_covariate_to_keep,])).dot(Q.T)
            else:
                self._C = np.corrcoef(self._expression)
            self._C /= np.trace(self._C)
        
        if design_matrix is not None:
            self._expression = limma(exprs=self._expression, pheno=pd.DataFrame(design_matrix), covariate_formula = "~ -1")

    def _compute_giraffe(self):
        """
        Giraffe optimization.
        :return: R : matrix g x tf, partial effects between tanscription factor and gene.
                 TFA : matrix tf x n, transcrption factor activity.
        """

        variables_to_adjust = 0
        if self._adjusting is not None:
            variables_to_adjust = self._adjusting.shape[1]

        giraffe_model = Model(np.random.random((self._motif.shape[1], self._expression.shape[1])), self._motif, variables_to_adjust)
        optim = torch.optim.Adam(giraffe_model.parameters(), lr=self._lr)  # We run Adam to optimize f(R, TFA)
        if self._l > 0:
            optim = ProxAdam(giraffe_model.parameters(), lr=self._lr, lambda_ = self._l)

        last_loss = float('inf')
        converged = False
        for i in range(self._max_iter):
            pred = giraffe_model(
                torch.Tensor(self._expression),
                torch.Tensor(self._ppi),
                torch.Tensor(self._C),
                self._lam,
                self._balance_fn,
                self._adjusting,
                self._save_computation
            )  # Compute f(R, TFA)
            loss = F.mse_loss(pred, torch.norm(
                torch.Tensor(np.zeros((3, 3))))).sqrt()  # Minimization problem: loss = ||f(R, TFA)||
            optim.zero_grad()  # Reset gradients
            loss.backward()
            optim.step()  # Adam step
            if i >= self._min_iter and abs(last_loss - loss.detach().numpy()) / last_loss < 1e-3:
                converged=True
                break
            last_loss = loss.detach().numpy()

        if not converged:
            logging.warn('GIRAFFE did not converge and was stopped after ' + str(self._max_iter) + ' iterations. Try increasing the max_iter parameter.')
        R = giraffe_model.R.detach().numpy()
        TFA = torch.abs(giraffe_model.TFA.detach()).numpy()
        return R, TFA

    def get_regulation(self):
        return self._R

    def get_tfa(self):
        return self._TFA
        # TFA_giraffe = np.zeros(self._TFA.shape)
        # for i in range(self._TFA.shape[1]):
        #    TFA_giraffe[:, i] = LinearRegression(fit_intercept=False, positive=True).fit(self._R, self._expression[:, i]).coef_
        # return TFA_giraffe


class Model(nn.Module):
    def __init__(
            self,
            tfa_prior,
            motif,
            variables_to_adjust
    ):
        super().__init__()
        self.TFA = nn.Parameter(torch.Tensor(tfa_prior))
        self.R = nn.Parameter(torch.Tensor(motif))
        self.variables_to_adjust = variables_to_adjust
        self.coefs = nn.Parameter(torch.Tensor(torch.ones((motif.shape[0], self.variables_to_adjust))))

    def forward(
            self,
            Y,
            PPI,
            C,
            lam,
            balance_fn,
            adjusting,
            save_computation
    ):
        if adjusting is None:
            L1 = torch.norm(Y - torch.matmul(self.R, torch.abs(self.TFA))) ** 2
            L2 = torch.norm(torch.matmul(torch.t(self.R), self.R) - PPI) ** 2
            L3 = torch.norm(torch.Tensor(torch.zeros((2, 2))))
            if not save_computation:
                L3 = torch.norm(torch.matmul(self.R, torch.t(self.R)) - C) ** 2
            L4 = torch.norm(torch.matmul(torch.abs(self.TFA), torch.t(torch.abs(self.TFA))) - PPI) ** 2
            L5 = torch.norm(self.R) ** 2
            weights = self._get_weights(lam, balance_fn, L1, L2, L3, L4, L5)
            return weights[0] * L1 + weights[1] * L2 + weights[2] * L3 + weights[3] * L4 + weights[4] * L5
        else :
            L1 = torch.norm(Y - torch.matmul(torch.hstack([self.R, self.coefs]), torch.vstack([torch.abs(self.TFA), torch.t(adjusting)]))) ** 2
            L2 = torch.norm(torch.matmul(torch.t(self.R), self.R) - PPI) ** 2
            L3 = torch.norm(torch.Tensor(torch.zeros((2, 2))))
            if not save_computation:
                L3 = torch.norm(torch.matmul(self.R, torch.t(self.R)) - C) ** 2
            L4 = torch.norm(torch.matmul(torch.vstack([torch.abs(self.TFA), torch.t(adjusting)]), torch.t(torch.vstack([torch.abs(self.TFA), torch.t(adjusting)]))) - torch.hstack([torch.vstack([PPI, torch.Tensor(np.ones((self.variables_to_adjust, self.R.shape[1])))]), torch.Tensor(torch.ones((self.R.shape[1] + self.variables_to_adjust, self.variables_to_adjust)))])) ** 2
            L5 = torch.norm(torch.hstack([self.R, self.coefs])) ** 2
            weights = self._get_weights(lam, balance_fn, L1, L2, L3, 0, torch.Tensor([0]))
            return weights[0] * L1 + 3 * weights[1] * L2 + weights[2] * L3 + weights[3] * L4 + weights[4] * L5

    def _get_weights(self, lam, balance_fn, L1, L2, L3, L4, L5):
        weights = [1, 1, 1, 1, 1]
        if lam is not None:
            weights = lam
        elif balance_fn is not None:
            weights = balance_fn(L1, L2, L3, L5)
        else:
            sum = L1.item() + L2.item() + L3.item() + L5.item()
            weights = [1 - L1.item() / sum, 1 - L2.item() / sum, 1 - L3.item() / sum, 1, 1]
        return weights

class ProxAdam(Adam ):
    def __init__(self, params, lr=required, lambda_=0):

        kwargs = dict(lr=lr)
        super().__init__(params, **kwargs)
        self._lambda = lambda_
        proxs = [self.soft_thresholding]
        if len(proxs) != len(self.param_groups):
            raise ValueError("Invalid length of argument proxs: {} instead of {}".format(len(proxs), len(self.param_groups)))

        for group, prox in zip(self.param_groups, list(proxs)):
            group.setdefault('prox', prox)

    def step(self):
        # perform a gradient step
        super().step()

        for group in self.param_groups:
            prox = group['prox']

            # apply the proximal operator to each parameter in a group
            for p in group['params']:
                p.data = prox(p.data)

    def soft_thresholding(self, x):
        return torch.sign(x) * torch.nn.functional.relu(torch.abs(x) - self._lambda)
