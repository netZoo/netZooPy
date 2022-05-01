import numpy as np
def onlineCoexpression(si,n,mi,std,cov):
        """ 
        Description:
            onlineCoexpression computes the correlation matrix of n
            samples deprived of sample si, using the correlation matrix
            of n samples. ~4-8x faster when number of genes ~ number of
            observations and 12x-35x when number of samples is much
            larger than the variables.
            Particularly intersting in iterative computation of several
            coexpression matrix with large samlples or when computing large 
            matrices that require consequent GPU/CPU memory.

        Inputs:
            si : k samples to remove as a k by genes vector
            n  : number of all the samples
            mi : mean of all the samples
            std: std of all the samples
            cov: covariance matrix of all the samples 

        Outputs:
            onCoex : co-expression matrix deprived of samples si

         Authors: 
            Marouen Ben Guebila, Daniel Morgan
        """
        # First we compute the new mean online
        newm   = (1/(n-1)) * ( (mi*n)- si)
        # Then we compute the new std online using the orthogonality trick
        newstd = np.sqrt((np.square(std) - (1/n) * np.square(si - newm)) * ((n-1)/(n-2)))
        # Then we compute the new covariance online
        onCov= (1/(n-2)) * ( (cov*(n-1)) - ( (n/(n-1)) * (np.matmul((si-mi).T,(si-mi))) ) )
        # Finally, we derive the new coexpression online
        onCoex= onCov / np.matmul(newstd.T,newstd)
        # We set the diagonal explicitly to avoid numerical stability
        np.fill_diagonal(onCoex, 1)

        return onCoex
