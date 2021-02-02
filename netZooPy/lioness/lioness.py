from __future__ import print_function
import os, os.path,sys
import numpy as np
import pandas as pd
from .timer import Timer
sys.path.insert(1,'../panda')
from netZooPy.panda.panda import Panda
from joblib.externals.loky import set_loky_pickler
from joblib import parallel_backend
from joblib import Parallel, delayed
from joblib import wrap_non_picklable_objects

class Lioness(Panda):
    """
    Description:
        Using LIONESS to infer single-sample gene regulatory networks.
        1. Reading in PANDA network and preprocessed middle data
        2. Computing coexpression network
        3. Normalizing coexpression network
        4. Running PANDA algorithm
        5. Writing out LIONESS networks

    Inputs:
        Panda: Panda object.

    Methods:
        __init__              : Intialize instance of Lioness class and load data.
        __lioness_loop        : The LIONESS algorithm.
        save_lioness_results  : Saves LIONESS network.

    Example:
        from netZooPy.lioness.lioness import Lioness
        To run the Lioness algorithm for single sample networks, first run PANDA using the keep_expression_matrix flag, then use Lioness as follows:
        panda_obj = Panda('../../tests/ToyData/ToyExpressionData.txt', '../../tests/ToyData/ToyMotifData.txt', '../../tests/ToyData/ToyPPIData.txt', remove_missing=False, keep_expression_matrix=True)
        lioness_obj = Lioness(panda_obj)

        Save Lioness results:
        lioness_obj.save_lioness_results('Toy_Lioness.txt')
        Return a network plot for one of the Lioness single sample networks:
        plot = AnalyzeLioness(lioness_obj)
        plot.top_network_plot(column= 0, top=100, file='top_100_genes.png')

        Example lioness output:
        Sample1 Sample2 Sample3 Sample4
        -------------------------------
        -0.667452814003	-1.70433776179	-0.158129613892	-0.655795512803
        -0.843366539284	-0.733709815256	-0.84849895139	-0.915217389738
        3.23445386464	2.68888472802	3.35809757371	3.05297381396
        2.39500370135	1.84608635425	2.80179804094	2.67540878165
        -0.117475863987	0.494923925853	0.0518448588965	-0.0584810456421

        TF, Gene and Motif order is identical to the panda output file.

    Authors: 
        cychen, davidvi, dcolinmorgan

    Reference:
        Kuijjer, Marieke Lydia, et al. "Estimating sample-specific regulatory networks." Iscience 14 (2019): 226-240.
    """
    def __init__(self, obj, computing='cpu', precision='double',ncores=1,start=1, end=None, save_dir='lioness_output', save_fmt='npy'):
        """
        Description:
            Initialize instance of Lioness class and load data.

        Inputs:
            obj             : PANDA object, generated with keep_expression_matrix=True.
            obj.motif_matrix: TF DNA motif binding data in tf-by-gene format.
                              If set to None, Lioness will be performed on gene coexpression network.
            computing       : 'cpu' uses Central Processing Unit (CPU) to run PANDA
                              'gpu' use the Graphical Processing Unit (GPU) to run PANDA
            precision       : 'double' computes the regulatory network in double precision (15 decimal digits).
                              'single' computes the regulatory network in single precision (7 decimal digits) which is fastaer, requires half the memory but less accurate.
            start           : Index of first sample to compute the network.
            end             : Index of last sample to compute the network.
            save_dir        : Directory to save the networks.
            save_fmt        : Save format.
                              '.npy': (Default) Numpy file.
                              '.txt': Text file.
                              '.mat': MATLAB file.
        """
        # Load data
        with Timer("Loading input data ..."):
            self.export_panda_results = obj.export_panda_results
            self.expression_matrix = obj.expression_matrix
            self.motif_matrix = obj.motif_matrix
            self.ppi_matrix = obj.ppi_matrix
            self.correlation_matrix=obj.correlation_matrix
            if precision=='single':
                self.correlation_matrix=np.float32(self.correlation_matrix)
                self.motif_matrix=np.float32(self.motif_matrix)
                self.ppi_matrix=np.float32(self.ppi_matrix)
            self.computing=computing
            self.n_cores=int(ncores)
            if hasattr(obj,'panda_network'):
                self.network = obj.panda_network
            elif hasattr(obj,'puma_network'):
                self.network = obj.puma_network
            else:
                print('Cannot find panda or puma network in object')
                raise AttributeError('Cannot find panda or puma network in object')
            del obj

        # Get sample range to iterate
        self.n_conditions = self.expression_matrix.shape[1]
        self.indexes = range(self.n_conditions)[start-1:end]  # sample indexes to include

        # Create the output folder if not exists
        self.save_dir = save_dir
        self.save_fmt = save_fmt
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        # Run LIONESS
        if int(self.n_conditions) >= int(self.n_cores) and self.computing=='cpu':
            self.total_lioness_network=Parallel(n_jobs=self.n_cores)(self.__par_lioness_loop(i) for i in (self.indexes))
            
        elif self.computing=='gpu':
            for i in self.indexes:
                self.total_lioness_network = self.__lioness_loop(i)
        #        # self.export_lioness_results = pd.DataFrame(self.total_lioness_network)

        # create result data frame
        self.export_lioness_results = pd.DataFrame(self.total_lioness_network)
        self.save_lioness_results()
    
    # if int(self.n_conditions) >= int(self.n_cores) and self.computing=='cpu':  #### tried to make this function robust but seems it cannot be done
    #     @delayed
    #     @wrap_non_picklable_objects
        
    def __lioness_loop(self,i):
        """
        Description:
            Initialize instance of Lioness class and load data.

        Outputs:
            self.total_lioness_network: An edge-by-sample matrix containing sample-specific networks.
        """
        # for i in self.indexes:
        print("Running LIONESS for sample %d:" % (i+1))
        idx = [x for x in range(self.n_conditions) if x != i]  # all samples except i
        with Timer("Computing coexpression network:"):
            if self.computing=='gpu':
                import cupy as cp
                correlation_matrix = cp.corrcoef(self.expression_matrix[:, idx])
                if cp.isnan(correlation_matrix).any():
                    cp.fill_diagonal(correlation_matrix, 1)
                    correlation_matrix = cp.nan_to_num(correlation_matrix)
                correlation_matrix=cp.asnumpy(correlation_matrix)
            else:
                correlation_matrix = np.corrcoef(self.expression_matrix[:, idx])
                if np.isnan(correlation_matrix).any():
                    np.fill_diagonal(correlation_matrix, 1)
                    correlation_matrix = np.nan_to_num(correlation_matrix)

        with Timer("Normalizing networks:"):
            correlation_matrix_orig = correlation_matrix # save matrix before normalization
            correlation_matrix = self._normalize_network(correlation_matrix)

        with Timer("Inferring LIONESS network:"):
            if self.motif_matrix is not None:
                del correlation_matrix_orig
                subset_panda_network = self.panda_loop(correlation_matrix, np.copy(self.motif_matrix), np.copy(self.ppi_matrix),self.computing)
            else:
                del correlation_matrix
                subset_panda_network = correlation_matrix_orig

        lioness_network = self.n_conditions * (self.network - subset_panda_network) + subset_panda_network

        with Timer("Saving LIONESS network %d to %s using %s format:" % (i+1, self.save_dir, self.save_fmt)):
            path = os.path.join(self.save_dir, "lioness.%d.%s" % (i+1, self.save_fmt))
            if self.save_fmt == 'txt':
                np.savetxt(path, lioness_network)
            elif self.save_fmt == 'npy':
                np.save(path, lioness_network)
            elif self.save_fmt == 'mat':
                from scipy.io import savemat
                savemat(path, {'PredNet': lioness_network})
            else:
                print("Unknown format %s! Use npy format instead." % self.save_fmt)
                np.save(path, lioness_network)
        if self.computing=='gpu' and i == 0:
            self.total_lioness_network = np.fromstring(np.transpose(lioness_network).tostring(),dtype=lioness_network.dtype)
        elif self.computing=='gpu' and i!=0:
            self.total_lioness_network=np.column_stack((self.total_lioness_network ,np.fromstring(np.transpose(lioness_network).tostring(),dtype=lioness_network.dtype)))

        return self.total_lioness_network
    
    @delayed
    @wrap_non_picklable_objects
    def __par_lioness_loop(self,i):
        """
        Description:
            Initialize instance of Lioness class and load data.

        Outputs:
            self.total_lioness_network: An edge-by-sample matrix containing sample-specific networks.
        """
        # for i in self.indexes:
        print("Running LIONESS for sample %d:" % (i+1))
        idx = [x for x in range(self.n_conditions) if x != i]  # all samples except i
        with Timer("Computing coexpression network:"):
            if self.computing=='gpu':
                import cupy as cp
                correlation_matrix = cp.corrcoef(self.expression_matrix[:, idx])
                if cp.isnan(correlation_matrix).any():
                    cp.fill_diagonal(correlation_matrix, 1)
                    correlation_matrix = cp.nan_to_num(correlation_matrix)
                correlation_matrix=cp.asnumpy(correlation_matrix)
            else:
                correlation_matrix = np.corrcoef(self.expression_matrix[:, idx])
                if np.isnan(correlation_matrix).any():
                    np.fill_diagonal(correlation_matrix, 1)
                    correlation_matrix = np.nan_to_num(correlation_matrix)

        with Timer("Normalizing networks:"):
            correlation_matrix_orig = correlation_matrix # save matrix before normalization
            correlation_matrix = self._normalize_network(correlation_matrix)

        with Timer("Inferring LIONESS network:"):
            if self.motif_matrix is not None:
                del correlation_matrix_orig
                subset_panda_network = self.panda_loop(correlation_matrix, np.copy(self.motif_matrix), np.copy(self.ppi_matrix),self.computing)
            else:
                del correlation_matrix
                subset_panda_network = correlation_matrix_orig

        lioness_network = self.n_conditions * (self.network - subset_panda_network) + subset_panda_network

        with Timer("Saving LIONESS network %d to %s using %s format:" % (i+1, self.save_dir, self.save_fmt)):
            path = os.path.join(self.save_dir, "lioness.%d.%s" % (i+1, self.save_fmt))
            if self.save_fmt == 'txt':
                np.savetxt(path, lioness_network)
            elif self.save_fmt == 'npy':
                np.save(path, lioness_network)
            elif self.save_fmt == 'mat':
                from scipy.io import savemat
                savemat(path, {'PredNet': lioness_network})
            else:
                print("Unknown format %s! Use npy format instead." % self.save_fmt)
                np.save(path, lioness_network)
        # if i == 0:
            # self.total_lioness_network = np.fromstring(np.transpose(lioness_network).tostring(),dtype=lioness_network.dtype)
        # else:
        #    self.total_lioness_network=np.column_stack((self.total_lioness_network ,np.fromstring(np.transpose(lioness_network).tostring(),dtype=lioness_network.dtype)))

        self.total_lioness_network = np.fromstring(np.transpose(lioness_network).tostring(),dtype=lioness_network.dtype)
        return self.total_lioness_network

    def save_lioness_results(self, file='lioness'):
        """
        Description:
            Saves LIONESS network.

        Outputs:
            file: Path to save the network.
        """
        #self.lioness_network.to_csv(file, index=False, header=False, sep="\t")
        if self.save_fmt=='npy':
            np.save(file,np.transpose(self.total_lioness_network))
        else:
            self.np.savetxt(file+'.txt', np.transpose(self.total_lioness_network), delimiter="\t",header="")
        return None



