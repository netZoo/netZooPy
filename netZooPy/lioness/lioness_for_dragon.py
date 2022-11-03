from __future__ import print_function
import os.path as path
import os, sys
import numpy as np
import pandas as pd
import numpy.ma as ma

sys.path.insert(1, "../panda")
from netZooPy.dragon import *
from .timer import Timer

class LionessDragon():
    """
    Description
    ----------
        This function uses LIONESS to infer single-sample partial correlation networks.
        Currently, all networks are written as ROWS to a single output file.
        Note this differs from LIONESS, which writes networks as columns.
        The purpose of the file format is to enable reading in one network at a time for downstream analyses.

    Methods
    ----------
        __init__             : Initialize instance of LionessDragon class and load data.
        __lioness_loop       : Run the LIONESS algorithm, saving individual networks one-at-a-time as computed


    Example
    ----------
        
    Authors
    ----------
        katehoffshutta
    """

    _indexes = []

    def __init__(self,layer1,layer2,output_dir="dragon-lioness-output"):

        """
        Description
        ----------
            Initialize instance of LionessDragon class and load data.

        Parameters
        ----------
            layer1 : str
                A file path to an n x (p1+1) matrix of omics layer 1, cleaned and ready for analysis.
                One column of the matrix must be labeled "id" and index the inner_join between
                omics layers. 

            layer2 : str
                A file path to an n x (p2+1) matrix of omics layer 2, cleaned and ready for analysis
                One column of the matrix must be labeled "id" and index the inner_join between
                omics layers.

            output_dir : str
                A folder for storing the output networks. default: "lioness-dragon-output"
        """

        # assign output directory
        self._outdir = output_dir

        # load data
        print("[LIONESS-DRAGON] Loading input data ...")

        self._layer_1 = pd.read_csv(layer1,sep="\t", header=0,index_col=0)
        self._layer_2 = pd.read_csv(layer2,sep="\t",header=0,index_col=0)

        # merge to ensure ordering matches
        self._all_data = pd.merge(self._layer_1,self._layer_2,on="id",how="inner", suffixes=('_layer1', '_layer2'))
        self._indexes = range(self._all_data.shape[0])

        print("[LIONESS-DRAGON] Fitting overall DRAGON network ...")
        # run the first round of DRAGON
        all_data = self._all_data

        # split merged data back to layers
        meth_data = all_data.filter(regex='layer1')
        exp_data = all_data.filter(regex='layer2')

        # run DRAGON and store in self._network
        lambdas, lambdas_landscape = estimate_penalty_parameters_dragon(meth_data,exp_data)
        self._network = get_partial_correlation_dragon(meth_data,exp_data,lambdas)

    def lioness_loop(self):
        """
        Description
        ----------

            Run LIONESS with DRAGON.
            Write each network to an individual file.

        Outputs
        ----------
           In output directory, writes a sample-specific network adjacency matrix for
           each sample
        """
        
        if not os.path.exists(self._outdir):
            os.makedirs(self._outdir)

        for i in self._indexes:
            outfile=self._outdir + "/lioness-dragon-" + str(i) + ".csv"
            outfile=open(outfile,'w')

            print("[LIONESS-DRAGON] Running LIONESS-DRAGON for sample %d:" % (i + 1))
            idx = [x for x in self._indexes if x != i] 
            with Timer("[LIONESS-DRAGON] Running DRAGON to fit partial correlation network:"):

                # subset leave-one-out data
                all_data = self._all_data.iloc[idx]

                # split merged data back to layers
                meth_data = all_data.filter(regex='layer1')
                exp_data = all_data.filter(regex='layer2')

                # calculate penalty parameters
                lambdas, lambdas_landscape = estimate_penalty_parameters_dragon(meth_data,exp_data)
                
                # calculate partial correlations
                sub_lioness_network = get_partial_correlation_dragon(meth_data,exp_data,lambdas)

                # apply LIONESS formula to get individual network
                lioness_network = len(self._indexes) * (self._network - sub_lioness_network) + sub_lioness_network
    
                np.savetxt(outfile,lioness_network,delimiter=",", header="")
                outfile.close()
            
        return

