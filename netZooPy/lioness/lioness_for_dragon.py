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
        The adjacency matrix of each network is written as an individual file to a specified output directory.

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

    def __init__(self,layer1,layer2,output_dir="dragon-lioness-output",merge_col="id",ext1="_layer1",ext2="_layer2",delim=","):

        """
        Description
        ----------
            Initialize instance of LionessDragon class, load data, and fit the overall network.

        Parameters
        ----------
            layer1 : str
                A file path to an n x (p1+1) matrix of omics layer 1, cleaned and ready for analysis.
                One column of the matrix must index the inner_join between
                omics layers. The name of this column is specified with the merge_col argument. Default is "id".

            layer2 : str
                A file path to an n x (p2+1) matrix of omics layer 2, cleaned and ready for analysis
                One column of the matrix must index the inner_join between
                omics layers. The name of this column is specified with the merge_col argument. Default is "id". 

            output_dir : str
                A folder for storing the output networks. default: "lioness-dragon-output"

            merge_col : str
                Name of the index column for joining the two omics layers. default: "id"

            ext1 : str
                Extension to add to omics labels from layer 1. Default: "_layer1"

            ext2 : str
                Extension to add to omics labels from layer 2. Default: "_layer2"

            delim : str
                Delimiter for input files. Default: ","
        """

        # assign output directory
        self._outdir = output_dir

        # assign extensions
        self._ext1 = ext1
        self._ext2 = ext2
        self._merge_col = merge_col

        # load data
        print("[LIONESS-DRAGON] Loading input data ...")

        self._layer_1 = pd.read_csv(layer1,sep=delim, header=0,index_col=0)

        self._layer_1 = self._layer_1.add_suffix(ext1)
        self._layer_1 = self._layer_1.rename(index=str, columns={self._merge_col+ext1:self._merge_col})
        print(self._layer_1.index)

        self._layer_2 = pd.read_csv(layer2,sep=delim,header=0,index_col=0)

        self._layer_2 = self._layer_2.add_suffix(ext2)
        self._layer_2 = self._layer_2.rename(index=str, columns={self._merge_col+ext1:self._merge_col})
        print(self._layer_2.index)

        self._all_data = pd.merge(self._layer_1,self._layer_2,on = self._merge_col, how="inner") 
        print(self._all_data.index)


        self._indexes = range(self._all_data.shape[0])
        self._cutoff = len(self._indexes)


        self._identifiers = self._all_data.index
        self._lambdas = [0,0]

        self._identifiers = self._all_data.index
        self._lambdas = [0,0]
        print("[LIONESS-DRAGON] Fitting overall DRAGON network ...")
        # run the first round of DRAGON
        all_data = self._all_data

        print("[LIONESS-DRAGON] Splitting data back to separate layers ...")

        # split merged data back to layers
        data_layer1 = self._all_data.filter(regex=ext1,axis=1)
        data_layer2 = self._all_data.filter(regex=ext2,axis=1)
       
        # run DRAGON and store in self._network
        lambdas, lambdas_landscape = estimate_penalty_parameters_dragon(data_layer1,data_layer2)
        self._network = get_partial_correlation_dragon(data_layer1,data_layer2,lambdas)
        self._lambdas = lambdas
        print("[LIONESS-DRAGON] Finished fitting overall DRAGON network ...")

    def set_cutoff(self,cutoff=0):
        self._cutoff = cutoff
        
    def lioness_loop(self,reestimate_lambda=False):

        """
        Description
        ----------

            Run LIONESS with DRAGON.
            Write each network to an individual file.

        Parameters
        ----------
            reestimate_lambda : bool
                If false (default), estimate a single shrinkage parameter lambda
                on the full sample and apply this parameter in estimating each leave-one-out
                network. If true, reestimate a separate lambda within each leave-one-out
                network.

        Outputs
        ----------
           In output directory, writes a sample-specific network adjacency matrix for
           each sample
    
        """
        
        if not os.path.exists(self._outdir):
            os.makedirs(self._outdir)

        print("[LIONESS-DRAGON] Preparing to run a total of " + str(len(self._indexes)) + " networks")
        print("[LIONESS-DRAGON] reestimate_lambda parameter set to: "+str(reestimate_lambda))
        for i in self._indexes[0:self._cutoff]:
            outfile=self._outdir + "/lioness-dragon-" + str(self._identifiers[i]) + ".csv"
            outfile=open(outfile,'w')

            print("[LIONESS-DRAGON] Running LIONESS-DRAGON for sample %d:" % (i + 1))
            idx = [x for x in self._indexes if x != i] 
            with Timer("[LIONESS-DRAGON] Running DRAGON to fit partial correlation network:"):

                # subset leave-one-out data
                all_data = self._all_data.iloc[idx]
               
                # split merged data back to layers

                data_layer1 = all_data.filter(regex=self._ext1)
                data_layer2 = all_data.filter(regex=self._ext2)

                if reestimate_lambda:
                    lambdas, lambdas_landscape = estimate_penalty_parameters_dragon(data_layer1,data_layer2)
                else:
                    lambdas = self._lambdas
                
                # calculate partial correlations
                sub_lioness_network = get_partial_correlation_dragon(data_layer1,data_layer2,lambdas)
                
                # apply LIONESS formula to get individual network
                lioness_network = len(self._indexes) * (self._network - sub_lioness_network) + sub_lioness_network
        
                lioness_df = pd.DataFrame(lioness_network,columns = data_layer1.keys().append(data_layer2.keys()))
                lioness_df.to_csv(outfile)
                outfile.close()
            
        return
