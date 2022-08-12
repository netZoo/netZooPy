from __future__ import print_function

import sys
sys.path.insert(1,'../panda')
from netZooPy.panda.panda import Panda
from .lioness import Lioness
from netZooPy.panda.analyze_panda import AnalyzePanda
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

class AnalyzeLioness(Lioness):
    """
        Plots LIONESS network.

    Parameters
    ------------
        lioness_data: object
            lioness object.
    
    """
    def __init__(self, lioness_data):
        """
        
            Intialize instance of AnalyzeLioness class and load variables.
        """
        self.export_panda_results = lioness_data.export_panda_results
        self.lioness_results = lioness_data.export_lioness_results
        return None
    def top_network_plot(self, index = 0, top = 100, file = 'lioness_top_100.png'):
        """
            Network of top genes.

        Parameters
        -----------
            index: int (defaults to 0)
                Index of sample to plot.
            top   : int (defaults to 100)
                Top number of genes to plot.
            file  : str
                File to save the network plot.
        """
        # we added 2 to the index to account for TF and Gene columns
        self.export_panda_results[['force']] = self.lioness_results.iloc[:,index+2]
        plot = AnalyzePanda(self)
        plot.top_network_plot(top, file)
        return None
