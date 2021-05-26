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
    Description:
        Plots LIONESS network.

    Inputs:
        Panda: LIONESS object.
    
    Methods:
        __init__            : Intialize instance of AnalyzePanda class.
        top_network_plot    : Selects top genes.
    """
    def __init__(self, lioness_data):
        """
        Description:
            Intialize instance of AnalyzeLioness class and load variables.

        Inputs:
            lioness_data : LIONESS object 
        """
        self.export_panda_results = lioness_data.export_panda_results
        self.lioness_results = lioness_data.export_lioness_results
        return None
    def top_network_plot(self, index = 0, top = 100, file = 'lioness_top_100.png'):
        """
        Description:
            Selects top genes.

        Inputs:
            column: Index of sample to plot.
            top   : Top number of genes to plot.
            file  : File to save the network plot.
        """
        self.export_panda_results[['force']] = self.lioness_results.iloc[index,:]
        plot = AnalyzePanda(self)
        plot.top_network_plot(top, file)
        return None
