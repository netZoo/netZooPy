from __future__ import print_function

from pypanda import Panda
from pypanda import Lioness
from pypanda import AnalyzePanda

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

class AnalyzeLioness(Lioness):
    '''Network plot for Lioness data.'''
    def __init__(self, lioness_data):
        '''Load variables from lioness.'''
        self.export_panda_results = lioness_data.export_panda_results
        self.lioness_results = lioness_data.export_lioness_results
        return None
    def top_network_plot(self, column = 0, top = 100, file = 'lioness_top_100.png'):
        '''Select top genes.'''
        self.export_panda_results[['force']] = self.lioness_results.ix[:,column]
        plot = AnalyzePanda(self)
        plot.top_network_plot(top, file)
        return None
