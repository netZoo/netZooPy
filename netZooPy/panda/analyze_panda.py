from __future__ import print_function
from .panda import Panda
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

class AnalyzePanda(Panda):
    """
    Description:
        Plots PANDA network.

    Inputs:
        Panda: PANDA object.
    
    Methods:
        __init__            : Intialize instance of AnalyzePanda class.
        top_network_plot    : Selects top genes.
        __shape_plot_network: Creates plot.
        __create_plot       : Runs plot.
    """
    def __init__(self, panda_data):
        """
        Description:
            Intialize instance of AnalyzePanda class.

        Inputs:
            panda_data : PANDA object 
        """
        if not hasattr(panda_data,'export_panda_results'):
            print()
            raise AttributeError("Panda object does not contain the export_panda_results attribute.\n"+
                "Run Panda with the flag release_memory=False")

        '''Load variables from PANDA.'''
        self.panda_results = pd.DataFrame(panda_data.export_panda_results, columns=['tf','gene','motif','force'])
        #self.panda_results = panda_data.export_panda_results
        return None

    def top_network_plot(self, top = 100, file = 'panda_top_100.png'):
        """
        Description:
            Selects top genes.

        Inputs:
            top : Top number of genes to plot.
            file: File to save the network plot.
        """
        subset_panda_results = self.panda_results.sort_values(by=['force'], ascending=False)
        subset_panda_results = subset_panda_results[subset_panda_results.tf != subset_panda_results.gene]
        subset_panda_results = subset_panda_results[0:top]
        self.__shape_plot_network(subset_panda_results = subset_panda_results, file = file)
        return None

    def __shape_plot_network(self, subset_panda_results, file = 'panda.png'):
        """
        Description:
            Creates plot.

        Inputs:
            subset_panda_results : Reduced PANDA network to the top genes.
            file                 : File to save the network plot.
        """
        #reshape data for networkx
        unique_genes = list(set(list(subset_panda_results['tf'])+list(subset_panda_results['gene'])))
        unique_genes = pd.DataFrame(unique_genes)
        unique_genes.columns = ['name']
        unique_genes['index'] = unique_genes.index
        subset_panda_results = subset_panda_results.merge(unique_genes, how='inner', left_on='tf', right_on='name')
        subset_panda_results = subset_panda_results.rename(columns = {'index': 'tf_index'})
        subset_panda_results = subset_panda_results.drop(['name'], 1)
        subset_panda_results = subset_panda_results.merge(unique_genes, how='inner', left_on='gene', right_on='name')
        subset_panda_results = subset_panda_results.rename(columns = {'index': 'gene_index'})
        subset_panda_results = subset_panda_results.drop(['name'], 1)
        links = subset_panda_results[['tf_index', 'gene_index', 'force']]
        self.__create_plot(unique_genes = unique_genes, links = links, file = file)
        return None

    def __create_plot(self, unique_genes, links, file = 'panda.png'):
        """
        Description:
            Runs the plot.

        Inputs:
            unique_genes : Unique list of PANDA genes.
            links        : Edgdes of the subset PANDA network to the top genes.
            file         : File to save the network plot.
        """
        #plot
        g = nx.Graph()
        g.clear()
        plt.clf()
        g.add_nodes_from(unique_genes['index'])
        edges = []
        for i in range(0, len(links)):
            edges = edges + [(links.iloc[i]['tf_index'], links.iloc[i]['gene_index'], links.iloc[i]['force']/200)]
        g.add_weighted_edges_from(edges)
        labels = {}
        for i in range(0, len(unique_genes)):
            labels[unique_genes.iloc[i]['index']] = unique_genes.iloc[i]['name']
        pos = nx.spring_layout(g)
        nx.draw_networkx(g, pos, labels=labels,
                 node_size=20, font_size=3,
                 alpha=0.3, width = 0.5, linewidths =0.5)
        plt.axis('off')
        plt.savefig(file, dpi=300)
        return None
