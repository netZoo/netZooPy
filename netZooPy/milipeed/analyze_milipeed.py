from __future__ import print_function

import sys, os
sys.path.insert(1,'../panda')
from netZooPy.panda.panda import Panda
from .milipeed import Milipeed
from netZooPy.panda.analyze_panda import AnalyzePanda

import numpy as np
import collections
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import netZooPy
import graphviz
import mpld3
from datetime import datetime, date

class AnalyzeMilipeed(Milipeed):
    '''GLM MILIPEED links discriminated by age, sex, BMI, FEV and PY.'''
    def __init__(self, metadata='/udd/redmo/analyses/MILIPEED/subj_metadata.txt',gene_subset="/udd/redmo/analyses/MILIPEED/COPD_GWAS_genes.csv",outdir='mili_analysis'):
        '''Load variables from Milipeed.'''
        metadata = pd.read_csv(metadata,sep='\t',header=0,index_col=0)
        # metadata.columns=['ID','age','sex','BMI','FEV','PY']
        ##can include categorizations of variables as well, such as defining COPD/not and binning ages, hoever, then you will want to 
        # metadata['COPD'] = np.where(metadata['FEV']>.7, 'COPD', 'NO')
        # metadata['age_range'] = np.where(metadata['age']<50, 'fourties', np.where(metadata['age']<60,'fifties',np.where(metadata['age']<70,'sixties',np.where(metadata['age']<80,'seventies','eighties'))))
        # metadata.set_index('ID', inplace=True)
        metadata.reindex(subjects)

        appended_data = []
        path=self.save_dir
        traces= os.listdir(path)
        for j,trace in enumerate(traces):
            filepath = os.path.join(path, trace)
            if self.save_fmt == 'npy':
                data=pd.melt(pd.DataFrame(np.load(filepath)))
            else:
                data=pd.melt(pd.DataFrame(np.loadtxt(filepath)))
            if j==0:
                append_data=pd.DataFrame(data.value)
            else:
                append_data=pd.concat([append_data,pd.DataFrame(data.value)],axis=1)

        self.population=pd.concat([metadata,append_data.T],axis=1) ## hopefully this works, untested on small jupyter        
        for covariate in self.metadata.columns: ##all other covariates to numeric
            self.population[covariate]=pd.to_numeric(self.population[covariate])
        # self.population['PY']=pd.to_numeric(self.population['PY'])
        self.population=self.population.round({'age':0})
        if self.population['age'] is not object: ##convert 
            self.population['age']=self.population['age'].astype(object)
        self.date="{:%d.%m.%Y}".format(datetime.now())

        self.milipeed_analysis= self.__analysis_loop()

    def LM(self,gene):
        self.population[gene]=pd.to_numeric(self.population[gene])
        # fmla = (gene+"~ age+ PY+ FEV+sex")
        fmla = (gene + "~"+self.metadata.columns) ##restrict metadata input to those for use in GLM
        model = sm.formula.glm(fmla, family=sm.families.Gaussian(),data=self.population[gene]).fit()
        self.results=pd.DataFrame(results_summary_to_dataframe(model,gene))  
        return self.results

    def results_summary_to_dataframe(results,gene):
        pvals = results.tvalues
        coeff = results.params
        # conf_lower = results.conf_int()[0]
        # conf_higher = results.conf_int()[1]
        self.results_df = pd.DataFrame({gene+"pvals":pvals,gene+"coeff":coeff,})
        self.results_df = self.results_df[[gene+"coeff",gene+"pvals"]]
        return self.results_df

    def __analysis_loop(self):
        # gene=self.population.columns[8]
        # results=LM(self.population,gene)
        # results.T.to_csv(('/udd/redmo/analyses/MILIPEED/MILI_'+set+'_indv_'+ccc+".txt"),sep='\t')
        for gene in self.population.columns[len(metadata):self.population.shape[1]]: ### columns are links now if above append after T worked
            self.results=LM(self.population,gene)   ## ^^ check if len(metadata) == 8
            self.results.T.to_csv(('/udd/redmo/analyses/MILIPEED/MILI_analysis'+self.date+".txt"),sep='\t')
            return self.results.T

## still working on

    def top_network_plot(self, column = 0, top = 100, file = 'milipeed_top_100.png'):
        '''Select top genes.'''
        self.export_panda_results[['force']] = self.milipeed_results.iloc[:,column]
        plot = AnalyzePanda(self)
        plot.top_network_plot(top, file)
        return None

