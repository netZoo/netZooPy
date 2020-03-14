from __future__ import print_function

import sys, os, glob,re 
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
from multiprocessing import Process
from datetime import datetime, date
import statsmodels.api as sm

class AnalyzeMilipeed(Milipeed):
    '''GLM MILIPEED links discriminated by age, sex, BMI, FEV and PY.'''
    def __init__(self,input_path,gene_subset,mili_nets='/udd/redmo/analyses/MILIPEED/mili_subj.txt',links_file='/udd/redmo/analyses/MILIPEED/milipeed_links.txt',meta='/udd/redmo/analyses/MILIPEED/subj_metadata.txt',outdir='.'):
    # def __init__(self,input_path,gene_subset,omili_nets,links_file,meta,utdir='.',):
        '''Load variables from Milipeed.'''
        self.metadata = pd.read_csv(meta,sep='\t',header=0,index_col=0)
        subjmeta = pd.read_csv(mili_nets,sep='\t',names=['subj'],index_col=0)
        self.metadata=self.metadata.merge(subjmeta,left_on=self.metadata.index,right_on=subjmeta.index)
        total_links=pd.read_csv(links_file,sep='\t',names=['TF','gene'])

        # metadata.columns=['ID','age','sex','BMI','FEV','PY']
        ##can include categorizations of variables as well, such as defining COPD/not and binning ages, hoever, then you will want to 
        # metadata['COPD'] = np.where(metadata['FEV']>.7, 'COPD', 'NO')
        # metadata['age_range'] = np.where(metadata['age']<50, 'fourties', np.where(metadata['age']<60,'fifties',np.where(metadata['age']<70,'sixties',np.where(metadata['age']<80,'seventies','eighties'))))
        # metadata.set_index('ID', inplace=True)
        # metadata.reindex(subjects) ##if theyve been sorted since extraction

        
        # path=self.save_dir
        self.path=input_path
        self.outdir=outdir
        traces= os.listdir(input_path)
        # if traces[1].endswith('npy'):
        #     append_data = pd.DataFrame()
        #     for j,trace in enumerate(traces):
        #         filepath = os.path.join(input_path, trace)
                
        #         data=pd.melt(pd.DataFrame(np.load(filepath,mmap_mode='r')))
        #         append_data=append_data.append(pd.DataFrame(data.value),ignore_index=True)
        #     self.population=pd.concat([self.metadata,append_data.T],axis=1,sort=True) ## hopefully this works, untested on small jupyter        
        #     # for covariate in self.metadata.columns: ##all other covariates to numeric
        #     #     self.population[covariate]=pd.to_numeric(self.population[covariate])
        #     # self.population['PY']=pd.to_numeric(self.population['PY'])
        #     self.population=self.population.round({'age':0})
        #     if self.population['age'] is not object: ##convert 
        #         self.population['age']=self.population['age'].astype(object)
        #     self.date="{:%d.%m.%Y}".format(datetime.now())

        #     self.milipeed_analysis= self.__analysis_loop()

        # elif gene_subset is not None:
        append_data = pd.DataFrame()
        gene_sub=pd.read_csv(gene_subset,sep='\t',names=['gene'])
        for j,trace in enumerate(traces):
            filepath = os.path.join(input_path, trace)
            data=pd.DataFrame(pd.read_csv(filepath,sep='\t',header=None,index_col=None))
            data.index=total_links['gene']
            subnet=data.merge(gene_sub,left_on=data.index,right_on='gene')

            append_data=pd.concat([append_data,pd.DataFrame(subnet[0])],sort=True,axis=1)
            # self.append_data=self.append_data.T
            # tmp=pd.DataFrame(pd.read_csv(links_file,header=None,index_col=None, skiprows=i*2500, nrows=2500))                
            # tmp=pd.read_csv(links_file,sep='\t',header=None,dtype=str,index_col=None, skiprows=i*2500, nrows=2500)
        tmp=total_links.merge(gene_sub)
        tmp['TF']=tmp['TF'].str.replace('-','')
        tmp['gene']=tmp['gene'].str.replace('-','')
        append_data=append_data.T
        append_data.index=self.metadata.index
        append_data.columns=tmp['TF']+"_"+tmp['gene']
        self.population=self.metadata.merge(append_data,left_index=True,right_index=True)
        self.population=self.population.round({'age':0})
        # if self.population['age'] is not object: ##convert 
            # self.population['age']=self.population['age'].astype(object)
        self.date="{:%d.%m.%Y}".format(datetime.now())
        del append_data, tmp
        # self.milipeed_analysis= runInParallel(self.__analysis_loop(i),)
        self.milipeed_analysis=self.__analysis_loop()


            
            # links_file.merge

        # else:

        #     traces= glob.glob(input_path+'/*.txt')
        #     # links=pd.melt(pd.DataFrame(pd.read_csv(filepath,sep='\t',header=0,index_col=0)))
        #     # meta='/udd/redmo/analyses/MILIPEED/subj_metadata.txt'

        #     # path='/udd/redmo/analyses/SPIDER/lioness_output2'
        #     # links=9734704
        #     # traces= os.listdir(path)
        #     # filepath = os.path.join(path, trace)
        #     for i in range(np.int(np.round(links/len(total_links),0))): ## test 2500 links at a time, could parallelize
        #         append_data = pd.DataFrame()
        #         traces.sort(key=lambda f: int(re.sub('\D', '', f)))
        #         for j,trace in enumerate(traces):
        #             filepath = os.path.join(input_path, trace)
        #             data=pd.DataFrame(pd.read_csv(filepath,sep='\t',header=None,index_col=None, skiprows=i*2500, nrows=2500))
        #             append_data=pd.concat([append_data,pd.DataFrame(data)],sort=True,axis=1)
        #         # self.append_data=self.append_data.T
        #         # tmp=pd.DataFrame(pd.read_csv(links_file,header=None,index_col=None, skiprows=i*2500, nrows=2500))                
        #         tmp=pd.read_csv(links_file,sep='\t',header=None,dtype=str,index_col=None, skiprows=i*2500, nrows=2500)
        #         tmp[0]=tmp[0].str.replace('-','')
        #         tmp[1]=tmp[1].str.replace('-','')
        #         append_data=append_data.T
        #         append_data.index=self.metadata.index
        #         append_data.columns=tmp[0]+"_"+tmp[1]
        #         self.population=self.metadata.merge(append_data,left_index=True,right_index=True)
        #         self.population=self.population.round({'age':0})
        #         if self.population['age'] is not object: ##convert 
        #             self.population['age']=self.population['age'].astype(object)
        #         self.date="{:%d.%m.%Y}".format(datetime.now())
        #         del append_data, tmp
        #         # self.milipeed_analysis= runInParallel(self.__analysis_loop(i),)
        #         self.milipeed_analysis=self.__analysis_loop()
               

    def iLiM(self,gene):
        self.population[gene]=pd.to_numeric(self.population[gene])
        # fmla = (gene+"~ age+ PY+ FEV+sex")
        fmla = (str(gene) + "~"+self.metadata.columns[1]+"+"+self.metadata.columns[2]+"+"+self.metadata.columns[3]+"+"+self.metadata.columns[4]+"+"+self.metadata.columns[5]) ##restrict metadata input to those for use in GLM
        ## ^^ make this fill automagically
        # model = sm.formula.glm(fmla, family=sm.families.Gaussian(),data=self.population[gene]).fit()
        model = sm.formula.glm(fmla, family=sm.families.Gaussian(),data=self.population[['age','sex','BMI','FEV1','packyears',gene]]).fit()

        self.results=pd.DataFrame(self.results_summary_to_dataframe(model,gene))  
        return self.results

    def results_summary_to_dataframe(self,results,gene):
        pvals = results.tvalues
        coeff = results.params
        # conf_lower = results.conf_int()[0]
        # conf_higher = results.conf_int()[1]
        self.results_df = pd.DataFrame({gene+"pvals":pvals,gene+"coeff":coeff,})
        self.results_df = self.results_df[[gene+"coeff",gene+"pvals"]]
        return self.results_df

    def __analysis_loop(self):
        # count=1
        # gene=self.population.columns[8]
        # results=LM(self.population,gene)
        # results.T.to_csv(('/udd/redmo/analyses/MILIPEED/MILI_'+set+'_indv_'+ccc+".txt"),sep='\t')
        for count,gene in enumerate(self.population.columns[(self.metadata.shape[1]+1):self.population.shape[1]]): ### columns are links now if above append after T worked
            self.results=self.iLiM(gene)   ## ^^ check if len(metadata) == 8
            self.results.T.to_csv(os.path.join(self.outdir,self.date+".txt"),sep='\t',mode='a')
            if (count/100).is_integer():
                # print("determining diff links for:"+ gene+", no.:"+count)
                print(count)
            # count=count+1
        return self.results.T

    def runInParallel(*fns):
        proc = []
        for fn in fns:
            p = Process(target=fn)
            p.start()
            proc.append(p)
        for p in proc:
            p.join()

    def importAnalysis():
        analysis=pd.read_csv('/udd/redmo/analyses/SPIDER/lioness_output2mili_analysis08.03.2020.txt',sep='\t',skiprows=range(3,2500,3),header=0,index_col=0)


## still working on

    def top_network_plot(self, column = 0, top = 100, file = 'milipeed_top_100.png'):
        '''Select top genes.'''
        self.export_panda_results[['force']] = self.milipeed_results.iloc[:,column]
        plot = AnalyzePanda(self)
        plot.top_network_plot(top, file)
        return None

