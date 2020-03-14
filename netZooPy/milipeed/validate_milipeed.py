from __future__ import print_function

import sys, os
sys.path.insert(1,'../panda')
from netZooPy.panda.panda import Panda
from .milipeed import Milipeed

import re, netZooPy, graphviz, glob, os, collections 
import numpy as np
from IPython.display import Image
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
from sklearn import metrics
import plotly.tools as tls
from sklearn.preprocessing import label_binarize
from netZooPy.panda.panda import Panda
import subprocess
# import matplotlib.backends.backend_pdf

class ValidateMilipeed(Milipeed):
    '''GLM MILIPEED links discriminated by age, sex, BMI, FEV and PY.'''
    def __init__(self,chip_meta,outdir='mili_validation'):
        '''Load variables from Milipeed.'''
        self.chip = self.restrict_ChIP()
        self.motif = self.format_motif(self.chip)
        self.validation=self.run_validation(self.motif)
        self.plots= self.plot_validation(self.validation)

        self.table = self.plot_validation()

        self.plot_heatmap()

    def restrict_ChIP(self,chip_meta='/udd/redmo/data/MotifPipeline/ENCODE/A549_hg19/metadata.tsv',selection='optimal IDR thresholded peaks',assembly='hg19'):
        tfdb=pd.read_csv(chip_meta,sep='\t',header=0)#,sep='\t',names=['motif','TF'])
        sub=tfdb[['File accession','Experiment target','Output type','Assembly']]
        sub['gene']=(sub['Experiment target'].str.split('-', expand=True).rename(columns=lambda x: f"string_{x+1}"))['string_1']
        sub=sub[sub['Output type']==selection]
        sub=sub[sub['Assembly']==assembly]
        del sub['Experiment target'],sub['Output type'],sub['Assembly']
        sub.to_csv(('A549_hg19/meta2IDR.txt'),sep='\t',header=False,index=False)

    def format_motif(self,motif_path='/udd/rekrg/EpiPANDA/FIMO_results/ScanBedResults/',out_path='/udd/redmo/data//MotifPipeline/hg19_refseq_100kb_tr/'):
        traces= os.listdir(motif_path)
        for j,trace in enumerate(traces):
            filepath = os.path.join(motif_path, trace)
            data=pd.read_csv(filepath,sep='\t',names=['loc','pwma','pval','gene','dist'])
            data['chr']=(data['loc'].str.split(':', expand=True).rename(columns=lambda x: f"string_{x+1}"))['string_1']
            data['tmp']=(data['loc'].str.split(':', expand=True).rename(columns=lambda x: f"string_{x+1}"))['string_2']
            data['start']=(data['tmp'].str.split('-', expand=True).rename(columns=lambda x: f"string_{x+1}"))['string_1']
            data['end']=(data['tmp'].str.split('-', expand=True).rename(columns=lambda x: f"string_{x+1}"))['string_2']
            del data['tmp'], data['loc']
            data=data[['chr','start','end','gene','pwm','pval','dist']]
            tt=trace.split('_', 1)[0]
            z=tfdb.TF[tt==tfdb.motif]
            data['name']='name'+data.index.map(str)
            if z.empty is False:
                zz=tfdb.TF[float(str(z)[0:3])]
                data.to_csv((out_path+zz+'.txt'),sep='\t',header=False,index=False)

    # def convert_assembly():
    #     ~/tools/liftOver ~/data/MotifPipeline/ENCODE/wgbsin/ENCFF005TID.txt ~/hg38ToHg19.over.chain ~/data/MotifPipeline/ENCODE/wgbsin/ENCFF005TID_hg19.txt unmatched.txt

    def run_validation(self,outdir,TFdir='/udd/redmo/data/MotifPipeline/ENCODE/A549_hg19',motifdir='/udd/redmo/data/MotifPipeline/hg19_refseq_100kb_tr/',bsfile='/udd/redmo/data/MotifPipeline/ENCODE/wgbsin/ENCFF005TID_hg19.txt'):
        valoutdir = outdir+'/miliVal_outdir';
        subprocess.check_call(["/udd/redmo/netZooPy/netZooPy/milipeed/validate_milipeed.sh TFdir motifdir bsfile valoutdir"],shell=True)
        return valoutdir

    def plot_validation(self):
        table=[]
        val_score = outdir+'/miliVal_outdir';
        pdf = matplotlib.backends.backend_pdf.PdfPages(val_score+"/val_output.pdf")
        traces= os.listdir(val_score)
        for j,trace in enumerate(traces):
            filepath = os.path.join(val_score, trace)
            data=pd.read_csv(filepath,sep='\t',names=["chrMotif", "MSS", "MES",'pwm','pval','gene','chrWGBS','WSS','WES',"wgbs",'chrChIP','CSS','CES',"ChIPTF",'idk']) 

        # Refit data to 0-1    
            data.ChIPTF=data.ChIPTF.replace('.','0')
            data.ChIPTF=pd.to_numeric(data.ChIPTF)
            data.ChIPTF=data.ChIPTF/1000
            data=data[data.wgbs<=100]
            data['weight']=((data.pwm-min(data.pwm))/(max(data.pwm)-min(data.pwm)))

        ## Call ChIP Binding (extreme/safe binarize)
            data3=data.copy()
            df=data3.groupby(data3.gene).size().reset_index(name='counts')
            data3 = data3.groupby(data3.gene).agg({"chrMotif":'first',"MSS":'first', "MES":'first','weight':'first',"wgbs":sum,"ChIPTF":'first'})
            data3=data3.merge(df,on='gene')
            data2=data3.copy()
            data3.ChIPTF[data3.ChIPTF>0.0001]=1 ## start with 1000 confidence
            data3.ChIPTF[data3.ChIPTF<0.0001]=0
            data.wgbs=1-(data.wgbs/100)
            data3.wgbs=1-(data3.wgbs/100)
        #     data3['density']=data3['counts']/((data3['End']-data3['Start'])/2)
            data3.wgbs[data3.wgbs<0]=0
            data2.wgbs=1-(data2.wgbs/100)
            data2.wgbs[data2.wgbs<0]=0
            
            # if min(data3.ChIPTF)==0:
            f= plt.figure(figsize=(14, 4))
        ## Plot all three value distributions fit to same 0-1 scale
            plt.subplot(1, 5, 1)
            colors = ['wgbs','motif','ChIP']
            plt.hist([data3.wgbs,data3.weight,data3.ChIPTF],alpha=.5,log=True,stacked=False,label=colors,color=['tab:blue','deeppink','k'])
            plt.legend(loc="upper center")
            plt.xlabel('Link Weight')
            plt.ylabel('Frequency')
            plt.title('A549-'+trace+' links')

        ## Plot scatter of WGBS and motif to ChIP
            plt.subplot(1, 5, 2)
            sns.scatterplot(data.wgbs,data.ChIPTF,label='wgbs',alpha=.5)
            sns.scatterplot(data.weight,data.ChIPTF,label='motif',alpha=.5,color='deeppink').set_title("ChIP vs Motif & wgbs")

        ## Plot all ROC curves
        ###PWM
            fpr, tpr, thresholds = metrics.roc_curve(data3.ChIPTF, data3.weight)
            roc_auc=metrics.auc(fpr, tpr)
            plt.subplot(1, 5, 3)
            plt.plot(fpr, tpr,
                     label='Motif (area = {0:0.2f})'
                           ''.format(roc_auc),
                     color='deeppink', linestyle=':', linewidth=4)
        ###WGBS
            fpr, tpr, thresholds = metrics.roc_curve(data3.ChIPTF, data3.wgbs)
            roc_auc=metrics.auc(fpr, tpr)
            plt.plot(fpr, tpr,
                     label='WGBS (area = {0:0.2f})'
                           ''.format(roc_auc),
                     color='tab:blue', linestyle=':', linewidth=4)

            plt.plot([0, 1], [0, 1], 'k--')
            plt.xlabel('False Positive Rate')
            plt.ylabel('True Positive Rate')
            plt.title('AUROC')
            plt.legend(loc="lower right")
            plt.subplots_adjust(wspace=.5)

        ## Overlap of link calls
            fig=plt.subplot(1, 5, 5)
            labels = ['TPTP','FPFP','FNFN','TNTN','TPFN','FNTP','FPTN','TNFP']
            pwm = [
                len(data[(data.weight>.5)&(data.ChIPTF>.5)&data.wgbs>.5]),
                len(data[(data.weight>.5)&(data.ChIPTF<.5)&data.wgbs>.5]),    
                len(data[(data.weight<.5)&(data.ChIPTF>.5)&data.wgbs<.5]),
                len(data[(data.weight<.5)&(data.ChIPTF<.5)&data.wgbs<.5]),
                len(data[(data.weight>.5)&(data.ChIPTF>.5)&data.wgbs<.5]),
                len(data[(data.weight<.5)&(data.ChIPTF>.5)&data.wgbs>.5]),
                len(data[(data.weight>.5)&(data.ChIPTF<.5)&data.wgbs<.5]), 
                len(data[(data.weight<.5)&(data.ChIPTF<.5)&data.wgbs>.5])]

            x = np.arange(len(labels))
            width = 0.35
            PWM=plt.bar(x ,np.round(np.log(pwm),decimals=4),width,label='scores',color='tab:blue')
            plt.xticks(x,labels,rotation=90)
            plt.ylabel('log count')
            plt.title('matched links \n motif then wgbs')

            plt.subplot(1, 5, 4)
            colors = ['motif']
            cc=np.unique(data.gene, return_counts=True)[1]
            plt.hist(cc,alpha=.5,log=True,stacked=False,label=colors,color=['deeppink'])
            plt.legend(loc="upper right")
            plt.xlabel('motif')
            plt.ylabel('Frequency')
            plt.title('freq wbgs event \n binding in same motif')
            # plt.show()
            Col=trace.split('/')[6]
            Col1=Col.split('_')[0]
            Col2=Col.split('_')[1] 
            
            Col3=roc_auc
            Col4=roc_auc2
            column = Col1, Col2, Col3, Col4
            table.append(column)
            
            plt.show()
        return table
        # for fig in range(1, plt.gcf().number + 1):
        #     pdf.savefig( fig )
        # pdf.close()

        def plot_heatmap(self):
            aurocs=pd.DataFrame(self.table)
            aurocs.columns=['cell','TF','pwm','wgbs']
            aurocs['diff']=aurocs['wgbs']-aurocs['pwm']
            del aurocs['wgbs'], aurocs['pwm']
            heat=aurocs.pivot_table(index=['TF'], columns='cell')
            plt.figure(figsize=(12, 20))
            # plt.imshow(heat, cmap='hot', interpolation='nearest')
            rdgn = sns.diverging_palette(h_neg=130, h_pos=10, s=99, l=55, sep=3, as_cmap=True)
            ax= sns.heatmap(heat, center=0,cmap=rdgn)
            # plt.show
            plt.savefig(self.outdir,dpi=300)


