import re, netZooPy, graphviz, glob, os, collections 
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
from sklearn import metrics
from sklearn.preprocessing import label_binarize
from netZooPy.panda.panda import Panda


def LM(DEE,gene):
    DEE[gene]=pd.to_numeric(DEE[gene])
    # DEE['age']=DEE['age'].astype(object)
    fmla = (gene+"~ age+ PY+ FEV+sex")
    model = sm.formula.glm(fmla, family=sm.families.Gaussian(),data=DEE).fit()
    results=pd.DataFrame(results_summary_to_dataframe(model,gene))  
    return results

def results_summary_to_dataframe(results,gene):
    pvals = results.tvalues
    coeff = results.params
    conf_lower = results.conf_int()[0]
    conf_higher = results.conf_int()[1]
    results_df = pd.DataFrame({gene+"pvals":pvals,gene+"coeff":coeff,})
    results_df = results_df[[gene+"coeff",gene+"pvals"]]
    return results_df

## compute distinct links per age per set of MILIPEED LTCOPD nets
set='milipeed'
control = pd.read_csv('/udd/redmo/analyses/MILIPEED/control.txt',sep=',',header=0 ,usecols=['number','age.x','FEV1FVC.x','sex.x','BMI.x','packyears.x','sapphire'])
control.columns=['iden','no','age','sex','bmi','FEV','PY']
control['ID']=(control['iden'].str.split('-', expand=True).rename(columns=lambda x: f"string_{x+1}"))['string_2']
control.sort_values(by='ID',ascending=False)
del control['iden']

case = pd.read_csv('/udd/redmo/analyses/MILIPEED/case.txt',sep=',',header=0 ,usecols=['number','age.x','FEV1FVC.x','sex.x','BMI.x','packyears.x','sapphire'])
case.columns=['iden','no','age','sex','bmi','FEV','PY']
case['ID']=(case['iden'].str.split('-', expand=True).rename(columns=lambda x: f"string_{x+1}"))['string_2']
case.sort_values(by='ID',ascending=False)
del case['iden']
EE = control.append(pd.DataFrame(data = case), ignore_index=True)
# EE
EE['COPD'] = np.where(EE['FEV']>.7, 'COPD', 'NO')
EE['age_range'] = np.where(EE['age']<50, 'fourties', np.where(EE['age']<60,
        'fifties',np.where(EE['age']<70,'sixties',np.where(EE['age']<80,'seventies','eighties'))))
del case, control

names1 = [os.path.basename(x) for x in glob.glob('/udd/redmo/analyses/MILIPEED/' + set +'_control/*.pairs')]
names2 = [os.path.basename(x) for x in glob.glob('/udd/redmo/analyses/MILIPEED/' + set +'_case/*.pairs')]
names=names1+names2
j=[i.split('-', 1)[1] for i in names]
names=[i.split('_', 1)[0] for i in j]
EE.set_index('ID', inplace=True)

EE=EE.reindex(names)

r=pd.read_csv('/udd/redmo/analyses/MILIPEED/links2.txt',usecols=[0,1],sep='\t',header=0)
# r.columns=['TF','gene']
control = pd.read_csv('/udd/redmo/analyses/MILIPEED/'+set+'_control/results.txt',sep='\t',header=0)
case = pd.read_csv('/udd/redmo/analyses/MILIPEED/'+set+'_case/results.txt',sep='\t',header=0)
JJ = pd.concat([r.TF,r.gene,control, case], axis=1)
del case, control

COPD_GWAS_genes = pd.read_csv("/udd/redmo/analyses/MILIPEED/COPD_GWAS_genes.csv",usecols=[0])
gwasTF=JJ[JJ['TF'].isin(COPD_GWAS_genes.closest_gene)]
gwasgene=JJ[JJ['gene'].isin(COPD_GWAS_genes.closest_gene)]
gwasONLY=gwasTF[gwasTF['gene'].isin(COPD_GWAS_genes.closest_gene)]
csv=gwasgene
links=gwasgene.TF+'_'+gwasgene.gene
# str.replace(links,'-','')
# links.replace('-','')
DDD=csv.transpose()
DDD.columns=links
DDD=DDD.drop(DDD.index[0:2])
DEE = np.concatenate([EE,DDD], axis=1)
DEE=pd.DataFrame(DEE)
DEE.columns=(EE.join(DDD)).columns ## remove 0:78 when full population can be run
DEE.index=EE.index
DEE.age=pd.to_numeric(DEE.age)
DEE=DEE.round({'age':0})
DEE['FEV']=pd.to_numeric(DEE['FEV'])
DEE['PY']=pd.to_numeric(DEE['PY'])
DEE['age']=DEE['age'].astype(object)
from datetime import datetime, date
ccc="{:%d.%m.%Y}".format(datetime.now())
gene=DEE.columns[8]
results=LM(DEE,gene)
results.T.to_csv(('/udd/redmo/analyses/MILIPEED/MILI_'+set+'_indv_'+ccc+".txt"),sep='\t')
for gene in DEE.columns[8:DEE.shape[1]]:
    results=LM(DEE,gene)
    results.T.to_csv(('/udd/redmo/analyses/MILIPEED/MILI_'+set+'_indv_'+ccc+".txt"),sep='\t',mode='a',header=False)
    print(gene)