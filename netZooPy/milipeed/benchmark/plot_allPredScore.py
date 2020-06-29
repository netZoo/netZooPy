import re, glob, os, collections, scipy
import numpy as np
# import dask
from matplotlib import colors as mcolors
from IPython.display import Image
# import networkx as nx
import matplotlib.pyplot as plt
from datetime import datetime, date
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
from sklearn import metrics
# import sklearn.svm as svm
from scipy import stats

def plot_predScore(directory):
	allbox = pd.DataFrame()
	traces= glob.glob(directory+'/meltbox.txt')
	for trace in traces:
		buffer=os.path.basename(trace).split('_')[2]
        metlbox=pd.read_csv(trace,sep='\t')
        meltbox['buffer']=buffer
        allbox=pd.concat([allbox,meltbox],axis=1)
		#after collctive plotting

		# allbox=pd.concat([meltbox1010,meltbox55,meltbox00])
		allbox['cell_buff']=allbox['data - cell line']+'_'+allbox['buffer'].astype(str)

	plt.figure(figsize=(12, 5))
	plt.xticks(rotation=30)
	# col_list = [A549]
	cells=['A549','GM12878', 'HeLa', 'HepG2', 'K562','SKNSH']
	tests=['pwm','me','wg']
	for cell in cells:
	    for test in tests:
	        allbox2=allbox[allbox['cell_buff'].str.contains(pat=cell)]
	        allbox2=allbox2[allbox2['cell_buff'].str.contains(pat=test)]

	        g=sns.boxplot(x='cell_buff',y='AUROC', data=allbox2)
	        g=sns.swarmplot(x='cell_buff',y='AUROC', data=allbox2,
	                      size=2, color=".3", linewidth=0)
	        g.set(ylim=(0, 1))

	        plt.savefig(outdir"sthml_"+cell+test+"_allbox_buff.png",dpi=300,bbox_inches = "tight")
	        plt.close()