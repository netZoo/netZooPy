
import glob, os
import numpy as np
from matplotlib import colors as mcolors
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from pathlib import Path
import warnings

warnings.filterwarnings('ignore')
warnings.simplefilter('ignore')
# import matplotlib.backends.backend_pdf
table=[]
tab=[]

def predScore(indir='data/MotifPipeline/sthlm_motif1010',outdir='data/MotifPipeline/test/sthlm_motif1010/',cell=None,TF=None):
    traces= glob.glob(indir+'/*')
    traces = list(filter(lambda file: os.stat(file).st_size > 0, traces))
    WWW=3
    
    Path(outdir).mkdir(parents=True, exist_ok=True)
    # pdf = matplotlib.backends.backend_pdf.PdfPages("output.pdf")
    # X=np.random.randint(0,high=len(traces))

    for jac,trace in enumerate(traces):
        if cell is not None and TF is not None:
            indices = [i for i, s in enumerate(traces) if cell+'_'+TF in s]
            trace=traces[indices[0]]
        elif cell is not None and TF is None:
            indices = [i for i, s in enumerate(traces) if cell in s]
            trace=traces[indices[0]]
        elif cell is None and TF is not None:
            indices = [i for i, s in enumerate(traces) if TF in s]
            trace=traces[indices[0]]
        data=pd.read_csv(trace,sep='\t')
        if (data.shape[1])>16:
            data=pd.read_csv(trace,sep='\t',usecols=[0,1,2,4,10,11,12,13,17,18],names=["chr", "start", "end",'weight',"hits1",'W1','array','ChIPTF','gene'])
        elif(data.shape[1])<12:
            data=pd.read_csv(trace,sep='\t',usecols=[0,1,2,4,10,11,12,13,14],names=["chr", "start", "end",'weight',"hits1",'W1','array','ChIPTF','gene'])
        else:
            data=pd.read_csv(trace,sep='\t',usecols=[0,1,2,3,7,8,9,10,14],names=["chr", "start", "end",'weight',"hits1",'W1','array','location','ChIPTF'])

        Col1=os.path.basename(trace).split('_')[0] #cell
        Col2=os.path.basename(trace).split('_')[1] #TF
        
        table2=[]
        tbl=[]
        tmpTBL2=[]
        tmpTBL=[]
        table3=[]
        tblC=[]
        tmpTBLC=[]
        data=data.fillna(0)
        data.ChIPTF=data.ChIPTF.replace('.',0)
        data.ChIPTF[(data.ChIPTF==Col2)]=1
        data=data[(data.ChIPTF==0)|(data.ChIPTF==1)]
        data.ChIPTF=pd.to_numeric(data.ChIPTF)
        if np.sum(data.ChIPTF)>5:
    #     if not data.empty:
    #         data=data[np.abs(data.c2-data.c3)<=1000]
            try:
                data.weight=(data.weight-data.weight.min())/(data.weight.max()-data.weight.min())

            #     data=data[data.ChIPTF!='.']

                data['wgbs']=data.W1/100
                data.wgbs=1-data.wgbs
                data.array=1-data.array
            #     trace=traces[zzzz]
                data['cell']=Col1#os.path.basename(trace).split('_')[0]
                data['TF']=Col2#os.path.basename(trace).split('_')[1]
                tab.append(data)
                data4=data.copy()
            ##examining mearray repeats indicates multiple wgbs events 
                data2=data.copy()
                data2=data2.drop_duplicates()
                df2=data.groupby([data2.chr,data2.start,data2.end,data2.weight,data2.array]).size().reset_index(name='counts')
            ##examining wgbs repeats indicates multiple mearray events 
                data3=data.copy()
            #     del data3['ChIPTF'],data3['gene']
            #     data3=data3[data3['weight']!=data3['wgbs']]
                data3=data3.drop_duplicates()
                df3=data.groupby([data3.chr,data3.start,data3.end,data3.weight,data3.wgbs]).size().reset_index(name='counts')

                data4=data.copy()
            #     del data3['ChIPTF'],data3['gene']
            #     data3=data3[data3['weight']!=data3['wgbs']]
                data4=data4.drop_duplicates()
                df4=data.groupby([data4.chr,data4.start,data4.end,data4.array,data4.wgbs]).size().reset_index(name='counts')


                df=data.groupby([data.chr,data.start,data.end]).size().reset_index(name='counts')
            #     data = data.groupby([data.chr,data.start,data.end]).agg({'weight':'mean',"wgbs":'mean',"ChIPTF":'mean',"array":'min'})
                data = data.groupby([data.chr,data.start,data.end]).agg({'weight':'mean',"wgbs":'mean',"ChIPTF":'max',"array":'mean','hits1':'mean','W1':'mean'})

                data=data.merge(df,on=['chr','start','end'])

                ##make random bimodal distribution for comparison
                down=np.floor(len(data.array)/2)
                up=np.ceil(len(data.array)/2)
                x=1-(np.random.pareto(10, down.astype(int))+ 0)
                y=(np.random.pareto(10, up.astype(int))+ 0)
                z=np.concatenate((x,y),axis=None)
                z[z>1]=1
                z[z<0]=0

                f= plt.figure(figsize=(8, 8))

            ## Plot all three value distributions fit to same 0-1 scale
                plt.subplot(2, 2, 1)
        #         data.ChIPTF[data.ChIPTF!=0]=1

                colors = ['motif','wgbs','array','ChIP']
                plt.hist([pd.to_numeric(data['weight']),pd.to_numeric(data['wgbs']),pd.to_numeric(data['array']),pd.to_numeric(data['ChIPTF'])],alpha=.5,log=True,stacked=False,label=colors,color=['k','tab:blue','tab:orange','deeppink'])
                plt.legend(loc="best")#,bbox_to_anchor=(1,1))
                plt.xlabel('Link Weight')
                plt.ylabel('Frequency')
                plt.title('A. Freq of '+trace.split('/')[3]+'\n link weight')

            ## Plot scatter of WGBS and motif to ChIP
            #     plt.subplot(1, 5, 2)
            #     sns.scatterplot(data.wgbs,data.ChIPTF,label='wgbs',alpha=.5)
            #     sns.scatterplot(data.weight,data.ChIPTF,label='motif',alpha=.5,color='deeppink').set_title("ChIP vs Motif & wgbs")

            ## Plot all ROC curves
            ###PWM
        #         data.ChIPTF[data.ChIPTF!=0]=1

                plt.subplot(2, 2, 3)
                plt.plot([0, 1], [0, 1], 'k--')
                fpr, tpr, thresholds = metrics.roc_curve(data.ChIPTF, data.weight)
                roc_auc=metrics.auc(fpr, tpr)
                plt.plot(fpr, tpr,
                         label='Motif (area = {0:0.2f})'
                               ''.format(roc_auc),
                         color='k', linestyle=':', linewidth=4)

            #     rpwma=np.random.permutation(data.weight)
        #         fpr, tpr, thresholds = metrics.roc_curve(data.ChIPTF,data.shW)
        #         rand_pwm_auroc=metrics.auc(fpr, tpr)
        #         plt.plot(fpr, tpr, 
        #                  label='shuf Motif ({0:0.2f})'
        #                        ''.format(rand_pwm_auroc),
        #                  color='xkcd:light grey', linestyle='-', linewidth=2)
            ###WGBS
                fpr, tpr, thresholds = metrics.roc_curve(data.ChIPTF, data.wgbs)
                roc_auc2=metrics.auc(fpr, tpr)
                plt.plot(fpr, tpr,
                         label='WGBS (area = {0:0.2f})'
                               ''.format(roc_auc2),
                         color='tab:blue', linestyle=':', linewidth=4)

        #     #     rwgbsa=np.random.permutation(data.wgbs)
        #         fpr, tpr, thresholds = metrics.roc_curve(data.ChIPTF, data.shWG)
        #         rand_wgbs_auroc=metrics.auc(fpr, tpr)
        #         plt.plot(fpr, tpr, 
        #                  label='shuf WGBS ({0:0.2f})'
        #                        ''.format(rand_wgbs_auroc),
        #                  color='xkcd:light blue', linestyle='-', linewidth=2)

                fpr, tpr, thresholds = metrics.roc_curve(data.ChIPTF, data.array)
                roc_auc5=metrics.auc(fpr, tpr)
                plt.plot(fpr, tpr,
                         label='MeArray (area = {0:0.2f})'
                               ''.format(roc_auc5),
                         color='tab:orange', linestyle=':', linewidth=4)

            #     rarray=np.random.permutation(data.array)
        #         fpr, tpr, thresholds = metrics.roc_curve(data.ChIPTF, data.shA)
        #         rand_arr_auroc=metrics.auc(fpr, tpr)
        #         plt.plot(fpr, tpr, 
        #                  label='shuf MeArray ({0:0.2f})'
        #                        ''.format(rand_arr_auroc),
        #                  color='xkcd:light orange', linestyle='-', linewidth=2)
        #         
                dataxxx=data[((data['hits1']))>=10]
                fpr, tpr, thresholds = metrics.roc_curve(dataxxx.ChIPTF, dataxxx.wgbs)
                rand_arr_auroc=metrics.auc(fpr, tpr)
                plt.plot(fpr, tpr, 
                         label='threshold WGBS ({0:0.2f})'
                               ''.format(rand_arr_auroc),
                         color='xkcd:red', linestyle='-', linewidth=2)
                fpr, tpr, thresholds = metrics.roc_curve(dataxxx.ChIPTF, dataxxx.array)
                rand_arr_auroc=metrics.auc(fpr, tpr)
                plt.plot(fpr, tpr, 
                         label='thresold array ({0:0.2f})'
                               ''.format(rand_arr_auroc),
                         color='xkcd:green', linestyle='-', linewidth=2)

                fpr, tpr, thresholds = metrics.roc_curve(data.ChIPTF, z)
                rra=metrics.auc(fpr, tpr)
                plt.plot(fpr, tpr, 
                         label='random ({0:0.2f})'
                               ''.format(rra),
                         color='pink', linestyle='-', linewidth=2)

                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.title('C. AUROC')
                plt.legend(loc="best")#,bbox_to_anchor=(1.75,1), title="E.")
                plt.subplots_adjust(wspace=.5)

                plt.subplot(2, 2, 4)
                precision, recall, thresholds = metrics.precision_recall_curve(data.ChIPTF, data.weight)
                roc_auc3=metrics.average_precision_score(data.ChIPTF, data.weight)
                plt.plot(recall,precision, 
                         label='Motif (avg Prec = {0:0.2f})'
                               ''.format(roc_auc3),
                         color='k', linestyle=':', linewidth=4)
        #         rpwmp=np.random.permutation(data.weight)
        #         precision, recall, thresholds = metrics.precision_recall_curve(data.ChIPTF, data.shW)
        #         rand_pwm_aupr=metrics.average_precision_score(data.ChIPTF, data.shW)
        #         plt.plot(recall,precision, 
        #                  label='shuf Motif ({0:0.2f})'
        #                        ''.format(rand_pwm_aupr),
        #                  color='xkcd:light grey', linestyle='-', linewidth=2)

                precision, recall, thresholds = metrics.precision_recall_curve(data.ChIPTF, data.wgbs)
                roc_auc4=metrics.average_precision_score(data.ChIPTF, data.wgbs)
                plt.plot(recall,precision, 
                         label='WGBS (avg Prec = {0:0.2f})'
                               ''.format(roc_auc4),
                         color='tab:blue', linestyle=':', linewidth=4)
                rwgp=np.random.permutation(data.wgbs)
        #         precision, recall, thresholds = metrics.precision_recall_curve(data.ChIPTF, data.shWG)
        #         rand_wgbs_aupr=metrics.average_precision_score(data.ChIPTF, data.shWG)
        #         plt.plot(recall,precision, 
        #                  label='shuf WGBS ({0:0.2f})'
        #                        ''.format(rand_wgbs_aupr),
        #                  color='xkcd:light blue', linestyle='-', linewidth=2)

                precision, recall, thresholds = metrics.precision_recall_curve(data.ChIPTF, data.array)
                roc_auc6=metrics.average_precision_score(data.ChIPTF, data.array)
                plt.plot(recall,precision, 
                         label='MeArray (avg Prec = {0:0.2f})'
                               ''.format(roc_auc6),
                         color='tab:orange', linestyle=':', linewidth=4)
        #         rmep=np.random.permutation(data.array)
        #         precision, recall, thresholds = metrics.precision_recall_curve(data.ChIPTF, data.shA)
        #         rand_me_aupr=metrics.average_precision_score(data.ChIPTF, data.shA)
        #         plt.plot(recall,precision, 
        #                  label='shuf MeArray ({0:0.2f})'
        #                        ''.format(rand_me_aupr),
        #                  color='xkcd:light orange', linestyle='-', linewidth=2)
                precision, recall, thresholds = metrics.precision_recall_curve(dataxxx.ChIPTF, dataxxx.wgbs)
                rand_me_aupr=metrics.average_precision_score(dataxxx.ChIPTF, dataxxx.wgbs)
                plt.plot(recall,precision, 
                         label='threshold WGBS ({0:0.2f})'
                               ''.format(rand_me_aupr),
                         color='xkcd:red', linestyle='-', linewidth=2)

                precision, recall, thresholds = metrics.precision_recall_curve(dataxxx.ChIPTF, dataxxx.array)
                rand_me_aupr=metrics.average_precision_score(dataxxx.ChIPTF, dataxxx.array)
                plt.plot(recall,precision, 
                         label='threshold array ({0:0.2f})'
                               ''.format(rand_me_aupr),
                         color='xkcd:green', linestyle='-', linewidth=2)

        #         randrand=np.random.rand(data.array)
                precision, recall, thresholds = metrics.precision_recall_curve(data.ChIPTF, z)
                rr=metrics.average_precision_score(data.ChIPTF, z)
                plt.plot(recall,precision, 
                         label='random ({0:0.2f})'
                               ''.format(rr),
                         color='pink', linestyle='-', linewidth=2)

                plt.plot([1, 0], [0, 1], 'k--')
                plt.xlabel('Recall or TPR')
                plt.ylabel('Precision')
                plt.title('D. AUPR')
                plt.legend(loc="best")#,bbox_to_anchor=(-1,0))
                plt.subplot(2, 2, 4)
                for jj,start in enumerate(data.start):
                    zz=df3[df3['start']==start]
                    if zz.count!=1:
                        yy=[np.abs(x - y) for i,x in enumerate(zz.wgbs) for j,y in enumerate(zz.wgbs) if i != j]
            #     for jj,start in enumerate(data33.start):
                    zzz=df2[df2['start']==start]
                    if zzz.count!=1:
                        yyy=[np.abs(x - y) for i,x in enumerate(zzz.array) for j,y in enumerate(zzz.array) if i != j]

                    zzzz=df4[df4['start']==start]
                    if zzzz.count!=1:
                        yyyy=[np.abs(x - y) for i,x in enumerate(zzzz.wgbs) for j,y in enumerate(zzzz.array) if i != j]

                    table2.append(np.mean(yy))
                    tmpTBL2.append(len(zz))#/np.math.factorial(len(yy)))
                    tbl.append(np.mean(yyy))
                    tmpTBL.append(len(zzz))#/np.math.factorial(len(yyy)))
                    tblC.append(np.mean(yyyy))
                    tmpTBLC.append(len(zzzz))
                colors = ['wgbs X motif','array X motif']
        #         plt.scatter(table2,tmpTBL2,alpha=1,label='wgbs',color=['w'])
                plt.subplot(2, 2, 2)
        #         plt.scatter(table2,tmpTBL2,alpha=.5,label='wgbs',color=['tab:blue'])
        #         plt.scatter(tbl,tmpTBL,alpha=1,label='array',color=['tab:orange'])
                plt.tight_layout()
        #         plt.figure(figsize=(5,5))
                ar = pd.concat([pd.DataFrame(table2),pd.DataFrame(tmpTBL2)],axis=1,sort=False)
                wg = pd.concat([pd.DataFrame(tbl),pd.DataFrame(tmpTBL)],axis=1,sort=False)
                chi = pd.concat([pd.DataFrame(tblC),pd.DataFrame(tmpTBLC)],axis=1,sort=False)

                ar['type']='array'
                wg['type']='wgbs'
                chi['type']='ChIPTF'
                ar['length']=data.end-data.start
                wg['length']=data.end-data.start
                chi['length']=data.end-data.start
                result = pd.concat([ar, wg, chi])
                result.columns=['A','B','type','length']
                sns.violinplot(x="B", y="A", hue="type",
                    inner="quart",#split=True, #palette="muted",
                    palette={"array": 'tab:orange', "wgbs": 'tab:blue', 'ChIPTF':'deeppink'},
                    data=result,scale="area",
                    height=12, aspect=1,cut=0)
                sns.despine(left=True)


                plt.legend(loc="best")#,bbox_to_anchor=(-1,0))
                plt.xlabel('motif hits')
                plt.ylabel('avg abs pair-diff')
                number=str(np.mean(data.end-data.start))
                plt.title('B. mean absolute pairwise \n difference by CGs/probes \n per motif region ('+number+'bp)')
                plt.savefig(outdir+Col1+"_"+Col2+'.png')
        #         j=np.where(np.array(data.ChIPTF).reshape(1,-1)>.5,1,0)
        #         z=np.where(np.array(data.weight).reshape(1,-1)>.5,1,0)
        #         mcm = metrics.confusion_matrix(j.ravel(),z.ravel())
        #         tn = mcm[0, 0]
        #         tp = mcm[1, 1]
        #         fn = mcm[1, 0]
        #         fp = mcm[0, 1]
        #         mcc1= np.divide((tp* tn - fp * fn), np.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))

        #     #     jj=np.where(np.array(data.ChIPTF).reshape(1,-1)>.5,1,0)
        #         zz=np.where(np.array(data.wgbs).reshape(1,-1)>.5,1,0)
        #         mcm2 = metrics.confusion_matrix(j.ravel(),zz.ravel())
        #         tn2 = mcm2[0, 0]
        #         tp2 = mcm2[1, 1]
        #         fn2 = mcm2[1, 0]
        #         fp2 = mcm2[0, 1]
        #         mcc2= np.divide((tp2* tn2 - fp2 * fn2), np.sqrt((tp2+fp2)*(tp2+fn2)*(tn2+fp2)*(tn2+fn2)))

            #     jjj=np.where(np.array(data.array).reshape(1,-1)>.5,1,0)
        #         zzz=np.where(np.array(data.array).reshape(1,-1)>.5,1,0)
        #         mcm3 = metrics.confusion_matrix(j.ravel(),zzz.ravel())
        #         tn3 = mcm3[0, 0]
        #         tp3 = mcm3[1, 1]
        #         fn3 = mcm3[1, 0]
        #         fp3 = mcm3[0, 1]
        #         mcc3= np.divide((tp3* tn3 - fp3 * fn3), np.sqrt((tp3+fp3)*(tp3+fn3)*(tn3+fp3)*(tn3+fn3)))



            #     Col=trace.split('/')[6]
            #     Col1=Col.split('_')[0]
            #     Col2=Col.split('_')[1] 

                Col3=roc_auc #motif auroc
                Col4=roc_auc2 #wgbs auroc
                Col5=roc_auc3 #motif aupr
                Col6=roc_auc4 #wgbs aupr
                Col10=roc_auc5 #mearray auroc
                Col11=roc_auc6 #mearray aupr
        #         Col7=mcc1 #motif mcc
        #         Col8=mcc2 #wgbs mcc
        #         Col12=mcc3 #mearray mcc
        #         Col9=np.nanmean(table2)
        #         Col13=np.nanmean(tbl)
        #         Col17=np.nanmean(tblC)
        #         Col14=np.nanmean(tmpTBL2)
        #         Col15=np.nanmean(tmpTBL)
        #         Col18=np.nanmean(tmpTBLC)
                Col16=data.size
                Col20=np.max(result.B)
                Col21=np.mean(result.length)
        #         column = Col1, Col2, Col3, Col4, Col10, Col5, Col6, Col11,Col7,Col8,Col12,Col9,Col13,Col14,Col15,Col16
                column = Col1, Col2, Col3, Col4, Col10, Col5, Col6,Col11,Col16,Col20,Col21

    #             table.append(column)

                np.transpose(pd.DataFrame((column))).to_csv(outdir+'sthlm_PRE_overall.txt',mode='a')#,header=['cell','tf','mauroc','wauroc','meauroc','maupr','waupr','meaupr','size','max_hit','mo_length'])
            except ValueError:
                pass
                print('no match for '+Col1+"_"+Col2)
                np.transpose(pd.DataFrame((column))).to_csv(outdir+'sthlm_PRE_overall.txt',mode='a')#,header=['cell','tf','mauroc','wauroc','meauroc','maupr','waupr','meaupr','size','max_hit','mo_length'])
        if cell is not None:
            plt.show()
            break
    #     plt.subplots_adjust(wspace=.4, hspace=.5)
    #     plt.show()
    # dataaaa=pd.concat(tab)
    # dataaaa.to_csv('data/MotifPipeline/encode_sthlm_df_TF2.txt')
    #         pd.DataFrame(table).to_csv('data/MotifPipeline/sthlm_PRE_overall_TF3.txt',mode='a',header=['cell','tf','mauroc','wauroc','meauroc','maupr','waupr','meaupr',
    # #                                                                               'mmcc','wmcc','memcc',
    #                                                                               'size','max_hit','mo_length'])