aurocs=pd.read_table('data/MotifPipeline/sthlm_PRE_overall_TF5eee.txt',sep=',',usecols=[1,2,3,4,5,6,7,8],names=['cell','TF','mo_auroc','wgbs_auroc','me_auroc','mo_aupr','wgbs_aupr','me_aupr'])
aurocs=aurocs[aurocs.cell!='0']
aurocs.drop(['mo_aupr','wgbs_aupr','me_aupr'],axis=1,inplace=True)


###PLOT
# measure="wg_auroc"
measure='auroc'

##sum across cells and sort
heat=aurocs.pivot_table(index=['TF'], columns='cell')
# heat=heat[[measure]]
heat['mean']=np.nanmean(heat,axis=1)
heat['count']=heat.isnull().sum(axis=1)
heat['count']=(21-heat['count'])
# heat=heat.sort_values(by=['count'],ascending=True)
# heat=heat.sort_values(by=['mean'],ascending=False)
heat['weight']=heat['mean']*(heat['count'])
heat=heat.sort_values(by=['weight'],ascending=False)

heat=heat.dropna(axis=1, how='all')
del heat['count']
del heat['mean']
del heat['weight']
box=pd.DataFrame(heat.to_records())
mask = box.isnull()
box.columns=['TF','me-A549','me-GM12878','me-H1','me-HeLa','me-HepG2','me-K562','me-SKNSH',
             'pwm-A549','pwm-GM12878','pwm-H1','pwm-HeLa','pwm-HepG2','pwm-K562','pwm-SKNSH',
             'gw-A549','wg-GM12878','wg-H1','wg-HeLa','wg-HepG2','wg-K562','wg-SKNSH']

box=box[['TF','pwm-A549','pwm-GM12878','pwm-H1','pwm-HeLa','pwm-HepG2','pwm-K562','pwm-SKNSH',
          'me-A549','me-GM12878','me-H1','me-HeLa','me-HepG2','me-K562','me-SKNSH',
             'gw-A549','wg-GM12878','wg-H1','wg-HeLa','wg-HepG2','wg-K562','wg-SKNSH']]

# box.columns=['TF','A549','GM12878','H1','HeLa','HepG2','K562','SKNSH']
meltbox=pd.melt(box,id_vars=['TF'])
# del box.TF
meltbox.columns=['TF','data - cell line','AUROC']
plt.figure(figsize=(12, 5))
plt.xticks(rotation=30)
sns.boxplot(x='data - cell line',y='AUROC', data=meltbox)
sns.swarmplot(x='data - cell line',y='AUROC', data=meltbox,
              size=2, color=".3", linewidth=0)

plt.savefig("netZooPy/tests/milipeed/ToyFigs/sthml_2Diff"+measure+"_box5eee.png",dpi=300,bbox_inches = "tight")
plt.show
# Tweak the visual presentation
# ax.xaxis.grid(True)
# ax.set(ylabel="")
sns.despine(trim=True, left=True)
# plt.show()