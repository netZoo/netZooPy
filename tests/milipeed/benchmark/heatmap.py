# table2=table
# table.to_csv('data/MotifPipeline/PRE_overall.txt')
aurocs=pd.read_table('data/MotifPipeline/sthlm_PRE_overall_TF5eee.txt',sep=',',usecols=[1,2,3,4,5,6,7,8],names=['cell','TF','mo_auroc','wgbs_auroc','me_auroc','mo_aupr','wgbs_aupr','me_aupr'])
aurocs=aurocs[aurocs.cell!='0']
# aurocs=pd.DataFrame(table)
# aurocs['wg_auroc']=aurocs['wgbs_auroc']-aurocs['mo_auroc']
# aurocs['me_auroc']=aurocs['me_auroc']-aurocs['mo_auroc']
# aurocs['wg_aupr']=aurocs['wgbs_aupr']-aurocs['mo_aupr']
# aurocs['me_aupr']=aurocs['me_aupr']-aurocs['mo_aupr']
# aurocs.drop(['mo_auroc','wgbs_auroc','me_auroc','mo_aupr','wgbs_aupr','me_aupr'],axis=1,inplace=True)
# aurocs=aurocs[aurocs['wgbs_auroc_imp']>.2]
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
# heat=heat[heat.index.str.contains(pat='CEB')]

# heat=heat[heat['count']==5]
# heat=heat[heat['mean']<0.1]
# del heat.iloc[:,0]
heat=heat.dropna(axis=1, how='all')
del heat['count']
del heat['mean']
del heat['weight']
# heat['sum']=heat.isnull().sum(axis=0)
# heat=heat.sort_values(by=['sum'],ascending=True)
# del heat['sum']
heat=pd.DataFrame(heat.to_records())

heat.columns=['TF','me-A549','me-GM12878','me-H1','me-HeLa','me-HepG2','me-K562','me-SKNSH',
             'pwm-A549','pwm-GM12878','pwm-H1','pwm-HeLa','pwm-HepG2','pwm-K562','pwm-SKNSH',
             'gw-A549','wg-GM12878','wg-H1','wg-HeLa','wg-HepG2','wg-K562','wg-SKNSH']

heat=heat[['TF','pwm-A549','pwm-GM12878','pwm-H1','pwm-HeLa','pwm-HepG2','pwm-K562','pwm-SKNSH',
          'me-A549','me-GM12878','me-H1','me-HeLa','me-HepG2','me-K562','me-SKNSH',
             'gw-A549','wg-GM12878','wg-H1','wg-HeLa','wg-HepG2','wg-K562','wg-SKNSH']]
heat=heat.set_index('TF')

plt.figure(figsize=(12, 30))
# grid_kws = {"height_ratios": (60,12), "hspace": .3}

# f, (ax, cbar_ax) = plt.subplots(2)#, gridspec_kw=grid_kws)

ax = sns.heatmap(heat, linewidth=.1,cmap="YlGnBu", cbar=False)#cbar_ax=cbar_ax,cbar_kws={"orientation": "horizontal"})
mappable = ax.get_children()[0]
plt.colorbar(mappable, ax = [ax],orientation = 'horizontal',pad=.04) #.02 with 12x60, .03 for 12x40
plt.xticks(rotation=30)
# plt.imshow(heat, cmap='hot', interpolation='nearest')
# rdgn = sns.diverging_palette(h_neg=200, h_pos=20, s=99, l=55, sep=3, as_cmap=True)
# ax= sns.heatmap(heat, center=0,cmap=rdgn,linewidth=.1)
# ax= sns.heatmap(heat, linewidth=.1,cmap="YlGnBu")

plt.savefig("netZooPy/tests/milipeed/ToyFigs/sthml_2Diff"+measure+"_heat5eee.png",dpi=300,bbox_inches = "tight")
plt.show