import pandas as pd
import numpy as np
import matplotlib.pylab as plt


cycle='2020072918'
#stormid='isaias09l'
stormid='nine09l'

basedir='/work/noaa/hwrf/save/hskim/Isaias/scripts/footprint'

M=pd.read_pickle(basedir+'/output/HMON/'+cycle+'/'+stormid.upper()+'_OHC.pkl')
W=pd.read_pickle(basedir+'/output/HWRF/'+cycle+'/'+stormid.upper()+'_OHC.pkl')
N=pd.read_pickle(basedir+'/output/GOFS3p1/'+cycle+'/ISAIAS09L_OHC.pkl')

ii=M.YMDH.unique()

wo,mo,no=[],[],[]
wstd,mstd,nstd=[],[],[]
for k, ic in enumerate(ii):
    mo.append(M.loc[M.YMDH==ic].mean(axis=0)[4])
    mstd.append(M.loc[M.YMDH==ic].std(axis=0)[4])

    wo.append(W.loc[W.YMDH==ic].mean(axis=0)[4])
    wstd.append(W.loc[W.YMDH==ic].std(axis=0)[4])

    no.append(N.loc[N.YMDH==ic].mean(axis=0)[4])
    nstd.append(N.loc[N.YMDH==ic].std(axis=0)[4])

mo=np.array(mo)
mstd=np.array(mstd)
wo=np.array(wo)
wstd=np.array(wstd)
no=np.array(no)
nstd=np.array(nstd)

plt.figure(figsize=(8,5))
plt.fill_between(ii,no-nstd,no+nstd,color='k',alpha=0.2)
plt.plot(ii,no,'-ok',markersize=4,label='GOFS3.1')
plt.fill_between(ii,mo-mstd,mo+mstd,color='g',alpha=0.2)
plt.plot(ii,mo,'-og',markersize=4,label='HMON')
plt.fill_between(ii,wo-wstd,wo+wstd,color='b',alpha=0.2)
plt.plot(ii,wo,'-ob',markersize=4,label='HWRF')
plt.ylim([-20,100])
plt.axhline(y=0,color='gray',linewidth=0.7)
plt.legend()
plt.title(stormid.upper()+'\n mean OHC [kJcm$^{-3}$]'+' for '+cycle,fontsize=12)
plt.ylabel('Ocean Heat Content',fontsize=11)
plt.savefig('output/cf_meanOHC_34kt_'+stormid+'_'+cycle+'.png',bbox_inches='tight')

plt.show()
