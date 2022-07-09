""" estimate footprint values (at the 34-kt radius) for the ocean skill metrics -
    OHC 
    Z26
    MLD_temp_crit
    MLT_temp_crit
    MLD_dens_crit
    MLT_dens_crit
    T100
    PEA

Libraries for the metrics:
  from Upper_ocean_metrics import OHC_from_profile, 
                                  MLD_temp_crit, 
                                  MLD_dens_crit,
                                  T100, 
                                  Potential_energy_anomaly100

--------
history:
--------
           edited by hsk 5/5/2021 """

################################################################################
#%% Frist, user inputs for a storm of interest, tcID, and cycle
storm='isaias'
tcid='09l'
iplot = False
i2pickl = False

# Load Python modules and libraries
import xarray as xr
import netCDF4
import numpy as np
from datetime import datetime, timedelta
import seawater as sw

from utils4HWRF import readTrack6hrly, Rwinds6hr, btRwinds
from geoutils import haversine
import os, sys, glob
import pandas as pd
import matplotlib.pylab as plt

from Upper_ocean_metrics import OHC_from_profile, MLD_temp_crit, MLD_dens_crit,\
                                T100, Potential_energy_anomaly100

from utils import coast180

#########################################################
# Find BT coincide with a cycle of interest
# Bdeck
abdeck='/work/noaa/hwrf/noscrub/input/abdeck'
btk=os.path.join(abdeck,'btk')+'/bal'+tcid[:2]+'2020.dat'
bdn,bln,blt,bpmn,bvmx=readTrack6hrly(btk)
bln=-1*bln+360

#%% 3-hourly 1/12-degree global ocean anal GOFS3.1
# NCODA url and download dataset
url_anal = 'http://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0/ts3z'
anal = xr.open_dataset(url_anal,decode_times=False)

analx,analy = np.meshgrid(anal.lon.data,anal.lat.data)
adum=np.ones(analx.shape)*np.nan
ttm = anal.time
tm = netCDF4.num2date(ttm[:],ttm.units)

#---------- building dataframe for each footprint
names=['YMDH','lon','lat','SST','SSS','OHC','Z26','MLDTc','MLTTc','MLDdc','MLTdc','PEA']
#DF=pd.DataFrame()
#c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13=[],[],[],[],[],[],[],[],[],[],[],[],[]
for N,D in enumerate(bdn):

    #dfc=pd.DataFrame()

    dR=haversine(analx,analy,bln[N],blt[N])/1000.

    ymdh=datetime.strftime(D,'%Y%m%d%H')

    # find a set of radii for the 34-, 50- and 64-kt wind, and apply the 34-kt wind radius as for a footprint
    brads=btRwinds(btk,ymdh)    
    if brads[0] != 0:
      dR[dR>brads[0]]=np.nan

# Find the lat and lon index in the GOFS grid
      ii=np.where(~np.isnan(dR))
      i0=np.unique(ii[0])
      i1=np.unique(ii[1])
      ni0=len(i0)
      ni1=len(i1)

# Footprint (x,y)
      lns,lts=np.meshgrid(analx[0,i1],analy[i0,0])
    
      one=dR/dR
      oneq=np.nan*np.ones([ni0,ni1])
      for r,R in enumerate(i0):
         for c,C in enumerate(i1):
           oneq[r,c]=one[R,C]

#------  each forecast lead time @6 hourly
# Find the time index 
      oktimem = np.where(tm == D)[0]
      time_anal = tm[oktimem[0]]

################################################################################
#%% estimate (OHC.Z26,MLD_temp_cr,MLT_temp_cr,MLD_dens_cr,MLT_dens_cr,T100,PEA) 
#   at each grid cell on a footprint

      # mixed depth by a temperature criteria
      dtemp = 0.2
      ref_depth = 10 # meters

      # mixed depth by a density criteria
      drho = 0.125
      ref_depth = 10 # meters

      Ta=np.squeeze(anal.water_temp[oktimem[0],:,i0,i1].data)
      Sa=np.squeeze(anal.salinity[oktimem[0],:,i0,i1].data)
      Za=np.tile(anal.depth.data,[ni0,ni1,1])
      Za=Za.transpose(2,0,1)
  
      Da=sw.dens(Sa,Ta,Za)

      c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13=[],[],[],[],[],[],[],[],[],[],[],[],[]
      knt=0
      for j in range(ni0):
         for i in range(ni1):
             c1.append(ymdh)
             c2.append(lns.flatten()[knt])
             c3.append(lts.flatten()[knt])

             ti=np.squeeze(Ta[:,j,i])
             if np.isnan(ti[0]):   # only within the radius and sea points
                #[cycle,fhr,lon,lat,SST,SSS,OHC,Z26,MLDtc,MMLTdc,MLDdc,MLTdc,T100,PEA]
                c4.append(ti[0])
                c5.append(ti[0])
                c6.append(ti[0])
                c7.append(ti[0])
                c8.append(ti[0])
                c9.append(ti[0])
                c10.append(ti[0])
                c11.append(ti[0])
                c12.append(ti[0])
                c13.append(ti[0])
             else:
                si=np.squeeze(Sa[:,j,i])
                zi=np.squeeze(Za[:,j,i])
                di=np.squeeze(Da[:,j,i])

                if ti[0]<=26:
                   o=0
                   z26=0
                else:
                  o,z26=OHC_from_profile(zi,ti,di)

                dc,tc = MLD_temp_crit(dtemp,ref_depth,zi,ti)

                ddc,dtc = MLD_dens_crit(drho,ref_depth,zi,ti,di)

                temp100 = T100(zi,ti)
                pea = Potential_energy_anomaly100(zi,di)
              
                c4.append(ti[0])
                c5.append(si[0])
                c6.append(o)
                c7.append(z26)
                c8.append(dc)
                c9.append(tc)
                c10.append(ddc)
                c11.append(dtc)
                c12.append(temp100)
                c13.append(pea)
  
         knt += 1
 
      if iplot:
         SST=np.reshape(c4,[ni0,ni1])
         SSS=np.reshape(c5,[ni0,ni1])
         OHC=np.reshape(c6,[ni0,ni1])
         mldtc=np.reshape(c8,[ni0,ni1])
         mlTtc=np.reshape(c9,[ni0,ni1])
         mlddc=np.reshape(c10,[ni0,ni1])
         mlTdc=np.reshape(c11,[ni0,ni1])
         t100=np.reshape(c12,[ni0,ni1])
         pea=np.reshape(c13,[ni0,ni1])

         cx,cy=coast180()
         cx=cx+360
     
         plt.figure(figsize=(13,10))
         plt.subplot(331)
         plt.contourf(lns,lts,SST*oneq,cmap='Spectral_r'); plt.colorbar()
         plt.plot(bln,blt,'-ko',markersize=3); plt.plot(cx,cy,'gray',alpha=0.6)
         plt.axis([lns.min(),lns.max(),lts.min(),lts.max()])
         plt.title('(A) SST [C]')
         plt.ylabel('latitude [N]')
         plt.subplot(332)
         plt.contourf(lns,lts,SSS*oneq,cmap='Spectral_r'); plt.colorbar()
         plt.plot(bln,blt,'-ko',markersize=3); plt.plot(cx,cy,'gray',alpha=0.6)
         plt.axis([lns.min(),lns.max(),lts.min(),lts.max()])
         plt.title('(B) SSS [psu]')
         plt.subplot(333)
         plt.contourf(lns,lts,OHC*oneq,cmap='Spectral_r'); plt.colorbar()
         plt.title('(C) OHC [kJ/cm$^2$]')
         plt.plot(bln,blt,'-ko',markersize=3); plt.plot(cx,cy,'gray',alpha=0.6)
         plt.axis([lns.min(),lns.max(),lts.min(),lts.max()])
         plt.subplot(334)
         plt.contourf(lns,lts,mldtc*oneq,cmap='Spectral_r'); plt.colorbar()
         plt.ylabel('latitude [N]')
         plt.plot(bln,blt,'-ko',markersize=3); plt.plot(cx,cy,'gray',alpha=0.6)
         plt.axis([lns.min(),lns.max(),lts.min(),lts.max()])
         plt.title('(D) MLD$_{Tc}$ [m]')
         plt.subplot(335)
         plt.contourf(lns,lts,mlddc*oneq,cmap='Spectral_r'); plt.colorbar()
         plt.plot(bln,blt,'-ko',markersize=3); plt.plot(cx,cy,'gray',alpha=0.6)
         plt.axis([lns.min(),lns.max(),lts.min(),lts.max()])
         plt.title('(E) MLD$_{dc}$ [m]')
         plt.subplot(336)
         plt.contourf(lns,lts,mlTtc*oneq,cmap='Spectral_r'); plt.colorbar()
         plt.plot(bln,blt,'-ko',markersize=3); plt.plot(cx,cy,'gray',alpha=0.6)
         plt.axis([lns.min(),lns.max(),lts.min(),lts.max()])
         plt.title('(E) MLT$_{tc}$ [C]')
         plt.subplot(337)
         plt.contourf(lns,lts,mlTdc*oneq,cmap='Spectral_r'); plt.colorbar()
         plt.plot(bln,blt,'-ko',markersize=3); plt.plot(cx,cy,'gray',alpha=0.6)
         plt.xlabel('lontitude [E]'); plt.ylabel('latitude [N]')
         plt.axis([lns.min(),lns.max(),lts.min(),lts.max()])
         plt.title('(F) MLT$_{tc}$ [C]')
         plt.subplot(338)
         plt.contourf(lns,lts,t100*oneq,cmap='Spectral_r'); plt.colorbar()
         plt.plot(bln,blt,'-ko',markersize=3); plt.plot(cx,cy,'gray',alpha=0.6)
         plt.xlabel('lontitude [E]')
         plt.axis([lns.min(),lns.max(),lts.min(),lts.max()])
         plt.title('(G) T100 [C]')
         plt.subplot(339)
         plt.contourf(lns,lts,pea*oneq,cmap='Spectral_r'); plt.plot(bln,blt,'-ko',markersize=3)
         plt.colorbar()
         plt.axis([lns.min(),lns.max(),lts.min(),lts.max()])
         plt.plot(bln,blt,'-ko',markersize=3); plt.plot(cx,cy,'gray',alpha=0.6)
         plt.xlabel('lontitude [E]')
         plt.title('(H) PEA [J/m$^3$]')
         plt.suptitle('GOFS3.1: \n'+storm.upper()+' Analysis Time '+datetime.strftime(D,"%HZ %Y/%m/%d"),fontsize=13)

         plt.savefig('GOFS3p1_'+storm+tcid+'_'+ymdh+'.png',bbox_inches='tight')

         plt.close('all')

    print('(OHC,Z26,MLD_temp_cr,MLT_temp_cr,MLD_dens_cr,MLT_dens_cr,T100,PEA) \
    estimates from GOFS3.1 on a footprint done for forcast hour ',str(N*6))

# --- option to save
if i2pickl:
    dict={'YMDH':c1,     'lead [hr]':c2,    'lon [E]':c3,  'lat [N]':c4,    'SST [C]':c5,\
      'SSS [PSU]':c6, 'OHC [kJ/cm^3]':c7,'Z26 [m]':c8,  'MLDtc [m]':c9,  'MLTtc [C]':c10,\
      'MLDdc [m]':c11,'MLTdc [C]':c12,   'T100 [C]':c13,'PEA [J/m^3]':c14}

    # DataFrame
    dfc = pd.DataFrame.from_dict(dict, orient='columns')

######################################################################################
# Save data to dataframe


################################################################################
