"""

 utils4HWRF.py
 --------
    contains functions based on python for generating graphics:
    0. atcfRmax6hr
    1. readBT6hrly
    2. readTrack
    3. readTrack6hrly
    4. astring_to_strings
    5. mapdomain(basin) - sets of (min,max) in longitude and latitude
                          for mapping.
    6. Colors_SaffirSimpson  - returns a set of colors 
			corresponding the Saffir-Simpson scale.
    7. fromAdeck - read adeck for operational guidance
    8. Rwinds6hr(bt,fhr) - returns an average radius for a set of 3 wind thresholds (34,50,64)-kt
    9. grb2subset2nc - returns a (set) of variables of interest from the atm grb2 file.
    *. btRwinds- read Bdeck and returns a set of mean radius for the wind threshold.

    modified 5/4/2021
    modified at 9/11/2019
    by Hyun-Sook Kim 9/13/2015
---------------------------------------------------------------

"""

import numpy as np
import numpy.ma as ma
import struct
import sys

#import matplotlib.mlab

from datetime import datetime, timedelta

# - -------------------------------------------------------------
def fromAdeck(opm,afile,YMDH):
    '''
    --- read adeck and extract forecast for a cycle YMDH from
        a model of interest(opm)
    '''
    import re
    
    with open(afile,'r') as infile:
       atxt=[];
       for n in infile.readlines():
          if re.search(opm, n) and re.search(YMDH, n):
            atxt.append(n)
        
    d0 =[ datetime.strptime(s.split(', ')[2],'%Y%m%d%H') for s in atxt ]
    dt = [ timedelta(int(s.split(', ')[5])) for s in atxt ]

    dts = np.array([ d0[n] + dt[n]/24 for n in range(0,len(dt)) ])
    lns = np.array([ float(s.split(', ')[7][:-1])/10. for s in atxt ])
    lts = np.array([ float(s.split(', ')[6][:-1])/10. for s in atxt ])
    pmin = np.array([ float(s.split(', ')[9]) for s in atxt ])
    vmax = np.array([ float(s.split(', ')[8]) for s in atxt ])

    rx=[]
    for n in range(13,17,1):
      rx.append([float(s.split(', ')[n]) for s in atxt])

    r34=np.mean(rx,axis=0)
    r34=np.array(r34)

    return (dts, lns, lts, pmin, vmax, r34)

def atcfRmax6hr(bt):
    '''
    --- read atcf file and returns a set of variables at 6 hourly including Rmax ---
    '''
    # - find indices for 6 hourly intervals 
    with open(bt,'r') as infile:
        atxt=[]
        for n in infile.readlines():
           s=n[:59]
           st=int(s.split(', ')[5])
           if st%6 == 0:
              atxt.append(n)
    
    fhrs=np.array([ int(s.split(', ')[5]) for s in atxt])
    sfr=np.unique(fhrs)

    # - build a file with the radius of the largest wind threshold 
    #   along (34, 50, 64) kt. 
    #
    afull=[]
    for n in range(0,len(sfr),1):
       ti=(fhrs==sfr[n]).astype(int)
       h=[x for x,y in enumerate(ti) if y==1]
       afull.append(atxt[h[-1]])
    afull=np.array(afull)

    # - select values
    d0 =[ datetime.strptime(s.split(', ')[2],'%Y%m%d%H') for s in afull ]
    dt = [ timedelta(int(s.split(', ')[5])) for s in afull ]

    dts = np.array([ d0[n] + dt[n]/24 for n in range(0,len(dt)) ])
    lns = np.array([ float(s.split(', ')[7][:-1])/10. for s in afull ])
    lts = np.array([ float(s.split(', ')[6][:-1])/10. for s in afull ])
    pmin = np.array([ float(s.split(', ')[9]) for s in afull ])
    vmax = np.array([ float(s.split(', ')[8]) for s in afull ])

    rx=[]
    for n in range(13,17,1):
      rx.append([float(s.split(', ')[n]) for s in afull])
   
    r34=np.mean(rx,axis=0)
    r34=np.array(r34)
    
    return (dts, lns, lts, pmin, vmax, r34)

# -------------------------------------------
def btRwinds(bt,YMDH,units=True):
    '''
    --- read btk file and returns a set of variables corresponding 
        to Vxx (34-kt, 50-kt and 64-kt wind thresholds at fhr ---
        
    --- units= True -> rout in [km]  (default)
              False -> rout in [nm]
    
    --- rout in 3 columns [35-kt,50-kt,64-kt] at maximum.

    '''
    with open(bt,'r') as infile:
        atxt=[]
        for n in infile.readlines():
           s=n[:59]
           st=int(s.split(', ')[5])
           if st%6 == 0:
              atxt.append(n)

    atxt = np.array(atxt)

    valid = np.array([ int(s.split(', ')[2]) for s in atxt])
    ii = [ x for x,y in enumerate(valid) if y==int(YMDH) ]

    atxt = atxt[ii]

    rout=np.zeros(3)
    for k in range(0,len(ii),1):
       rx=[]
       line=atxt[k].split(', ')
       rx = [ float(line[n]) for n in range(13,17,1)]
       if units:
          rout[k]=np.mean(rx)*1.852
       else:
          rout[k]=np.mean(rx)

    return(rout)

    # - find indices for 6 hourly intervals
def Rwinds6hr(bt,fhr,units=True):
    '''
    --- read atcf file and returns a set of variables corresponding
        to Vxx (34-kt, 50-kt and 64-kt wind thresholds at fhr ---

    --- units= True -> rout in [km]  (default)
              False -> rout in [nm]

    --- rout in 3 columns [35-kt,50-kt,64-kt] at maximum.

    '''
    with open(bt,'r') as infile:
        atxt=[]
        for n in infile.readlines():
           s=n[:59]
           st=int(s.split(', ')[5])
           if st%6 == 0:
              atxt.append(n)

    atxt = np.array(atxt)

    fhrs = np.array([ int(s.split(', ')[5]) for s in atxt])
    ii = [ x for x,y in enumerate(fhrs) if y==fhr ]

    atxt = atxt[ii]

    rout=np.zeros(3)
    for k in range(0,len(ii),1):
       rx=[]
       line=atxt[k].split(', ')
       rx = [ float(line[n]) for n in range(13,17,1)]
       if units:
          rout[k]=np.mean(rx)*1.852
       else:
          rout[k]=np.mean(rx)

    return(rout)

# -------------------------------------------

def readTrack(bt):
    '''
    --- read best track or atcf file ---
    '''
    with open(bt,'r') as infile:
        atmp=[]
        for n in infile.readlines():
                atmp.append(n[:59])

    atxt=np.unique(atmp)

    d0 =[ datetime.strptime(s.split(', ')[2],'%Y%m%d%H') for s in atxt ]
    dt = [ timedelta(int(s.split(', ')[5])) for s in atxt ]

    dts =[ d0[n] + dt[n]/24 for n in range(0,len(dt)) ]
    lns = [ float(s.split(', ')[7][:-1])/10. for s in atxt ]
    lts = [ float(s.split(', ')[6][:-1])/10. for s in atxt ]
    pmin = [ float(s.split(', ')[9]) for s in atxt ]
    vmax = [ float(s.split(', ')[8]) for s in atxt ]

    return (dts, lns, lts, pmin, vmax)

def readTrack6hrly(bt):
    '''
    --- read atcf file ---
    '''
    with open(bt,'r') as infile:
        atmp=[]
        for n in infile.readlines():
                atmp.append(n[:59])
    atxt=np.unique(atmp)

    idx=[]
    knt=0
    for s in atxt:
       if int(s.split(', ')[5])%6 == 0:
          idx.append(knt)
       knt += 1
    atxt=atxt[idx]
    if atxt[0]=='\n':
       atxt[1:]

    d0 =[ datetime.strptime(s.split(', ')[2],'%Y%m%d%H') for s in atxt ]
    dt = [ timedelta(int(s.split(', ')[5])) for s in atxt ]

    dts = np.array([ d0[n] + dt[n]/24 for n in range(0,len(dt)) ])
    lns = np.array([ float(s.split(', ')[7][:-1])/10. for s in atxt ])
    lts = np.array([ float(s.split(', ')[6][:-1])/10. for s in atxt ])
    #-- cross the dateline:
    knt=0
    for s in atxt:
       ns=s.split(', ')[7][-1]
       ts=s.split(', ')[6][-1]
       if ns=='W':
         lns[knt]=(lns[knt])*(1)
       if ts=='S':
         lts[knt]=-1*lts[knt]
       knt+=1
     
    lns=np.array(lns)
    lts=np.array(lts)
 
    pmin = np.array([ float(s.split(', ')[9]) for s in atxt ])
    vmax = np.array([ float(s.split(', ')[8]) for s in atxt ])

    return (dts, lns, lts, pmin, vmax)

def readBT6hrly(bt):
    '''
    --- read Best Track file ---
    '''
    with open(bt,'r') as infile:
        atmp=[]
        for n in infile.readlines():
                atmp.append(n[:59])
    atxt=np.unique(atmp)

    idx=[]
    knt=0
    for s in atxt:
       if int(s.split(', ')[5])%6 == 0:
          idx.append(knt)
       knt += 1
    atxt=atxt[idx]
     
    d0 =[ datetime.strptime(s.split(', ')[2],'%Y%m%d%H') for s in atxt ]
    dt = [ timedelta(int(s.split(', ')[5])) for s in atxt ]

    dts = np.array([ d0[n] + dt[n]/24 for n in range(0,len(dt)) ])
    lns = np.array([ float(s.split(', ')[7][:-1])/10. for s in atxt ])
    lts = np.array([ float(s.split(', ')[6][:-1])/10. for s in atxt ])
    pmin = np.array([ float(s.split(', ')[9]) for s in atxt ])
    vmax = np.array([ float(s.split(', ')[8]) for s in atxt ])

    return (dts, lns, lts, pmin, vmax)

def astring_to_strings(arg):
    switcher = {
        'l': "al",
        'e': "ep",
        'w': "wp",
        'c': "cp",
        'b': "io",
        's': "sp",
    }
    return switcher.get(arg, "nothing")

def mapdomain(basin):
    return{
        'al': [-100,-10,5,45],
        'ep': [-180,-80,5,40],
        'cp': [-210,-140,5,40],
        'wp': [100,180,0,40],
        'io': [30,100,5,35],
        'sn': [30,100,-35,-5],
        'sp': [100,200,-35,-5],
    }.get(basin, 0) # 0 for default of basin not found

# ##################################################
#     returns (knot,hPa, color)
# --------------------------------------------------
def Colors_SaffirSimpson(arg):
    return {
        'td': [34, 1020, 'lawngreen'],
        'ts': [64, 1020, 'yellow'],      # (Vmax,Pmin) upper bound
        'c1': [95,  979, 'coral'],         
        'c2': [112, 964, 'red'],
        'c3': [136, 944, 'magenta'],
        'c4': [137, 919, 'darkmagenta'],
        'c5': [137, 919, 'darkmagenta'],
    }.get(arg, 0)        # 0 for detault of arg not found
#    return {
#        'td': [34, 1020, 'cyan'],
#        'ts': [63, 1020, 'green'],
#        'c1': [83,  980, 'yellow'],
#        'c2': [95,  965, 'orange'],
#        'c3': [113, 945, 'red'],
#        'c4': [135, 920, 'magenta'],
#    }.get(arg, 0)        # 0 for detault of arg not found

def plot_Saffir_Simpson_WPscale(arg1,linethick,alpha):
    """ arg1 = vmax or pmin
    """
    from utils4HWRF import Colors_SaffirSimpson
    import matplotlib.pyplot as plt

    #   Saffir-Simpson: [knot, hPa, 'color']
    ss0=Colors_SaffirSimpson('td')
    ss1=Colors_SaffirSimpson('ts')
    ss2=Colors_SaffirSimpson('c1')
    ss3=Colors_SaffirSimpson('c2')
    ss4=Colors_SaffirSimpson('c3')
    ss5=Colors_SaffirSimpson('c4')

    plt.figure(plt.gcf().number)
    if ( arg1[0].lower()=='v' ):
       plt.axhline(ss0[0],linestyle='-',color=ss0[-1],linewidth=linethick,label='TD',alpha=alpha)
       plt.axhline(ss1[0],linestyle='-',color=ss1[-1],linewidth=linethick,label='TS',alpha=alpha)
       plt.axhline(ss2[0],linestyle='-',color=ss2[-1],linewidth=linethick,label='C1',alpha=alpha)
       plt.axhline(ss3[0],linestyle='-',color=ss3[-1],linewidth=linethick,label='C2',alpha=alpha)
       plt.axhline(ss4[0],linestyle='-',color=ss4[-1],linewidth=linethick,label='C3',alpha=alpha)
       plt.axhline(ss5[0],linestyle='-',color=ss5[-1],linewidth=linethick,label='C4',alpha=alpha)
    if ( arg1[0].lower()=='p' ):
       plt.axhline(ss2[1],linestyle='-',color=ss2[-1],linewidth=linethick,alpha=alpha)
       plt.axhline(ss3[1],linestyle='-',color=ss3[-1],linewidth=linethick,alpha=alpha)
       plt.axhline(ss4[1],linestyle='-',color=ss4[-1],linewidth=linethick,alpha=alpha)
       plt.axhline(ss5[1],linestyle='-',color=ss5[-1],linewidth=linethick,alpha=alpha)

    return

def shade_Saffir_Simpson_WPscale(arg1,linethick,alpha):
    """ arg1 = vmax or pmin
    """
    from utils4HWRF import Colors_SaffirSimpson
    import matplotlib.pyplot as plt

    #   Saffir-Simpson: [knot, hPa, 'color']
    ss0=Colors_SaffirSimpson('td')
    ss1=Colors_SaffirSimpson('ts')
    ss2=Colors_SaffirSimpson('c1')
    ss3=Colors_SaffirSimpson('c2')
    ss4=Colors_SaffirSimpson('c3')
    ss5=Colors_SaffirSimpson('c4')

    plt.figure(plt.gcf().number)
    if ( arg1[0].lower()=='v' ):
       plt.axhline(ss0[0],linestyle='-',color=ss0[-1],linewidth=linethick,label='TD',alpha=alpha)
       plt.axhline(ss1[0],linestyle='-',color=ss1[-1],linewidth=linethick,label='TS',alpha=alpha)
       plt.axhline(ss2[0],linestyle='-',color=ss2[-1],linewidth=linethick,label='C1',alpha=alpha)
       plt.axhline(ss3[0],linestyle='-',color=ss3[-1],linewidth=linethick,label='C2',alpha=alpha)
       plt.axhline(ss4[0],linestyle='-',color=ss4[-1],linewidth=linethick,label='C3',alpha=alpha)
       plt.axhline(ss5[0],linestyle='-',color=ss5[-1],linewidth=linethick,label='C4',alpha=alpha)
    if ( arg1[0].lower()=='p' ):
       plt.axhline(ss2[1],linestyle='-',color=ss2[-1],linewidth=linethick,alpha=alpha)
       plt.axhline(ss3[1],linestyle='-',color=ss3[-1],linewidth=linethick,alpha=alpha)
       plt.axhline(ss4[1],linestyle='-',color=ss4[-1],linewidth=linethick,alpha=alpha)
       plt.axhline(ss5[1],linestyle='-',color=ss5[-1],linewidth=linethick,alpha=alpha)

    return

def SaffirSimpsonColor_Vmax(Vkt):
    """ return colors """
    if Vkt <34:
       return (Colors_SaffirSimpson('td')[-1])
    elif (Vkt>=34.0 and Vkt<=63.0):
       return (Colors_SaffirSimpson('ts')[-1])
    elif (Vkt>63.0 and Vkt<=82.0):
       return (Colors_SaffirSimpson('c1')[-1])
    elif (Vkt>82.0 and Vkt<=95.0):
       return (Colors_SaffirSimpson('c2')[-1])
    elif (Vkt>96.0 and Vkt<=112.0):
       return (Colors_SaffirSimpson('c3')[-1])
    elif (Vkt>113.0 and Vkt<=136.0):
       return (Colors_SaffirSimpson('c4')[-1])
    else:
       # (Vkt>136.0):
       return (Colors_SaffirSimpson('c5')[-1])

def SaffirSimpson_track(lons,lats,Vkt,alpha,linewidth):
   
    from utils4HWRF import SaffirSimpsonColor_Vmax
    import matplotlib.pylab as plt
 
    lnc=[]
    ltc=[] 
    for k in range(len(lons)-1):
       lnc.append([lons[k],lons[k+1]])
       ltc.append([lats[k],lats[k+1]])

    for N in range(len(lons)-1):
       plt.plot(lnc[N],ltc[N],color=SaffirSimpsonColor_Vmax(Vkt[N]),linewidth=linewidth,alpha=alpha)


def grb2subset2nc(grb2,vars,WErange,SNrange):
    """ Subsetting HWRF/HMOM parent domain
    using -new_grid 
    and save it to netCDF. -hsk 2019 """

    w2path='/gpfs/hps/usrx/local/nceplibs/grib_util.v1.1.1/exec'
    #vars=':(WTMP|LAND):'
    cwd=os.getcwd()
    tmpdir=os.path.join(cwd,'tmp')

    if not os.path.isdir(tmpdir):
      p=Path(tmpdir)
      p.mkdir(parents=True)

#    if os.path.exists(os.path.join(tmpdir,ncout)):
#      os.remove(os.path.join(tmpdir,ncout))

    ncout=os.path.join(os.path.join(cwd,'tmp'),'prs.'+grbf[-9:-5]+'.nc')
    cmd=os.path.join(w2path,'wgrib2 ')+grbf+' -new_grid latlon '+WErange+' '+SNrange+' hsk.grb2'
    print('... cmd=',cmd)
    os.system (cmd)
    cmd=os.path.join(w2path,'wgrib2 ')+'hsk.grb2 -netcdf '+ncout+' -match '+'"'+vars+'"'
    os.system (cmd)
    return(ncout)

#--- modified Rankine Axisymmetric Vortex profile:
def rankin_vortex(r,Rmax,q):
    """ Wood and White (JAS, 2011)
        inputs:
            r [km]   : distance from a TC center
            Rmax [km]: the maximum radius
            q        : arbitry power (higher integer - high gradient around Rmax;
                                      a fraction - flat profile: low gradient)
    by HSK Feb 2020 """
    rho=r/Rmax
    return ( r/((1+rho**(2**q))**1/q) )


