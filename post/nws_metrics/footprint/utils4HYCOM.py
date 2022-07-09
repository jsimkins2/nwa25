"""

 utils4hycom.py
 --------
    contains functions based on python for generating graphics:
    1. readVar - read a HYCOM archive .[ab] file,
    2. readBin - read HYCOM binary [ab] file (level data),
    3. readgrids - read regional.grid.[ab]
    4. readdepth - read regional.depth.[ab]
    5. HYCOMdate2normal
    6. datetime2matlabdn(dt):
    7. find_ijs(in2d, target):
    8. readBinz - read HYCOM 3z level [ab], calling parse_z

    x. arakawa
    
    by Hyun-Sook Kim 9/13/2015
---------------------------------------------------------------
"""
import numpy as np
import numpy.ma as ma
import struct
import sys

from datetime import datetime, timedelta

def parse_b(inFile,fileType='archive'):
    
    hycom={}  #define empty dictionary

    inFile=inFile+'.b'
    lines=[line.rstrip() for line in open(inFile)]

    if fileType=='relax' or fileType=='depth':
        idm=int([line.split() for line in lines if 'i/jdm' in line][0][2])                
        jdm=int([line.split() for line in lines if 'i/jdm' in line][0][3][:-1])          
    elif fileType=='forcing':
        idm=int([line.split() for line in lines if 'i/jdm' in line][0][2])
        jdm=int([line.split() for line in lines if 'i/jdm' in line][0][3])
    elif fileType=='archive' or fileType=='grid':
        idm=int([line.split() for line in lines if 'longitudinal' in line][0][0])
        jdm=int([line.split() for line in lines if 'latitudinal' in line][0][0])
        
    hycom['idm']=idm
    hycom['jdm']=jdm

    count=0
    for line in lines:
        count+=1
        if (fileType=='archive' and line[0:5] == 'field') or \
            (fileType == 'relax' and line[0:5]=='i/jdm') or \
            (fileType == 'depth' and line[0:5]=='i/jdm') or \
            (fileType == 'forcing' and line[0:5]=='i/jdm') or \
	    (fileType == 'grid' and line[0:5]=='plon:'):
            break
        
    if (fileType=='grid'):
       count=count-1
    if (fileType=='depth'):
       count=count+3
    lines=lines[count:]

    if fileType=='relax' or fileType=='forcing':
        vars=[line.split()[0][:-1] for line in lines]
    if fileType=='depth':
        vars=[line.split()[1] for line in lines]
    else:
        vars=[line.split()[0] for line in lines]

    count=0
    for var in vars:
        count += 1
        if hycom.__contains__(var):
            hycom[var]=hycom[var] +[count] 
        else:
            hycom[var]=[count] 
        
    return hycom

def readVar(inFile,fileType,fieldName,pntidx=None):

    hycom=parse_b(inFile,fileType)

    ijdm=hycom['idm']*hycom['jdm']
    npad=4096-(ijdm%4096)
    fld2=ma.array([],fill_value=1e30)

    inFile=inFile+'.a'

    fid=open(inFile,'rb')
    Index=hycom[fieldName]
 
    for lyr in range(0,len(Index),1):
        fid.seek((hycom[fieldName][lyr-1]-1)*4*(npad+ijdm),0)
        fld=fid.read(ijdm*4)
        fld=struct.unpack('>'+str(ijdm)+'f',fld)
        fld=np.array(fld)

        if pntidx is not None:
            fld=fld[pntidx]
            if fld2.size == 0:
                fld2=fld.copy()
            else:
                fld2=np.vstack((fld2,fld))
        else:
            fld=ma.reshape(fld,(hycom['jdm'],hycom['idm']))
            if fld2.size == 0:
                fld2=fld.copy()
            else:
                fld2=ma.dstack((fld2,fld))

    fid.close()
    fld2=ma.masked_greater(fld2,1e10)
    return fld2

def parse_l(inFile,fileType='archive'):

    vlist={}

    inFile=inFile+'.b'
    lines=[line.rstrip() for line in open(inFile)]

    idm=int([line.split() for line in lines if 'longitudinal' in line][0][0])
    jdm=int([line.split() for line in lines if 'latitudinal' in line][0][0])
    vlist['idm']=idm
    vlist['jdm']=jdm
   
    count=0
    for line in lines:
        count+=1
        if (fileType=='archive' and line[0:5] == 'field'): 
           break

    lines=lines[count:]
    vars=[line.split()[0] for line in lines]

    count=0
    for var in vars:
        count += 1
        if vlist.__contains__(var):
            vlist[var]=vlist[var] +[count]
        else:
            vlist[var]=[count]

    #print vlist

    return vlist
def parse_z(inFile,nlevels,fileType='3z'):

    vlist={}

    inFile=inFile+'.b'
    lines=[line.rstrip() for line in open(inFile)]

    idm=int([line.split() for line in lines if 'longitudinal' in line][0][0])
    jdm=int([line.split() for line in lines if 'latitudinal' in line][0][0])
    vlist['idm']=idm
    vlist['jdm']=jdm

    count=0
    for line in lines:
        count+=1
        if (fileType=='3z' and line[0:5] == 'field'):
           break

    lines=lines[count:]
    vars=[line.split()[2] for line in lines]
   
    count=0
    for var in vars:
       count += 1
       if vlist.__contains__(var):
          vlist[var]=vlist[var] +[count]
       else:
          vlist[var]=[count]

    return vlist

def readBinz(zbinFile,fileType,fieldName,pntidx=None):

    myvar=parse_z(zbinFile,fileType)

    #print("myvar=",myvar)

    ijdm = myvar['idm'] * myvar['jdm']
    npad = 4096 - (ijdm%4096)
    fld2 = ma.array([],fill_value=1.2676506002282294e+30)

    aFile = zbinFile+'.a'
    fid   = open(aFile,'rb')
    
    #print myvar

    Indx=myvar[fieldName]

    #print 'index and fieldName=',Indx,fieldName

    for knt in range(0,len(Indx),1):
        #print 'readBin: lyr=',lyr
        lyr = Indx[knt]
        fid.seek((lyr-1)*4*(npad+ijdm),0)
        fld = fid.read(ijdm*4)
        fld = struct.unpack('>'+str(ijdm)+'f',fld)
        fld = np.array(fld)

        if pntidx is not None:
           fld = fld[pntidx]
           if fld2.size == 0:
                fld2=fld.copy()
           else:
                fld2=np.vstack((fld2,fld))
        else:
           fld = ma.reshape(fld,(myvar['jdm'],myvar['idm']))
           if fld2.size == 0:
                fld2=fld.copy()
           else:
                fld2=ma.dstack((fld2,fld))

    fid.close()
    fld2=ma.masked_greater(fld2,1e10)
    return fld2

def readBin(binFile,fileType,fieldName,pntidx=None):

    myvar=parse_l(binFile,fileType)

    ijdm = myvar['idm'] * myvar['jdm']
    npad = 4096 - (ijdm%4096)
    fld2 = ma.array([],fill_value=1.2676506002282294e+30)

    aFile = binFile+'.a'
    fid   = open(aFile,'rb')

    Indx=myvar[fieldName]
   
    for knt in range(0,len(Indx),1):
       #print 'readBin: lyr=',lyr
       lyr = Indx[knt]
       fid.seek((lyr-1)*4*(npad+ijdm),0)
       fld = fid.read(ijdm*4)
       fld = struct.unpack('>'+str(ijdm)+'f',fld)
       fld = np.array(fld)
       if pntidx is not None:
          fld = fld[pntidx]
          if fld2.size == 0:
            fld2=fld.copy()
          else:
            fld2=np.vstack((fld2,fld))
       else:
          fld = ma.reshape(fld,(myvar['jdm'],myvar['idm']))
          if fld2.size == 0:
            fld2=fld.copy()
          else:
            fld2=ma.dstack((fld2,fld))

    fid.close()
    fld2=ma.masked_greater(fld2,1e10)
    return fld2

def readgrids(rFile,fieldName,layers,pntidx=None):

    lonlat=parse_b(rFile,'grid')

    ijdm = lonlat['idm'] * lonlat['jdm']
    npad = 4096 - (ijdm%4096)
    fld2 = ma.array([],fill_value=1e30)

    aFile = rFile+'.a'
    fid   = open(aFile,'rb')
    
    for lyr in layers:
       fid.seek((lonlat[fieldName][lyr-1]-1)*4*(npad+ijdm),0)
       fld = fid.read(ijdm*4)
       fld = struct.unpack('>'+str(ijdm)+'f',fld)
       fld = np.array(fld)

       if pntidx is not None:
         fld = fld[pntidx]
         if fld2.size == 0:
            fld2=fld.copy()
         else:
            fld2=np.vstack((fld2,fld))
       else:
         fld = ma.reshape(fld,(lonlat['jdm'],lonlat['idm']))
         if fld2.size == 0:
            fld2=fld.copy()
         else:
            fld2=ma.dstack((fld2,fld))

    fid.close()
    fld2=ma.masked_greater(fld2,1e10)
    return fld2

def readdepth(dFile,pntidx=None):

    fieldName='depth'
    lonlat=parse_b(dFile,fieldName)

    ijdm = lonlat['idm'] * lonlat['jdm']
    npad = 4096 - (ijdm%4096)
    fld2 = ma.array([],fill_value=1e30)

    aFile = dFile+'.a'
    fid   = open(aFile,'rb')

    lyr=0
    fid.seek((lonlat[fieldName][lyr-1]-1)*4*(npad+ijdm),0)
    fld = fid.read(ijdm*4)
    fld = struct.unpack('>'+str(ijdm)+'f',fld)
    fld = np.array(fld)

    if pntidx is not None:
       fld = fld[pntidx]
       if fld2.size == 0:
          fld2=fld.copy()
       else:
          fld2=np.vstack((fld2,fld))
    else:
       fld = ma.reshape(fld,(lonlat['jdm'],lonlat['idm']))
       if fld2.size == 0:
          fld2=fld.copy()
       else:
          fld2=ma.dstack((fld2,fld))

    fid.close()
    fld2=ma.masked_greater(fld2,1e10)
    return fld2

def HYCOMday2normal(hytime):
    """
    --- converts HYCOM model days to gregorian date ---
    """
    tzero=datetime(1901,1,1,0,0)
    myday=tzero+timedelta(hytime-1)
    # myday.ctime()

    return myday

def get_hycomtime(bfile,fileType):
    """
    --- read b file and return a time series in gregorian
    
    Caveat: it applies to only forcing.*.b and archv.*.b.
    """

    bfile=bfile+'.b'
    lines=[line.rstrip() for line in open(bfile)]

    count=0
    for line in lines:
        count+=1
        if (fileType=='archive' and line[0:5] == 'field') or \
            (fileType == 'forcing' and line[0:5]=='i/jdm'): 
            break

    lines=lines[count:]

    htimes=[HYCOMday2normal(float(line.split()[3])) for line in lines]
    
    return(htimes)
    
def datetime2matlabdn(dt):
    mdn = dt + timedelta(days = 366) 
    frac_seconds = (dt-datetime(dt.year,dt.month,dt.day,0,0,0)).seconds / (24.0 * 60.0 * 60.0)
    frac_microseconds = dt.microsecond / (24.0 * 60.0 * 60.0 * 1000000.0)
    return mdn.toordinal() + frac_seconds + frac_microseconds

def find_ijs(in2d, target):
    for row, i in (in2d):
       try:	
           column = i.index(target)
       except ValueError:
           continue
       return row, column
    return -1

def arakawa (u,v):
    u=0.5*(u+np.roll(u,-1,axis=1))
    v=0.5*(v+np.roll(v,-1,axis=0))
    return u,v

def fixed3zlayers():
   """ returns a set of depths for z-levels for HYCOM 3z.[ab]
   """

   z=[1,3,5,7.5,10,15,20,25,30,40,50,60,70,80,100,120,140,160,180,200,220,240,250,275,300,350]
   return([-1*a for a in z])

def EkmanPumping(lns,lts,txs,tys):
# Compute Ekman Pumping
#   input: Wind stress on (Ulon,Ulat) 
#          Grids for a scalar variable (lns,lts)
#
#   output: vertical velocity induced by wind stress [m/s]
#
#   history
#   -------
#       Oct 2019, hsk
#

   #from geo4HYCOM import arakawa
   #import numpy as np

   #dx=haversine(hlon[0,1:],hlat[0,1:],hlon[0,:-1],hlat[0,:-1])
   dx=np.gradient(lns,axis=1)*60*1.852*1000 	# [m]
   dy=np.gradient(lts,axis=0)*60*1.852*1000

   tx,ty=arakawa(txs,tys)
   dtx=np.gradient(tx,axis=1)
   dty=np.gradient(ty,axis=0)

   sx=dty/dx
   sy=dtx/dy

   we=(sx-sy)/1024

   fc=2*7.295E-5*np.cos(np.mean(lts)*np.pi/180)

   return(we/fc)


def readrestart(basin,inFile,fieldName,pntindx=None):
   """ reads a restart.[a/b] and extract 2D field 
       for a var of interest(fieldName).
   
       note: (u,v,temp,saln,dp) are of 2 time steps.
       When one of these inquired, the values will come from the 2nd time-step.
       -----------------------------------------------------------------------
   by HSK July 2020 """

   hycom=parse_r(inFile)
   out=hycombasin(basin)
   idm=out[0]
   jdm=out[1]
   kdm=out[2]

   #dkm=41
   ijdm=idm*jdm
   npad=4096-(ijdm%4096)
   fld2=ma.array([],fill_value=1e30)

   inFile=inFile+'.a'
   print('processing',fieldName,' from ',inFile)
   
   # delete values in (u,v,temp,saln,dp) for the 1st time-step
   del hycom['u'][:kdm]
   del hycom['v'][:kdm]
   del hycom['dp'][:kdm]
   del hycom['temp'][:kdm]
   del hycom['saln'][:kdm]

   fid=open(inFile,'rb')
   for k,indx in enumerate(hycom[fieldName]):  # for (u,v,dp,temp,saln)
      fid.seek((indx-1)*4*(npad+ijdm),0)
      fld=fid.read(ijdm*4)
      fld=struct.unpack('>'+str(ijdm)+'f',fld)
      fld=np.array(fld)

      if pntindx is not None:
            fld=fld[pntidx]
            if fld2.size == 0:
                fld2=fld.copy()
            else:
                fld2=np.vstack((fld2,fld))
      else:
            fld=ma.reshape(fld,jdm,idm)
            if fld2.size == 0:
                fld2=fld.copy()
            else:
                fld2=ma.dstack((fld2,fld))
   fid.close()
   fld2=ma.masked_greater(fld2,1e10)
   return fld2

def parse_r(inFile,fileType='restart'):
   """ for reading a restart.[ab] """
   hycom={}
   inFile=inFile+'.b'
   lines=[line.rstrip() for line in open(inFile)]

   # --- remove 2 lines
   count=0
   for line in lines:
     if (line[:9]=="u       :"):
        break
     count+=1

   lines=lines[count:]
   vars=[line.split()[0] for line in lines]
      
   count=0
   for var in vars:
      count += 1
      if hycom.__contains__(var):
          hycom[var]=hycom[var] +[count]
      else:
          hycom[var]=[count]

   return hycom
    
def hycombasin(basin):
   """ extract idm,jdm,kdm,iref,jref for a regional domain dimensions 
	and a lower left corner index (iref,jref) in the global 1/12-degree HYCOM.

        Hyun-Sook Kim July 15, 2020."""

   Mybasin={'name': ['idm', 'jdm', 'kdm', 'iref', 'jref'],
            'hat10': [1135,633,41,2346,1518],
            'hep20': [1226,476,41, 1318,1568],
            'hwp30': [ 989,480,41,  311,1568],
            'hin40': [ 712,331,41, 4113,1568],
            'hsn50': [1215,567,41, 3968, 909],
            'hsp60': [1313,566,41,   11, 909],
            'hcp70': [ 876,476,41, 1074,1568]}

   return(Mybasin[basin])


