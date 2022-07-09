#
""" NOTES:
      - requires Python 2.4 or greater
      - elements of the lists must be hashable
      - order of the original lists is not preserved
"""

import numpy as np
import os
from scipy.ndimage import map_coordinates
import csv

def unique(a):
    """ return the list with duplicate elements removed """
    return list(set(a))

def intersect(a, b):
    """ return the intersection of two lists """
    return list(set(a) & set(b))

def union(a, b):
    """ return the union of two lists """
    return list(set(a) | set(b))

def interp3d(x, y, z, v, xi, yi, zi, **kwargs):
    """Sample a 3D array "v" with pixel corner locations at "x","y","z" at the
    points in "xi", "yi", "zi" using linear interpolation. Additional kwargs
    are passed on to ``scipy.ndimage.map_coordinates``."""

    def index_coords(corner_locs, interp_locs):
        index = np.arange(len(corner_locs))
        if np.all(np.diff(corner_locs) < 0):
            corner_locs, index = corner_locs[::-1], index[::-1]
        return np.interp(interp_locs, corner_locs, index)

    orig_shape = np.asarray(xi).shape
    xi, yi, zi = np.atleast_1d(xi, yi, zi)
    for arr in [xi, yi, zi]:
        arr.shape = -1

    output = np.empty(xi.shape, dtype=float)
    coords = [index_coords(*item) for item in zip([x, y, z], [xi, yi, zi])]

    map_coordinates(v, coords, order=1, output=output, **kwargs)

    return output.reshape(orig_shape)

def closest(points,x):
    return sorted(points, key=lambda p: abs(p[0] - x)) [:2]


def NDBCstations():
    tfile='/work/noaa/hwrf/save/hskim/PostPy/parm/NDBCstations.txt'
    infile = open(tfile,'r')
    lines=infile.readlines()[7:]
    infile.close()

    nlat=[]
    for line in lines[:-1]:
        s=float(line[29:31])+float(line[33:35])/60.+float(line[37:39])/3600.
        nlat.append(s)
        del (s)
	    
    nlon=[]
    for line in lines[:-1]:
        m=(-1.)*(float(line[46:49])+float(line[51:53])/60.+float(line[55:58])/3600.)
        nlon.append(m)
        del(m)

    return nlon, nlat

def _mkdirs(path):
    if not os.path.exists(path):
      l=[]
      p = "/"
      l = path.split("/")
      i = 1
      while i < len(l):
        p = p + l[i] + "/"
        i = i + 1
        if not os.path.exists(p):
            os.makedirs(p,mode=0o755, exist_ok=False)

def _mkdir_recursive(self, path):
    sub_path = os.path.dirname(path)
    if not os.path.exists(sub_path):
        self._mkdir_recursive(sub_path)
    if not os.path.exists(path):
        os.makedirs(path)

def coast180():
   inf=open('/work/noaa/hwrf/save/hskim/PostPy/fixdata/coastlines_180x180.dat','r')
   c1 = []
   c2 = []
   hsk=np.genfromtxt(inf)
   c1=hsk[:,0]
   c2=hsk[:,1]
   return (c1,c2)

def coast360():
   inf=open('/work/noaa/hwrf/save/hskim/PostPy/fixdata/coastlines_180x180.dat','r')
   c1 = []
   c2 = []

   hsk=np.genfromtxt(inf)
   c1=hsk[:,0]+360.
   c2=hsk[:,1]
   return (c1,c2)

def medfilt (x, k):
   """Apply a length-k median filter to a 1D array x.
   Boundaries are extended by repeating endpoints.
   """
   assert k % 2 == 1, "Median filter length must be odd."
   assert x.ndim == 1, "Input must be one-dimensional."
   k2 = (k - 1) // 2
   y = np.zeros ((len (x), k), dtype=x.dtype)
   y[:,k2] = x
   for i in range (k2):
       j = k2 - i
       y[j:,i] = x[:-j]
       y[:j,i] = x[0]
       y[:-j,-(i+1)] = x[j:]
       y[-j:,-(i+1)] = x[-1]
   return np.median (y, axis=1)

def rmse(predictions,targets):
    return np.sqrt(((predictions - targets)**2).mean())

def r_NDBC(infile,var,PER_HOUR=True):
   """ read in xxxxh201[x].txt of a NDBC buoy ascii file,
       and returns "var" in the entire time periods of the data contains

   header:
    0   1  2  3  4  5    6    7    8      9     10  11    12    13    14    15    16   17
   ----------------------------------------------------------------------------------------
   #YY  MM DD hh mm WDIR WSPD GST  WVHT   DPD   APD MWD   PRES  ATMP  WTMP  DEWP  VIS  TIDE
   #yr  mo dy hr mn degT m/s  m/s     m   sec   sec degT   hPa  degC  degC  degC   mi    ft

       -hsk 9/2019
   """
   from datetime import datetime, timedelta
   #import datetime
   import numpy as np
   from datetime_misc import datetime2matlabdn

   dic={'YY':0,'MM':1,'DD':2,'hh':3,'mm':4,'WDIR':5,'WSPD':6,'GST':7,'WVHT':8,'DPD':9,'APD':10,
       'MWD':11,'PRES':12,'ATMP':13,'WTMP':14,'DEWP':15,'VIS':16,'TIDE':17}

   #indx=dic[var.upper()]
   #print(indx)

   b=np.genfromtxt(infile)
   dns=[]
   dts=[]
   for  k,l in enumerate(b):
      dtmp=datetime(int(l[0]),int(l[1]),int(l[2]),int(l[3]),int(l[4]),0)
      ntmp=datetime2matlabdn(dtmp)
      if PER_HOUR:
         if ntmp%1 >= 0.5:
            dtmp=dtmp-timedelta(int(l[4])/24/60)+timedelta(1/24)
         else:
            dtmp=dtmp-timedelta(int(l[4])/24/60)

      dns.append(datetime2matlabdn(dtmp))
      dts.append(dtmp)

   ovar=b[:,dic[var.upper()]]
   ovar[ovar==999.0]=np.nan
   ovar[ovar==99.0]=np.nan
   return(dns,dts,ovar)

def xy4NDBC(IDofI):   
   """ read in list.ids and returns a set of (lon,lat) position 
       for the station(id).

       -hsk 9/2019
   """

   print('xy4NDBC: ... find (lon,lat) for '+IDofI)
   intxt='//work/noaa/hwrf/save/hskim/PostPy/parm/NDBCstations.txt'
   with open(intxt) as f:
     for k,line in enumerate(f):
        if k>=6:
           each=line.split()
           lat=float(each[2])+float(each[3])/60.0+float(each[4])/3600.0
           lon=float(each[6])+float(each[7])/60.0+float(each[8])/3600.0
           if each[0]==IDofI.upper():
              return(each[0],-1*np.round(lon*1000)/1000,np.round(lat*1000)/1000)

def rotate_vec(x,y,angle):
   """ angle [deg] """

   import numpy as np

   d2rad=angle*np.pi/180

   xo=x*np.cos(d2rad)+y*np.sin(d2rad)
   yo=(-1)*x.np.sin(d2rad)+y*np.cos(d2rad)

   return (xo,yo)

#---- new color for Spectral_r
def hsk_Spectral_r():
   from matplotlib import cm
   from matplotlib.colors import ListedColormap, LinearSegmentedColormap
   spectral_r = cm.get_cmap('Spectral_r', 256)
   newcolors = spectral_r(np.linspace(0, 1, 256))
   pink = np.array([248/256, 24/256, 148/256, 1])
   darkred1=np.array([150/256,0,150/256,1])
   darkred2=np.array([100/256,0,100/256,1])
   newcolors[-15:, :] = darkred2
   newcolors[-25:-15,:]=darkred1
   return(ListedColormap(newcolors))

def radial_profile(data, center):
   y, x = np.indices((data.shape))
   r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
   r = r.astype(np.int)

   tbin = np.bincount(r.ravel(), data.ravel())
   nr = np.bincount(r.ravel())
   radialprofile = tbin / nr
   return radialprofile

