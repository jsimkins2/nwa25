
#%% Function that converts geographic coordinates to GOFS3.1 coordinates
def Geo_coord_to_GOFS_coord(long,latg):
    import numpy as np
    long = np.asarray(long)
    if np.ndim(long) > 0:
        lonm = [360 + ln if ln<0 else ln for ln in long]
    else:
        lonm = [360 + long if long<0 else ii][0]
    latm = latg
    return lonm, latm

################################################################################
#%% First, we need a data set. We will use the publicly available data Output
# from global ocean model GOFS3.1
# at (longitude,latitude) = (-70,24) on July 3 2020 at 00 UTC

url_model = 'http://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0/ts3z'
latg = 24
long = -70
date_ini = '2020/07/30/12'

import xarray as xr
import netCDF4
import numpy as np
from datetime import datetime, timedelta
import seawater as sw

from Upper_ocean_metrics import OHC_from_profile, MLD_temp_crit, MLD_dens_crit,\
                                T100, Potential_energy_anomaly100

# Download dataset
model = xr.open_dataset(url_model,decode_times=False)

lat_model = np.asarray(model.lat[:])
lon_model = np.asarray(model.lon[:])
depth_model = np.asarray(model.depth[:])
ttm = model.time
tm = netCDF4.num2date(ttm[:],ttm.units)

# Find the time index
tini = datetime.strptime(date_ini,'%Y/%m/%d/%H')
oktimem = np.where(tm == tini)[0]
time_model = tm[oktimem[0]]

# We need to convert geographic coordinates to GOFS3.1 coordinates
lonG,latG = Geo_coord_to_GOFS_coord(long,latg)
# Find the lat and lon index in the GOFS grid
# these are just coordinate transformations to locate the correct lat/lon grid cell we want to use
oklonm = np.interp(lonG,lon_model,np.arange(len(lon_model))).astype('int')
oklatm = np.interp(latG,lat_model,np.arange(len(lat_model))).astype('int')

# Retrieve temperature and salinity
temp_model = np.asarray(model.variables['water_temp'][oktimem[0],:,oklatm,oklonm])
salt_model = np.asarray(model.variables['salinity'][oktimem[0],:,oklatm,oklonm])
# Get density
dens_model = sw.dens(salt_model,temp_model,depth_model)

################################################################################
#%% Example1: Ocean heat content from the global ocean model GOFS3.1
# at (longitude,latitude) = (-70,24) on July 3 2020 at 00 UTC

OHC = OHC_from_profile(depth_model,temp_model,dens_model)
print('The OHC from GOFS3.1 at (lon,lat) = (-71,24) on July 3 2020 at 00 UTC is ',\
      OHC, ' kJ/m^3 \n')

################################################################################
#%% Example2: Mixed layer depth and Mixed layer temperature based on the
# temperature criteria: T - T10 <= 0.2 (de Boyer Montégut, Clément, et al.
# "Mixed layer depth over the global ocean: An examination of profile data and
# a profile‐based climatology."
#Journal of Geophysical Research: Oceans 109.C12 (2004).),
# using the global ocean model GOFS3.1
# at (longitude,latitude) = (-70,24) on July 3 2020 at 00 UTC

dtemp = 0.2
ref_depth = 10 # meters
MLD_temp_crit, MLT_temp_crit = MLD_temp_crit(dtemp,ref_depth,depth_model,\
                                             temp_model)
print('The mixed layer temperature (based on a temperatere criteria) from GOFS3.1 at (lon,lat) = (-71,24) on July 3 2020 at 00 UTC is ',\
      MLT_temp_crit, ' degrees C \n')

################################################################################
#%% Example3: Mixed layer depth and Mixed layer temperature based on the
# density criteria: rho10 - rho <= 0.125 kg/m^3 (de Boyer Montégut, Clément,
# et al. "Mixed layer depth over the global ocean: An examination of profile
# data and a profile‐based climatology." Journal of Geophysical Research: Oceans
# 109.C12 (2004)),
# using the global ocean model GOFS3.1
# at (longitude,latitude) = (-70,24) on July 3 2020 at 00 UTC

drho = 0.125
ref_depth = 10 # meters
MLD_dens_crit, MLT_dens_crit = MLD_dens_crit(drho,ref_depth,depth_model,\
                                             temp_model,dens_model)
print('The mixed layer temperature (based on a density criteria) from GOFS3.1 at (lon,lat) = (-71,24) on July 3 2020 at 00 UTC is ',\
      MLT_dens_crit, ' degrees C \n')

################################################################################
#%% Example4: Depth average temperature in the top 100 meters
# using the global ocean model GOFS3.1
# at (longitude,latitude) = (-70,24) on July 3 2020 at 00 UTC

temp100 = T100(depth_model,temp_model)
print('The depth average temperature in the top 100 meters from GOFS3.1 at (lon,lat) = (-71,24) on July 3 2020 at 00 UTC is ',\
      temp100, ' degrees C \n')
################################################################################
#%% Example5: Potential energy anomaly (J/m^3) in the top 100 meters of the water column
# using the global ocean model GOFS3.1
# at (longitude,latitude) = (-70,24) on July 3 2020 at 00 UTC

PEA = Potential_energy_anomaly100(depth_model,dens_model)
print('The potential energy anomaly in the top 100 meters from GOFS3.1 at (lon,lat) = (-71,24) on July 3 2020 at 00 UTC is ',\
      PEA, ' J/m^3 \n')

################################################################################
