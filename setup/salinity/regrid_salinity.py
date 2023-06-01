#!/usr/bin/env python
# coding: utf-8

# # Regrid Chassignet & Xu Salinity File to Gridtools NWA25

# In[1]:


import xesmf as xe
import xarray as xr
import numpy as np


# In[2]:


old_ds = xr.open_dataset('sss_mom6.nc', decode_times='False')
old_ds


# In[3]:


old_ds = old_ds.rename({"Longitude": "lon", "Latitude": "lat"})
old_ds


# In[4]:


new_ds = xr.open_dataset('ocean_mask.nc')
new_ds


# In[5]:


new_ds = new_ds.rename({"x": "lon", "y": "lat"})
new_ds


# In[6]:


lon_centers = new_ds.lon.values
lon_corners = 0.5 * (lon_centers[:-1, :-1] + lon_centers[1:, 1:])
nrows, ncols = lon_centers.shape
lon_b = np.zeros((nrows+1, ncols+1))
lon_b[1:-1, 1:-1] = lon_corners

# Estimate the values on the left and right columns
lon_b[1:-1, 0] = 2 * lon_centers[:-1, 0] - lon_b[1:-1, 1]
lon_b[1:-1, -1] = 2 * lon_centers[:-1, -1] - lon_b[1:-1, -2]

# Estimate the values on the top and bottom rows
lon_b[0, 1:-1] = lon_centers[0, :-1]
lon_b[-1, 1:-1] = lon_centers[-1, :-1]

lon_b[0, 0] = 2 * lon_centers[0, 0] - lon_b[1, 1]
lon_b[0, -1] = 2 * lon_centers[0, -1] - lon_b[1, -2]
lon_b[-1, 0] = 2 * lon_centers[-1, 0] - lon_b[-2, 1]
lon_b[-1, -1] = 2 * lon_centers[-1, -1] - lon_b[-2, -2]


# In[7]:


lat_centers = new_ds.lat.values
lat_corners = 0.5 * (lat_centers[:-1, :-1] + lat_centers[1:, 1:])
nrows, ncols = lat_centers.shape
lat_b = np.zeros((nrows+1, ncols+1))
lat_b[1:-1, 1:-1] = lat_corners

# Estimate the values on the left and right columns
lat_b[1:-1, 0] = 2 * lat_centers[:-1, 0] - lat_b[1:-1, 1]
lat_b[1:-1, -1] = 2 * lat_centers[:-1, -1] - lat_b[1:-1, -2]

# Estimate the values on the top and bottom rows
lat_b[0, 1:-1] = lat_centers[0, :-1]
lat_b[-1, 1:-1] = lat_centers[-1, :-1]

lat_b[0, 0] = 2 * lat_centers[0, 0] - lat_b[1, 1]
lat_b[0, -1] = 2 * lat_centers[0, -1] - lat_b[1, -2]
lat_b[-1, 0] = 2 * lat_centers[-1, 0] - lat_b[-2, 1]
lat_b[-1, -1] = 2 * lat_centers[-1, -1] - lat_b[-2, -2]


# In[204]:


nxp = np.append(new_ds['nx'].values, len(new_ds['nx'].values))
nyp = np.append(new_ds['ny'].values, len(new_ds['ny'].values))


# In[18]:


new_ds = xr.open_dataset('ocean_mask.nc')
lat = new_ds.y.data
lon = new_ds.x.data
# add nxp and nyp dimensions for the lat/lon corners to latch onto
new_ds = new_ds.expand_dims({'nyp':(len(new_ds.ny) + 1)})
new_ds = new_ds.expand_dims({'nxp':(len(new_ds.nx) + 1)})
new_ds


# In[19]:


new_ds = new_ds.drop_vars('mask')
new_ds = new_ds.drop_vars('x')
new_ds = new_ds.drop_vars('y')
new_ds['lat'] = (('ny', 'nx'), lat)
new_ds['lon'] = (('ny', 'nx'), lon)
new_ds


# In[20]:


new_ds = new_ds.assign(lon_b=xr.DataArray(lon_b, dims=("nyp", "nxp")))
new_ds = new_ds.assign(lat_b=xr.DataArray(lat_b, dims=("nyp", "nxp")))


# In[21]:


new_ds


# In[ ]:


regridder = xe.Regridder(old_ds, new_ds, "bilinear")
dr_out = regridder(old_ds)

# rename ny nx
dr_out = dr_out.rename({'ny':'Y'},{'nx':'X'})

# assign coordinates for the Y and X dimensions
dr_out = dr_out.assign_coords(Y=b.Y)
dr_out = dr_out.assign_coords(X=b.X)

# assign dattype as double
dr_out.Y = dr_out.Y.astype('double')
dr_out.X = dr_out.X.astype('double')

# write to netcdf
dr_out.to_netcdf('gtnwa25_sss.nc')

# In[ ]:




