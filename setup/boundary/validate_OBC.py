#!/usr/bin/env python

import os
import numpy as np
import xarray as xr

obc_dir = '/home/cermak/workdir/configs/nwa25/OBC'
f1 = os.path.join(obc_dir, 'indiv1/uv_001_1996.nc')
f2 = os.path.join(obc_dir, 'indiv.ans/uv_001_1996.nc')

ds1 = xr.open_dataset(f1)
ds2 = xr.open_dataset(f2)

allvars = list(ds1.variables)
coords = list(ds1.coords)
onlyvars = set(allvars) - set(coords)

notEqual = {
        'count': 0,
        'coord': [],
        'variable': [],
}

for i in coords:
    if not(np.all(ds1[i] == ds2[i])):
        notEqual['coord'].append(i)
        notEqual['count'] = notEqual['count'] + 1

for i in onlyvars:
    if not(np.all(ds1[i] == ds2[i])):
        notEqual['variable'].append(i)
        notEqual['count'] = notEqual['count'] + 1

if notEqual['count'] > 0:
    print("Some items do not match:")
    if len(notEqual['coord']) > 0:
        print("  Coordinates:\n    %s" % (", ".join(notEqual['coord'])))
    if len(notEqual['variable']) > 0:
        print("  Variables:\n    %s" % (", ".join(notEqual['variable'])))
else:
    print("NetCDF files are equivalent.")
