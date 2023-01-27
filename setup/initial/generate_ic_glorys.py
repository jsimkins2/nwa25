import numpy as np
import xarray
import xesmf

# https://github.com/raphaeldussin/HCtFlood
import sys
sys.path.append('/home/Andrew.C.Ross/git/HCtFlood')
from HCtFlood import kara as flood
sys.path.append('../boundary')
from boundary import rotate_uv

import depths


def write_initial(glorys_file, vgrid_file, grid_file, output_file):
    vgrid = xarray.open_dataarray(vgrid_file)
    z = depths.vgrid_to_layers(vgrid)
    ztarget = xarray.DataArray(
        z,
        name='zl',
        dims=['zl'], 
        coords={'zl': z}, 
    )
    glorys = (
        xarray.open_dataset(glorys_file)
        [['thetao', 'so', 'zos', 'uo', 'vo']]
        .rename({'longitude': 'lon', 'latitude': 'lat'})
    )
    # Round time down to midnight
    #glorys.time.data[0] = numpy.datetime64(f'{year}-01-01T00:00:00.000000000')
    #from IPython import embed; embed()
    year = str(glorys.time.data[0])[0:4]
    glorys.time.data[0] = np.datetime64(f'{year}-01-01T00:00:00.000000000')

    # Interpolate GLORYS vertically onto target grid.
    # Depths below bottom of GLORYS are filled by extrapolating the deepest available value.
    revert = glorys.interp(depth=ztarget, kwargs={'fill_value': 'extrapolate'}).ffill('zl', limit=None)

    # Flood temperature and salinity over land. 
    flooded = xarray.merge((
        flood.flood_kara(revert[v], zdim='zl') for v in ['thetao', 'so', 'uo', 'vo']
    ))

    # flood zos separately to avoid the extra z=0 added by flood_kara.
    flooded['zos'] = flood.flood_kara(revert['zos']).isel(z=0).drop('z')

    # Horizontally interpolate the vertically interpolated and flooded data onto the MOM grid. 
    target_grid = xarray.open_dataset(grid_file)
    target_grid['x'] -= 360.0
    target_t = (
        target_grid
        [['x', 'y']]
        .isel(nxp=slice(1, None, 2), nyp=slice(1, None, 2))
        .rename({'y': 'lat', 'x': 'lon', 'nxp': 'xh', 'nyp': 'yh'})
    )
    # Interpolate u and v onto supergrid to make rotation possible
    target_uv = (
        target_grid
        [['x', 'y']]
        .rename({'y': 'lat', 'x': 'lon'})
    )
    
    regrid_kws = dict(method='nearest_s2d', reuse_weights=True, periodic=False)

    glorys_to_t = xesmf.Regridder(glorys, target_t, filename='regrid_glorys_tracers.nc', **regrid_kws)
    glorys_to_uv = xesmf.Regridder(glorys, target_uv, filename='regrid_glorys_uv.nc', **regrid_kws)

    interped_t = glorys_to_t(flooded[['thetao', 'so', 'zos']]).drop(['lon', 'lat'])

    # Interpolate u and v, rotate, then extract individual u and v points
    interped_uv = glorys_to_uv(flooded[['uo', 'vo']]).drop(['lon', 'lat'])
    urot, vrot = rotate_uv(interped_uv['uo'], interped_uv['vo'], target_grid['angle_dx'])
    uo = urot.isel(nxp=slice(0, None, 2), nyp=slice(1, None, 2)).rename({'nxp': 'xq', 'nyp': 'yh'})
    uo.name = 'uo'
    vo = vrot.isel(nxp=slice(1, None, 2), nyp=slice(0, None, 2)).rename({'nxp': 'xh', 'nyp': 'yq'})
    vo.name = 'vo'
    
    interped = (
        xarray.merge((interped_t, uo, vo))
        .transpose('time', 'zl', 'yh', 'yq', 'xh', 'xq')
    )

    # Rename to match MOM expectations.
    interped = interped.rename({
        'thetao': 'temp',
        'so': 'salt',
        'zos': 'ssh',
        'uo': 'u',
        'vo': 'v'
    })

    # Fix output metadata, including removing all _FillValues.
    all_vars = list(interped.data_vars.keys()) + list(interped.coords.keys())
    encodings = {v: {'_FillValue': None} for v in all_vars}
    encodings['time'].update({'dtype':'float64', 'calendar': 'gregorian'})
    interped['zl'].attrs = {
        'long_name': 'Layer pseudo-depth, -z*',
         'units': 'meter',
         'cartesian_axis': 'Z',
         'positive': 'down'
    }

    interped.to_netcdf(
        output_file,
        format='NETCDF4',
        engine='netcdf4',
        encoding=encodings,
        unlimited_dims='time'
    )


def main():
    glorys_file = '/Volumes/A1/workdir/james/glorys/1996/GLORYS_REANALYSIS_1996-01-01.nc'
    vgrid_file = '/Volumes/A1/workdir/james/gt_nwa25_input/gridInfo/vgrid_75_2m.nc'
    grid_file = '/Volumes/A1/workdir/james/gt_nwa25_input/gridInfo/ocean_hgrid.nc'
    output_file = '/Volumes/A1/workdir/james/gt_nwa25_input/initial/glorys_ic_75z.nc'
    write_initial(glorys_file, vgrid_file, grid_file, output_file)


if __name__ == '__main__':
    main()
