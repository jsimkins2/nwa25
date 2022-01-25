# this code has been adapted from code created by Andrew Ross - https://github.com/andrew-c-ross/nwa-shared

import numpy as np
import xarray
import xesmf
import sys
from HCtFlood import kara as flood


def vgrid_to_interfaces(vgrid, max_depth=6500.0):
    """Convert layer thicknesses to interface depths.
    Args:
        vgrid: array of layer thicknesses.
        max_depth: maximum depth of the model. The lowest interface depth will be set to this.
    Returns:
        Array of interface depths.     
    """
    if isinstance(vgrid, xarray.DataArray):
        vgrid = vgrid.data
    zi = np.concatenate([[0], np.cumsum(vgrid)])
    zi[-1] = max_depth
    return zi


def vgrid_to_layers(vgrid, max_depth=6500.0):
    """Convert layer thicknesses to depths of layer midpoints.
    Args:
        vgrid: array of layer thicknesses.
        max_depth: maximum depth of the model. The lowest interface depth will be set to this.
    Returns:
        Array of layer depths.     
    """
    if isinstance(vgrid, xarray.DataArray):
        vgrid = vgrid.data
    ints = vgrid_to_interfaces(vgrid, max_depth=max_depth)
    z = (ints + np.roll(ints, shift=1)) / 2
    layers = z[1:]
    return layers


def rotate_uv(u, v, angle, in_degrees=True):
    """Rotate velocities from earth-relative to model-relative.
    Args:
        u: west-east component of velocity.
        v: south-north component of velocity.
        angle: angle of rotation from true north to model north.
        in_degrees (bool): typically angle is in radians, but set this to True if it is in degrees.
    Returns:
        Model-relative west-east and south-north components of velocity.
    """
    if in_degrees:
        angle = np.radians(angle)
    urot = np.cos(angle) * u + np.sin(angle) * v
    vrot = -np.sin(angle) * u + np.cos(angle) * v
    return urot, vrot


def interpolate_flood_tracers(ds, target_grid):
    """Interpolate and flood data at tracer points (temperature, salinity, free surface).
    Args:
        ds (xarray.Dataset): Dataset with variables temp, salt, and ssh.
        target_grid (xarray.Dataset): Model supergrid with variables x, y and coords nxp, nyp.
    Returns:
        xarray.Dataset: Dataset flooded and interpolated to MOM tracer grid. 
    """
    # Flood temperature and salinity over land.
    flooded = xarray.merge((
        flood.flood_kara(ds[v], zdim='zl') for v in ['temp', 'salt']
    ))
    
    # Flood ssh separately to avoid extra z=0
    flooded['ssh'] = flood.flood_kara(ds['ssh']).isel(z=0).drop('z')
    
    # Interpolate
    target_points = (
        target_grid
        [['x', 'y']]
        .isel(nxp=slice(1, None, 2), nyp=slice(1, None, 2))
        .rename({'y': 'lat', 'x': 'lon', 'nxp': 'xh', 'nyp': 'yh'})
    )
    soda_to_mom = xesmf.Regridder(
        flooded, 
        target_points, 
        method='bilinear', 
        filename='regrid_soda_tracers.nc',
        reuse_weights=True,
        periodic=True
    )
    interped = soda_to_mom(flooded).drop(['lon', 'lat'])
    return interped


def interpolate_flood_velocity(ds, target_grid):
    """Interpolate and flood velocity data.
    Args:
        ds (xarray.Dataset): Dataset with variables u and v.
        target_grid (xarray.Dataset): Model supergrid with variables x, y and coords nxp, nyp.
    Returns:
        xarray.Dataset: Dataset flooded and interpolated to MOM velocity grid. 
    """
    # Flood over land.
    flooded = xarray.merge((
        flood.flood_kara(ds[v], zdim='zl') for v in ['u', 'v']
    ))

    # Interpolate u and v onto supergrid to make rotation possible
    target_uv = (
        target_grid
        [['x', 'y']]
        .rename({'y': 'lat', 'x': 'lon'})
    )
    soda_to_uv = xesmf.Regridder(
        ds, target_uv, 
        filename='regrid_soda_uv.nc',
        method='nearest_s2d',
        reuse_weights=True,
        periodic=True
    )
    interped_uv = soda_to_uv(flooded[['u', 'v']]).drop(['lon', 'lat'])
    urot, vrot = rotate_uv(interped_uv['u'], interped_uv['v'], target_grid['angle_dx'])
    # Subset onto u and v points.
    uo = urot.isel(nxp=slice(0, None, 2), nyp=slice(1, None, 2)).rename({'nxp': 'xq', 'nyp': 'yh'})
    uo.name = 'u'
    vo = vrot.isel(nxp=slice(1, None, 2), nyp=slice(0, None, 2)).rename({'nxp': 'xh', 'nyp': 'yq'})
    vo.name = 'v'
    
    interped = (
        xarray.merge((uo, vo))
        .transpose('time', 'zl', 'yh', 'yq', 'xh', 'xq')
    )

    return interped


def write_initial(soda_file, vgrid_file, grid_file, start_date, output_file):
    """Interpolate initial conditions for MOM from a SODA file and write to a new file.
    Args:
        soda_file (str): Path to SODA file to use for initial conditions.
        vgrid_file (str): Path to vertical grid to interpolate data to.
        grid_file (str): Path to horizontal grid file (ocean_hgrid.nc) to interpolate data to.
        start_date (np.datetime64): Overwrite the SODA datetime with this datetime. Useful if model start date and SODA 5-day dates do not match.
        output_file (str): Write resulting initial conditions to this file.
    """
    vgrid = xarray.open_dataarray(vgrid_file)
    z = vgrid_to_layers(vgrid)
    ztarget = xarray.DataArray(
        z,
        name='zl',
        dims=['zl'], 
        coords={'zl': z}
    )

    soda = (
        xarray.open_dataset(soda_file)
        .rename({'st_ocean': 'z'})
        [['temp', 'salt', 'ssh', 'u', 'v']]
    )

    # Interpolate SODA vertically onto target grid.
    # Depths below bottom of SODA are filled by extrapolating the deepest available value.
    revert = soda.interp(z=ztarget, kwargs={'fill_value': 'extrapolate'}).ffill('zl', limit=None)

    # Split SODA into data on tracer and velocity points
    tracers = revert[['temp', 'salt', 'ssh']].rename({'xt_ocean': 'lon', 'yt_ocean': 'lat'})
    velocity = revert[['u', 'v']].rename({'xu_ocean': 'lon', 'yu_ocean': 'lat'})

    # Horizontally interpolated the vertically interpolated
    # and flooded data onto the MOM grid.
    grid = xarray.open_dataset(grid_file)

    interped = xarray.merge((
        interpolate_flood_tracers(tracers, grid),
        interpolate_flood_velocity(velocity, grid)
    ))

    # Overwrite the SODA file time with the intended model start date.
    interped['time'] = (('time', ), [start_date])

    # Fix output metadata, including removing all _FillValues.
    all_vars = list(interped.data_vars.keys()) + list(interped.coords.keys())
    encodings = {v: {'_FillValue': None} for v in all_vars}
    encodings['time'].update({'dtype':'float64', 'calendar': 'gregorian'})
    interped['zl'].attrs = {
        'units': 'meter',
        'cartesian_axis': 'Z',
        'positive': 'down'
    }

    interped.to_netcdf(
        output_file,
        format='NETCDF3_64BIT',
        engine='netcdf4',
        encoding=encodings,
        unlimited_dims='time'
    )


def main():
    # Use SODA data centered on 1992-12-30.
    # Model start date is 1993-01-01.
    # https://dsrs.atmos.umd.edu/DATA/soda3.12.2/REGRIDED/ocean/soda3.12.2_5dy_ocean_reg_1993_01_04.nc
    soda_file = '/home/james/SODA/5day/soda3.3.1_5dy_ocean_reg_2010_01_05.nc'
    start_date = np.datetime64('2010-01-05T00:00:00')
    # Used in filename below, don't change
    start_str = np.datetime_as_string(start_date, unit='D')
    
    # Save the ICs here:
    output_file = f'/home/james/initCond/nwa25/soda_ic_75z_{start_str}.nc'

    # Model vertical grid:
    vgrid_file = '/home/james/gridInfo/nwa25/vgrid_75_2m.nc'

    # Model horizontal grid:
    grid_file = '/home/james/gridInfo/nwa25/ocean_hgrid.nc'

    write_initial(soda_file, vgrid_file, grid_file, start_date, output_file)


if __name__ == '__main__':
    main()
