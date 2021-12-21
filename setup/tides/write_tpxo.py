import datetime as dt
import numpy as np
from os import path
import pandas as pd
import xarray
from boundary import Segment

# xarray gives a lot of unnecessary warnings
import warnings
warnings.filterwarnings('ignore')


def write_tpxo(constituents, tpxo_dir, segments, horizontal_subset):
    tpxo_h = (
        xarray.open_dataset(path.join(tpxo_dir, 'h_tpxo9.v5a.nc'))
        .rename({'lon_z': 'lon', 'lat_z': 'lat', 'nc': 'constituent'})
        .isel(constituent=constituents, **horizontal_subset)
    )
    h = tpxo_h['ha'] * np.exp(-1j * np.radians(tpxo_h['hp']))
    tpxo_h['hRe'] = np.real(h)
    tpxo_h['hIm'] = np.imag(h)
    tpxo_u = (
        xarray.open_dataset(path.join(tpxo_dir, 'u_tpxo9.v5a.nc'))
        .rename({'lon_u': 'lon', 'lat_u': 'lat', 'nc': 'constituent'})
        .isel(constituent=constituents, **horizontal_subset)
    )
    tpxo_u['ua'] *= 0.01  # convert to m/s
    u = tpxo_u['ua'] * np.exp(-1j * np.radians(tpxo_u['up']))
    tpxo_u['uRe'] = np.real(u)
    tpxo_u['uIm'] = np.imag(u)
    tpxo_v = (
        xarray.open_dataset(path.join(tpxo_dir, 'u_tpxo9.v5a.nc'))
        .rename({'lon_v': 'lon', 'lat_v': 'lat', 'nc': 'constituent'})
        .isel(constituent=constituents, **horizontal_subset)
    )
    tpxo_v['va'] *= 0.01  # convert to m/s
    v = tpxo_v['va'] * np.exp(-1j * np.radians(tpxo_v['vp']))
    tpxo_v['vRe'] = np.real(v)
    tpxo_v['vIm'] = np.imag(v)
    # Tidal amplitudes are currently constant over time.
    # Seem to need a time dimension to have it read by MOM.
    # But also, this would later allow nodal modulation
    # or other long-term variations to be added.
    # the date should begin 1 month before first day of simulation
    times = xarray.DataArray(
        pd.date_range('2009-12-01', periods=1),
        dims=['time']
    )
    for seg in segments:
        print(seg)
        seg.regrid_tidal_elevation(
            tpxo_h[['lon', 'lat', 'hRe']],
            tpxo_h[['lon', 'lat', 'hIm']],
            times,
            flood=True
        )
        seg.regrid_tidal_velocity(
            tpxo_u[['lon', 'lat', 'uRe']],
            tpxo_u[['lon', 'lat', 'uIm']],
            tpxo_v[['lon', 'lat', 'vRe']],
            tpxo_v[['lon', 'lat', 'vIm']],
            times,
            flood=True
        )


def main():
    """
    constituents in TPXO9:
    con =
    "m2  ",  0
    "s2  ",  1
    "n2  ",  2
    "k2  ",  3
    "k1  ",  4
    "o1  ",  5
    "p1  ",  6
    "q1  ",  7
    "mm  ",  8
    "mf  ",  9
    "m4  ", 10
    "mn4 ", 11
    "ms4 ", 12
    "2n2 ", 13
    "s1  "  14
    """
    # Use these constituents:
    # (MOM6 does not currently support m4, mn4, ms4, 2n2, or s1.)
    constituents = range(0, 10)

    # Subset TPXO9 to a region roughly around the NWA25 domain
    # for computational efficiency.
    horizontal_subset = dict(ny=slice(500,1100), nx=slice(1500,None))

    # Path to tpxo9.v1 data
    tpxo_dir = '/home/james/tpxo9v5a/DATA/'

    # Save output here
    output_dir = '/home/james/tpxo9v5a/DATA/'

    # Setup NWA boundaries
    hgrid = xarray.open_dataset('/home/james/gridInfo/nwa25/ocean_hgrid.nc')
    # segments taken from nwa25_obc.py
    segments = [
        Segment(1, 'south', hgrid, output_dir=output_dir),
        Segment(2, 'north', hgrid, output_dir=output_dir),
        Segment(3, 'east', hgrid, output_dir=output_dir)
    ]

    write_tpxo(constituents, tpxo_dir, segments, horizontal_subset)


    
if __name__ == '__main__':
    main()
