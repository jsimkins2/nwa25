# this code has been adapted from code created by Andrew Ross - https://github.com/andrew-c-ross/nwa-shared
from glob import glob
from subprocess import run
import os 
import xarray
from boundary import Segment

# xarray gives a lot of unnecessary warnings
import warnings
warnings.filterwarnings('ignore')


def write_year(year, glorys_dir, segments, is_first_year=False):
    glorys = (
        xarray.open_mfdataset(sorted(glob(os.path.join(glorys_dir, f'{year}/GLORYS_REANALYSIS_{year}-*.nc'))))
        .rename({'latitude': 'lat', 'longitude': 'lon', 'depth': 'z'})
    )

    # Floor first time down to midnight so that it matches initial conditions
    if is_first_year:
       tnew = xarray.concat((glorys['time'][0].dt.floor('1d'), glorys['time'][1:]), dim='time')
       glorys['time'] = ('time', tnew)

    for seg in segments:
        seg.regrid_velocity(glorys['uo'], glorys['vo'], suffix=year, flood=False)
        for tr in ['thetao', 'so']:
            seg.regrid_tracer(glorys[tr], suffix=year, flood=False)
        seg.regrid_tracer(glorys['zos'], suffix=year, flood=False)
        

# this is an xarray based way to concatenate the obc yearly files into one file (per variable of output)
# the former way to do this was based on ncrcat from NCO tools
def ncrcat_rename(nsegments, ncrcat_outdir,output_dir, delete_old_files=False):
    rename_dict={'thetao':'temp', 'so':'salt', 'zos':'zeta', 'uv':'uv'}
    for var in ['thetao', 'so', 'zos', 'uv']:
        for seg in range(1, nsegments+1):
            comb = xarray.open_mfdataset(f'{output_dir}{var}_00{seg}_*')
            if var!='uv':
                comb=comb.rename({f'{var}_segment_00{seg}':f'{rename_dict[var]}_segment_00{seg}'})
                if var!='zos':
                    comb=comb.rename({f'dz_{var}_segment_00{seg}':f'dz_{rename_dict[var]}_segment_00{seg}'})
            # Fix output metadata, including removing all _FillValues.
            all_vars = list(comb.data_vars.keys()) + list(comb.coords.keys())
            encodings = {v: {'_FillValue': None} for v in all_vars}
            encodings['time'].update({'dtype':'float64', 'calendar': 'gregorian'})
            comb.to_netcdf(f'{ncrcat_outdir}{rename_dict[var]}_00{seg}.nc',
                           encoding=encodings,
                           unlimited_dims='time',
                           format='NETCDF3_64BIT')
            print(f'concatenated and renamed {rename_dict[var]}_00{seg}.nc')
            del(comb)
            if delete_old_files==True:
                os.remove(f'{output_dir}{var}_00{seg}_*')


def main():
    first_year = 2012
    glorys_dir = '/Volumes/A1/workdir/james/glorys/'
    output_dir = '/Volumes/A1/workdir/james/nwa25_input/boundary/boundary_final/'
    ncrcat_outdir = '/Volumes/A1/workdir/james/nwa25_input/boundary/boundary_final/'
    hgrid = xarray.open_dataset('/home/james/gridInfo/nwa25/ocean_hgrid.nc')
    segments = [
        Segment(1, 'south', hgrid, output_dir=output_dir),
        Segment(2, 'north', hgrid, output_dir=output_dir),
        Segment(3, 'east', hgrid, output_dir=output_dir)
        ]

    for y in range(2012,2013):
        print(y)
        write_year(y, glorys_dir, segments,is_first_year=y==first_year)
    
    ncrcat_rename(len(segments), ncrcat_outdir, output_dir)




if __name__ == '__main__':
    main()
