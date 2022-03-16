#add the two precipitation schemes together
import numpy as np
import xarray as xr
from numba import vectorize
import argparse


def trr_from_crr_and_lsrr(crr,lsrr):
    """Add together ERA5 convective rain rate and ERA5 large scale rain rate to give us a total rain rate
    Args:
        crr (xarray.DataArray): convective rain rate
        lsrr (xarray.DataAray): large scale rain rate
    """

    trr = crr.drop(args["crrvar"])
    trr['trr'] = crr[args["crrvar"]] + lsrr[args["lsrrvar"]]
    trr['trr'].attrs = {'units': 'kg m-2 s-1','long_name': 'Total Rainfall Rate'}

    return trr


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="aggregate ERA5 convective + ERA5 large-scale rainfall rates"
    )
    parser.add_argument(
        "--crr",
        type=str,
        required=True,
        help="name of convective rainfall rate file",
    )
    parser.add_argument(
        "--lsrr",
        type=str,
        required=True,
        help="name of large-scale rainfall rate file",
    )
    parser.add_argument(
        "--crrvar",
        type=str,
        required=False,
        default='crr',
        help="name of convective rainfall rate variable",
    )
    parser.add_argument(
        "--lsrrvar",
        type=str,
        required=False,
        default='lsrr',
        help="name of large-scale rainfall rate variable",
    )
    parser.add_argument(
        "-o",
        "--fileout",
        type=str,
        required=True,
        help="name of produced total rainfall rate file",
    )
    parser.add_argument(
        "-c",
        "--chunks",
        type=int,
        required=False,
        default=10,
        help="size of chunk in time",
    )

    args = vars(parser.parse_args())

    crr = xr.open_dataset(args["crr"], chunks={'time': args["chunks"]})
    lsrr = xr.open_dataset(args["lsrr"], chunks={'time': args["chunks"]})

    trr = trr_from_crr_and_lsrr(crr,lsrr)
    trr.trr.encoding = {k: v for k, v in crr.crr.encoding.items() if k in {'_FillValue', 'missing_value', 'dtype'}}
    #trr.trr.encoding.update({'add_offset': None, 'scale_factor': None})
    trr.to_netcdf(args["fileout"], format="NETCDF4_CLASSIC")



###  python /glade/scratch/jsimkins/ERA5/scripts/precip_gen.py --crr /glade/scratch/jsimkins/ERA5/original/ERA5_convective_rain_rate_1997.nc --lsrr /glade/scratch/jsimkins/ERA5/original/ERA5_large_scale_rain_rate_1997.nc -o /glade/scratch/jsimkins/ERA5/original/ERA5_total_rain_rate_1997.nc
