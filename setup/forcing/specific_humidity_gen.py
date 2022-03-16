# original code from Raf Dussin : https://github.com/raphaeldussin/newERAtools/blob/master/humidity.py
import numpy as np
import xarray as xr
from numba import vectorize
import argparse


def huss_from_dtas_and_sp(dtas, sp):
    """Compute specific humidity from dewpoint and sea level pressure
    Args:
        dtas (xarray.DataArray): dewpoint at 2m
        sp (xarray.DataAray): surface pressure
    """
    huss = xr.apply_ufunc(huss_from_dtas_and_sp_raw, dtas, sp, dask='parallelized',
                        output_dtypes=[dtas.dtype])
    return huss

def huss_from_dtas_and_sp_raw(dtas, sp):
    """
    compute specific humidity from dewpoint and sea level pressure (ufunc)
    Equation taken from: https://confluence.ecmwf.int/pages/viewpage.action?pageId=171411214
    Parameters
    ----------
    dtas : float/np.ndarray
        dew point at 2m
    sp : float/np.ndarray
        surface pressure
    Returns
    -------
    huss : float/np.ndarray
        specific humidity at 2m
    """
    Rdry=287.0597
    Rvap=461.5250
    a1=611.21
    a3=17.502
    a4=32.19
    T0=273.16
    
    # Calculation of E saturation water vapour from Teten's formula
    E=a1**(a3*(dtas-T0)/(dtas-a4))

    # Calculation of saturation specific humidity at 2m qsat  (equal to huss)

    huss=(Rdry/Rvap)*E/(sp-((1-Rdry/Rvap)*E))
   
    return huss


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="compute specific humidity file from era5 dewpoint temperature at 2m and surface pressure (dtas and sp)"
    )
    parser.add_argument(
        "--dtas",
        type=str,
        required=True,
        help="name of dewpoint file",
    )
    parser.add_argument(
        "--sp",
        type=str,
        required=True,
        help="name of pressure file",
    )
    parser.add_argument(
        "--dtasvar",
        type=str,
        required=False,
        default='d2m',
        help="name of dewpoint variable",
    )
    parser.add_argument(
        "--spvar",
        type=str,
        required=False,
        default='sp',
        help="name of surface pressure variable",
    )
    parser.add_argument(
        "-o",
        "--fileout",
        type=str,
        required=True,
        help="name of produced specific humidity file",
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

    dtas = xr.open_dataset(args["dtas"], chunks={'time': args["chunks"]})[args["dtasvar"]]
    sp = xr.open_dataset(args["sp"], chunks={'time': args["chunks"]})[args["spvar"]]

    huss = huss_from_dtas_and_sp(dtas, sp)

    huss.to_dataset(name='huss').to_netcdf(args["fileout"])


