import numpy as np
import xarray as xr
import argparse

# functions borrowed from metpy
def mixing_ratio(partial_press, total_press, molecular_weight_ratio=0.622):
    return (molecular_weight_ratio * partial_press
                / (total_press - partial_press))

def specific_humidity_from_mixing_ratio(mr):
    return mr / (1 + mr)


def saturation_vapor_pressure(temperature):
    sat_pressure_0c = 6.112e2 # Pa
    return sat_pressure_0c * np.exp(17.67 * (temperature - 273.15) # K -> C
                                        / (temperature - 29.65))   # K -> C


def saturation_mixing_ratio(total_press, temperature):
    return mixing_ratio(saturation_vapor_pressure(temperature), total_press)


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

    sphum = sphum.to_dataset()
    tdew = xr.open_dataset(args["dtas"], chunks={'time': args["chunks"]})[args["dtasvar"]]
    pair = xr.open_dataset(args["sp"], chunks={'time': args["chunks"]})[args["spvar"]]

    smr = saturation_mixing_ratio(pair, tdew)
    sphum = specific_humidity_from_mixing_ratio(smr)
    sphum.name = 'huss'
    sphum = sphum.to_dataset()


    sphum.to_dataset(name='huss').to_netcdf(args["fileout"])


