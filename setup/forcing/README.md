# ERA5 Forcing Data

As of March 2022, the ESMG lab has decided to use ERA5 for our forcing data. This higher resolution (both spatially and temporally, when compared to JRA55-DO) data should serve our 1/25th of a degree resolution domain well. This dataset requires extra preparation before it can be used within MOM6. 

Data Link: https://confluence.ecmwf.int/display/CKB/How+to+download+ERA5

Note that the ESMG lab has this dataset locally downloaded on Antares at `/Volumes/P5/ERA5/`. 

## Preparation Steps

1) **Download data** (or navigate to Antares if in ESMG lab)
2) **Subset data** to domain of interest: [subset_calc_era5](https://github.com/jsimkins2/nwa25/blob/main/setup/forcing/subset_calc_era5.ipynb)
    2a) In the above script, [total rain rate](https://github.com/jsimkins2/nwa25/blob/main/setup/forcing/total_rainfall_rate_gen.py) and [specific humidity](https://github.com/jsimkins2/nwa25/blob/main/setup/forcing/specific_humidity_gen.py) are both calculated since these variables are not provided by ERA5.
3) **Pad each data file**: [pad_era5](https://github.com/jsimkins2/nwa25/blob/main/setup/forcing/pad_era5.ipynb). Select the xarray method or the CDO method from the script. MOM6 performs a temporal interpolation and needs extra time values on either side of the dataset. Note that for Climate Data Operators Mergetime to be work properly, 2 timeslots must be added to either side of the data file to keep the dimensions intact.
4) **Flood each data file**: [flood_era5](https://github.com/jsimkins2/nwa25/blob/main/setup/forcing/flood_era5.ipynb). Forcing data over land cells should be removed and then interpolated based on data over ocean cells. This ensures that data along the coasts don't get confused between variables right on the edge of ocean/land. 

### Notes on ERA5

- One must compile MOM6 with updated FMS coupler surface_flux file to account for variables such as wind and temperature having different heights. This allows one to edit `input.nml` and have MOM6 perform an interpolation for variables at differnet heights. 
- When first attempting to use ERA5 files, MOM6 would crash after a few minutes and claim `ABORT MPT ERROR: MPI_COMM_WORLD rank 178 has terminated without calling MPI_Finalize() aborting job`. This was occurring because the latitudes in the ERA5 data are stored from 90N to 0. Flipping the latitudes removed this error. The scripts above take care of this issue. 
- ERA5 radiation variables come in J/m-2 instead of W/m-2. The flooding routine above converts the units to W/m-2.


# JRA55-DO Forcing Data
In runs before March 2022, JRA55-DO served as our forcing dataset. 

Data Link: https://climate.mri-jma.go.jp/pub/ocean/JRA55-do/

### Steps

1) Download Data using instructions from link above
2) Subset data
3) Pad data
4) Flood data
