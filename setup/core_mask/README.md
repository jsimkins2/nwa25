# Core Masking

MOM6 allows a user to mask out cores that are over 100% land cells. We can mask out these cores since the data will be NaN'd out following computation anyway. This does not improve the overall speed of a model run, but it increases efficiency and thus shrinks core hours per day of simulation. In our NWA25 case, our grid contains roughly 40% land. We were able to mask out ~32% of cores as land cores. Thus, our core hours per day of simulation cost was brought down from ~60 core hours / day of simulation to ~41 core hours / day of simulation.

## Software

[FRE-NCtools](https://github.com/NOAA-GFDL/FRE-NCtools)

The `check_mask` function will create a mask table text file that can be used by MOM6. This simple text file provides MOM6 with coordinates of cores that are 100% land. for example, if we start with a core layout of 20 x 36 and core 10,12 is 100% land, 10,12 will be added as a line in the mask table text file.

### Check_mask Usage

`module load netcdf`
`/t5/workdir/cermak/local/bin/check_mask --grid_file ocean_mosaic.nc --ocean_topog ocean_topog.nc --layout 20,36`

### Mask Table

In this example, the `check_mask` function above returned the file `mask_table.228.20x36`. The naming of this file indicates that 228 cores of our 720 total cores used (20x36 layout) are land. Thus, only 492 cores have ocean. Now we need to tell MOM6 to use this file. Mask tables are defaultly read from the INPUT directory of a MOM6 simulation. Place the file in this directory. Then, we must update our `MOM_layout` (and `SIS_layout` if you're running ice too). 

- Add these lines to your `MOM_layout` and `SIS_layout`
`LAYOUT = 20,36`
`MASKTABLE = "mask_table.228.20x36" ! 20*36-228 = 492 PEs`

### Requesting Specific Number of Cores

The original layout we began with in this example was 20x36 (720 total cores). Only 492 cores matter, however. So, we need to tell the batch scheduler that we only want to request 492 cores. On Cheyenne, this can be done like so...

`mpiexec_mpt -np 492 /glade/work/jsimkins/runs-MOM6/nwa25/MOM6`

------------------------------------------
# Expanding IO Layout

FMS allows a user to expand the IO layout for output files. When the layout is expanded, FMS will use a requested number of  write nodes as opposed to the default value of 1 write node. For example, if you request a 2,2 `IO_layout`, then 4 write nodes will be used and you will have 4 output files for every output set requested in `diag_table`. In this case, the horizontal grid is split into 4 quadrants. Note that this increased number of write nodes does not currently enhance the speed of a simulation. FMS is working on this capability, but as of 2/10/2022 it does not enhance the speed of a simluation; it only provides a method to distribute the data for large output files. 

If you wish to use an IO layout, the following lines must be added to `MOM_layout` and `SIS_layout`

`IO_LAYOUT = 2,2`

In this example, 4 output files will be written (each will be 1 horizontal quadrant). They will output in the form of "output_file_0000", "output_file_0001", etc.

## Combining Output

If you enabled an IO layout, the files can be spatially combined. For smaller files, the Python package `xarray` can be used to combine these files together spatially.

#### mppnccombine-fast

https://github.com/coecms/mppnccombine-fast

Vague Installation Recipe:

https://github.com/jsimkins2/nwa25/blob/main/setup/core_mask/install_mppnccombine-fast.md

This software stems from the FRE-NCtools mppnccombine function detailed below with the major difference being that mppnccombine-fast actually works. This allows a user to stitch together the output files by leveraging HDF5 and parallel processing. 

`mpirun -n 10 /t5/workdir/cermak/local/bin/mppnccombine-fast --output combined_output output_file_0000 output_file_0001 output_file_0002 output_file_0003`

 
#### Xarray

- Xarray `open_mfdataset`
```
import xarray as xr
ds = xr.open_mfdataset("./output_file_000*", combine='by_coords')
```

- Xarray `open_dataset`
```
import xarray as xr
ds0=xr_open_dataset("./output_file_0000")
ds1=xr_open_dataset("./output_file_0001")
ds2=xr_open_dataset("./output_file_0002")
ds3=xr_open_dataset("./output_file_0003")
ds_combined=xr.combine_by_coords([ds0,ds1,ds2,ds3],combine_attrs="drop")
```

#### FRE-NCtools mppnccombine

For larger files, [FRE-NCtools](https://github.com/NOAA-GFDL/FRE-NCtools) has a function called `mppnccombine` that can stitch files together. It works like so;

- mppnccombine
`/t5/workdir/cermak/local/bin/mppnccombine /t0/scratch/james/output-MOM6/nwa25/chl_update/20100105.ocean_5day.nc.000*`

However, this function quickly runs into memory issues. Even when using a batch scheduler and enabling parallel processing, this function only worked past the memory issues when the individual files were ~3gb.

These functions were used this to split files in half by the time dimension:

```
ncks -d time,0,10 output_file_0_10_0000 output_file_0000
ncks -d time,11,20 output_file_11_20_0000 output_file_0000
```

When the file size is small enough, mppnccombine can actually work. 

