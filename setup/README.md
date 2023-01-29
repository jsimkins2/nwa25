# Conda Environment Instructions

These scripts have been run successfully using a Python 3.7 conda environment.
Here is a recipe for creating this. 

1) `conda create -n xesmf_env python=3.7`
2) `conda activate xesmf_env`
3) `conda install -c conda-forge xesmf esmpy=8.2.0`
4) `conda install -c conda-forge dask netcdf4`

Please refer to [CONDA](CONDA.md) link for expanded information for installation of
esmf, esmpy and xesmf for python.

## Notes

- Note that `esmpy=8.2.0` must be [installed in the same instance](https://github.com/JiaweiZhuang/xESMF/issues/47#issuecomment-665516640) of `xesmf` installation.  The packages must be built with MPI (either `mpich` or `openmpi`)

## MPI

To ensure the installed conda packages are installed with MPI, look for `mpich` or `openmpi` in the `Build`
column.  Here is some examples using `conda search esmpy`:

```
# Name                       Version           Build              Channel
esmpy                          8.4.0 mpi_mpich_py311hf216de5_101  conda-forge
esmpy                          8.4.0 mpi_openmpi_py311hcdf3ef6_1  conda-forge
esmpy                          8.4.0 nompi_py311h8e2db7d_1  conda-forge
```

In this example, the first two packages are built with MPI since ('mpich` and `openmpi`) appear in 
the build title.  The last package does not have MPI support and will not work if you attempt to
utilize operations that require MPI to work.

# Flooding

If you require flooding of OBCs, then you also need to install:
 * https://github.com/raphaeldussin/HCtFlood

# Site specific notes

## antares

- An example installed environment for this can be found at `/home/james/anaconda3/envs/xesmf`
