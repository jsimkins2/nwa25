# Conda Environment Instructions

## Python 3.7.x

1) `conda create -n xesmf_env python=3.7`
2) `conda activate xesmf_env`
3) `conda install -c conda-forge xesmf=0.3.0 esmpy=8.2.0 bottleneck=1.3.5`
4) `conda install -c conda-forge dask=2021.10.0 netcdf4`

## Python 3.9.x

NOTE: This method takes a long time for the resolver to find all the
right package combinations.

`$ conda create -n xesmf_env_test -c conda-forge xesmf=0.3.0 esmpy=8.2.0 bottleneck=1.3.5 dask=2021.10.0 netcdf4`

Last checked on Jan 26, 2023.

## Python 3.10.x

1) `conda create -n xesmf_env python=3.10.8`
2) `conda activate xesmf_env`
3) `conda install -c conda-forge xesmf=0.3.0 esmpy=8.2.0 bottleneck=1.3.5`
4) `conda install -c conda-forge dask=2021.10.0 netcdf4`

Last checked on Jan 26, 2023.

## Python 3.11.x

There are reported speed improvements to this version of python.

## Jupyter

5) `conda install -c conda-forge jupyter jupyterlab numba nodejs`

## Flooding

If flooding is required:

6) Install HCtFlood

```
git clone https://github.com/raphaeldussin/HCtFlood
cd src/HCtFlood
python -m pip install -e .
```

If HCtFlood requires more iterations to converge on a solution
here is the needed modification:
```
diff --git a/HCtFlood/kara.py b/HCtFlood/kara.py
index 539050b..00201f0 100644
--- a/HCtFlood/kara.py
+++ b/HCtFlood/kara.py
@@ -106,7 +106,7 @@ def flood_kara_ma(masked_array, spval=1e+15):


 @njit
-def flood_kara_raw(field, mask, nmax=1000):
+def flood_kara_raw(field, mask, nmax=2000):
     """Extrapolate land values onto land using the kara method
     (https://doi.org/10.1175/JPO2984.1)
```

# Notes

- ESMF packages must have mpi support, see below.
- You must use `xesmf=0.3.0` to be able to create and reuse weight files.
- There is also a problem with a later version of dask, recommend `dask=2021.10.0`

After the conda packages are installed, run: `conda list | grep mpi`, the following
packages should be similar in appearance:

```
esmf                      8.2.0           mpi_mpich_h4975321_100    conda-forge
esmpy                     8.2.0           mpi_mpich_py37h7352969_101    conda-forge
fftw                      3.3.10          nompi_h77c792f_102    conda-forge
hdf5                      1.12.1          mpi_mpich_h9c45103_0    conda-forge
libnetcdf                 4.8.1           mpi_mpich_h319fa22_1    conda-forge
mpi                       1.0                       mpich    conda-forge
mpi4py                    3.1.3            py37h52370cb_1    conda-forge
mpich                     4.0.2              h846660c_100    conda-forge
netcdf-fortran            4.5.4           mpi_mpich_h1364a43_0    conda-forge
netcdf4                   1.5.8           nompi_py37hf784469_101    conda-forge
```

It is very important that esmf, esmpy, libnetcdf, hdf5 and netcdf-fortran have
`mpi_mpich` within the build name (3rd column) of the package listing.  
If they show up as `nompi` then ESMF will not work.
