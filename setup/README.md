# Conda Environment Instructions

These scripts have been run successfully on Antares using a Python 3.7 conda environment located here - `/home/james/anaconda3/envs/xesmf/bin`. Here is a recipe for creating this. 

1) `conda create -n xesmf_env python=3.7`
2) `conda activate xesmf_env`
3) `conda install -c conda-forge xesmf esmpy=8.2.0`
4) `conda install -c conda-forge dask netCDF4`



Notes:

- Note that `esmpy=8.2.0` must be [installed in the same instance](https://github.com/JiaweiZhuang/xESMF/issues/47#issuecomment-665516640) of `xesmf` installation.
- If you're running on antares, my environment for this can be found at `/home/james/anaconda3/envs/xesmf`
