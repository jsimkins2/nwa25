Must use netcdf library that has been built with openmpi

Triton Install

`git clone https://github.com/coecms/mppnccombine-fast`
`module load netcdf`
`module load cmake`
`cd src/mppnccombine-fast`
`export HDF5_ROOT=/opt/ohpc/pub/libs/gnu8/openmpi3/hdf5/1.10.5`
`export NETCDF_ROOT=/opt/ohpc/pub/libs/gnu8/openmpi3/netcdf/4.7.1`
`make`
`make PREFIX=/home/cermak/workdir/local install`


Theoretical Cheyenne Install (never got working due to path issues for HDF)

`module load netcdf-mpi`
`module load cmake`
`cd src/mppnccombine-fast`
`export HDF5_ROOT=/glade/u/apps/ch/opt/hdf5/1.10.1/gnu/7.1.0`
`export NETCDF_ROOT=/glade/u/apps/ch/opt/netcdf-mpi/4.8.1/mpt/2.22/gnu/11.1.0/`
`make`
`make PREFIX=/glade/work/jsimkins/mppnccombine-fast/ install`

