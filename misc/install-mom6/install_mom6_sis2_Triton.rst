Installing MOM6+SIS2 on Rutger's Triton
============================================================

1) Clone the MOM6 Examples Github Repository

.. code-block :: bash
  
    git clone --recursive https://github.com/NOAA-GFDL/MOM6-examples.git MOM6-examples

2) Update all submodules within cloned repository

.. code-block :: bash

   cd MOM6-examples
   git submodule init
   git submodule update --recursive

3) Double check that MOM6 submodules are up to date

.. code-block :: bash

   cd src/MOM6
   git submodule init
   git submodule update

4) Create a build directory within ~/MOM6-examples/

.. code-block :: bash

    cd ../../
    mkdir build

5A) Edit the `triton-gnu.mk` flags for netCDF4 use

.. code-block :: bash

    vim ~/MOM6-examples/src/mkmf/templates/triton-gnu.mk

- Ensure the FFLAGS and CFLAGS are using netCDF4. (i.e., `FFLAGS += -DGFORTRAN -Duse_netCDF4` and `CFLAGS := -D__IFC -Duse_netCDF4`
- using netCDF3 will result in restart errors during model execution


5B) Update `surface_flux.F90` in `src/coupler/`

Copy code from updated `surface_flux.F90` [script](https://raw.githubusercontent.com/NOAA-GFDL/FMScoupler/a81da6c81f40b64acf0ecfa55f46ba5b5b9d4f1f/shared/surface_flux.F90) and overwrite the file.

.. code-block :: bash

   vim /home/james/MOM6-examples/src/coupler/surface_flux.F90
   
6) Create bash file in MOM6-examples

.. code-block :: bash

    cd ~/MOM6-examples/
    vim build_mom6_sis2.bash

7) Copy the following code to build_mom6_sis2.bash

.. code-block :: bash

    #!/bin/bash
    
    mkdir -p /home/james/MOM6-examples/build/gnu/ice_ocean_SIS2/repro/
    
    module load netcdf-fortran

    mkdir -p build/gnu/shared/repro/
    (cd build/gnu/shared/repro/; rm -f path_names; \
    ../../../../src/mkmf/bin/list_paths -l ../../../../src/FMS; \
    ../../../../src/mkmf/bin/mkmf -t ../../../../src/mkmf/templates/triton-gnu.mk -p libfms.a -c "-Duse_libMPI -Duse_netCDF" path_names)

    (cd build/gnu/shared/repro/; source ../../env; make NETCDF=4 REPRO=1 libfms.a -j)

    mkdir -p build/gnu/ice_ocean_SIS2/repro/
    (cd build/gnu/ice_ocean_SIS2/repro/; rm -f path_names; \
    ../../../../src/mkmf/bin/list_paths -l ./ ../../../../src/MOM6/config_src/{infra/FMS1,memory/dynamic_symmetric,drivers/FMS_cap,external} ../../../../src/MOM6/src/{*,*/*}/ ../../../../src/{atmos_null,coupler,land_null,ice_param,icebergs,SIS2,FMS/coupler,FMS/include}/)
    (cd build/gnu/ice_ocean_SIS2/repro/; \
    ../../../../src/mkmf/bin/mkmf -t ../../../../src/mkmf/templates/triton-gnu.mk -o '-I../../shared/repro' -p MOM6 -l '-L../../shared/repro -lfms' -c '-Duse_AM3_physics -D_USE_LEGACY_LAND_' path_names )

    (cd build/gnu/ice_ocean_SIS2/repro/; source ../../env; make NETCDF=4 REPRO=1 MOM6 -j)



8) Run build_mom6_sis2.bash

.. code-block :: bash
    
    chmod u+x build_mom6_sis2.bash
    ./build_mom6_sis2.bash


