Installing MOM6+SIS2 on NCAR's Cheyenne with Updated FMScoupler
============================================================
Note that you'll have to change all instances of jsimkins to your cheyenne username

1) Navigate to work directory on glade

.. code-block :: bash

    cd /glade/work/jsimkins

2) Clone the MOM6 Examples Github Repository

.. code-block :: bash
  
    git clone --recursive https://github.com/NOAA-GFDL/MOM6-examples.git MOM6-examples

3) Update all submodules within cloned repository

.. code-block :: bash

   cd MOM6-examples
   git submodule init
   git submodule update --recursive

4) Double check that MOM6 submodules are up to date

.. code-block :: bash

   cd src/MOM6
   git submodule init
   git submodule update

5) Update `surface_flux.F90` in `src/coupler/`

Copy code from updated `surface_flux.F90` [script](https://raw.githubusercontent.com/NOAA-GFDL/FMScoupler/a81da6c81f40b64acf0ecfa55f46ba5b5b9d4f1f/shared/surface_flux.F90) and overwrite the file.

.. code-block :: bash

   vim /glade/work/jsimkins/MOM6-examples/src/coupler/surface_flux.F90
   

6) Create a build directory within /glade/work/jsimkins/MOM6-examples/

.. code-block :: bash

    cd /glade/work/jsimkins/MOM6-examples/
    mkdir build

7) Clone NCAR MOM6-Cases to /glade/work/jsimkins/

.. code-block :: bash

    cd ../
    git clone --recursive https://github.com/NCAR/MOM6-cases.git

8) Copy NCAR MOM6-Cases Cheyenne.mk files to MOM6-examples

.. code-block :: bash

    cp /glade/work/jsimkins/MOM6-cases/src/mkmf/templates/cheyenne* /glade/work/jsimkins/MOM6-examples/src/mkmf/templates/

9) Create bash file in MOM6-examples

.. code-block :: bash

    cd /glade/work/jsimkins/MOM6-examples/
    vim build_mom6_sis2.bash

10) Copy the following code to build_mom6_sis2.bash

.. code-block :: bash

    #!/bin/bash
    mkdir -p /glade/work/jsimkins/MOM6-examples/build/intel/ice_ocean_SIS2/repro/

    module load ncarenv
    module load intel
    module load netcdf
    module load ncarcompilers
    module load mpt/2.19

    mkdir -p build/intel/shared/repro/
    (cd build/intel/shared/repro/; rm -f path_names; \
    ../../../../src/mkmf/bin/list_paths -l ../../../../src/FMS; \
    ../../../../src/mkmf/bin/mkmf -t ../../../../src/mkmf/templates/cheyenne-intel.mk -p libfms.a -c "-Duse_libMPI -Duse_netCDF" path_names)

    (cd build/intel/shared/repro/; source ../../env; make NETCDF=4 REPRO=1 libfms.a -j)

    mkdir -p build/intel/ice_ocean_SIS2/repro/
    (cd build/intel/ice_ocean_SIS2/repro/; rm -f path_names; \
    ../../../../src/mkmf/bin/list_paths -l ./ ../../../../code/ ../../../../src/MOM6/config_src/{infra/FMS1,memory/dynamic_symmetric,drivers/FMS_cap,external} ../../../../src/MOM6/src/{*,*/*}/ ../../../../src/{atmos_null,coupler,land_null,ice_param,icebergs,SIS2,FMS/coupler,FMS/include}/)
    (cd build/intel/ice_ocean_SIS2/repro/; \
    ../../../../src/mkmf/bin/mkmf -t ../../../../src/mkmf/templates/cheyenne-intel.mk -o '-I../../shared/repro' -p MOM6 -l '-L../../shared/repro -lfms' -c '-Duse_AM3_physics -D_USE_LEGACY_LAND_' path_names )

    (cd build/intel/ice_ocean_SIS2/repro/; source ../../env; make NETCDF=4 REPRO=1 MOM6 -j)

11) Make build_mom6_sis2.bash executable

.. code-block :: bash

    chmod +x build_mom6_sis2.bash

12) Run build_mom6_sis2.bash

.. code-block :: bash

    ./build_mom6_sis2.bash


