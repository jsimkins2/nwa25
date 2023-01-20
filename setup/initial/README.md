# Initial Conditions 

- `generate_ic_glorys.py` is the most up to date script and is used for this set of runs.
- The final IC file is saved as NETCDF4, **however, MOM6 seems to demand a NETCDF3 file**. So, after you output the file, please convert this to NETCDF3. Saving the xarray object to NETCDF3 eats up more memory and can cause your systme to crash...hence the extra step
- Check the units of your angle_dx variable. You may need to edit rotate_uv function if your angle_dx is in radians. This assumes that your angle_dx is in degrees. The python gridtools library outputs angle_dx in degrees.

### Adding IC file to MOM_input

```
! === module MOM_state_initialization ===
THICKNESS_CONFIG = "coord"      ! default = "uniform"
                                ! A string that determines how the initial layer thicknesses are specified for a
                                ! new run:
                                !     file - read interface heights from the file specified
                                !     thickness_file - read thicknesses from the file specified
                                !       by (THICKNESS_FILE).
                                !     coord - determined by ALE coordinate.
                                !     uniform - uniform thickness layers evenly distributed
                                !       between the surface and MAXIMUM_DEPTH.
                                !     list - read a list of positive interface depths.
                                !     DOME - use a slope and channel configuration for the
                                !       DOME sill-overflow test case.
                                !     ISOMIP - use a configuration for the
                                !       ISOMIP test case.
                                !     benchmark - use the benchmark test case thicknesses.
                                !     Neverworld - use the Neverworld test case thicknesses.
                                !     search - search a density profile for the interface
                                !       densities. This is not yet implemented.
                                !     circle_obcs - the circle_obcs test case is used.
                                !     DOME2D - 2D version of DOME initialization.
                                !     adjustment2d - 2D lock exchange thickness ICs.
                                !     sloshing - sloshing gravity thickness ICs.
                                !     seamount - no motion test with seamount ICs.
                                !     dumbbell - sloshing channel ICs.
                                !     soliton - Equatorial Rossby soliton.
                                !     rossby_front - a mixed layer front in thermal wind balance.
                                !     USER - call a user modified routine.
TS_CONFIG = "file"              !
                                ! A string that determines how the initial tempertures and salinities are
                                ! specified for a new run:
                                !     file - read velocities from the file specified
                                !       by (TS_FILE).
                                !     fit - find the temperatures that are consistent with
                                !       the layer densities and salinity S_REF.
                                !     TS_profile - use temperature and salinity profiles
                                !       (read from TS_FILE) to set layer densities.
                                !     benchmark - use the benchmark test case T & S.
                                !     linear - linear in logical layer space.
                                !     DOME2D - 2D DOME initialization.
                                !     ISOMIP - ISOMIP initialization.
                                !     adjustment2d - 2d lock exchange T/S ICs.
                                !     sloshing - sloshing mode T/S ICs.
                                !     seamount - no motion test with seamount ICs.
                                !     dumbbell - sloshing channel ICs.
                                !     rossby_front - a mixed layer front in thermal wind balance.
                                !     SCM_CVMix_tests - used in the SCM CVMix tests.
                                !     USER - call a user modified routine.
TS_FILE = "ic/soda_ic_75z_2010-01-05.nc" !
                                ! The initial condition file for temperature.
TEMP_IC_VAR = "temp"            ! default = "PTEMP"
                                ! The initial condition variable for potential temperature.
SALT_IC_VAR = "salt"            ! default = "SALT"
                                ! The initial condition variable for salinity.
VELOCITY_CONFIG = "file"        ! default = "zero"
                                ! A string that determines how the initial velocities are specified for a new
                                ! run:
                                !     file - read velocities from the file specified
                                !       by (VELOCITY_FILE).
                                !     zero - the fluid is initially at rest.
                                !     uniform - the flow is uniform (determined by
                                !       parameters INITIAL_U_CONST and INITIAL_V_CONST).
                                !     rossby_front - a mixed layer front in thermal wind balance.
                                !     soliton - Equatorial Rossby soliton.
                                !     USER - call a user modified routine.
VELOCITY_FILE = "ic/soda_ic_75z_2010-01-05.nc" !
                                ! The name of the velocity initial condition file.
                                ```
