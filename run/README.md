# Successful Simulation Configurations


This folder contains successful simulation configurations for NWA25 runs on various machines. A description of operations is as follows.

### Cheyenne

The original testing configuration for NWA25. This set of configs was adapated from Chassignet & Xu 2017 NA12 configuration. 3 years of simulation were performed before we ran out of time on Cheyenne. 

### Sim1.0

A slight adaptation of Cheyenne parameters with input from Mehmet Ilicak. The following parameters were important settings for proper Gulf Stream Path.

  LAPLACIAN = False                                   
  SMAGORINSKY_KH = False                                
  SMAG_LAP_CONST = 2.0       
  SMAGORINSKY_AH = True                                                         
  SMAG_BI_CONST = 0.06                                                           
  BIHARMONIC = True                                                               
  AH = 0.0                                              
  AH_VEL_SCALE = 0.0         
  MIXEDLAYER_RESTRAT = True      !
  FOX_KEMPER_ML_RESTRAT_COEF = 1.0
  MLE_USE_PBL_MLD = True          
  MLE_FRONT_LENGTH = 5000.0       

This simulation is currently being run on Gaea and so far 5 years of simulated time have been run, with more years ongoing.

### SIMULATION 1.2

Updated OBC generated files to reflect updates to xesmf and xarray. Python file for OBC generation has been updated in setup.

### Simulation 1.3

DT_THERM=1800
DT=300
dt_cpld=900
dt_atmos=900

Lowering DT_THERM timestep to 1800.

Per Andrew Ross:

To summarize, NWA12 behaves normally with DT_THERM = 1200 or 1800. If I try to increase DT_THERM to 3600, however, the Gulf Stream rapidly dissipates with a bizarre wave traveling to the south, and then reestablishes as a very weak current that doesn't separate at all. At the same time there appears to be extensive mixing in the Gulf Stream region, evident as cooling above 1000 m and warming below 1000 m all the way down to 5000 m. The initial dissipation starts immediately after switching to DT_THERM = 3600, and the total collapse and reestablishment happens within the first few months. Afterwards the weak boundary current remains roughly stable for years. I've attached some figures showing sea level as this is happening.

Previously I've run different experiments to try to identify the issue. I completely turned off tides and the open boundaries but still encountered the problem. The only control that I've found is the biharmonic viscosity (either Smagorinsky or using the velocity scale); with a low viscosity, the Gulf Stream dissipates much more rapidly, while with a high viscosity, the dissipation still happens but it occurs a lot slower.


### Simulation 1.4

DT_THERM=1200
DT=300
dt_cpld=600
dt_atmos=600

Lowering DT_THERM timestep to 1200
Lowering dt_cpld timestep to 600
Lowering dt_atmos timestep to 600


### Simulation 1.5

DT_THERM=1200
DT=300
dt_cpld=600
dt_atmos=600

Arakawa lamb coriolis scheme instead of sadourny energy


### Simulation 1.6

DT_THERM=1200
DT=300
dt_cpld=600
dt_atmos=600

Sadourney Enstropy coriolis scheme instead of sadourny energy


### Simulation 1.7

DT_THERM=900
DT=300
dt_cpld=600
dt_atmos=600

# Simulation 1.8

DT_THERM=300
DT=300
dt_cpld=600
dt_atmos=600

# Simulation 2.0

Alistair's recommendations to Enrique to attempt to fix the Gulf Stream issue (protruding too far south in previous runs and northern boundary pushing cold water far too south than it should be). This conversations occurred at the end of December, 2022.

DT_THERM == DT_CPLD == 1200
DT = 300
DT_ATMOS can be set tot he frequency of our forcing (hourly?) or to same as DT_THERM 
KV = 1 e. -4
SMAG_BI = 0.06
LAPLACIAN = TRUE
SMAGORINSKY_KH = TRUE
SMAGORINSKY_LAPLACIAN_CONSTANT = 0.15

RSLN_SCALED_KH, _KHTN, _KHTL   should all be set to =FALSE


# Simulation 2.1

Same as Simulation 2.0, but BBL_EFIC = 0.1 (per Dujuan's recommendations) AND turning salinity off.


# Simulation 2.2

Simulation 2.0 but with Closed Boundary Conditions
