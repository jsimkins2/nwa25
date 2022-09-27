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


