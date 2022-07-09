import numpy as np

################################################################################
def footprint_OHC(depth,temp,dens):
    # This function Calculates the ocean heat content from a temperature and
    # density profile (Leipper, Dale F., and Douglas Volgenau. "Hurricane heat
    # potential of the Gulf of Mexico". Journal of Physical Oceanography 2.3
    #(1972): 218-224).
    # Inputs: 1D vectors depth, temperature and density
    # Output: Ocean heat content of the water column in kJ/cm^3

    cp = 3985 #Heat capacity in J/(kg K)
    depth=np.abs(depth)
    dT=temp-26
    ii=np.argwhere(dT<0)
    dT[ii[:,0],ii[:,1],ii[:,2]]=0.
    j1=dT.shape[1]
    j2=dT.shape[2]

    oout=np.ones([j1,j2])*np.nan
    for k in range(j1):
        dt=np.squeeze(dT[:,k,:])
        dt[np.isnan(dt)]=0

        dz=np.squeeze(depth[:,k,:])
        dd=np.squeeze(dens[:,k,:])
        rho=np.nanmean(dd,axis=0)
        o=cp*rho*np.trapz(dt.T,dz.T)
    
        oout[k,:]=o

    ohc=oout*1e-7   # [J/m^3] to [kJ/cm^3]
    return ohc

################################################################################
def OHC_from_profile(depth,temp,dens):
    # This function Calculates the ocean heat content from a temperature and
    # density profile (Leipper, Dale F., and Douglas Volgenau. "Hurricane heat
    # potential of the Gulf of Mexico". Journal of Physical Oceanography 2.3
    #(1972): 218-224).
    # Inputs: 1D vectors depth, temperature and density
    # Output: Ocean heat content of the water column in kJ/cm^3

    cp = 3985 #Heat capacity in J/(kg K)
    ok26 = temp >= 26
    depth = np.abs(depth)

    if len(depth[ok26]) != 0:
        if np.nanmax(depth[ok26])<2:
            OHC = np.nan
        else:
            rho0 = np.nanmean(dens[ok26])
            OHC = np.abs(cp * rho0 * np.trapz(temp[ok26]-26,depth[ok26]))
            OHC = OHC * 10**(-7) # in kJ/cm^2
    else:
        OHC = np.nan
    Z26=np.nanmax(depth[ok26])

    return OHC, Z26

################################################################################
def MLD_temp_crit(dtemp,ref_depth,depth,temp):
    # This function calculates the mixed layer depth and Mixed layer temperature
    # based on a temperature criteria: T - T_at_ref_depth <= dtemp
    # Inputs
    # dtemp: delta temperature from the mixed layer depth definition used
    # ref_depth: Reference depth from the mixed layer depth definition used
    # depth and temp: 1D vectors depth and temperature
    # Output
    # MLD and MLT: mixed layer depth and Mixed layer temperature

    ok_ref_depth = np.where(depth >= ref_depth)[0][0]
    temp_ref_depth = temp[ok_ref_depth]
    delta_T = temp_ref_depth - temp
    ok_mld_temp = np.where(delta_T <= dtemp)[0]

    if ok_mld_temp.size == 0:
        MLD = np.nan
        MLT = np.nan
    else:
        MLD = depth[ok_mld_temp[-1]]
        MLT = np.nanmean(temp[ok_mld_temp])

    return MLD, MLT

################################################################################
def MLD_dens_crit(drho,ref_depth,depth,temp,dens):
    # This function calculates the mixed layer depth and Mixed layer temperature
    # based on a density criteria: rho_at_ref_depth - rho <= drho
    # Inputs
    # drho: delta density from the mixed layer depth definition used
    # ref_depth: Reference depth from the mixed layer depth definition used
    # depth, temp and dens: 1D vectors depth, temperature and density
    # Output
    # MLD and MLT: mixed layer depth and Mixed layer temperature

    ok_ref_depth = np.where(depth >= ref_depth)[0][0]
    rho_ref_depth = dens[ok_ref_depth]
    delta_rho = -(rho_ref_depth - dens)
    ok_mld_rho = np.where(delta_rho <= drho)[0]

    if ok_mld_rho.size == 0:
        MLD = np.nan
        MLT = np.nan
    else:
        MLD = depth[ok_mld_rho[-1]]
        MLT = np.nanmean(temp[ok_mld_rho])

    return MLD, MLT

################################################################################
def T100(depth,temp):
    # This function calculates the depth average temperature in the top 100
    # meters
    # Inputs
    # depth, temp: 1D vectors depth and temperature
    # Output
    # T100: depth average temperature in the top 100 meters

    okd = np.abs(depth) <= 100
    if len(np.where(np.isnan(temp[okd]))[0])>10:
        T100 = np.nan
    else:
        T100 = np.nanmean(temp[okd])
    return T100

################################################################################
def Potential_energy_anomaly100(depth,dens):
    # This function calculates the potential energy anomaly
    # (Simpson J, Brown J, Matthews J, Allen G (1990) Tidal straining, density
    # currents and stirring in the control of estuarine stratification.
    # Estuaries 13(2):125–132), in the top 100 meters
    # Inputs
    # depth, dens: 1D vectors depth and density
    # Output
    # PEA: potential energy anomaly in J/m^3

    g = 9.8 #m/s
    dindex = np.fliplr(np.where(np.asarray(np.abs(depth)) <= 100))[0]
    if len(dindex) == 0:
        PEA = np.nan
    else:
        zz = np.asarray(np.abs(depth[dindex]))
        denss = np.asarray(dens[dindex])
        ok = np.isfinite(denss)
        z = zz[ok]
        densi = denss[ok]
        if len(z)==0 or len(densi)==0 or np.min(zz) > 10 or np.max(zz) < 30:
            PEA = np.nan
        else:
            if z[-1] - z[0] > 0:
                # So PEA is < 0
                # sign = -1
                # Adding 0 to sigma integral is normalized
                z = np.append(0,z)
            else:
                # So PEA is < 0
                # sign = 1
                # Adding 0 to sigma integral is normalized
                z = np.flipud(z)
                z = np.append(0,z)
                densit = np.flipud(densi)

            # adding density at depth = 0
            densitt = np.interp(z,z[1:],densit)
            density = np.flipud(densitt)

            # defining sigma
            max_depth = np.nanmax(zz[ok])
            sigma = -1*z/max_depth
            sigma = np.flipud(sigma)

            rhomean = np.trapz(density,sigma,axis=0)
            drho = rhomean - density
            torque = drho * sigma
            PEA = g * max_depth * np.trapz(torque,sigma,axis=0)

    return PEA
