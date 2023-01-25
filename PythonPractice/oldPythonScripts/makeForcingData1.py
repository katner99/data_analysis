#!/usr/bin/env python
import sys
import shutil
import os.path
import datetime
import numpy as np
from netCDF4 import Dataset
import math
import logging
import glob
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
# from file forcing_classes_functions.py
from forcing_classes_functions import   forceVariable # class to contain data on variables
from forcing_classes_functions import   fx,fy,NorthPoleWindVectorCoherance,Interpolate, Interpolate_perDay # functions for processing data
from forcing_classes_functions import   DegreesToKelvin, KelvinToDegrees, PCtoFraction # functions for converting data

# ---------------------- 00100 UTILITIES --------------------------------------
# Get command line arguments
if len(sys.argv) != 5:
    sys.exit("Stopped - Incorrect arguements. Use: makeForcingData.py <type:ncep/jra55> <year> <daily/6hour> <varNum>")

type = sys.argv[1]              # reanalysis name, NCEP, JRA55, ERA5
year = sys.argv[2]              # year to run through preprocessing
daily = sys.argv[3]             # how many values given for each day
varNum = int(sys.argv[4])       # variable to preprocess


# set up logging to file 
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename='forcing.'+str(year)+"_"+str(varNum)+"_"+str(type)+"_"+str(daily)+".log",
                    filemode='w')
# define a Handler which writes INFO messages or higher
console = logging.StreamHandler()
console.setLevel(logging.INFO)
# set a format which is simpler for console use
formatter = logging.Formatter('%(levelname)-5s %(message)s')
# tell the handler to use this format
console.setFormatter(formatter)
# add the handler to the root logger
logging.getLogger('').addHandler(console)
log = logging.getLogger("Run")

# Get command line arguments
if len(sys.argv) != 5 and (type != "ncep" or type != "jra55" or type!= "era5"):
    log.error("Command line arguments not correct")
    log.error("Use: makeForcingData.py <type:ncep/jra55> <year> <daily/6hour> <varNum>")
    sys.exit("Stopped")

log.info("Processing for type: "+str(type))
log.info("Processing for year: "+str(year))

# ------------------------- 00200 Define data --------------------------------------
origVars = []
# using a class (forceVariable) to store data and relevant functions
    #    self.fileName = fileName
    #    self.nc_id = nc_id
    #    self.varId = varId
    #    self.inVarId = inVarId
    #    self.limits = limits
    #    self.inBulk = inBulk
    #    self.factor = factor
    #    self.directConversion = directConversion
    #    self.units = units
    #    self.conversionFunc = conversionFunc
    #    self.type_tag = type_tag
    #    self.type_var = type_var
    #    self.type_units = type_units
    #    self.type_factor = type_factor
    #    self.type_time_res = type_time_res
    #    self.type_steps_from_midnight = type_steps_from_midnight


# ------------------------- 00210 NCEP variables --------------------------------------   
# air temperature and pressure required, but also included in humidity calculation with the third file
if varNum == -1: # ncep processing
    origVars.append(forceVariable("air.2m.gauss."+str(year)+".nc", -1, "air", "air", (-100,100), 1, True, False, "degC", KelvinToDegrees,
                                    "anl_surf.011_tmp.reg_tl319.", "var11", "degK", 1, 6, 0))
    origVars.append(forceVariable("pres.sfc.gauss."+str(year)+".nc", -1, "pres", "pres", (10000,200000), 1, True, True, "kg/m^2/s", None,
                                    "anl_surf.001_pres.reg_tl319.","var1", "kg/m^2/s", 1, 6, 0 ))
    origVars.append(forceVariable("shum.2m.gauss."+str(year)+".nc", -1, "humidity", "shum", (0,10), 1, True, False, "kg/kg", None,
                                    "anl_surf.051_spfh.reg_tl319.","var51", "kg/kg", 1, 6, 0))

    # origVars.append(forceVariable("shum.2m.gauss."+str(year)+".nc", -1, "humidity", "shum", (0,10), 1, True, True, "kg/kg", None,
    #                                 "anl_surf.051_spfh.reg_tl319.","var51", "kg/kg", 1, 6, 0))

    # # # # keep these two together to calc the wind speed
    origVars.append(forceVariable("uwnd.10m.gauss."+str(year)+".nc", -1, "wspd", "uwnd", (0,100), 1, True, False, "m/s", None,
                                    "anl_surf.033_ugrd.reg_tl319.","var33", "m/s", 1, 6, 0 ))
    origVars.append(forceVariable("vwnd.10m.gauss."+str(year)+".nc", -1, "wspd", "vwnd", (0,100), 1, True, False, "m/s" , None,
                                    "anl_surf.034_vgrd.reg_tl319.","var34", "m/s", 1, 6, 0 ))

    origVars.append(forceVariable("uwnd.10m.gauss."+str(year)+".nc", -1, "uwnd", "uwnd",(-100,100), 1, True, False, "m/s", None,
                                    "anl_surf.033_ugrd.reg_tl319.","var33", "m/s", 1, 6, 0 ))
    origVars.append(forceVariable("vwnd.10m.gauss."+str(year)+".nc", -1, "vwnd", "vwnd",(-100,100), 1, True, False, "m/s" , None,
                                    "anl_surf.034_vgrd.reg_tl319.","var34", "m/s", 1, 6, 0 ))

    origVars.append(forceVariable("tcdc.eatm.gauss."+str(year)+".nc", -1, "tcdc", "tcdc", (0,1), 1, True, False, "", PCtoFraction,
                                    "fcst_surf.071_tcdc.reg_tl319.","var71", "", 1, 3, 1 ))
    origVars.append(forceVariable("prate.sfc.gauss."+str(year)+".nc", -1, "prate", "prate",(0,100), 1, True, True, "kg/m^2/s", None,
                                    "fcst_phy2m.061_tprat.reg_tl319.","var61", "kg/m^2/s", 1.15741E-5, 3, 1 )) #mm/day -> kg/m^2/s

# # # two components for surface stress. Again these need to be kept sequential
    origVars.append(forceVariable("uflx.sfc.gauss."+str(year)+".nc", -1, "uflx", "uflx", (-100,100), 1, False, False, "", None,
                                    "fcst_phy2m.124_uflx.reg_tl319.","var124", "N/m^2", 1, 3, 1 ))
    origVars.append(forceVariable("vflx.sfc.gauss."+str(year)+".nc", -1, "vflx", "vflx", (-100,100), 1, False, False, "", None,
                                    "fcst_phy2m.125_vflx.reg_tl319.","var125", "N/m^2", 1, 3, 1 ))


# add extra forcing files for CORE bulk formulation (the above is for CLIO, low res was based on this)
# origVars.append(forceVariable("snow.sfc.gauss."+str(year)+".nc", -1, "snow", "snow", (0,1000), 1, True, True, "kg/m^2/s", None,
#                                 "fcst_phy2m.064_srweq.reg_tl319.","var64", "kg/m^2/s", 1.15741E-5, 3, 1 ))
# origVars.append(forceVariable("swrad.sfc.gauss."+str(year)+".nc", -1, "srad", "srad", (0,3000), 1, True, True, "W/m^2", None,
#                                 "fcst_phy2m.161_csdsf.reg_tl319.","var161", "kg/m^2/s", 1, 3, 1 ))
# origVars.append(forceVariable("lwrad.sfc.gauss."+str(year)+".nc", -1, "lrad", "lrad", (0,1000), 1, True, True, "W/m^2", None,
#                                 "fcst_phy2m.163_csdlf.reg_tl319.","var163", "kg/m^2/s", 1, 3, 1 ))
# origVars.append(forceVariable("swrad.sfc.gauss."+str(year)+".nc", -1, "srad", "srad",(0,3000), 1, True, True, "W/m^2", None,
#                                 "fcst_phy2m.204_dswrf.reg_tl319.","var204", "kg/m^2/s", 1, 3, 1 ))
# origVars.append(forceVariable("lwrad.sfc.gauss."+str(year)+".nc", -1, "lrad", "lrad", (0,1000), 1, True, True, "W/m^2", None,
#                                 "fcst_phy2m.205_dlwrf.reg_tl319.","var205", "kg/m^2/s", 1, 3, 1 ))

    
    
# ------------------------- 00210 JRA55 variables --------------------------------------  
# self, fileName, nc_id, varId, inVarId, limits, factor, inBulk, directConversion, units, conversionFunc, 
#                   type_tag, type_var, type_units, type_factor, type_time_res,type_steps_from_midnight):

# tas - air temperature
if varNum == 0:
    origVars.append(forceVariable("air.jraso."+str(year)+".nc", -1, "air", "air",(0,400), 1, True, True, "degK", None,
                                "tas_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_", "tas", "degK", 1, 3, 0))
# huss - specific humidty
if varNum == 1:
    origVars.append(forceVariable("shum.jraso."+str(year)+".nc", -1, "shum", "shum", (0,1), 1, True, True, "kg/kg", None,
                                 "huss_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_","huss", "kg/kg", 1, 3, 0))

# psl - pressure at sea level
if varNum == 2:
    origVars.append(forceVariable("pres.jraso."+str(year)+".nc", -1, "pres", "pres", (10000,200000), 1, True, True, "kg/m^2/s", None,
                                "psl_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_","psl", "kg/m^2/s", 1, 3, 0 ))

# rsds - short wave rad
if varNum == 3:
    origVars.append(forceVariable("rsds.jraso."+str(year)+".nc", -1, "rsds", "rsds",  (0,3000), 1, True, True, "W/m^2", None,
                                "rsds_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_","rsds", "W/m^2", 1, 3, 0 ))
# rlds - long wave rad
if varNum == 4:
    origVars.append(forceVariable("rlds.jraso."+str(year)+".nc", -1, "rlds", "rlds",  (0,1000), 1, True, True, "W/m^2", None,
                                "rlds_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_","rlds", "W/m^2", 1, 3, 0 ))

# prra - precipitate
if varNum == 5:
    origVars.append(forceVariable("prra.jraso."+str(year)+".nc", -1, "prra", "prra", (0,100), 1, True, True, "kg/m^2/s", None,
                                "prra_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_","prra", "kg/m^2/s", 1, 3, 0 )) 

# prsn - precipitate (snow)
if varNum == 6:
    origVars.append(forceVariable("snow.jraso."+str(year)+".nc", -1, "prsn", "prsn", (0,100), 1, True, True, "kg/m^2/s", None,
                                "prsn_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_","prsn", "kg/m^2/s", 1, 3, 0 )) 

# uas # vas
if varNum == 7:
    origVars.append(forceVariable("uflx.jraso."+str(year)+".nc", -1, "uflx", "uflx", (-100,100), 1, False, True, "m/s", None,
                                "uas_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_","uas", "m/s", 1, 3, 0 )) 
    origVars.append(forceVariable("vflx.jraso."+str(year)+".nc", -1, "vflx", "vflx", (-100,100), 1, False, True, "m/s", None,
                                "vas_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_","vas", "m/s", 1, 3, 0 )) 


# ------------------------- 00220 ERA5 variables --------------------------------------  
# self, fileName, nc_id, varId, inVarId, limits, factor, inBulk, directConversion, units, conversionFunc, 
#                   type_tag, type_var, type_units, type_factor, type_time_res,type_steps_from_midnight):

# t2m = temperature
if varNum == 8:
    origVars.append(forceVariable("air.era."+str(year)+".nc", -1, "air", "air",(0,400), 1, True, True, "degK", None,
                                "2m_temperature_ERA5_", "t2m", "degK", 1, 1, 0))
# specific humidty, to be calculated needs the dewpoint temperature and the surface pressure:
# d2m = dewpoint temperature
if varNum == 9:
    origVars.append(forceVariable("shum.era."+str(year)+".nc", -1, "humidity", "shum", (0,1), 1, True, True, "kg/kg", None,
                                 "2m_dewpoint_temperature_ERA5_","d2m", "degK", 1, 1, 0))

# sp - surface pressure
if varNum == 10:
    origVars.append(forceVariable("pres.era."+str(year)+".nc", -1, "pres", "pres", (10000,200000), 1, True, True, "kg/m^2/s", None,
                                "surface_pressure_ERA5_","sp", "kg/m^2/s", 1, 1, 0 ))

# msdwswrf - short wave rad
if varNum == 11:
    origVars.append(forceVariable("rsds.era."+str(year)+".nc", -1, "rsds", "rsds",  (0,3000), 1, True, True, "W/m^2", None,
                                "mean_surface_downward_short_wave_radiation_flux_ERA5_","msdwswrf", "W/m^2", 1, 1, 0 ))
# msdwlwrf - long wave rad
if varNum == 12:
    origVars.append(forceVariable("rlds.era."+str(year)+".nc", -1, "rlds", "rlds",  (0,1000), 1, True, True, "W/m^2", None,
                                "mean_surface_downward_long_wave_radiation_flux_ERA5_","msdwlwrf", "W/m^2", 1, 1, 0 ))

# tp - precipitate
if varNum == 13:
    origVars.append(forceVariable("prra.era."+str(year)+".nc", -1, "prra", "prra", (0,100), 1, True, True, "kg/m^2/s", None,
                                "total_precipitation_ERA5_","tp", "kg/m^2/s", 1, 1, 0 )) 

# sf - precipitate (snow)
if varNum == 14:
    origVars.append(forceVariable("snow.era."+str(year)+".nc", -1, "prsn", "prsn", (0,100), 1, True, True, "kg/m^2/s", None,
                                "snowfall_ERA5_","sf", "kg/m^2/s", 1, 1, 0 )) 

# u10 # v10 
if varNum == 15:
    origVars.append(forceVariable("uflx.era."+str(year)+".nc", -1, "uflx", "uflx", (-100,100), 1, False, True, "m/s", None,
                                "10m_u_component_of_wind_ERA5_","u10", "m/s", 1, 1, 0 )) 
    origVars.append(forceVariable("vflx.era."+str(year)+".nc", -1, "vflx", "vflx", (-100,100), 1, False, True, "m/s", None,
                                "10m_v_component_of_wind_ERA5_","v10", "m/s", 1, 1, 0 )) 



# append origVars value with daily or time interval
if type == 'jra55' or type == 'era5':
    for i in range (0,len(origVars)):
        origVars[i].fileName = origVars[i].fileName[:-3]+"."+str(daily)+".nc"

# # # the skt file is skin/surface temperature and only needed for a forcing file
# # # called ncep_bulk_force, one has not been generated since 2008
# # # origVars.append(forceVariable("skt.sfc.gauss."+str(year)+".nc", -1, "-", (-100,100), 1, False, False, "", None ))

createMod = False
createOnlyMod = False # set to True if don't want to recreate uflx and vflx files

# add extra forcing files for CORE bulk formulation (the above is for CLIO, low res was based on this)
# Modulus of momentum fluxes - calc after u/v flux interpolations, at end when do kelvin file
# origVars.append(forceVariable(""+str(year)+".nc", -1, "tdif", (-100,100), 1, False, False, "", None,
#                                 "","", "N/m^2", 1, 3, 1 ))



# ----------------- 00230 pre-process JRA55 data --------------------------------------
# JRA55 files are month-by-month and 6 hourly so pre-processing steps need to be carried out to make add tags
# the idea is we populate a file with the same name as the ncep files for the interpolation, an equivalent file that can 
# be compared.
if type == 'jra55':
    # time:units = "hours since 1800-01-01 00:00:0.0" ;
    originDate = datetime.datetime(1800,1,1,0,0,0)
    for v in range(0,len(origVars)):
        if origVars[v].type_tag != "":

            # create equivalent file to the NCEP input
            equivFile = origVars[v].fileName
            # if file exists, skip
            if os.path.isfile(equivFile):
                print(equivFile, " exists")
            else:
                # get list of jra files
                if len(sorted(glob.glob(origVars[v].type_tag+str(int(year)-1)+"*.nc"))) > 0:
                    typePrevMonth = [sorted(glob.glob(origVars[v].type_tag+str(int(year)-1)+"*.nc"))[-1]]
                else:
                    typePrevMonth = [sorted(glob.glob(origVars[v].type_tag+year+"*.nc"))[0]]
 
                typeFiles = sorted(glob.glob(origVars[v].type_tag+year+"*.nc"))
                               
                if len(sorted(glob.glob(origVars[v].type_tag+str(int(year)+1)+"*.nc"))) > 0:
                    typeNextMonth = [sorted(glob.glob(origVars[v].type_tag+str(int(year)+1)+"*.nc"))[0]]
                else:
                    typeNextMonth = [sorted(glob.glob(origVars[v].type_tag+year+"*.nc"))[-1]]

                typePrevMonth.extend(typeFiles)
                typePrevMonth.extend(typeNextMonth)
                typeFiles = typePrevMonth
                print(typeFiles)

                # open first to get dimensions
                nc_type_id=Dataset(typeFiles[0],"r")
                type_lat = nc_type_id.variables['lat'][:] 
                type_lon = nc_type_id.variables['lon'][:] 
                type_time = nc_type_id.variables['time'][:] 

                nc_eqiv_id = Dataset(equivFile, 'w', format=nc_type_id.data_model)
                # nc_eqiv_id = Dataset(equivFile, 'w', format='NETCDF3_CLASSIC')

                log.info("Creating equivalent to NCEP for jra55: "+equivFile+ " "+str(nc_eqiv_id.data_model))
                nc_eqiv_id.createDimension("lat", len(type_lat)) 
                nc_eqiv_id.createDimension("lon", len(type_lon)) 
                nc_eqiv_id.createDimension("time", None)
                nc_eqiv_lat = nc_eqiv_id.createVariable("lat", nc_type_id.variables["lat"].dtype, ("lat"))
                nc_eqiv_lon = nc_eqiv_id.createVariable("lon", nc_type_id.variables["lon"].dtype, ("lon"))
                nc_eqiv_time = nc_eqiv_id.createVariable("time", nc_type_id.variables["time"].dtype, ("time"))
                # varName = origVars[v].fileName.split('.')[0]
                varName = origVars[v].varId
                nc_eqiv_var = nc_eqiv_id.createVariable(varName, nc_type_id.variables[origVars[v].type_var].dtype, ("time", "lat", "lon"))

                # do a check that type_lat is the right way up for what we want
                latSwap = 0
                if type_lat[-1] > type_lat[0]:
                    log.info("Swapping latitudes as upside down")
                    type_lat = type_lat * -1
                    latSwap = 1
                nc_eqiv_id.variables['lat'][:] = type_lat

                nc_eqiv_id.variables['lon'][:] = type_lon
                nc_eqiv_var.setncattr("units", origVars[v].type_units)
                nc_eqiv_var.setncattr("desc", "type55 data in NCEP equivalent formatting")

                if daily == 'daily':
                    log.info("Creating Daily Averages")
                    #  create daily averaged version of data
                    step=int(24/origVars[v].type_time_res)
                else:
                    log.info("Creating 6-hour Averages")
                    step=int(6/origVars[v].type_time_res)

                print(step)

                timePoint = 0

                halfStep = int(step/2)
                midOff = origVars[v].type_steps_from_midnight
                for f in range(1,len(typeFiles)-1): # take end ones off as prev/next years
                    # read in file
                    # print("b: ",typeFiles[f])
                    nc_type_id=Dataset(typeFiles[f],"r")
                    type_var = nc_type_id.variables[origVars[v].type_var][:].data * origVars[v].type_factor

                    if len(type_var.shape) > 3:
                        type_var = type_var[:,0,:,:]
                    # read in previous file 
                    nc_type_prev_id=Dataset(typeFiles[f-1],"r")
                    type_var_prev = nc_type_prev_id.variables[origVars[v].type_var][:].data * origVars[v].type_factor
                    # print(type_var_prev.shape,origVars[v].type_var)
                    if len(type_var_prev.shape) > 3:
                        type_var_prev = type_var_prev[:,0,:,:]
                    # type_var_prev_day = type_var_prev[-3:-1,:,:]
                    type_var_prev_day = type_var_prev[(-1-halfStep-midOff):-1,:,:]

                    # # read in next file 
                    # nc_type_next_id=Dataset(typeFiles[f+1],"r")
                    # type_var_next = nc_type_next_id.variables[origVars[v].type_var][:] 
                    # type_var_next_day = type_var_next[0:2,:,:,:]

                    # add slices to end
                    # type_var = np.concatenate((type_var_prev_day,type_var,type_var_next_day),axis=0)
                    type_var = np.concatenate((type_var_prev_day,type_var),axis=0)

                    # simplified mean
                    for t in range(0,type_var.shape[0]-step,step):
                        # print("d: ",t,timePoint,step,typeFiles[f],type_var.shape[0],halfStep)
                        type_mean = np.mean( type_var[t:t+step,:,:], axis=0 )

                        # for ii in range(t,t+step):
                        #     print(ii,t,step,type_var[ii,99,99], type_mean[99,99])
                        # input("------")

                        if latSwap == 0:
                            nc_eqiv_id.variables[varName][timePoint,:,:] = type_mean
                        else:
                            nc_eqiv_id.variables[varName][timePoint,:,:] = np.flip(type_mean,axis=0)

                        timePoint = timePoint+1

                    # for t in range(halfStep,type_var.shape[0]-halfStep,step):
                        # print("d: ",t,timePoint,step,typeFiles[f],type_var.shape[0],halfStep)
                        # # type_midday_mean = np.mean(np.stack((type_var[t-2,:,:],type_var[t+2,:,:]),axis=0),axis=0)
                        # print(t-halfStep,t+halfStep,t-(halfStep-1))
                        # # print(type_var.shape)

    
                        # type_midday_mean = np.mean(np.stack((type_var[t-halfStep,:,:],type_var[t+halfStep,:,:]),axis=0),axis=0)
                        # print(type_var[t-halfStep,200,500],type_var[t+halfStep,200,500])
                        # type_midday_mean = type_midday_mean.reshape((1, type_midday_mean.shape[0], type_midday_mean.shape[1]))

                        # type_daily = np.mean(np.concatenate( (type_var[t-(halfStep):t+(halfStep),:,:], type_midday_mean) ,axis=0),axis=0)
                        # # type_daily = type_midday_mean

                        # print(type_var[t-(halfStep):t+(halfStep),200,500])

                        # for ii in range(t-halfStep,t+halfStep):
                        #     print(ii,t,halfStep,type_var[ii,200,500], type_midday_mean[0,200,500], type_daily[200,500])
                        # input("------")
                        
                        # if latSwap == 0:
                        #     nc_eqiv_id.variables[varName][timePoint,:,:] = type_daily
                        # else:
                        #     nc_eqiv_id.variables[varName][timePoint,:,:] = np.flip(type_daily,axis=0)
                        # timePoint = timePoint+1

                # create time dimension
                startOfYear = datetime.datetime(int(year),1,1,0,0,0)
                timePoints = []

                if daily == 'daily':
                    for d in range(0,timePoint):
                        dayOfYear = startOfYear + datetime.timedelta(days=d)
                        hoursSince = (dayOfYear-originDate).days * 24
                        timePoints.append(hoursSince)
                else:
                    for t in range(0,timePoint):
                        dayOfYear = startOfYear + datetime.timedelta(seconds=6*t*60*60)
                        hoursSince = (dayOfYear-originDate).days * 24
                        timePoints.append(hoursSince)


                nc_eqiv_id.variables['time'][:] = timePoints


                # print("e: ",typeFiles[f])

                    # # The original NCEP data times were midnight so average with this at centre
                    # for t in range(0,type_var.shape[0],4):
                    #     # takes midday values at ends of the day and take averages
                    #     t_mid = t+2
                    #     t_mid_prev = t-2
                    #     if t_mid_prev > 0 and t_mid < type_var.shape[0]:
                    #             type_midday_mean = np.mean(np.stack((type_var[t_mid_prev,0,:,:],type_var[t_mid,0,:,:]),axis=0),axis=0)
                    #     else:
                    #         if t_mid_prev < 0:
                    #             type_midday_mean = np.mean(np.stack((type_var[t_mid_prev,0,:,:],type_var[t_mid,0,:,:]),axis=0),axis=0)

                    #     # takes midnight values at ends of the day and takes averages
                    #     if t+4 < type_var.shape[0]:
                    #         type_midnight_mean = np.mean(np.stack((type_var[t,0,:,:],type_var[t+4,0,:,:]),axis=0),axis=0)
                    #     else:
                    #         if f<len(typeFiles):
                    #             # take first value from next file
                    #             type_midnight_mean = np.mean(np.stack((type_var[t,0,:,:],type_var_next[0,0,:,:]),axis=0),axis=0)
                    #         else:
                    #             # instead of reading next year in, just take one midnight value
                    #             type_midnight_mean = type_var[t,0,:,:]



                    #     type_midnight_mean = type_midnight_mean.reshape((1, type_midnight_mean.shape[0], type_midnight_mean.shape[1]))
                    #     type_daily = np.mean(np.concatenate( (type_var[t+1:t+4,0,:,:], type_midnight_mean) ,axis=0),axis=0)
                    
                #         nc_eqiv_id.variables[origVars[v].varId][day,:,:] = type_daily
                #         day = day+1

                # input("")
                nc_eqiv_id.close()


# ----------------- 00240 pre-process ERA5 data --------------------------------------
# ERA5 data is yearly with hourly intervals
if type == 'era5':
    # time:units = "hours since 1900-01-01 00:00:0.0" ;
    originDate = datetime.datetime(1900,1,1,0,0,0)
    #go through each variable
    for v in range(0,len(origVars)):
        # if the variable name isn't empty
        if origVars[v].type_tag != "":
            # create equivalent file to the NCEP input
            equivFile = origVars[v].fileName
            # if file exists, skip
            if os.path.isfile(equivFile):
                print(equivFile, " exists")
            else:
                # get list of type files
                # uses previous year to get dimensions (?)
                #if len(sorted(glob.glob(origVars[v].type_tag+str(int(year)-1)+"*.nc"))) > 0:
                #    eraPrevYear = [sorted(glob.glob(origVars[v].type_tag+str(int(year)-1)+"*.nc"))[-1]]
                #else:
                #    eraPrevYear = [sorted(glob.glob(origVars[v].type_tag+year+"*.nc"))[0]]
 
                #eraFiles = sorted(glob.glob(origVars[v].type_tag+year+"*.nc"))
                               
                #if len(sorted(glob.glob(origVars[v].type_tag+str(int(year)+1)+"*.nc"))) > 0:
                #    eraNextMonth = [sorted(glob.glob(origVars[v].type_tag+str(int(year)+1)+"*.nc"))[0]]
                #else:
                #    eraNextMonth = [sorted(glob.glob(origVars[v].type_tag+year+"*.nc"))[-1]]

                #eraPrevMonth.extend(eraFiles)
                #eraPrevMonth.extend(eraNextMonth)
                #eraFiles = eraPrevMonth
                #print(eraFiles)

                # open first to get dimensions
                nc_type_id=Dataset(origVars[v].type_tag+year+".nc","r")
                era_lat = nc_type_id.variables['latitude'][:] 
                era_lon = nc_type_id.variables['longitude'][:] 
                era_time = nc_type_id.variables['time'][:] 

                nc_eqiv_id = Dataset(equivFile, 'w', format=nc_type_id.data_model)
                # nc_eqiv_id = Dataset(equivFile, 'w', format='NETCDF3_CLASSIC')

                log.info("Creating equivalent to NCEP for ERA5: "+equivFile+ " "+str(nc_eqiv_id.data_model))
                nc_eqiv_id.createDimension("lat", len(era_lat)) 
                nc_eqiv_id.createDimension("lon", len(era_lon)) 
                nc_eqiv_id.createDimension("time", None)
                nc_eqiv_lat = nc_eqiv_id.createVariable("lat", nc_type_id.variables["latitude"].dtype, ("lat"))
                nc_eqiv_lon = nc_eqiv_id.createVariable("lon", nc_type_id.variables["longitude"].dtype, ("lon"))
                nc_eqiv_time = nc_eqiv_id.createVariable("time", nc_type_id.variables["time"].dtype, ("time"))
                # varName = origVars[v].fileName.split('.')[0]
                
                varName = origVars[v].varId
                nc_eqiv_var = nc_eqiv_id.createVariable(varName, nc_type_id.variables[origVars[v].type_var].dtype, ("time", "lat", "lon"))
                    
                # do a check that type_lat is the right way up for what we want
                latSwap = 0
                if era_lat[-1] > era_lat[0]:
                    log.info("Swapping latitudes as upside down")
                    era_lat = era_lat * -1
                    latSwap = 1
                
                nc_eqiv_id.variables['lat'][:] = era_lat
                nc_eqiv_id.variables['lon'][:] = era_lon
                nc_eqiv_var.setncattr("units", origVars[v].type_units)
                nc_eqiv_var.setncattr("desc", "ERA5 data in NCEP equivalent formatting")

                if daily == 'daily':
                    log.info("Creating Daily Averages")
                    #  create daily averaged version of data
                    step=int(24/origVars[v].type_time_res)
                else:
                    log.info("Creating 6-hour Averages")
                    step=int(6/origVars[v].type_time_res)

                print(step)

                timePoint = 0
                halfStep = int(step/2)
                midOff = origVars[v].type_steps_from_midnight

                # create time dimension
                startOfYear = datetime.datetime(int(year),1,1,0,0,0)
                timePoints = []

                if daily == 'daily':
                    for d in range(0,timePoint):
                        dayOfYear = startOfYear + datetime.timedelta(days=d)
                        hoursSince = (dayOfYear-originDate).days * 24
                        timePoints.append(hoursSince)
                else:
                    for t in range(0,timePoint):
                        dayOfYear = startOfYear + datetime.timedelta(seconds=6*t*60*60)
                        hoursSince = (dayOfYear-originDate).days * 24
                        timePoints.append(hoursSince)


                nc_eqiv_id.variables['time'][:] = timePoints


                # print("e: ",typeFiles[f])

                    # # The original NCEP data times were midnight so average with this at centre
                    # for t in range(0,type_var.shape[0],4):
                    #     # takes midday values at ends of the day and take averages
                    #     t_mid = t+2
                    #     t_mid_prev = t-2
                    #     if t_mid_prev > 0 and t_mid < type_var.shape[0]:
                    #             type_midday_mean = np.mean(np.stack((type_var[t_mid_prev,0,:,:],type_var[t_mid,0,:,:]),axis=0),axis=0)
                    #     else:
                    #         if t_mid_prev < 0:
                    #             type_midday_mean = np.mean(np.stack((type_var[t_mid_prev,0,:,:],type_var[t_mid,0,:,:]),axis=0),axis=0)

                    #     # takes midnight values at ends of the day and takes averages
                    #     if t+4 < type_var.shape[0]:
                    #         type_midnight_mean = np.mean(np.stack((type_var[t,0,:,:],type_var[t+4,0,:,:]),axis=0),axis=0)
                    #     else:
                    #         if f<len(typeFiles):
                    #             # take first value from next file
                    #             type_midnight_mean = np.mean(np.stack((type_var[t,0,:,:],type_var_next[0,0,:,:]),axis=0),axis=0)
                    #         else:
                    #             # instead of reading next year in, just take one midnight value
                    #             type_midnight_mean = type_var[t,0,:,:]



                    #     type_midnight_mean = type_midnight_mean.reshape((1, type_midnight_mean.shape[0], type_midnight_mean.shape[1]))
                    #     type_daily = np.mean(np.concatenate( (type_var[t+1:t+4,0,:,:], type_midnight_mean) ,axis=0),axis=0)
                    
                #         nc_eqiv_id.variables[origVars[v].varId][day,:,:] = type_daily
                #         day = day+1

                # input("")
                nc_eqiv_id.close()
   

# ----------------- 00300 OPEN orig FILES --------------------------------------

log.info("Running script to convert orig files to Forcing files for NEMO models")
log.info("Required files: ")

for v in range(0,len(origVars)):
    if os.path.isfile(origVars[v].fileName):
        log.info(str(v)+": "+origVars[v].fileName+" found")
    else:
        log.error(origVars[v].fileName+" does not exist")

log.info("opening input orig files")

for v in range(0,len(origVars)):
    origVars[v].nc_id = Dataset(origVars[v].fileName,"r")
    log.info(str(v)+": "+origVars[v].fileName+" opened. Type: "+origVars[v].nc_id.data_model)

# ----------------- 00400 OPEN MESH_MASK FILES --------------------------------------
nc_mesh_id=Dataset("meshmask.nc","r")
log.info("Mesh opened. Type: "+nc_mesh_id.data_model)

# positions of output grid pixels
nc_lat = nc_mesh_id.variables['nav_lat'][:] 
nc_lon = nc_mesh_id.variables['nav_lon'][:]
nc_lon_orig = nc_lon.copy() # copy taken as adjusted for 0-360 for interpolation and then original written out

# lamda (latidude) and phi (longitude) of points on ORCA grid. nav_lon/lat is equivalent to glamt and gphit
# postions therefore defined at u-, v- and f-points
glamu = nc_mesh_id.variables['glamu'][:] 
gphiu = nc_mesh_id.variables['gphiu'][:] 
glamv = nc_mesh_id.variables['glamv'][:] 
gphiv = nc_mesh_id.variables['gphiv'][:] 
glamf = nc_mesh_id.variables['glamf'][:] 
gphif = nc_mesh_id.variables['gphif'][:] 

# convert to 2D files
glamu = glamu[0,:,:]
gphiu = gphiu[0,:,:]
glamv = glamv[0,:,:]
gphiv = gphiv[0,:,:]
glamf = glamf[0,:,:]
gphif = gphif[0,:,:]

# change longitude to be 0 to 360, not -180 to 180, ouput will be -180 to 180
# copy maintained for output purposes
for iy, ix in np.ndindex(nc_lon.shape):
    if nc_lon[iy,ix] < 0:
        nc_lon[iy,ix] = 360 + nc_lon[iy,ix]
    if glamu[iy,ix] < 0:
        glamu[iy,ix] = 360 + glamu[iy,ix]
    if glamv[iy,ix] < 0:
        glamv[iy,ix] = 360 + glamv[iy,ix]
    if glamf[iy,ix] < 0:
        glamf[iy,ix] = 360 + glamf[iy,ix]

#  get dtypes so can close mesh file
nav_dtype = nc_mesh_id.variables["nav_lon"].dtype
nc_mesh_id.close()

# ----------------- 00400 CREATE OUTPUT FILES (BULK) --------------------------------------
if type == 'jra55' or 'era5':
    outputFile = "bulk_"+str(year)+"_"+str(varNum)+"_"+str(type)+"_"+str(daily)+".nc"
else:
    outputFile = "bulk_"+str(year)+".nc"
nc_out_id = Dataset(outputFile, 'w', format='NETCDF4_CLASSIC')
log.info("Creating: "+outputFile)

out_xDim = nc_lat.shape[1]
out_yDim = nc_lat.shape[0]
log.info("Dimensions: "+str(out_yDim)+","+str(out_xDim))

# Set the available dimensions determined from meshmask file
nc_out_id.createDimension("x", out_xDim) 
nc_out_id.createDimension("y", out_yDim) 
nc_out_id.createDimension("deptht", 1)
nc_out_id.createDimension("time_counter", None)

# Set the nav-lon and lat variables
nc_v_nav_lon = nc_out_id.createVariable("nav_lon", nav_dtype, ("y","x")) #nc_mesh_id.variables["nav_lon"].dtype, ("y","x"))
nc_v_nav_lat = nc_out_id.createVariable("nav_lat", nav_dtype, ("y","x")) #nc_mesh_id.variables["nav_lat"].dtype, ("y","x"))
nc_out_id.variables['nav_lon'][:] = nc_lon_orig
nc_out_id.variables['nav_lat'][:] = nc_lat

# the limits are duplicated from original orig forcing runs
nc_v_nav_lon.setncattr("units", "degrees_east")
nc_v_nav_lon.setncattr("valid_min", np.array(-179.7507,'f'))
nc_v_nav_lon.setncattr("valid_max", np.array(180.,'f'))
nc_v_nav_lon.setncattr("long_name", "Longitude")
nc_v_nav_lat.setncattr("units", "degrees_east")
nc_v_nav_lat.setncattr("valid_min", np.array(-78.19058,'f'))
nc_v_nav_lat.setncattr("valid_max", np.array(89.6139,'f'))
nc_v_nav_lat.setncattr("long_name", "Latitude")

# Set the depth
nc_v_d = nc_out_id.createVariable("deptht", nav_dtype, ("y","x")) #, nc_mesh_id.variables["nav_lev"].dtype, ("deptht"))
nc_v_d.setncattr("units", "m")
nc_v_d.setncattr("valid_min", np.array(4.999938,'f'))
nc_v_d.setncattr("valid_max", np.array(5250.227,'f'))
nc_v_d.setncattr("title", "deptht")
nc_v_d.setncattr("long_name", "Vertical T levels")
nc_out_id.variables['deptht'][:] = 4.999938

# Set the time variable - copy from first file in list
nc_v_time = nc_out_id.createVariable("time_counter", 'f', ('time_counter'))
nc_v_time.setncattr("units", "seconds since    0-01-01 00:00:00")
nc_v_time.setncattr("calendar", "gregorian")
nc_v_time.setncattr("title", "time")
nc_v_time.setncattr("long_name", "   0-JAN-01 00:00:00")

# this will duplicate the time variable and so will compensate for leap-years

if len(origVars) > 0:
    times = origVars[0].nc_id.variables['time'][:]
    # convert to hours since start of the year
    times_since_startOfYear = []
    times = times * (60*60) # convert to seconds from hours
    for t in range(0,len(times)):
        times_since_startOfYear.append(times[t] - times[0])

# set the sizes of the output dimensions (remember time may be 366 due to leap year)
    out_dims = [len(times),nc_out_id.dimensions["y"].size,nc_out_id.dimensions["x"].size]
    log.info("Output Dimensions [y,x]: "+str(out_dims))

# ----------------- 00400 CREATE OUTPUT FILES (WIND STRESS VECTOR FILES) --------------------------------------
if createOnlyMod == False and (varNum == 7 or varNum == -1):
    # there is no depth component to this file and the times are set as hours from an epoch, not seconds
    # from start of year as previously in the bulk file
    outputFile_sx = "taux_1d_"+str(year)+"_"+str(daily)+".nc"
    outputFile_sy = "tauy_1d_"+str(year)+"_"+str(daily)+".nc"
    nc_out_sx_id = Dataset(outputFile_sx, 'w', format='NETCDF4_CLASSIC')
    nc_out_sy_id = Dataset(outputFile_sy, 'w', format='NETCDF4_CLASSIC')
    log.info("Creating: "+outputFile_sx+" "+outputFile_sy )

    # Set the available dimensions
    nc_out_sx_id.createDimension("x", out_xDim)
    nc_out_sx_id.createDimension("y", out_yDim)
    nc_out_sx_id.createDimension("time_counter", None)
    nc_out_sy_id.createDimension("x", out_xDim)
    nc_out_sy_id.createDimension("y", out_yDim)
    nc_out_sy_id.createDimension("time_counter", None)

    # Set the nav-lon and lat variables
    nc_v_nav_sx_lon = nc_out_sx_id.createVariable("nav_lon", nav_dtype, ("y","x")) #, nc_mesh_id.variables["nav_lon"].dtype, ("y","x"))
    nc_v_nav_sx_lat = nc_out_sx_id.createVariable("nav_lat", nav_dtype, ("y","x")) #, nc_mesh_id.variables["nav_lat"].dtype, ("y","x"))
    nc_v_nav_sy_lon = nc_out_sy_id.createVariable("nav_lon", nav_dtype, ("y","x")) #, nc_mesh_id.variables["nav_lon"].dtype, ("y","x"))
    nc_v_nav_sy_lat = nc_out_sy_id.createVariable("nav_lat", nav_dtype, ("y","x")) #, nc_mesh_id.variables["nav_lat"].dtype, ("y","x"))
    nc_out_sx_id.variables['nav_lon'][:] = nc_lon_orig
    nc_out_sx_id.variables['nav_lat'][:] = nc_lat
    nc_out_sy_id.variables['nav_lon'][:] = nc_lon_orig
    nc_out_sy_id.variables['nav_lat'][:] = nc_lat

    nc_v_nav_sx_lon.setncattr("units", "degrees_east")
    nc_v_nav_sx_lon.setncattr("valid_min", np.array(-179.7507,'f'))
    nc_v_nav_sx_lon.setncattr("valid_max", np.array(180.,'f'))
    nc_v_nav_sx_lon.setncattr("long_name", "Longitude")
    nc_v_nav_sy_lon.setncattr("units", "degrees_east")
    nc_v_nav_sy_lon.setncattr("valid_min", np.array(-179.7507,'f'))
    nc_v_nav_sy_lon.setncattr("valid_max", np.array(180.,'f'))
    nc_v_nav_sy_lon.setncattr("long_name", "Longitude")

    nc_v_nav_sx_lat.setncattr("units", "degrees_east")
    nc_v_nav_sx_lat.setncattr("valid_min", np.array(-78.19058,'f'))
    nc_v_nav_sx_lat.setncattr("valid_max", np.array(89.6139,'f'))
    nc_v_nav_sx_lat.setncattr("long_name", "Latitude")
    nc_v_nav_sy_lat.setncattr("units", "degrees_east")
    nc_v_nav_sy_lat.setncattr("valid_min", np.array(-78.19058,'f'))
    nc_v_nav_sy_lat.setncattr("valid_max", np.array(89.6139,'f'))
    nc_v_nav_sy_lat.setncattr("long_name", "Latitude")

    # Set the time variable - copy from first file in list
    nc_v_sx_time = nc_out_sx_id.createVariable("time_counter", 'f', ('time_counter'))
    nc_v_sx_time.setncattr("units", "seconds since    0-01-01 00:00:00")
    nc_v_sx_time.setncattr("calendar", "gregorian")
    nc_v_sx_time.setncattr("title", "time")
    nc_v_sx_time.setncattr("long_name", "   0-JAN-01 00:00:00")
    nc_v_sy_time = nc_out_sy_id.createVariable("time_counter", 'f', ('time_counter'))
    nc_v_sy_time.setncattr("units", "seconds since    0-01-01 00:00:00")
    nc_v_sy_time.setncattr("calendar", "gregorian")
    nc_v_sy_time.setncattr("title", "time")
    nc_v_sy_time.setncattr("long_name", "   0-JAN-01 00:00:00")

# if len(origVars) > 0:
#     times = origVars[0].nc_id.variables['time'][:]

#     # comment out as incrementing time_counter causing mem problems on hi res data
#     # nc_out_sx_id.variables['time_counter'][:] = times
#     # nc_out_sy_id.variables['time_counter'][:] = times

#     # set the sizes of the output dimensions (remember time may be 366 due to leap year)
#     out_dims = [len(times),nc_out_sx_id.dimensions["y"].size,nc_out_sx_id.dimensions["x"].size]
#     log.info("Output Dimensions [y,x]: "+str(out_dims))

# ---------------------- 00400 PROCESS EACH VARAIBLE ---------------------------

v_iter =  iter(range(0,len(origVars))) # Use an iterator as need to skip some files

for v in v_iter: 
    log.info("Parameter: "+origVars[v].varId)

    if origVars[v].inBulk == True: 

        # ---------------------- 00500 FILL BULK FILE ---------------------------
        nc_vars = [var for var in origVars[v].nc_id.variables] 
        in_varName = origVars[v].inVarId # nc_vars[3] # the variable from the input file

        nc_v = nc_out_id.createVariable(origVars[v].varId, origVars[v].nc_id.variables[in_varName].dtype, ('time_counter', 'y', 'x'))
        
        if type != 'jra55' and type != 'era5':
            nc_v.setncattr("long_name", origVars[v].nc_id.variables[in_varName].getncattr("long_name"))
            # set the units explicitly if defined, copy if not
            if len(origVars[v].units) > 0:
                nc_v.setncattr("units", origVars[v].units)
            else:
                nc_v.setncattr("units", origVars[v].nc_id.variables[in_varName].getncattr("units"))
        else:
            nc_v.setncattr("units", origVars[v].units)
            nc_v.setncattr("long_name", in_varName)

        nc_v.setncattr("short_name", in_varName)
        nc_v.setncattr("axis", "TYX")
        nc_v.setncattr("interval_operation",  np.array(86400.,'f'))
        nc_v.setncattr("interval_write", np.array(86400.,'f'))
        nc_v.setncattr("associate", "time_counter nav_lat nav_lon")
        nc_v.setncattr("missing_value", np.array(1E20,'f'))

        # ----------------- 00500 GET INPUT DATA  --------------------------------------

        if (origVars[v].directConversion == False and origVars[v].conversionFunc == None):
            # the direct conversion flag is None (null) so that preprocessing is required on the input data
            log.info("Preprocessing: "+origVars[v].varId)
            data_in = origVars[v].nc_id.variables[in_varName][:]

            # ----------------- 00600 PRE-PROCESS WIND SPEED DATA  --------------------------------------
            if origVars[v].varId == "wspd" or origVars[v].varId == "uwnd" or origVars[v].varId == "vwnd":
                if origVars[v].varId == "wspd":
                    u_data_in = origVars[v].nc_id.variables['uwnd'][:]
                    v_data_in = origVars[v+1].nc_id.variables['vwnd'][:]
                    data_in = np.sqrt(np.square(u_data_in) + np.square(v_data_in))
                    next(v_iter, None) # increment the variable on one to avoid repeating, parms still same
                else:
                    u_v = -1
                    v_v = -1
                    for vv in range(0,len(origVars)):
                        if origVars[vv].varId == "uwnd":
                            u_v = vv
                        if origVars[vv].varId == "vwnd":
                            v_v = vv

                    # the wind speed is calculated from two adjacent files in the list
                    u_data_in = origVars[u_v].nc_id.variables['uwnd'][:]
                    v_data_in = origVars[v_v].nc_id.variables['vwnd'][:]
                    data_in = np.sqrt(np.square(u_data_in) + np.square(v_data_in))


                log.info("Correcting for coherance at North Pole")
                data_in_lon = origVars[v].nc_id.variables['lon'][:] 
                data_in_lat = origVars[v].nc_id.variables['lat'][:]
                u_data_in, v_data_in = NorthPoleWindVectorCoherance(u_data_in, v_data_in, data_in_lon, data_in_lat, False)

                if origVars[v].varId == "uwnd":
                    data_in = u_data_in

                if origVars[v].varId == "vwnd":
                    data_in = v_data_in

            # ----------------- 00700 PRE-PROCESS HUMIDITY DATA  --------------------------------------
            if origVars[v].varId == "humidity" and type != "era5":

                # Loop over time so as not to overload mem
                tLim = len(origVars[v].nc_id.dimensions['time'])
                data_in = np.zeros( (   len(origVars[v].nc_id.dimensions['time']),
                                        len(origVars[v].nc_id.dimensions['lat']), 
                                        len(origVars[v].nc_id.dimensions['lon']) ) )

                for t in range(0,tLim):
                    # need to calculate humidity as follows:
                    #   (sp hum*pressure)/(mol wt water/mol wt air +.378*sp hum)/dew point 
                    shum_data_in = origVars[v].nc_id.variables['shum'][t,:,:]
                    # read in again the air and pressure components
                    air_T_data_in = origVars[v-1].nc_id.variables['air'][t,:,:]
                    if origVars[v-1].nc_id.variables['air'].getncattr('units') =='degK':
                        air_T_data_in = air_T_data_in - 273.15 # convert to deg C if needed
                    pres_data_in = origVars[v-1].nc_id.variables['pres'][t,:,:]

                    # perform calculation on input data (i.e. pre interpolation to avoid error)
                    a = pres_data_in * shum_data_in
                    b = 0.622 + 0.377 * shum_data_in

                    # this was causing memory issues so splitting it in to stripes.
                    # c = 611.0 *  np.power(10, (7.5 * air_T_data_in / (air_T_data_in + 237.29)) )
                    # data_in[t,:,:] = a/b/c

                    c3 = 7.5 * air_T_data_in
                    c4 = air_T_data_in + 237.29
                    xlim = data_in.shape[2]
                    ylim = data_in.shape[1]
                    for ty in range(0,ylim):
                        c_temp = 611.0 *  np.power(10, c3[ty,:]/c4[ty,:] )
                        data_in[t,ty,:] = a[ty,:] / b[ty,:] / c_temp
            elif origVars[v].varId == "humidity" and type == "era5":
                # Loop over time so as not to overload mem
                tLim = len(origVars[v].nc_id.dimensions['time'])
                data_in = np.zeros( (   len(origVars[v].nc_id.dimensions['time']),
                                        len(origVars[v].nc_id.dimensions['lat']), 
                                        len(origVars[v].nc_id.dimensions['lon']) ) )

                for t in range(0,tLim):
                    # need to calculate humidity as follows: 
                    dew_temp = origVars[v].nc_id.variables['shum'][t,:,:]
                    sur_pres = origVars[v+1].nc_id.variables['pres'][t,:,:]
                    # calculate vapour pressure (e)
                    e = 0.6112*math.exp(17.67*(dew_temp-273.15)/(243.5 + (dew_temp-273.15)))
                        
                    # calculate mixing ratio (r)
                    r = (0.622*e*1000)/(sur_pres - (e*1000))
                    
                    # calculate specific humidity (q)
                    q = r/(1 + r)
        else:
        # ----------------- 00800 DATA REQUIRES NO PRE-PROCESSING  -------------------------------------
            data_in = origVars[v].nc_id.variables[in_varName][:]


        # lon and lat data from input file
        data_in_lon = origVars[v].nc_id.variables['lon'][:] 
        data_in_lat = origVars[v].nc_id.variables['lat'][:]

        # Set output array
        data_out = np.zeros(out_dims,dtype=float)

        # dimensions on input data 
        in_lon_dim = len(data_in_lon)
        in_lat_dim = len(data_in_lat)
        xLim = data_out.shape[2]
        yLim = data_out.shape[1]
        tLim = data_out.shape[0]
        
        # tLim = 8

        start = datetime.datetime.now() # time stamp if needed for optimisation
        lats = data_in_lat
        lons = data_in_lon

        # change longitude to be 0 to 360, not -180 to 180-
        for iy, ix in np.ndindex(nc_lon.shape):
            if nc_lon[iy,ix] < 0:
                nc_lon[iy,ix] = 360 + nc_lon[iy,ix]

        # ----------------- 00800 INTERPOLATE THE DATA ONTO THE NEW GRID --------------------------------------

        if type != "jra55" and type != "era5":
        # see Interpolate function for details
            log.info("Interpolating: "+origVars[v].varId)
            log.info("Times: "+str(data_in.shape[0]))
            data_out = Interpolate(out_dims, xLim, yLim, tLim, data_in, lons, lats, nc_lat, nc_lon,
                                    origVars[v].directConversion, origVars[v].conversionFunc, origVars[v].limits, type, "T")
            # write out data
            nc_out_id.variables[origVars[v].varId][:] = data_out
        else:
            log.info("Interpolating perDay: "+origVars[v].varId)
            
            for t in range(0,tLim): # t is lon
                data_out = Interpolate_perDay(out_dims, xLim, yLim, tLim, data_in[t,:,:], lons, lats, nc_lat, nc_lon,
                                            origVars[v].directConversion, origVars[v].conversionFunc, origVars[v].limits, type, "T")
                # write out data
                nc_out_id.variables[origVars[v].varId][t,:,:] = data_out
                # nc_out_id.variables['time_counter'][t] = np.array(times_since_startOfYear)[t]
                log.info("Time point: "+str(t)+" for "+origVars[v].varId)

            
    else:

# ----------------- 00900 THE WIND STRESS FILES  --------------------------------------
# this data is not to be put in the bulk file but seperate TAU (wind stress) files

        if origVars[v].varId == "uflx":
            # ----------------- 01000 WIND STRESS PRE-PROCESSING  --------------------------------------
            # the wind stress is calculated from two adjacent files in the list
            u_data_in = origVars[v].nc_id.variables['uflx'][:]
            v_data_in = origVars[v+1].nc_id.variables['vflx'][:]
            data_in_lon = origVars[v].nc_id.variables['lon'][:] 
            data_in_lat = origVars[v].nc_id.variables['lat'][:]
            next(v_iter, None) # increment the variable on to avoid repeating, parms still same

            log.info("Correcting for coherance at North Pole")
            u_data_in, v_data_in = NorthPoleWindVectorCoherance(u_data_in, v_data_in, data_in_lon, data_in_lat, False)

        # get name of variable to write
        nc_vars = [var for var in origVars[v].nc_id.variables] 
        in_varName_x = nc_vars[3] # the variable from the input file
        nc_vars = [var for var in origVars[v+1].nc_id.variables] 
        in_varName_y = nc_vars[3] # the variable from the input file

        # create variables
        nc_vx = nc_out_sx_id.createVariable(origVars[v].varId, origVars[v].nc_id.variables[in_varName_x].dtype, ('time_counter', 'y', 'x'))
        nc_vy = nc_out_sy_id.createVariable(origVars[v+1].varId, origVars[v+1].nc_id.variables[in_varName_y].dtype, ('time_counter', 'y', 'x'))
        if type != 'jra55' and type != 'era5':
            nc_vx.setncattr("long_name", origVars[v].nc_id.variables[in_varName_x].getncattr("long_name"))
            nc_vy.setncattr("long_name", origVars[v+1].nc_id.variables[in_varName_y].getncattr("long_name"))
        if len(origVars[v].units) > 0:
            nc_vx.setncattr("units", origVars[v].units)
            nc_vy.setncattr("units", origVars[v+1].units)
        else:
            nc_vx.setncattr("units", origVars[v].nc_id.variables[in_varName_x].getncattr("units"))
            nc_vy.setncattr("units", origVars[v+1].nc_id.variables[in_varName_y].getncattr("units"))
        nc_vx.setncattr("short_name", in_varName_x)
        nc_vx.setncattr("axis", "TYX")
        nc_vx.setncattr("interval_operation",  np.array(86400.,'f'))
        nc_vx.setncattr("interval_write", np.array(86400.,'f'))
        nc_vx.setncattr("associate", "time_counter nav_lat nav_lon")
        nc_vx.setncattr("missing_value", np.array(1E20,'f'))
        nc_vy.setncattr("short_name", in_varName_y)
        nc_vy.setncattr("axis", "TYX")
        nc_vy.setncattr("interval_operation",  np.array(86400.,'f'))
        nc_vy.setncattr("interval_write", np.array(86400.,'f'))
        nc_vy.setncattr("associate", "time_counter nav_lat nav_lon")
        nc_vy.setncattr("missing_value", np.array(1E20,'f'))

        # Set output array
        data_out = np.zeros(out_dims,dtype=float)

        # dimensions on input data 
        in_lon_dim = len(data_in_lon)
        in_lat_dim = len(data_in_lat)
        xLim = data_out.shape[2]
        yLim = data_out.shape[1]
        tLim = data_out.shape[0]

        # tLim = 3

        start = datetime.datetime.now() # time stamp if needed for optimisation
        lats = data_in_lat
        lons = data_in_lon

        # change longitude to be 0 to 360, not -180 to 180
        for iy, ix in np.ndindex(nc_lon.shape):
            if nc_lon[iy,ix] < 0:
                nc_lon[iy,ix] = 360 + nc_lon[iy,ix]

        # we need to set the arrays for the U and V input values interpolated to the repective 
        # u and v points on the output grid. Each output grid position will have a contribution
        # from the input U and V data.
        u_onto_u = np.zeros((out_dims[0],out_dims[1],out_dims[2]),dtype=float)
        u_onto_v = np.zeros((out_dims[0],out_dims[1],out_dims[2]),dtype=float)
        v_onto_v = np.zeros((out_dims[0],out_dims[1],out_dims[2]),dtype=float)
        v_onto_u = np.zeros((out_dims[0],out_dims[1],out_dims[2]),dtype=float)

       # ----------------- 01200 APPLYING COMPONENTS OF VECTOR TO GRID --------------------------------------       

        # split up the calc into strips to prevent memory problems
        xlim = u_onto_u.shape[2]
        ylim = u_onto_u.shape[1]
        x_u = np.zeros((ylim,xlim))
        y_u = np.zeros((ylim,xlim))
        m_u = np.zeros((ylim,xlim))
        x_v = np.zeros((ylim,xlim))
        y_v = np.zeros((ylim,xlim))
        m_v = np.zeros((ylim,xlim))
        x_fu = np.zeros((ylim,xlim))
        y_fu = np.zeros((ylim,xlim))
        m_fu = np.zeros((ylim,xlim))
        x_fv = np.zeros((ylim,xlim))
        y_fv = np.zeros((ylim,xlim))
        m_fv = np.zeros((ylim,xlim))
        gsin_u = np.zeros((ylim,xlim))
        gcos_u = np.zeros((ylim,xlim))
        gsin_v = np.zeros((ylim,xlim))
        gcos_v = np.zeros((ylim,xlim))
        u_out = np.zeros((ylim,xlim))
        v_out = np.zeros((ylim,xlim))

        for ty in range(0,ylim):
            print(ty)
            # the following is translated explicitly from the original IDL code (written 1996)
            # copies precede each section
            # ;     ... north pole direction & modulous (at u-point)
            #       zxnpu = 0. - fsxspp( glamu, gphiu )
            #       zynpu = 0. - fsyspp( glamu, gphiu )
            #       znnpu = zxnpu*zxnpu + zynpu*zynpu
            x_u[ty,:] = 0. - fx(glamu[ty,:],gphiu[ty,:])
            y_u[ty,:] = 0. - fy(glamu[ty,:],gphiu[ty,:])
            m_u[ty,:] = x_u[ty,:]*x_u[ty,:] + y_u[ty,:]*y_u[ty,:]

            # ;     ... north pole direction & modulous (at v-point)
            #       zxnpv = 0. - fsxspp( glamv, gphiv )
            #       zynpv = 0. - fsyspp( glamv, gphiv )
            #       znnpv = zxnpv*zxnpv + zynpv*zynpv
            x_v[ty,:] = 0. - fx(glamv[ty,:],gphiv[ty,:])
            y_v[ty,:] = 0. - fy(glamv[ty,:],gphiv[ty,:])
            m_v[ty,:] = x_v[ty,:]*x_v[ty,:] + y_v[ty,:]*y_v[ty,:]
    
            # the numpy and IDL array functions to a shift have different directions and axis numbers
            # the np.roll axis and values were selected to match the output from the original IDL code
            # ;     ... j-direction: f-point segment direction (u-point)
            #       zxffu=  fsxspp( glamf, gphif ) - fsxspp( shift(glamf, 0, 1), shift(gphif, 0, 1) )
            #       zyffu=  fsyspp( glamf, gphif ) - fsyspp( shift(glamf, 0, 1), shift(gphif, 0, 1) )
            #       zmnpfu= sqrt ( znnpu * ( zxffu*zxffu + zyffu*zyffu )  )
            # IDL and python axis ordering and directions are different, this combination provides the best
            # # results when compared to the input data
            sh = -1 
            if type == 'jra55':
                sh = 1
            ax = 0
            # ax = 0 is rolling elements along that axis, i
            # x_fu[ty,:] = fx(glamf[ty,:], gphif[ty,:]) - fx( np.roll(glamf[ty,:],sh,axis=ax), np.roll(gphif[ty,:],sh,axis=ax)  )
            # y_fu[ty,:] = fy(glamf[ty,:], gphif[ty,:]) - fy( np.roll(glamf[ty,:],sh,axis=ax), np.roll(gphif[ty,:],sh,axis=ax)  )
            if ty-sh < len(glamf[:,0]):
                x_fu[ty,:] = fx(glamf[ty,:], gphif[ty,:]) - fx( glamf[ty-sh,:], gphif[ty-sh,:] )
                y_fu[ty,:] = fy(glamf[ty,:], gphif[ty,:]) - fy( glamf[ty-sh,:], gphif[ty-sh,:] )
            else:
                x_fu[ty,:] = fx(glamf[ty,:], gphif[ty,:]) - fx( glamf[0,:], gphif[0,:] )
                y_fu[ty,:] = fy(glamf[ty,:], gphif[ty,:]) - fy( glamf[0,:], gphif[0,:] )

            m_fu[ty,:] = np.sqrt( m_u[ty,:] * (x_fu[ty,:]*x_fu[ty,:] + y_fu[ty,:]*y_fu[ty,:]))

            # # ;     ... i-direction: f-point segment direction (v-point)
            # #       zxffv=  fsxspp( glamf, gphif ) - fsxspp( shift(glamf, 1, 0), shift(gphif, 1, 0) )
            # #       zyffv=  fsyspp( glamf, gphif ) - fsyspp( shift(glamf, 1, 0), shift(gphif, 1, 0) )
            # #       zmnpfv= sqrt ( znnpv * ( zxffv*zxffv + zyffv*zyffv )  )
            sh = 1
            ax = 1 # as we're taking stripes, this is just the array
            ax = 0
            x_fv[ty,:] = fx(glamf[ty,:], gphif[ty,:]) - fx( np.roll(glamf[ty,:],sh,axis=ax), np.roll(gphif[ty,:],sh,axis=ax)  )
            y_fv[ty,:] = fy(glamf[ty,:], gphif[ty,:]) - fy( np.roll(glamf[ty,:],sh,axis=ax), np.roll(gphif[ty,:],sh,axis=ax)  )
            m_fv[ty,:] = np.sqrt( m_v[ty,:] * (x_fv[ty,:]*x_fv[ty,:] + y_fv[ty,:]*y_fv[ty,:]))
    
            # # ;     ... cosinus and sinus using scalar and vectorial products
            # #       gsinu = ( zxnpu*zyffu - zynpu*zxffu ) / zmnpfu
            # #       gcosu = ( zxnpu*zxffu + zynpu*zyffu ) / zmnpfu
            gsin_u[ty,:] = (x_u[ty,:] * y_fu[ty,:] - y_u[ty,:] * x_fu[ty,:]) / m_fu[ty,:]
            gcos_u[ty,:] = (x_u[ty,:] * x_fu[ty,:] + y_u[ty,:] * y_fu[ty,:]) / m_fu[ty,:]
        
            # # ;     ... cosinus and sinus using scalar and vectorial products
            # # ;     (caution, rotation of 90 degres)
            # #       gsinv = ( zxnpv*zxffv + zynpv*zyffv ) / zmnpfv
            # #       gcosv =-( zxnpv*zyffv - zynpv*zxffv ) / zmnpfv
            gsin_v[ty,:] =      (x_v[ty,:] * x_fv[ty,:] + y_v[ty,:] * y_fv[ty,:]) / m_fv[ty,:]
            gcos_v[ty,:] = -1.0*(x_v[ty,:] * y_fv[ty,:] - y_v[ty,:] * x_fv[ty,:]) / m_fv[ty,:]

            # #       ind = where(abs(glamf-shift(glamf, 0, 1)) LT 1.e-8)
            # #       gsinu(ind) = 0.d
            #       gcosu(ind) = 1.d

        indexs = np.where(m_fu == 0)
        gsin_u[indexs] = 0
        gcos_u[indexs] = 1

        #       ind = where(abs(gphif-shift(gphif, 1, 0)) LT 1.e-8)
        #       gsinv(ind) = 0.d
        #       gcosv(ind) = 1.d
        indexs = np.where(m_fv == 0)
        gsin_v[indexs] = 0
        gcos_v[indexs] = 1
        # ---------------------------


        # ----------------- 01100 INTERPOLATION OF WIND STRESS POINTS  --------------------------------------
        if type != "jra55" and type != "era5":
            log.info("U component at u-point")
            u_onto_u = Interpolate(out_dims, xLim, yLim, tLim, u_data_in, lons, lats, gphiu, glamu, 
                                    origVars[v].directConversion, origVars[v].conversionFunc, origVars[v].limits, type, "U")
            log.info("U component at v-point")
            u_onto_v = Interpolate(out_dims, xLim, yLim, tLim, u_data_in, lons, lats, gphiv, glamv, 
                                    origVars[v].directConversion, origVars[v].conversionFunc, origVars[v].limits, type, "V")
            log.info("V component at u-point")
            v_onto_u = Interpolate(out_dims, xLim, yLim, tLim, v_data_in, lons, lats, gphiu, glamu, 
                                    origVars[v+1].directConversion, origVars[v+1].conversionFunc, origVars[v+1].limits, type, "U")
            log.info("V component at v-point")
            v_onto_v = Interpolate(out_dims, xLim, yLim, tLim, v_data_in, lons, lats, gphiv, glamv, 
                                    origVars[v+1].directConversion, origVars[v+1].conversionFunc, origVars[v+1].limits, type, "V")
        else:
            for t in range(0,tLim): # t is lon
                # log.info("U component at u-point perDay "+str(t))
                # u_onto_u[t,:,:] = Interpolate_perDay(out_dims, xLim, yLim, tLim, u_data_in[t,:,:], lons, lats, gphiu, glamu, 
                #                         origVars[v].directConversion, origVars[v].conversionFunc, origVars[v].limits, type, "U")
                # log.info("U component at v-point perDay "+str(t))
                # u_onto_v[t,:,:] = Interpolate_perDay(out_dims, xLim, yLim, tLim, u_data_in[t,:,:], lons, lats, gphiv, glamv, 
                #                         origVars[v].directConversion, origVars[v].conversionFunc, origVars[v].limits, type, "V")
                # log.info("V component at u-point perDay "+str(t))
                # v_onto_u[t,:,:] = Interpolate_perDay(out_dims, xLim, yLim, tLim, v_data_in[t,:,:], lons, lats, gphiu, glamu, 
                #                         origVars[v+1].directConversion, origVars[v+1].conversionFunc, origVars[v+1].limits, type, "U")
                # log.info("V component at v-point perDay "+str(t))
                # v_onto_v[t,:,:] = Interpolate_perDay(out_dims, xLim, yLim, tLim, v_data_in[t,:,:], lons, lats, gphiv, glamv, 
                #                         origVars[v+1].directConversion, origVars[v+1].conversionFunc, origVars[v+1].limits, type, "V")

                u_onto_u = np.zeros((out_dims[1],out_dims[2]),dtype=float)
                u_onto_v = np.zeros((out_dims[1],out_dims[2]),dtype=float)
                v_onto_v = np.zeros((out_dims[1],out_dims[2]),dtype=float)
                v_onto_u = np.zeros((out_dims[1],out_dims[2]),dtype=float)


                log.info("U component at u-point perDay "+str(t))
                u_onto_u[:,:] = Interpolate_perDay(out_dims, xLim, yLim, tLim, u_data_in[t,:,:], lons, lats, gphiu, glamu, 
                                        origVars[v].directConversion, origVars[v].conversionFunc, origVars[v].limits, type, "U")
                log.info("U component at v-point perDay "+str(t))
                u_onto_v[:,:] = Interpolate_perDay(out_dims, xLim, yLim, tLim, u_data_in[t,:,:], lons, lats, gphiv, glamv, 
                                        origVars[v].directConversion, origVars[v].conversionFunc, origVars[v].limits, type, "V")
                log.info("V component at u-point perDay "+str(t))
                v_onto_u[:,:] = Interpolate_perDay(out_dims, xLim, yLim, tLim, v_data_in[t,:,:], lons, lats, gphiu, glamu, 
                                        origVars[v+1].directConversion, origVars[v+1].conversionFunc, origVars[v+1].limits, type, "U")
                log.info("V component at v-point perDay "+str(t))
                v_onto_v[:,:] = Interpolate_perDay(out_dims, xLim, yLim, tLim, v_data_in[t,:,:], lons, lats, gphiv, glamv, 
                                        origVars[v+1].directConversion, origVars[v+1].conversionFunc, origVars[v+1].limits, type, "V")


            #  for t in range(0,tLim): # t is lon

                for ty in range(0,ylim):
                # adjust as v in opp direction to original data
                    # u_out[ty,:]   = u_onto_u[t,ty,:]  * gcos_u[ty,:]  + v_onto_u[t,ty,:]  * gsin_u[ty,:]
                    # v_out[ty,:]   = v_onto_v[t,ty,:]  * gcos_v[ty,:]  - u_onto_v[t,ty,:]  * gsin_v[ty,:]
                    u_out[ty,:]   = u_onto_u[ty,:]  * gcos_u[ty,:]  + v_onto_u[ty,:]  * gsin_u[ty,:]
                    v_out[ty,:]   = v_onto_v[ty,:]  * gcos_v[ty,:]  - u_onto_v[ty,:]  * gsin_v[ty,:]

                nc_out_sx_id.variables[in_varName_x][t,:,:] = u_out
                if type == 'jra55':
                    nc_out_sy_id.variables[in_varName_y][t,:,:] = v_out
                else:
                    nc_out_sy_id.variables[in_varName_y][t,:,:] = -1*v_out

                nc_out_sx_id.variables['time_counter'][t] = times_since_startOfYear[t]
                nc_out_sy_id.variables['time_counter'][t] = times_since_startOfYear[t]

        log.info("Momentum transfer interpolated from input data")
        # # ----------------- 01200 APPLYING COMPONENTS OF VECTOR TO GRID --------------------------------------       

        # # split up the calc into strips to prevent memory problems
        # xlim = u_onto_u.shape[2]
        # ylim = u_onto_u.shape[1]
        # x_u = np.zeros((ylim,xlim))
        # y_u = np.zeros((ylim,xlim))
        # m_u = np.zeros((ylim,xlim))
        # x_v = np.zeros((ylim,xlim))
        # y_v = np.zeros((ylim,xlim))
        # m_v = np.zeros((ylim,xlim))
        # x_fu = np.zeros((ylim,xlim))
        # y_fu = np.zeros((ylim,xlim))
        # m_fu = np.zeros((ylim,xlim))
        # x_fv = np.zeros((ylim,xlim))
        # y_fv = np.zeros((ylim,xlim))
        # m_fv = np.zeros((ylim,xlim))
        # gsin_u = np.zeros((ylim,xlim))
        # gcos_u = np.zeros((ylim,xlim))
        # gsin_v = np.zeros((ylim,xlim))
        # gcos_v = np.zeros((ylim,xlim))
        # u_out = np.zeros((ylim,xlim))
        # v_out = np.zeros((ylim,xlim))

        # for ty in range(0,ylim):
        #     print(ty)
        #     # the following is translated explicitly from the original IDL code (written 1996)
        #     # copies precede each section
        #     # ;     ... north pole direction & modulous (at u-point)
        #     #       zxnpu = 0. - fsxspp( glamu, gphiu )
        #     #       zynpu = 0. - fsyspp( glamu, gphiu )
        #     #       znnpu = zxnpu*zxnpu + zynpu*zynpu
        #     x_u[ty,:] = 0. - fx(glamu[ty,:],gphiu[ty,:])
        #     y_u[ty,:] = 0. - fy(glamu[ty,:],gphiu[ty,:])
        #     m_u[ty,:] = x_u[ty,:]*x_u[ty,:] + y_u[ty,:]*y_u[ty,:]

        #     # ;     ... north pole direction & modulous (at v-point)
        #     #       zxnpv = 0. - fsxspp( glamv, gphiv )
        #     #       zynpv = 0. - fsyspp( glamv, gphiv )
        #     #       znnpv = zxnpv*zxnpv + zynpv*zynpv
        #     x_v[ty,:] = 0. - fx(glamv[ty,:],gphiv[ty,:])
        #     y_v[ty,:] = 0. - fy(glamv[ty,:],gphiv[ty,:])
        #     m_v[ty,:] = x_v[ty,:]*x_v[ty,:] + y_v[ty,:]*y_v[ty,:]
    
        #     # the numpy and IDL array functions to a shift have different directions and axis numbers
        #     # the np.roll axis and values were selected to match the output from the original IDL code
        #     # ;     ... j-direction: f-point segment direction (u-point)
        #     #       zxffu=  fsxspp( glamf, gphif ) - fsxspp( shift(glamf, 0, 1), shift(gphif, 0, 1) )
        #     #       zyffu=  fsyspp( glamf, gphif ) - fsyspp( shift(glamf, 0, 1), shift(gphif, 0, 1) )
        #     #       zmnpfu= sqrt ( znnpu * ( zxffu*zxffu + zyffu*zyffu )  )
        #     # IDL and python axis ordering and directions are different, this combination provides the best
        #     # # results when compared to the input data
        #     sh = -1 
        #     if type == 'jra55':
        #         sh = 1
        #     ax = 0
        #     # ax = 0 is rolling elements along that axis, i
        #     # x_fu[ty,:] = fx(glamf[ty,:], gphif[ty,:]) - fx( np.roll(glamf[ty,:],sh,axis=ax), np.roll(gphif[ty,:],sh,axis=ax)  )
        #     # y_fu[ty,:] = fy(glamf[ty,:], gphif[ty,:]) - fy( np.roll(glamf[ty,:],sh,axis=ax), np.roll(gphif[ty,:],sh,axis=ax)  )
        #     if ty-sh < len(glamf[:,0]):
        #         x_fu[ty,:] = fx(glamf[ty,:], gphif[ty,:]) - fx( glamf[ty-sh,:], gphif[ty-sh,:] )
        #         y_fu[ty,:] = fy(glamf[ty,:], gphif[ty,:]) - fy( glamf[ty-sh,:], gphif[ty-sh,:] )
        #     else:
        #         x_fu[ty,:] = fx(glamf[ty,:], gphif[ty,:]) - fx( glamf[0,:], gphif[0,:] )
        #         y_fu[ty,:] = fy(glamf[ty,:], gphif[ty,:]) - fy( glamf[0,:], gphif[0,:] )

        #     m_fu[ty,:] = np.sqrt( m_u[ty,:] * (x_fu[ty,:]*x_fu[ty,:] + y_fu[ty,:]*y_fu[ty,:]))

        #     # # ;     ... i-direction: f-point segment direction (v-point)
        #     # #       zxffv=  fsxspp( glamf, gphif ) - fsxspp( shift(glamf, 1, 0), shift(gphif, 1, 0) )
        #     # #       zyffv=  fsyspp( glamf, gphif ) - fsyspp( shift(glamf, 1, 0), shift(gphif, 1, 0) )
        #     # #       zmnpfv= sqrt ( znnpv * ( zxffv*zxffv + zyffv*zyffv )  )
        #     sh = 1
        #     ax = 1 # as we're taking stripes, this is just the array
        #     ax = 0
        #     x_fv[ty,:] = fx(glamf[ty,:], gphif[ty,:]) - fx( np.roll(glamf[ty,:],sh,axis=ax), np.roll(gphif[ty,:],sh,axis=ax)  )
        #     y_fv[ty,:] = fy(glamf[ty,:], gphif[ty,:]) - fy( np.roll(glamf[ty,:],sh,axis=ax), np.roll(gphif[ty,:],sh,axis=ax)  )
        #     m_fv[ty,:] = np.sqrt( m_v[ty,:] * (x_fv[ty,:]*x_fv[ty,:] + y_fv[ty,:]*y_fv[ty,:]))
    
        #     # # ;     ... cosinus and sinus using scalar and vectorial products
        #     # #       gsinu = ( zxnpu*zyffu - zynpu*zxffu ) / zmnpfu
        #     # #       gcosu = ( zxnpu*zxffu + zynpu*zyffu ) / zmnpfu
        #     gsin_u[ty,:] = (x_u[ty,:] * y_fu[ty,:] - y_u[ty,:] * x_fu[ty,:]) / m_fu[ty,:]
        #     gcos_u[ty,:] = (x_u[ty,:] * x_fu[ty,:] + y_u[ty,:] * y_fu[ty,:]) / m_fu[ty,:]
        
        #     # # ;     ... cosinus and sinus using scalar and vectorial products
        #     # # ;     (caution, rotation of 90 degres)
        #     # #       gsinv = ( zxnpv*zxffv + zynpv*zyffv ) / zmnpfv
        #     # #       gcosv =-( zxnpv*zyffv - zynpv*zxffv ) / zmnpfv
        #     gsin_v[ty,:] =      (x_v[ty,:] * x_fv[ty,:] + y_v[ty,:] * y_fv[ty,:]) / m_fv[ty,:]
        #     gcos_v[ty,:] = -1.0*(x_v[ty,:] * y_fv[ty,:] - y_v[ty,:] * x_fv[ty,:]) / m_fv[ty,:]

        #     # #       ind = where(abs(glamf-shift(glamf, 0, 1)) LT 1.e-8)
        #     # #       gsinu(ind) = 0.d
        #     #       gcosu(ind) = 1.d

        # indexs = np.where(m_fu == 0)
        # gsin_u[indexs] = 0
        # gcos_u[indexs] = 1

        # #       ind = where(abs(gphif-shift(gphif, 1, 0)) LT 1.e-8)
        # #       gsinv(ind) = 0.d
        # #       gcosv(ind) = 1.d
        # indexs = np.where(m_fv == 0)
        # gsin_v[indexs] = 0
        # gcos_v[indexs] = 1
        # # ---------------------------


# # commment out full array operation vesion
        #         # the following is translated explicitly from the original IDL code (written 1996)
    #         # copies precede each section
    #         # ;     ... north pole direction & modulous (at u-point)
    #         #       zxnpu = 0. - fsxspp( glamu, gphiu )
    #         #       zynpu = 0. - fsyspp( glamu, gphiu )
    #         #       znnpu = zxnpu*zxnpu + zynpu*zynpu
    #         x_u = 0. - fx(glamu,gphiu)
    #         y_u = 0. - fy(glamu,gphiu)
    #         m_u = x_u*x_u + y_u*y_u

    #         # ;     ... north pole direction & modulous (at v-point)
    #         #       zxnpv = 0. - fsxspp( glamv, gphiv )
    #         #       zynpv = 0. - fsyspp( glamv, gphiv )
    #         #       znnpv = zxnpv*zxnpv + zynpv*zynpv
    #         x_v = 0. - fx(glamv,gphiv)
    #         y_v = 0. - fy(glamv,gphiv)
    #         m_v = x_v*x_v + y_v*y_v
    
    #         # the numpy and IDL array functions to a shift have different directions and axis numbers
    #         # the np.roll axis and values were selected to match the output from the original IDL code
    #         # ;     ... j-direction: f-point segment direction (u-point)
    #         #       zxffu=  fsxspp( glamf, gphif ) - fsxspp( shift(glamf, 0, 1), shift(gphif, 0, 1) )
    #         #       zyffu=  fsyspp( glamf, gphif ) - fsyspp( shift(glamf, 0, 1), shift(gphif, 0, 1) )
    #         #       zmnpfu= sqrt ( znnpu * ( zxffu*zxffu + zyffu*zyffu )  )
    #         # IDL and python axis ordering and directions are different, this combination provides the best
    #         # results when compared to the input data
    #         sh = -1
    #         ax = 0
    #         x_fu = fx(glamf, gphif) - fx( np.roll(glamf,sh,axis=ax), np.roll(gphif,sh,axis=ax)  )
    #         y_fu = fy(glamf, gphif) - fy( np.roll(glamf,sh,axis=ax), np.roll(gphif,sh,axis=ax)  )
    #         m_fu = np.sqrt( m_u * (x_fu*x_fu + y_fu*y_fu))

    #         # ;     ... i-direction: f-point           segment direction (v-point)
    #         #       zxffv=  fsxspp( glamf, gphif ) - fsxspp( shift(glamf, 1, 0), shift(gphif, 1, 0) )
    #         #       zyffv=  fsyspp( glamf, gphif ) - fsyspp( shift(glamf, 1, 0), shift(gphif, 1, 0) )
    #         #       zmnpfv= sqrt ( znnpv * ( zxffv*zxffv + zyffv*zyffv )  )

    #         sh = 1
    #         ax = 1
    #         x_fv = fx(glamf, gphif) - fx( np.roll(glamf,sh,axis=ax), np.roll(gphif,sh,axis=ax)  )
    #         y_fv = fy(glamf, gphif) - fy( np.roll(glamf,sh,axis=ax), np.roll(gphif,sh,axis=ax)  )
    #         m_fv = np.sqrt( m_v * (x_fv*x_fv + y_fv*y_fv))
    
    #         # ;     ... cosinus and sinus using scalar and vectorial products
    #         #       gsinu = ( zxnpu*zyffu - zynpu*zxffu ) / zmnpfu
    #         #       gcosu = ( zxnpu*zxffu + zynpu*zyffu ) / zmnpfu
    #         gsin_u = (x_u * y_fu - y_u * x_fu) / m_fu
    #         gcos_u = (x_u * x_fu + y_u * y_fu) / m_fu
        
    #         # ;     ... cosinus and sinus using scalar and vectorial products
    #         # ;     (caution, rotation of 90 degres)
    #         #       gsinv = ( zxnpv*zxffv + zynpv*zyffv ) / zmnpfv
    #         #       gcosv =-( zxnpv*zyffv - zynpv*zxffv ) / zmnpfv
    #         gsin_v =      (x_v * x_fv + y_v * y_fv) / m_fv
    #         gcos_v = -1.0*(x_v * y_fv - y_v * x_fv) / m_fv

    #         #       ind = where(abs(glamf-shift(glamf, 0, 1)) LT 1.e-8)
    #         #       gsinu(ind) = 0.d
    #         #       gcosu(ind) = 1.d
    #         indexs = np.where(m_fu == 0)
    #         gsin_u[indexs] = 0
    #         gcos_u[indexs] = 1

    #         #       ind = where(abs(gphif-shift(gphif, 1, 0)) LT 1.e-8)
    #         #       gsinv(ind) = 0.d
    #         #       gcosv(ind) = 1.d
    #         indexs = np.where(m_fv == 0)
    #         gsin_v[indexs] = 0
    #         gcos_v[indexs] = 1
    
# ---------------------------

        if type != "jra55" and type != "era5":
            for t in range(0,tLim): # t is lon
                for ty in range(0,ylim):
                # adjust as v in opp direction to original data
                    u_out[ty,:]   = u_onto_u[t,ty,:]  * gcos_u[ty,:]  + v_onto_u[t,ty,:]  * gsin_u[ty,:]
                    v_out[ty,:]   = v_onto_v[t,ty,:]  * gcos_v[ty,:]  - u_onto_v[t,ty,:]  * gsin_v[ty,:]
                nc_out_sx_id.variables[in_varName_x][t,:,:] = u_out
        if type == 'jra55' or type == 'era5':
            nc_out_sy_id.variables[in_varName_y][t,:,:] = v_out
        else:
            nc_out_sy_id.variables[in_varName_y][t,:,:] = -1*v_out
            nc_out_sx_id.variables['time_counter'][t] = times_since_startOfYear[t]
            nc_out_sy_id.variables['time_counter'][t] = times_since_startOfYear[t]



        # u_out = u_onto_u * gcos_u + v_onto_u * gsin_u
        # v_out = v_onto_v * gcos_v - u_onto_v * gsin_v

        # # adjust as v in opp direction to original data
        # nc_out_sx_id.variables[in_varName_x][:] = u_out
        # nc_out_sy_id.variables[in_varName_y][:] = -1*v_out

        log.info("Momentum transfer calculated")

# ----------------- 01300 TIDY UP AND CLOSE FILES --------------------------------------     
    finish = datetime.datetime.now()
    duration = finish-start
    log.info("Duration: "+str(duration))



for v in range(0,len(origVars)):
    origVars[v].nc_id.close()
    log.info(str(v)+": "+origVars[v].fileName+" Closed.")

if createOnlyMod == False and (varNum == 7 or varNum == -1):
    log.info("Closing bulk file")
    nc_out_id.close()

    log.info("Closing u and v flux files")
    nc_out_sx_id.close()
    nc_out_sy_id.close()

    # ------------------01400 CREATE KELVIN VERSION OF BULK -------------------------------------------
    # take copy of bulk file and renams
    outputFile_kelvin = "bulk_kelvin_"+str(year)+".nc"
    shutil.copy(outputFile, outputFile_kelvin)

    log.info("Creating: "+outputFile_kelvin)

    orig_id = Dataset(outputFile_kelvin,"r+")
    log.info("Original bulk file opened. Type: "+orig_id.data_model)

    airData = orig_id.variables['air'][:]
    newAirData = DegreesToKelvin(airData)
    orig_id.variables['air'][:] = newAirData

    # change units
    orig_id.variables['air'].setncattr("units", "degK")

    log.info("Closing: "+outputFile_kelvin)
    orig_id.close()

    finish = datetime.datetime.now()
    duration = finish-start
    log.info("Duration: "+str(duration))


# ------------------01500 CREATE DTAU FILE -------------------------------------------

if createMod == True:
    log.info("Creating wind stress modulus file")

    outputFile_dtau = "tdau_"+str(year)+".nc"
    outputFile_sx = "taux_1d_"+str(year)+".nc"
    outputFile_sy = "tauy_1d_"+str(year)+".nc"  

    nc_out_dt_id = Dataset(outputFile_dtau, 'w')

    nc_out_sx_id = Dataset(outputFile_sx, 'r')
    nc_out_sy_id = Dataset(outputFile_sy, 'r')

    log.info("Opening: "+outputFile_sx+" "+outputFile_sy )

    # Set the available dimensions
    nc_out_dt_id.createDimension("x", out_xDim)
    nc_out_dt_id.createDimension("y", out_yDim)
    nc_out_dt_id.createDimension("time_counter", None)

    # Set the nav-lon and lat variables
    nc_v_nav_dt_lon = nc_out_dt_id.createVariable("nav_lon", nav_dtype, ("y","x")) #, nc_mesh_id.variables["nav_lon"].dtype, ("y","x"))
    nc_v_nav_dt_lat = nc_out_dt_id.createVariable("nav_lat", nav_dtype, ("y","x")) #, nc_mesh_id.variables["nav_lat"].dtype, ("y","x"))
    nc_out_dt_id.variables['nav_lon'][:] = nc_lon_orig
    nc_out_dt_id.variables['nav_lat'][:] = nc_lat

    nc_v_nav_dt_lon.setncattr("units", "degrees_east")
    nc_v_nav_dt_lon.setncattr("valid_min", np.array(-179.7507,'f'))
    nc_v_nav_dt_lon.setncattr("valid_max", np.array(180.,'f'))
    nc_v_nav_dt_lon.setncattr("long_name", "Longitude")
    nc_v_nav_dt_lat.setncattr("units", "degrees_east")
    nc_v_nav_dt_lat.setncattr("valid_min", np.array(-78.19058,'f'))
    nc_v_nav_dt_lat.setncattr("valid_max", np.array(89.6139,'f'))
    nc_v_nav_dt_lat.setncattr("long_name", "Latitude")

    # Set the time variable - copy from first file in list
    nc_v_dt_time = nc_out_dt_id.createVariable("time_counter", 'f', ('time_counter'))
    nc_v_dt_time.setncattr("units", "seconds since    0-01-01 00:00:00")
    nc_v_dt_time.setncattr("calendar", "gregorian")
    nc_v_dt_time.setncattr("title", "time")
    nc_v_dt_time.setncattr("long_name", "   0-JAN-01 00:00:00")

    # create variables
    nc_dt = nc_out_dt_id.createVariable('tdif', nc_out_sx_id.variables['uflx'].dtype, ('time_counter', 'y', 'x'))

    nc_dt.setncattr("long_name", "wind stress modulus")
    nc_dt.setncattr("units", "")
    nc_dt.setncattr("short_name", "wind stress modulus")
    nc_dt.setncattr("axis", "TYX")
    nc_dt.setncattr("interval_operation",  np.array(86400.,'f'))
    nc_dt.setncattr("interval_write", np.array(86400.,'f'))
    nc_dt.setncattr("associate", "time_counter nav_lat nav_lon")
    nc_dt.setncattr("missing_value", np.array(1E20,'f'))

    times = nc_out_sx_id.variables['time_counter'][:]

    out_dims = [len(times),nc_out_sx_id.dimensions["y"].size,nc_out_sx_id.dimensions["x"].size]
    log.info("Output Dimensions [y,x]: "+str(out_dims))

    # Set output array
    data_out = np.zeros(out_dims,dtype=float)
    print(data_out.shape)
    tLim = data_out.shape[0]

    for t in range(0,tLim):
        print(t)
        u_data = nc_out_sx_id.variables['uflx'][t,:,:]
        v_data = nc_out_sy_id.variables['vflx'][t,:,:]
        print('--')
        data_out[t,:,:] = np.sqrt(u_data*u_data + v_data*v_data)
        nc_out_dt_id.variables['tdif'][t,:,:] = data_out[t,:,:]
        # nc_out_dt_id.variables['tdif'][t,:,:] = np.mod(u_data,v_data)

    nc_out_dt_id.close()


# ------------------01500 COMPLETE-------------------------------------------