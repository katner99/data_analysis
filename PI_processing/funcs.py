import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import animation
from PIL import Image,ImageFilter
import numpy as np
import sys
import datetime
import netCDF4 as nc
from mitgcm_python.utils import mask_land_ice
from mitgcm_python.calculus import over_area
from mitgcm_python.file_io import read_binary

###################### OTHER FUNCTIONS ######################
# Fins the nearest value to a point
# INPUT:
# array = array you want to look through
# value = value you want to find
# OUTPUT:
# idx = returns the index with the closest value to the one chosen
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def interpolate_currents(u, v):
    u_interp = np.empty(np.shape(u))
    v_interp = np.empty(np.shape(v))
    u_interp[...,:-1] = 0.5*(u[...,:-1] + u[...,1:])
    u_interp[...,-1] = u[...,-1]
    v_interp[...,:-1,:] = 0.5*(v[...,:-1,:] + v[...,1:,:])
    v_interp[...,-1,:] = v[...,-1,:]
    return u_interp, v_interp

def read_netcdf(filepath, filename):
    id = nc.Dataset(filepath+filename, 'r')
    time = id.variables["time"][:]
    lat = id.variables["YC"][:]
    lon = id.variables["XC"][:]
    land_mask = id.variables["maskC"][:,:,:]
    ice_mask = id.variables["SIarea"][:,:,:]
    if var == "THETA":
        data = id.variables[var][:,11:21,:,:]
        cs = "coolwarm"
    if var == "SALT":
        data = id.variables[var][:,1,:,:]
        data[data == 0] = np.nan
        cs = "PRGn_r"
    if var == "SIfwmelt":
        data = id.variables[var][:,:,:]
        cs = "Blues_r" 
    return time, lat, lon, land_mask, ice_mask, data. cs
