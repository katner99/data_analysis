import sys
import os.path
import datetime
import numpy as np
from netCDF4 import Dataset, MFDataset
import math
import matplotlib.pyplot as plt

# load files
fileName = "surface_pressure_ERA5_1990.nc"
variable = "sp"
surface_pressure = Dataset(fileName)
s_p = surface_pressure.variables[variable][:24, :, :]

fileName = "2m_dewpoint_temperature_ERA5_1990.nc"
variable = "d2m"
total_moisture_mass = Dataset(fileName)
d_T = total_moisture_mass.variables[variable][:24, :, :]

variable = "time"
time = total_moisture_mass.variables[variable][:24]

print(time.shape)

e = np.zeros((d_T[0,:,:].shape))
r = np.zeros((d_T[0,:,:].shape))
q = np.zeros((d_T[0,:,:].shape))
nc_test_id = Dataset("specific_humidity.nc", 'w', format='NETCDF3_CLASSIC')
nc_test_id.createDimension("y", 721)
nc_test_id.createDimension("x", 1440)
nc_v = nc_test_id.createVariable("s_h", 'f', ('x', 'y'))
nc_test_id.variables['s_h'][:] = np.transpose(q)

for iy, ix in np.ndindex(e.shape):
    # calculate vapour pressure (e)
    e[iy,ix] = 0.6112*math.exp(17.67*(d_T[1,iy,ix]-273.15)/(243.5 + (d_T[1,iy,ix]-273.15)))
        
    # calculate mixing ratio (r)
    r[iy,ix] = (0.622*(e[iy,ix])*1000)/((s_p[1,iy,ix]) - (e[iy,ix]*1000))
    
    # calculate specific humidity (q)
    q[iy,ix] = r[iy,ix]/(1 + r[iy,ix])
print(np.min(q))
print(np.min(r))
print(np.min(e))
nc_test_id.variables['s_h'][:,:] = np.flip(np.transpose(q[:,:]),1)