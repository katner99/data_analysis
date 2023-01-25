#!/usr/bin/env python
import sys
import os.path
import datetime
import numpy as np
import glob
from netCDF4 import Dataset, MFDataset
import math
import matplotlib.pyplot as plt

# get list of files
fileName = sys.argv[1]
print(fileName)

variable = sys.argv[2] # use ALL for all
print(variable)

ncfile = Dataset(fileName)
temperature = ncfile.variables[variable][:]

temp = np.zeros((30,192,144))
nc_test_id = Dataset("dailytemp.nc", 'w', format='NETCDF3_CLASSIC')
nc_test_id.createDimension("x", 192)
nc_test_id.createDimension("y", 144)
nc_test_id.createDimension("t", None)
nc_v = nc_test_id.createVariable("dailytemperature", 'f', ('t', 'x', 'y'))
nc_test_id.variables['dailytemperature'][:] = temp

for day in range(30):
    temp[day,:,:] = np.transpose(np.mean(temperature[day*7+1:day*7+8,:,:],axis=0))
    nc_test_id.variables['dailytemperature'][day,:,:] = temp[day,:,:]
print(np.shape(temp))
print(np.shape(temperature))

