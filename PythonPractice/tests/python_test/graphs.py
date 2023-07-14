############################################################################
# Written By Katherine Turner
# Date: 25 October 2022
# Practice script to read sample model output data and produce a graph
############################################################################
import numpy as np 
import netCDF4 as nc
import matplotlib
matplotlib.use("Agg")
from matplotlib.axis import Axis
import matplotlib.pyplot as plt 

#load the path and data name to the script
path = "/users/katner33/PAS_test_output/PAS_test/output/189001/MITgcm/"
output1999 = nc.Dataset(path+"output.nc")

# we are going to create a grph of monthly temperatures over Milan in the year 1999
# from netCDF file load the temperature variable ant time variable
temperatureMilan = output1999.variables["EXFatemp"][:,45,9]
time1999 = output1999.variables["time"][:]

# check that these have been loaded correctly
print(len(temperatureMilan),len(time1999))

# plot a figure and we will save this to temp.png
fig = plt.figure()
#fig, ax = plt.subplots()

x = time1999
y = temperatureMilan
labels = ['January','February','March','April','May','June','July','August','September','October','November','December']

#ax = fig.add_subplot()
plt.plot(x, y, color='green', linewidth = 2)
plt.xticks(x, labels,rotation=45)
#plt.set_xticklabels()

plt.xlabel("time (Months)")
plt.ylabel("temperature (K)")
plt.title("Monthly temperature over Milan 1999")

fig.savefig("temp.png")
