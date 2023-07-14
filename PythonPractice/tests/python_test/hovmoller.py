###############################################################################
# Script written by Katherine Turner
# Date: 25/10/2022
# script to create a hovmoller diagram
###############################################################################

import netCDF4 as nc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

############# 00100 LOAD DATA ################

# load file
path = "/data/oceans_output/shelf/kaight/archer2_mitgcm/PAS_LENS001_O/output/"
hovmoller = nc.Dataset(path+"hovmoller.nc")

# load variables for time, height, and temperature
temp = hovmoller.variables["pine_island_bay_temp"][:]
time = hovmoller.variables["time"][:]
height = hovmoller.variables["Z"][:]

# check the dimensions of the variables downloaded
print(np.ndim(temp),len(time),len(height))

############## 00200 PREP DATA ###################
# create a 2D grid
[X, Y] = np.meshgrid(time, height)
Z = np.transpose(temp)

fig, ax = plt.subplots(1,1)

# plot contour lines
ax.contourf(X, Y, Z1, cmap = "Spectral")
ax.set_ylim(-1100,0)
ax.set_title('Pine Island Glacier Temperature Contour Plot')
ax.set_xlabel('Time')
ax.set_ylabel('Depth')
  
fig.savefig("hovmollerTemp.png")
