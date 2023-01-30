import netCDF4 as nc
import numpy as np
import sys
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


def main():
    year = str(2078)
    u = "EXFuwind"
    v = "EXFvwind"
    filepath = "/data/oceans_output/shelf/katner33/PIctrl_output/PAS_ctrl08/output/"+year+"01/MITgcm/"
    filename = "output.nc"
    id = nc.Dataset(filepath+filename, 'r')
    U = id.variables[u][0,0:-1:20,0:-1:20]
    V = id.variables[v][0,0:-1:20,0:-1:20]
    
    temp = id.variables["THETA"][0,11:21,:,:]
    temp_d = np.mean(temp, axis=0)
    print(np.shape(temp_d))
    time = id.variables["time"][:]
    lat = id.variables["YC"][:]
    lon = id.variables["XC"][:]
    land_mask = id.variables["maskC"][1,:,:]
    print(np.shape(land_mask))
    [X,Y] = np.meshgrid(lon[0:-1:20], lat[0:-1:20])
    [LON,LAT] = np.meshgrid(lon, lat)
    temp_d[land_mask == 0] = np.nan

    fig, ax = plt.subplots()
    axes = [ax, ax.twinx()]
    axes[0].contourf(LON, LAT, temp_d, cmap="coolwarm")
    axes[1].quiver(X, Y, U, V)
    plt.show()

    
if __name__ == '__main__':
    main() # run the program
