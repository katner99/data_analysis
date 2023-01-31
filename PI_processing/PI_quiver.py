import netCDF4 as nc
import numpy as np
import sys
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from mitgcm_python.plot_latlon import plot_vel

def main():
    year = str(sys.argv[1])
    var = str(sys.argv[2])

    # load the initial files
    filepath = "/data/oceans_output/shelf/katner33/PIctrl_output/PAS_ctrl08/output/"+year+"01/MITgcm/"
    filename = "output.nc"
    id = nc.Dataset(filepath+filename, 'r')

    if var == "wind":
        u = "EXFuwind"
        v = "EXFvwind"
        U = id.variables[u][0,0:-1:20,0:-1:20]
        V = id.variables[v][0,0:-1:20,0:-1:20]
    
        temp = id.variables["THETA"][0,11:21,:,:]
        temp_d = np.mean(temp, axis=0)
    
        time = id.variables["time"][:]
        lat = id.variables["YC"][:]
        lon = id.variables["XC"][:]
        land_mask = id.variables["maskC"][1,:,:]
       
        [X,Y] = np.meshgrid(lon[0:-1:20], lat[0:-1:20])
        [LON,LAT] = np.meshgrid(lon, lat)
        temp_d[land_mask == 0] = np.nan
        title = "Test"
        fig, ax = plt.subplots()

        cp = plt.contourf(LON, LAT, temp_d, cmap="coolwarm")
        cb = plt.colorbar(cp)
        quiv = plt.quiver(X, Y, U, V, color = "white")

        plt.show()

    if var == "currents":
        # load up the current speeds
        U = id.variables["UVEL"][0,:,:,:] # [time, Z, YC, XG]
        V = id.variables["VVEL"][0,:,:,:] # [time, Z, YG, XC]
        # load up the grid
        lat = id.variables["YC"][:]
        lon = id.variables["XC"][:]
        [LAT, LON] = np.meshgrid(lon, lat)
        # load up the corresponding masks
        umask = id.variables["hFacW"] # [Z, YC, XG]
        vmask = id.variables["hFacS"] # [Z, YG, XC]
        # apply the mask
        U[umask == 1] = np.nan
        V[vmask == 1] = np.nan
        # plot using Kaitlin's code
        plot_vel(U, V, [LON, LAT], vel_option="bottom")

if __name__ == '__main__':
    main() # run the program
