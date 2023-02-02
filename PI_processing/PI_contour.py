import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import sys
from plots import contour_plots


def main():
    filepath = "/data/oceans_output/shelf/katner33/PIctrl_output/PAS_ctrl08/output/208801/MITgcm/"
    filename = "output.nc"
    id = nc.Dataset(filepath+filename, 'r')
    temp_d = id.variables["THETA"][0,11:21,:,:]
    temp = np.mean(temp_d, axis=0)
    time = id.variables["time"][:]
    lat = id.variables["YC"][:]
    lon = id.variables["XC"][:]
    depth = id.variables["Depth"][:,:]
    ice = id.variables["maskC"][0,:,:]
    pl = contour_plots(lon, lat, temp, "Theta", depth, ice, year = "2088", cm = "coolwarm")
    plt.show(pl)

if __name__ == '__main__':
    main() # run the program
