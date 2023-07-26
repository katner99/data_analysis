import netCDF4 as nc
import numpy as np
import sys
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from mitgcm_python.plot_latlon import plot_vel
from funcs import interpolate_currents
from plots import quiver_plot, compare_quiver_plots
import time

def main():
    start_time = time.time()
    year = str(2099)
    var = str("wind")
    font_size = 20

    # load the initial files
    filepath1 = "/data/oceans_output/shelf/katner33/PIctrl_output/PAS_ctrl08/output/"+year+"01/MITgcm/"
    filepath2 = "/data/oceans_output/shelf/kaight/archer2_mitgcm/PAS_LENS001_O/output/"+year+"01/MITgcm/"
    filename = "output.nc"
    id1 = nc.Dataset(filepath1+filename, 'r')
    id2 = nc.Dataset(filepath2+filename, 'r')

    u = "EXFuwind"
    v = "EXFvwind"
    U = id1.variables[u][0,:,:]
    V = id1.variables[v][0,:,:]
    TAU1 = id1.variables["oceTAUX"][0,:,:]
    print("--- %s seconds ---" % (time.time() - start_time))
    tau_interp = interpolate_currents(TAU1, "zonal")
    print("--- %s seconds ---" % (time.time() - start_time))
    U2 = id2.variables[u][0,:,:]
    V2 = id2.variables[v][0,:,:]
    TAU2 = id2.variables["oceTAUX"][0,:,:]
    
    tau_interp2 = interpolate_currents(TAU2, "zonal")
    print("--- %s seconds ---" % (time.time() - start_time))
    lat = id1.variables["YC"][:]
    lon = id1.variables["XC"][:]
    depth = id1.variables["Depth"][:,:]
    ice_mask = id1.variables["maskC"][1,:,:]
    print("--- %s seconds ---" % (time.time() - start_time))
     
    quiver_plot(U, V, lon, lat, exp = "PI_ctrl", var = "wind_stress", depth = depth, ice_mask = ice_mask, year = year, interp = False, overlay = True, func = tau_interp, show = False, save = True)
    #compare_quiver_plots(U2, V2, U1, V1, tau_interp2, tau_interp1, "LENS", "PIctrl", lon, lat, var, depth, ice_mask, year, interp = False, show = False, save = True)

if __name__ == '__main__':
    main() # run the program
