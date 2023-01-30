import netCDF4 as nc
import numpy as np
import sys
from plots import animate_contour


def main():
    year = str(sys.argv[1])
    var = str(sys.argv[2])
    exp = str(sys.argv[3])
    filepath = "/data/oceans_output/shelf/katner33/PIctrl_output/PAS_ctrl08/output/"+year+"01/MITgcm/"
    filename = "output.nc"
    id = nc.Dataset(filepath+filename, 'r')
    if var == "THETA":
        data = id.variables[var][:,11:21,:,:]
        cs = "coolwarm"
    if var == "SALT":
        data = id.variables[var][:,1,:,:]
        data[data == 0] = np.nan
        print(np.nanmax(data))
        cs = "PRGn_r"
    if var == "SIfwmelt":
        data = id.variables[var][:,:,:]
        #data[data == 0] = np.nan
        cs = "Blues_r"    
    time = id.variables["time"][:]
    lat = id.variables["YC"][:]
    lon = id.variables["XC"][:]
    land_mask = id.variables["maskC"][1,:,:]
    ice_mask = id.variables["SIarea"][:,:,:]
    
    print(np.shape(data))
    #depth_averaged = np.mean(data, axis=1)
    animate_contour(lon, lat, data, year, exp, var, cs, True, False, land_mask, ice_mask, False, True)

    
if __name__ == '__main__':
    main() # run the program
