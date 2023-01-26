import netCDF4 as nc
import sys
from plots import animate_contour


def main():
    year = str(sys.argv[1])
    var = str(sys.argv[2])
    exp = str(sys.argv[3])
    filepath = "/data/oceans_output/shelf/katner33/PIctrl_output/PAS_ctrl/output/"+year+"01/MITgcm/"
    filename = "output.nc"
    id = nc.Dataset(filepath+filename, 'r')
    data = id.variables[var][:,1,:,:]
    time = id.variables["time"][:]
    lat = id.variables["YC"][:]
    lon = id.variables["XC"][:]
    land_mask = id.variables["maskC"][1,:,:]
    ice_mask = id.variables["SIarea"][:,:,:]
    
    animate_contour(lon, lat, data, year, exp, var, True, True, land_mask, ice_mask, True, False)

    
if __name__ == '__main__':
    main() # run the program