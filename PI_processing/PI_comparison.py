import netCDF4 as nc
import numpy as np
import sys
from plots import compare_contour_plots, compare_contour_plots_TSM


def main():
    # read in the year and variable
    year = str(sys.argv[1])
    var = str(sys.argv[2])
    month = sys.argv[3]
    exp = "PIctrl"

    # read in the data from two different ensembles
    filepath_ens08 = "/data/oceans_output/shelf/katner33/PIctrl_output/PAS_ctrl07/output/"+year+"01/MITgcm/"
    #filepath_ens08 = "/data/oceans_output/shelf/katner33/PIctrl_output/PAS_ctrl08/output/"+year+"01/MITgcm/"
    filepath_ens07 = "/data/oceans_output/shelf/kaight/archer2_mitgcm/PAS_LENS001_O/output/192001/MITgcm/"
    filename = "output.nc"
    id07 = nc.Dataset(filepath_ens07+filename, 'r')
    id08 = nc.Dataset(filepath_ens08+filename, 'r')

    # read in the remaining variables (these should be the same between the ensembles
    time = id07.variables["time"][:]
    lat = id07.variables["YC"][:]
    lon = id07.variables["XC"][:]
    ice_mask = id07.variables["maskC"][1,:,:]
    depth = id07.variables["Depth"][:,:]

    # read in the variables chosen and set the associated color schemes
    if var == "THETA":
        # create an average for 200m to 700m 
        data_4D_07 = id07.variables[var][month,11:21,:,:]
        data_4D_08 = id08.variables[var][month,11:21,:,:]
        data_07 = np.mean(data_4D_07, axis=0)
        data_08 = np.mean(data_4D_08, axis=0)
        cs = "coolwarm"
        compare_contour_plots(lon, lat, data_07, data_08, month, apply_mask = True, mask = land_mask, title1 = "ens_07", title2 = "ens_08", var = var, year = year, cm = cs, save = True, show = False)
    elif var == "SALT":
        data_07 = id07.variables[var][month,1,:,:]
        data_08 = id08.variables[var][month,1,:,:]
        data_07[data_07 == 0] = np.nan
        data_08[data_08 == 0] = np.nan
        cs = "PRGn_r"
        compare_contour_plots(lon, lat, data_07, data_08, month, apply_mask = True, mask = land_mask, title1 = "ens_07", title2 = "ens_08", var = var, year = year, cm = cs, save = True, show = False)
    elif var == "SIfwmelt":
        data_07 = id07.variables[var][month,:,:]
        data_08 = id08.variables[var][month,:,:]
        cs = "Blues_r"    
        compare_contour_plots(lon, lat, data_07, data_08, month, apply_mask = True, mask = land_mask, title1 = "ens_07", title2 = "ens_08", var = var, year = year, cm = cs, save = True, show = False)
    elif var == "TSM":  # create one beast of a plot with all Temperature, Salinity, and Melt (TSM)
        # temperature
        temp_4D_07 = id07.variables["THETA"][month,11:21,:,:]
        temp_4D_08 = id08.variables["THETA"][month,11:21,:,:]
        temp_07 = np.mean(temp_4D_07, axis=0)
        temp_08 = np.mean(temp_4D_08, axis=0)
        # salinity
        salt_07 = id07.variables["SALT"][month,1,:,:]
        salt_08 = id08.variables["SALT"][month,1,:,:]
        salt_07[salt_07 == 0] = np.nan
        salt_08[salt_08 == 0] = np.nan
        # melt
        melt_07 = id07.variables["SIfwmelt"][month,:,:]
        melt_08 = id08.variables["SIfwmelt"][month,:,:]
        compare_contour_plots_TSM("PIctrlVSLENS",lon, lat, temp_07, temp_08, salt_07, salt_08, melt_07, melt_08, month, ens1 = "LENS", ens2 = "PIctrl", depth = depth, ice_mask = ice_mask, year = "1920", save = True, show = True)
    
if __name__ == '__main__':
    main() # run the program
