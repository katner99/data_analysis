from directories_and_paths import *
from mitgcm_python.grid import Grid
from mitgcm_python.plot_latlon import latlon_plot, overlay_vectors
from mitgcm_python.utils import mask_3d
from mitgcm_python.plot_utils.latlon import overlay_vectors
from funcs import mask

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import numpy as np
import xarray as xr

def main():
    filepath_1 = output_path + "PAS_temp02/output/207401/MITgcm/output.nc"
    filepath_2 = output_path+ "PAS_ctrl02/output/207401/MITgcm/output.nc"

    title_beg = "A co "
    input_data_1 = xr.open_dataset(filepath_1)
    input_data_2 = xr.open_dataset(filepath_2)
    grid = Grid(grid_filepath)

    u_1 = np.average(mask(input_data_1.EXFuwind.values, grid), axis = 0)
    v_1 = np.average(mask(input_data_1.EXFvwind.values, grid), axis = 0)

    u_2 = np.average(mask(input_data_2.EXFuwind.values, grid), axis = 0)
    v_2 = np.average(mask(input_data_2.EXFvwind.values, grid), axis = 0)

    speed_1 = np.sqrt(u_1**2 + v_1**2)
    speed_2 = np.sqrt(u_2**2 + v_2**2)

    #speed_1 = input_data_1.EXFatemp.values[0,...]
    #speed_2 = input_data_2.EXFatemp.values[0,...]

    u = u_1 - u_2
    v = v_1 - v_2

    speed = speed_1 - speed_2
 
    fig, ax = latlon_plot(speed, grid, ctype='PiYG', include_shelf=True, title="Wind temp - ctrl", return_fig=True)

    overlay_vectors(ax, u, v, grid, chunk = 10, scale = 5)

    fig.savefig("wind_difference_TC.png")

    #fig.show()

    

if __name__ == '__main__':
    main() # run the program
