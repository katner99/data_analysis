import xarray as xr
from mitgcm_python.grid import Grid
from funcs import find_nearest
from directories_and_paths import grid_filepath

grid = Grid(grid_filepath)
grid_file = xr.open_dataset(grid_filepath, decode_times=False)

depth_range = [
    find_nearest(grid_file.Z.values, -200),
    find_nearest(grid_file.Z.values, -700),
]

lat_range_73 = find_nearest(grid_file.YC.values, -73)
lon_range_73 = [
    find_nearest(grid_file.XC.values, 252.8),
    find_nearest(grid_file.XC.values, 255),
]

# below are the definitions for the undercurrent slices
lat_range_U1 = [find_nearest(grid_file.YC.values, -73.5), find_nearest(grid_file.YC.values, -71.8)]
lon_range_U1 = find_nearest(grid_file.XC.values, 238)

lat_range_U2 = [find_nearest(grid_file.YC.values, -73), find_nearest(grid_file.YC.values, -69.5)]
lon_range_U2 = find_nearest(grid_file.XC.values, 250)  

lat_range_U3 = [find_nearest(grid_file.YC.values, -74), find_nearest(grid_file.YC.values, -70)]
lon_range_U3 = find_nearest(grid_file.XC.values, 243) 

lat_range_U4 = [find_nearest(grid_file.YC.values, -72), find_nearest(grid_file.YC.values, -70)]
lon_range_U4 = find_nearest(grid_file.XC.values, 255)


sv = 10 ** (-6)

lat_bin, lon_bin, time_bin = 192, 288, 365

days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
