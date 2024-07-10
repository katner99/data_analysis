import numpy as np
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from plots import read_mask 
from calcs import monthly_average
from directories_and_paths import output_path

exp = "WIND"
filepath_force_1 = f"{output_path}{exp}_files_temp/oceFWflx.nc"
filepath_force_2 = f"{output_path}LENS_files_temp/oceFWflx.nc"
filepath_ctrl = f"{output_path}CTRL_files_temp/oceFWflx.nc"

dataset_force_1 = xr.open_dataset(filepath_force_1, decode_times=False)
dataset_force_2 = xr.open_dataset(filepath_force_2, decode_times=False)
dataset_ctrl = xr.open_dataset(filepath_ctrl, decode_times=False)

data_force_1 = monthly_average(dataset_force_1["oceFWflx"].values[(1980 - 1920)*12:(2030 - 1920)*12,...]* 3600 * 24 * 365 / 1000)
data_force_2 = monthly_average(dataset_force_2["oceFWflx"].values[(1980 - 1920)*12:(2030 - 1920)*12,...]* 3600 * 24 * 365 / 1000)

data_ctrl = np.nanmean(dataset_ctrl["oceFWflx"].values* 3600 * 24 * 365 / 1000, axis = 0)
DATA = data_force_1 - data_force_2
print (np.shape(DATA))

set_up_data = xr.open_dataset(f"{output_path}average_CTRL_1920-1950.nc", decode_times=False)
set_up = read_mask(set_up_data)

fig,ax = plt.subplots()

count = 0
year = 1980

def animate(i):
    ax.clear()
    min_val = -2
    max_val = -min_val
    ax.contourf(set_up["X"], set_up["Y"], DATA[i,:,:], cmap="PRGn_r",extend="both",levels=np.linspace(min_val, max_val, 15))
    ax.contourf(set_up["X"], set_up["Y"],set_up["land_mask"],cmap=matplotlib.colors.ListedColormap(set_up["colors"]))
    ax.contour(
        set_up["X"], set_up["Y"], set_up["mask"], 2, cmap="Greys", linestyles="dashed"
    )
    ax.set_title(f"{exp} - LENS {year + i}")

interval = 2#in seconds     
ani = animation.FuncAnimation(fig,animate,frames=50, interval = 100)

ani.save(f"{exp}_freshwater_LENS.gif", fps = 2)
plt.show()