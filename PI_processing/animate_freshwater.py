import numpy as np
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from plots import read_mask 
from directories_and_paths import output_path

exp = "LENS"
filepath_force = f"{output_path}{exp}_files_temp/oceFWflx.nc"
filepath_ctrl = f"{output_path}CTRL_files_temp/oceFWflx.nc"

dataset_force = xr.open_dataset(filepath_force, decode_times=False)
dataset_ctrl = xr.open_dataset(filepath_ctrl, decode_times=False)

DATA = (dataset_force["oceFWflx"].values* 3600 * 24 * 365 / 1000) - (dataset_ctrl["oceFWflx"].values* 3600 * 24 * 365 / 1000)
print (np.shape(DATA))
set_up_data = xr.open_dataset(f"{output_path}average_CTRL_1920-1950.nc", decode_times=False)
set_up = read_mask(set_up_data)

fig,ax = plt.subplots()

count = 0
year = 1920

def animate(i):
    count = count + 1
    ax.clear()
    min_val = -5
    max_val = -min_val
    ax.contourf(set_up["X"], set_up["Y"], DATA[i,:,:], cmap="PRGn_r",extend="both",levels=np.linspace(min_val, max_val, 15))
    ax.contourf(set_up["X"], set_up["Y"],set_up["land_mask"],cmap=matplotlib.colors.ListedColormap(set_up["colors"]))
    ax.contour(
        set_up["X"], set_up["Y"], set_up["mask"], 2, cmap="Greys", linestyles="dashed"
    )
    if count == 12:
        count = 0
        year = year + 1
    ax.set_title(f"{exp} - CTRL {year}")


interval = 2#in seconds     
ani = animation.FuncAnimation(fig,animate,frames=180*12, interval = 100)

ani.save(f"{exp}_freshwater.gif", fps = 2)
plt.show()