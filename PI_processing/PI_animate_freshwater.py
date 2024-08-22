import numpy as np
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from tools.plots import read_mask 
from tools.calcs import monthly_average
from .directories_and_paths import output_path

def animate_comparison(var, exp_1, exp_2 = "CTRL", start = 1920, stop = 2101):
    filepaths = [f"{output_path}{exp}_files_temp/oceFWflx.nc" for exp in [exp_1, exp_2]]
    datasets = [xr.open_dataset(filepath, decode_times=False) for filepath in filepaths]
    data = [monthly_average(dataset[var].values[(start - 1920)*12:(stop - 1920)*12,...]* 3600 * 24 * 365 / 1000) for dataset in datasets]
    DATA = data[0] - data[1]

    set_up_data = xr.open_dataset(f"{output_path}average_CTRL_1920-1950.nc", decode_times=False)
    set_up = read_mask(set_up_data)

    fig,ax = plt.subplots()

    def animate(i):
        ax.clear()
        min_val = -2
        max_val = -min_val
        ax.contourf(set_up["X"], set_up["Y"], DATA[i,:,:], cmap="PRGn_r",extend="both",levels=np.linspace(min_val, max_val, 15))
        ax.contourf(set_up["X"], set_up["Y"],set_up["land_mask"],cmap=matplotlib.colors.ListedColormap(set_up["colors"]))
        ax.contour(
            set_up["X"], set_up["Y"], set_up["mask"], 2, cmap="Greys", linestyles="dashed"
        )
        ax.set_title(f"{exp_1} - {exp_2} {start + i}")
  
    ani = animation.FuncAnimation(fig,animate,frames=50, interval = 100)

    ani.save(f"{exp_1}vs{exp_2}_freshwater_{start}.gif", fps = 2)
    plt.show()

def main():
    # example use
    animate_comparison("oceFWflx", "LENS", exp_2 = "CTRL", start = 1920, stop = 2101)

