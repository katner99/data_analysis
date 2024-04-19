import sys

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import numpy as np
import xarray as xr

from funcs import find_nearest
from directories_and_paths import output_path
from mitgcm_python.grid import Grid

def main():
    exp = "TEMP"
    ensemble = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    loc = "U1"
    var = f"salt_{loc}"

    filepaths = [f"{output_path}{exp}_ens0{ens}_noOBC/ensemble_mean_{exp}.nc" for ens in ensemble]
    
    dataset = []
    for filepath in filepaths:
        try:
            open(filepath)
        except FileNotFoundError:
            sys.exit(f"Stopped - Could not find input file {filepath}")
    
        input_data = xr.open_dataset(filepath, decode_times = False)
        dataset.append(input_data[var].values)
    
if __name__ == '__main__':
    main() # run the program
