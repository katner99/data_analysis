import numpy as np
from scipy import stats
from directories_and_paths import output_path, lens_path

import xarray as xr

def calc_average(experiments, var):
    ensemble = [1,2,3,4,5,6,7,8,9]
    years = range(2070, 2100)

    expy_vals = []
    for exp in experiments:
        ensemble_data = []
        for ens in ensemble:
            time_mean = [] 
            for year in years:
                filepath = f"output/{year}01/MITgcm/output.nc"
                if exp == "LENS":
                    if ens < 6:
                        filepath = f"{lens_path}PAS_LENS00{ens}_noOBC/{filepath}"
                    else:
                        filepath = f"{output_path}PAS_LENS00{ens}_noOBC/{filepath}"
                else:
                    filepath = f"{output_path}{exp}_ens0{ens}_noOBC/{filepath}"

                try:
                    with xr.open_dataset(filepath) as data_member:
                        input_data = data_member[var].where(data_member.hFacC == 1).mean(dim = "time")
                        subset = input_data.sel(Z = slice(-200, -700)).mean(dim = "Z")                        
                        time_mean.append(subset)
                
                except Exception as e:
                    print(f"Error reading file: {filepath}, {e}")

            ensemble_data.append(xr.concat(time_mean, dim = "year").mean(dim = "year"))
            print(filepath)

        expy_vals.append(xr.concat(ensemble_data, dim = "ensemble")) 
    
    expy_vals = xr.concat(expy_vals, dim = "experiment")
    print(expy_vals)
    return expy_vals

def main():
    experiments = ["LENS", "WIND"]
    expy_vals = calc_average(experiments, "THETA")

    filename = f"{experiments[1]}_temp.nc"
    # save incase the next bit dies RIP
    expy_vals.to_netcdf(filename)

    theta_vals = expy_vals.THETA.values
    pval = np.empty(theta_vals[0,0,...].shape)

    for lat in range(len(expy_vals.YC.values)):
        for lon in range(len(expy_vals.XC.values)):
            values_1 = theta_vals[0, :, lat, lon]
            values_2 = theta_vals[1, :, lat, lon]
            test = stats.ttest_ind(values_1, values_2)
            pval[lat, lon] = test.pvalue


    trend = xr.Dataset(
        {
            'pvalue': (['lat', 'lon'], pval)
        },
        coords={'lat':expy_vals.YC.values, 'lon': expy_vals.XC.values}
    )

    filename = f"{experiments[1]}_temp.nc"
    # save incase the next bit dies RIP
    trend.to_netcdf(filename)
    

if __name__ == "__main__":
    main()

    

