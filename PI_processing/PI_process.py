from mitgcm_python.forcing import convert_pre_industrial_data
from mitgcm_python.file_io import read_binary

def main():
    process_var = True
    var_names = ['TREFHT', 'QBOT', 'PSL', 'UBOT', 'VBOT', 'PRECT', 'FLDS', 'FSDS']

    if process_var == True:
        for var in var_names:
            print(('Processing ' + var))
            convert_pre_industrial_data(var)

    year = str(sys.argv[1])
    option = "binary"
    single_plot = False
    exp = "PI"
    if option == "NCfile":
        filepath = "/data/oceans_input/raw_input_data/CESM/LENS/daily/TREFHT/"
        filename = "b.e11.B1850C5CN.f09_g16.005.cam.h1.TREFHT.04020101-04991231.nc"
        var = "TREFHT"
        id = nc.Dataset(filepath+filename, 'r')
        temp = id.variables[var][0:365,:,:]
        lat = id.variables["lat"][:]
        lon = id.variables["lon"][:]
                                
    if option == "binary":
        grid_sizes = [192, 288]
        dimensions = ('t','y','x')
        filepath = "/data/oceans_input/processed_input_data/CESM/PIctrl/"
        filename = filepath+"PIctrl_ens03_PRECT_"+str(year)
        data = read_binary(filename, grid_sizes, dimensions)
        time = range(0, 365)
        print(np.shape(data))
   
    
    if option == "mask":
        id = nc.Dataset("/data/oceans_output/shelf/kaight/archer2_mitgcm/PAS_LENS001_O/output/189001/MITgcm/output.nc", 'r')
        mask = id.variables["maskC"][0,1:384:2,1:576:2]
        theta = [[np.nan for i in range(cols)] for j in range(rows)]
        stop_lat = find_nearest(lat, -62.389)
        start_lat = find_nearest(lat, -75.638)
        lat = lat[start_lat:stop_lat]
        stop_lon = find_nearest(lon, 279.94)
        start_lon = find_nearest(lon, 220.05)
        lon = lon[start_lon:stop_lon]
        temp = temp[0,start_lat:stop_lat,start_lon:stop_lon]
    
    #for t in time:
    lon = range(0, 288)
    lat = range(0, 192)
    data = np.transpose(data[0,:,:])
    #print(lat, lon)
    contourplots(lon,lat,data,"temp"+year, year, exp)
    #make_timeseries_at_point(time, thetaold, thetanew, saltold[:,380,590], saltnew[:,380,590], "Potential surface temperature at 62°S 100°W in "+year, "Salinity at 62°S 100°W in "+year, year)
    #theta_over_area=make_timeseries_over_area(time,theta,"temperature over the area "+year, year)
    
   # filename = "/users/katner33/processed_data_test/TREFHT_1920"
   # grid_sizes = [192,288]
   # dimensions = ('t', 'y', 'x')
   # data = read_binary(filename, grid_sizes, dimensions)
   # print(data[1:10,1:10])

if __name__ == '__main__':
    main() # run the program
