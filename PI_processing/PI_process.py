from mitgcm_python.forcing import convert_pre_industrial_data
from mitgcm_python.file_io import read_binary
from mitgcm_python.grid import Grid

def main():
    ## check that enough arguments have been input by the user
    if len(sys.argv) != 2:
        sys.exit("Stopped - Incorrect number of arguements. Use python PI_process.py <option>")
    
    # check the option chosen is valid
    lemenu = ["open_binary", "process_var", "process_OBCS"]
    option = str(sys.argv[1])
    if option not in lemenu:
        sys.exit("Stopped - Invalid option. Please choose from my amazing menu selection of <oneVSone>, <the_monster>")
    
    # process netcdf files and write to binary files
    if option == lemenu[1]:
        var_names = ['TREFHT', 'QBOT', 'PSL', 'UBOT', 'VBOT', 'PRECT', 'FLDS', 'FSDS']
        for var in var_names:
            print(('Processing ' + var))
            convert_pre_industrial_data(var)

    # read and check the binary files
    if option == "binary":
        grid_filepath = "/data/oceans_output/shelf/kaight/archer2_mitgcm/PAS_LENS001_O/output/192001/MITgcm/output.nc"
        grid = Grid(grid_filepath)
        dimensions = ('t','y','x')
        filepath = "/data/oceans_input/processed_input_data/CESM/PIctrl/"
        filename = filepath+"PIctrl_ens03_PRECT_"+str(year)
        data = read_binary(filename, grid, dimensions)
        print(np.shape(data))
   
    
    if option == "OBCS":
        location = ["N", "W", "E"]
        variables = [""]
        for l in location:
            # set up the climatology years
            start_year = 2080
            end_year = 2100

            # set up the number of ensembles
            ensembles = range(1, 11)
            var = "VVEL"
            filepath = "/data/oceans_input/processed_input_data/CESM/PAS_obcs/LENS_obcs/LENS_ens"

            # load up the MITgcm grid from the first LENS run
            grid_filepath = "/data/oceans_output/shelf/kaight/archer2_mitgcm/PAS_LENS001_O/output/192001/MITgcm/output.nc"
            grid = Grid(grid_filepath)
            if location in ["E", "W"]:
                print("longitude")
                grid_sizes = [12, grid.nz, grid.ny]
                dimensions = ('t','z','x')
                data_ens = np.zeros([len(ensembles), grid.nz,  grid.ny, 12])
                data_year = np.zeros([21, grid.nz, grid.ny, 12])
                #grid_sizes = [12, grid.ny]
                #dimensions = ('t', 'x')
                #data_ens = np.zeros([len(ensembles), grid.ny, 12])
                #data_year = np.zeros([21, grid.ny, 12])
            elif location in ["N", "S"]:
                print("latitude")
                grid_sizes = [12, grid.nz, grid.nx]
                dimensions = ('t','z','y')
                data_ens = np.zeros([len(ensembles), 12, grid.nx, grid.nz])
                data_year = np.zeros([21, 12, grid.nx, grid.nz])
                #grid_sizes = [12, grid.nx]
                #dimensions = ('t','y')
                #data_ens = np.zeros([len(ensembles), 12, grid.nx])
                #data_year = np.zeros([21, 12, grid.nx])
            else:
                sys.exit()
    
                
            for year in range(start_year, end_year+1):
                for ens in ensembles: 
                    ens_str = str(ens).zfill(3)
                    filename = filepath + ens_str+"_"+var+"_"+location+"_"+str(year)

                    # read the variable
                    data = read_binary(filename, grid_sizes, dimensions)
                    #print(np.shape(data), np.shape(data_ens))
                    data_ens[ens-1,:, :, :] = data
                data_year[year-start_year, :, :, :] = np.mean(data_ens, axis = 0)
                print(year)
            data = np.mean(data_year, axis = 0)

            #ẁrite out the file
            out_filepath = "/data/oceans_input/processed_input_data/CESM/PI_CTRL_obcs/LENS_climatology_"+var+"_"+location+"_"+str(start_year)+"-"+str(end_year)
            write_binary(data, out_filepath)
        

if __name__ == '__main__':
    main() # run the program
