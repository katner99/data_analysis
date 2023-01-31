from mitgcm_python.forcing import convert_pre_industrial_data
from mitgcm_python.file_io import read_binary

def main():
    process_var = True
    var_names = ['TREFHT', 'QBOT', 'PSL', 'UBOT', 'VBOT', 'PRECT', 'FLDS', 'FSDS']

    if process_var == True:
        for var in var_names:
            print(('Processing ' + var))
            convert_pre_industrial_data(var)
    
   # filename = "/users/katner33/processed_data_test/TREFHT_1920"
   # grid_sizes = [192,288]
   # dimensions = ('t', 'y', 'x')
   # data = read_binary(filename, grid_sizes, dimensions)
   # print(data[1:10,1:10])

if __name__ == '__main__':
    main() # run the program
