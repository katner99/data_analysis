from mitgcm_python.ics_obcs import process_cesm_obcs_ts, process_cesm_obcs_non_ts
import logging

expt = 'PI_CTRL'
bdry_loc=['N', 'E', 'W']
out_dir = '/data/oceans_input/processed_input_data/CESM/PI_CTRL_obcs/'
file_year_start = 1920
file_year_end = 2100
ens = 1
overlap = 0
for i in range(4, 21):
    start_year = i * 100
    end_year = start_year + 99  
    [overlap, ens] = process_cesm_obcs_ts(expt, ens, bdry_loc, start_year, end_year, out_dir, overlap)
    #[overlap, ens] = process_cesm_obcs_non_ts(expt, ens, bdry_loc, start_year, end_year, out_dir, overlap)
    logging.info(ens, i*100)
    print(overlap, ens)
#process_cesm_obcs_ts (expt, ens, bdry_loc=['N', 'E', 'W'], start_year=None, end_year=None, out_dir='./'):