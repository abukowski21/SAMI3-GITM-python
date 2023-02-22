from datetime import datetime, timedelta
import pandas as pd, numpy as np, os, sys
import matplotlib.pyplot as plt
from tqdm import tqdm


path = "/home/axb170054/scratch/GITM-testing/test_folders/step_function_driving/SAMI3-stretch/"
t0 = datetime(2011,5,20)





processed_path = 'processed_files'


geo_grid_files = {'glat':'glatu.dat','glon':'glonu.dat','alt':'zaltu.dat','mlat':'blatu.dat','mlon':'blonu.dat','malt':'baltu.dat'}

data_files = {'edens':'deneu.dat', 'hplusdens':'deni1u.dat','oplusdens':'deni2u.dat', 'noplusdens':'deni3u.dat', 'o2plusdens':'deni4u.dat',
              'heplusdens':'deni5u.dat', 'n2plusdens':'deni6u.dat', 'nplusdens':'deni7u.dat','hdens':'denn1u.dat','odens':'denn2u.dat', 'nodens':'denn3u.dat',
              'o2dens':'denn4u.dat', 'hedens':'denn5u.dat', 'n2dens':'denn6u.dat', 'ndens':'denn7u.dat', 'vsi1u':'vsi1u.dat',
              'vsi2u':'vsi2u.dat','u1pu':'u1pu.dat','u3hu':'u3hu.dat','u1u':'u1u.dat','u2u':'u2u.dat'}

time_file = 'time.dat'




g_alt_grid_df = pd.DataFrame()

for f in geo_grid_files:
    file = open(os.path.join(path, geo_grid_files[f]), 'rb')
    raw = np.fromfile(file, dtype='float32')
    file.close()
    
    listy = raw[1:-1]
    g_alt_grid_df[f] = listy
    
    
grid_fname = os.path.join(path,'processed_files', 'GRID_FILE.csv')
try:
    if not os.path.isfile(grid_fname):
        g_alt_grid_df.to_csv(grid_fname)
except:
    os.mkdir(os.path.join(path,'processed_files'))
    g_alt_grid_df.to_csv(grid_fname)



###TIME STUFF

times = pd.read_fwf(os.path.join(path, time_file), names = ['istep', 'hour','minute','second', 'hrdelta'], infer_nrows=1150)
times.pop('istep');

times_list = []
for hr in times['hrdelta']:
    times_list.append(t0 + timedelta(hours = hr))


output_col_map = {'deneu.dat': 'edens', 'deni1u.dat': 'hplusdens', 'deni2u.dat': 'oplusdens', 'deni3u.dat': 'noplusdens', 'deni4u.dat': 'o2plusdens',
 'deni5u.dat': 'heplusdens', 'deni6u.dat': 'n2plusdens', 'deni7u.dat': 'nplusdens', 'denn1u.dat': 'hdens', 'denn2u.dat': 'odens', 'denn3u.dat': 'nodens',
 'denn4u.dat': 'o2dens', 'denn5u.dat': 'hedens', 'denn6u.dat': 'n2dens', 'denn7u.dat': 'ndens', 'vsi1u.dat':'vsi1u',
              'vsi2u.dat':'vsi2u','u1pu.dat':'u1pu','u3hu.dat':'u3hu','u1u.dat':'u1u','u2u.dat':'u2u'}



    
##GET DF_any column
def process_all_vars(file_dict, g_alt_grid_df, path, debug = False):

    fnames = []
    files = {}
    outlist = []
    
    for var in file_dict:
        fnames.append(os.path.join(path, data_files[var]))

    for i, varname in enumerate(file_dict.keys()):
        files[file_dict[varname]] = open(fnames[i], 'rb')

    for i, dtime in enumerate(tqdm(times_list)):
       # print(i, dtime) # either this or tqdm - i think this is better

        new_df_per_time = g_alt_grid_df

        for k in files.keys():
            
            raw = np.fromfile(files[k], dtype='float32', count = len(g_alt_grid_df)+2)
            listy = raw[1:-1]

            new_df_per_time[output_col_map[k]] = listy
            list_of_dtimes = [dtime for x in range(0,len(listy))]
            
        
        dtime_str = str(dtime.date())+'_'+str(dtime.time()).replace(':','-')+'.h5'
        #dtime_str = str(dtime)
        i_fname = os.path.join(path,'processed_files',dtime_str)
        #new_df_per_time.to_csv(i_fname, index = None, compression= 'gzip')
        
        #return new_df_per_time
        new_df_per_time.to_hdf(i_fname, key='df', mode='w')
        
        outlist.append(i_fname)
#         #break

    for k in files.keys():
            files[k].close()
    return outlist



b = process_all_vars(data_files,g_alt_grid_df,path)


print('done! files written are \n', b)