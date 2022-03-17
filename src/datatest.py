import obspy
from scipy import signal
from obspy import read
import numpy as np
import glob
import datatest_func as dt
import h5py
import obspy

#filelist = glob.glob('../data/*')
filelist = ('../data/N4.R58B..BHZ.2015.230')

for line in filelist:
    wv = read(line)
    wv_cp = wv.copy()

    dt.check_cont(wv_cp)
    

#save file as h5-format
    h5outfile = '../h5data/{}.h5'.format(line.split('/')[2])
    with h5py.File(h5outfile, 'w') as h5file:
         datapoints = wv_cp[0].data
         timecount = list(range(0,len(datapoints)))
         inputdata = [timecount, datapoints]
         h5file.create_dataset('time_data', data=inputdata)
         

#for line in filelist:
#    wv = read(line, headonly=True)
#    wv_cp = wv.copy()
#    print(wv_cp)

#    if len(wv_cp) > 1:
#       for i in range(len(wv_cp)):
#           wv_con += wv_cp[i]
#           print(wv_cp,'Data has been connected')
