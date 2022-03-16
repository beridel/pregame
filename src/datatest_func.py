import obspy
from scipy import signal
from obspy import read
import numpy as np
import glob
from obspy.core import UTCDateTime

# check if the contuinuous data is separated
def check_cont(wv):    
    if len(wv) != 1:
       for i in range(len(wv)):
           end_first = wv[i].stats.endtime
           start_second = wv[i+1].stats.starttime
           if end_first > start_second: #means overlapping 
               

    else:
       cont_wv = wv_cp

    return cont_wv


# make a waveform detrended linearly
def pre_linear(wv):
    wv.detrend('linear')
    wv.merge('method=0, fill_value=0')

    return wv
