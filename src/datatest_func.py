import obspy
from scipy import signal
from obspy import read
import numpy as np
import glob
from obspy.core import UTCDateTime

# check if the contuinuous data is separated
def check_cont(wv):    
    delta = wv[0].stats.delta
    if len(wv) != 1:
       for i in range(int(len(wv)-1)):
           t1_start = wv[i].stats.starttime
           t1_end = wv[i].stats.endtime
           t2_start = wv[i+1].stats.starttime
           t2_end = wv[i+1].stats.endtime

           dt = t1_end - t2_start
           
           if dt > 0: #means overlapping 
              seg_pre = wv[i].trim(t1_start, t1_end - dt)
              seg_aft = wv[i+1].trim(t2_start + dt, t2_end)
              
              


              #seg_mid = np.zeros((int(len(dt)/delta)))


    else:
       cont_wv = wv_cp

    return cont_wv


# make a waveform detrended linearly
def pre_linear(wv):
    wv.detrend('linear')
    wv.merge('method=0, fill_value=0')

    return wv
