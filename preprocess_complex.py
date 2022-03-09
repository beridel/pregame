#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 18 17:07:47 2022

@author: mouchonc
"""

import numpy as np
import obspy as obs
import datetime


for i in range(len(st)):
    print('Start processing of trace ',i+1)
    tr = st[i]
    samp = tr.stats.sampling_rate
    nbpts = samp * 86400
    merged = obs.Trace(data=np.zeros(86400 * int(samp),dtype=np.float32))
    start_trace = tr.stats.strattime.date
    end_trace = start_trace + datetime.timedelta(days=1) - 1. / samp
    merged.stats.sampling_rate = samp
    merged.stats.delta = 1 / samp
    merged.stats.starttime = start_trace
    merged.stats.station = tr.stats.station
    merged.stats.channel = tr.stats.channel
    merged.stats.network = tr.stats.network
    merged.stats.location = tr.stats.location
    merged.mask = np.zeros(merged.stats.npts, dtype=np.int32)
    
    ## Look at gaps
    gap_list = tr.get_gaps()
    for gap in gap_list:
        duration = gap[6]
        
    
    ## Check the start date hours
    sec_check_starttime = tr.stats.starttime.second == 0
    min_check_starttime = tr.stats.starttime.minute == 0
    hour_check_starttime = tr.stats.starttime.hour == 0
    if np.all((sec_check_starttime,min_check_starttime,hour_check_starttime)):
        print('Start time already at 00:00:00')
    else:
        print('Start time calibration ...')
        maskArray = np.ma.masked(tr, mask=)
    ## Check the end date hours
    sec_check_endtime = tr.stats.endtime.second == 59
    min_check_endtime = tr.stats.endtime.minute == 59
    hour_check_endtime = tr.stats.endtime.hour == 23
        
    if np.all((sec_check_endtime,min_check_endtime,hour_check_endtime)):
        print('End time already at 23:59:59')
    else:
        print('End time calibration ...')