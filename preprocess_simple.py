import numpy as np
import obspy as obs
import datetime

path = '/Users/mouchonc/Documents/GitHub/pregame/altered_data/'

# Load stream
st1 = obs.read(path + 'US.CBN.00.BH1.2015.230')
st2 = obs.read(path + 'US.CBN.00.BH2.2015.230')
st3 = obs.read(path + 'US.CBN.00.BHZ.2015.230')
st = st1 + st2 + st3
detrended = st.copy()
print(st)

 ## Detrend the data
detrended.detrend("linear")
st.plot()
detrended.plot()

## Merge the data
merged = detrended.copy()
merged.merge(method=1, fill_value='interpolate', interpolation_samples=-1)
merged.plot()

## Check the true day
for i in range(merged):
    tr = merged[i]
    start_day = tr.stats.starttime.date
    end_day = tr.stats.endtime.date

    test = start_day == end_day

    if test:
        true_day = start_day
    else:  ## Need to create a functionm for that
        # 1: create a list of datetime between t1 = tr.stats.starttime and t2 = tr.stats.endtime
        t1 = tr.stats.starttime
        t2 = tr.stats.endtime
        num = (tend - base) / 60
        date_list = [base + datetime.timedelta(minutes=x) for x in range(int(num))]
        date1 = []
        date2 = []
        # 2: Count the number of time you have day1 and day2 in the list
        for i in range(len(date_list)):
            date1.append(date_list[i].date == start_day)
            date2.append(date_list[i].date == end_day)
            if date1[i] == date2[i]:
                print("Error: dates are the same!")
        count1 = date1.count(True)
        count2 = date2.count(True)
        # 3: Evaluate which date is the most present in the list and define it as the true date
        if count1 == max(count1,count2):
            true_day = start_day
        elif count2 == max(count1,count2):
            true_day = end_day


## Trim the data
trimed = merged.copy()
delta =  23 * 3600 + 59 * 60 + 59.999999 # in sec
true_tstart = obs.UTCDateTime(true_day)
true_tend = obs.UTCDateTime(true_day) + delta
trimed.trim(true_tstart, true_tend, pad=True, fill_value=0)
trimed.plot()
