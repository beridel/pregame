import numpy as np
import obspy as obs
import datetime

st1 = obs.read('/Users/mouchonc/Documents/RESEARCH/Seismic_processing/altered_data/US.CBN.00.BH1.2015.230')
st2 = obs.read('/Users/mouchonc/Documents/RESEARCH/Seismic_processing/altered_data/US.CBN.00.BH2.2015.230')
st3 = obs.read('/Users/mouchonc/Documents/RESEARCH/Seismic_processing/altered_data/US.CBN.00.BHZ.2015.230')
st = st1 + st2 + st3
detrended = st.copy()
print(st)
detrended.detrend("linear")
st.plot()
detrended.plot()
merged = detrended.copy()
merged.merge(method=1, fill_value='interpolate', interpolation_samples=-1)
merged.plot()



# time = np.zeros((1,86400))
# print(tr.stats.starttime.time == datetime.time(0,0))
# start_day = tr.stats.starttime.date
# end_day = tr.stats.endtime.date
# st.plot()
# st1.plot()
# st2.plot()
# st3.plot()
#st = file.Stream
#ch = file.channel
 
