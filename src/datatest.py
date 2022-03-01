import obspy
from obspy import read
import numpy as np

wv = read('../data/*')
wv = wv.copy()
print(wv)
#wv.plot()

wv.detrend('linear')
wv.merge('method=0, fill_value=0')
