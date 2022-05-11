import numpy as np
import obspy as obs
import datetime

path = '/Users/mouchonc/Documents/GitHub/pregame/altered_data/'


def nextpow2(n):
    m_f = np.log2(n)
    m_i = np.ceil(m_f)
    return 2 ** m_i

def timeshift(raw_sig, dt):
    if dt == 0.:
        return raw_sig

    npts = len(raw_sig)
    nptsPow2 = nextpow2(int(npts))

    if nptsPow2 - npts < 202:
        nptsPow2 *= 2

    nptsPow2 = max(nptsPow2, 256)

    sig = np.zeros(int(nptsPow2))
    nbeg = 100
    sig[nbeg:nbeg+npts] = raw_sig

    a = raw_sig[0]
    b = raw_sig[-1]

    sig[:nbeg] = a
    sig[nbeg + npts:] = b

    sig_fft = np.fft.rfft(sig)
    freq = np.fft.fftfreq(len(sig))[:len(sig) // 2 + 1]

    return np.fft.irfft(sig_fft * np.exp(-2 * np.pi * 1j *
                                         freq * dt))[nbeg:nbeg + npts]

def synchronize(trace, Tref):
    sr = trace.stats.sampling_rate
    delta = 1. / sr

    nbeg = np.floor((trace.stats.starttime - Tref) * sr + 0.5)

    beginTimeSync = Tref + datetime.timedelta(seconds=nbeg * delta)
    dt_toremove = (trace.stats.starttime - beginTimeSync)

    trace.data = timeshift(trace.data, dt_toremove / delta)
    trace.stats.starttime = beginTimeSync

    return trace


# Load stream
st = obs.read(path + 'US.CBN.00.BH1.2015.230')
st += obs.read(path + 'US.CBN.00.BH2.2015.230')
st += obs.read(path + 'US.CBN.00.BHZ.2015.230')

 ## Detrend the data
detrended = st.copy()
detrended.detrend("linear")

## Merge the data
merged = detrended.copy()
merged.merge(method=1, fill_value='interpolate', interpolation_samples=-1)

## Trim the data
trimed = merged.copy()

## Print Stream
print(st)
print(detrended)
print(merged)

## Plot Stream
st.plot()
detrended.plot()
merged.plot()


#%% Check the true day
for i in range(len(merged)):
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
        num = (t2 - t1) / 60
        date_list = [t1 + datetime.timedelta(minutes=x) for x in range(int(num))]
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
print(merged)

delta =  23 * 3600 + 59 * 60 + 59.999999 # in sec
true_tstart = obs.UTCDateTime(true_day)
true_tend = obs.UTCDateTime(true_day + datetime.timedelta(days=1)) #+ delta
trimed.trim(true_tstart, true_tend, pad=True, fill_value=0)
trimed[2].trim(true_tstart, true_tend, pad=True, fill_value=0)
trimed.plot()

print(trimed)
