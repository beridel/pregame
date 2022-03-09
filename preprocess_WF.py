from .config import config as cfg
from .dataset import DayData
import os
import numpy as np
import obspy as obs
import multiprocessing as mp
from functools import partial
import datetime as dt
from scipy import stats
from obspy.clients.fdsn import Client
import traceback

import ipdb

def preprocess(station, channel, component, files, parallel=False, remove_ir=False, client='IRIS'):
    if not isinstance(files, list):
        files = [files]

    start_timer = dt.datetime.now()

    header = obs.Stream()
    for file in files:
        read_stream = obs.read(file, headonly=True)
        for trace in read_stream.traces:
            trace.file = file

        header.extend(read_stream)
    
    print('Finished reading headers for {}:{} in {}s'.format(station, channel, (dt.datetime.now() - start_timer).total_seconds()))

    header = header.select(station=station, channel=channel)
    if len(header.traces) == 0:
        print('No files to process for {sta}:{cmp}!'.format(sta=station, cmp=component))
        return

    header.sort(keys=['starttime'])
    n_traces = len(header.traces)

    start_times = np.zeros(n_traces)
    end_times = np.zeros(n_traces)
    for i in range(n_traces):
        trace = header.traces[int(i)]
        start_times[i] = trace.stats.starttime.timestamp
        end_times[i] = trace.stats.endtime.timestamp

    start_date = obs.UTCDateTime(np.min(start_times)).date
    end_date = obs.UTCDateTime(np.max(end_times)).date

    print('Preprocessing {}:{} from {} to {}...'.format(station, channel, start_date, end_date))

    day_data = DayData(obs.UTCDateTime(np.min(start_times)).date, station, component)
    _write_dirs(day_data, raw=True)
    
    start_timer = dt.datetime.now()
    dates = []
    date = start_date
    while date <= end_date:
        dates.append(date)
        date += dt.timedelta(days=1)
   
    client = Client(client)
    if parallel:
        n_proc = np.min([cfg.n_proc, mp.cpu_count() / 2]) # takes a lot of memory
        pool = mp.Pool(processes=int(cfg.n_proc))
        pool.map_async(partial(_preprocess_wrap, header=header, station=station, channel=channel, component=component, remove_ir=remove_ir, client=client), dates)

        pool.close()
        pool.join()
    
    else:
        for date in dates:
            _preprocess(date, header, station, channel, component, remove_ir, client)

    print('Done in {}s!'.format((dt.datetime.now() - start_timer).total_seconds()))


def _preprocess_wrap(date, header, station, channel, component, remove_ir, client):
    try:
        _preprocess(date, header, station, channel, component, remove_ir, client)
    except:
        print('%s: %s' % (date, traceback.format_exc()))


def _preprocess(date, header, station, channel, component, remove_ir, client):
    start_timer = dt.datetime.now()
    start_comp = np.array([trace.stats['starttime'] < obs.UTCDateTime(date + dt.timedelta(days=1)) for trace in header.traces], dtype=np.bool)
    end_comp = np.array([trace.stats['endtime'] > obs.UTCDateTime(date) for trace in header.traces], dtype=np.bool)

    header_day = header.copy()
    header_day.traces = [trace for trace, test in zip(header.traces, np.logical_and(start_comp, end_comp)) if test]

    raw_data = obs.Stream()
    load_files = []
    for trace in header_day.traces:
        if any([file == trace.file for file in load_files]):
            continue
        else:
            raw_data.extend(obs.read(trace.file))
            load_files.append(trace.file)
    
    day_data = DayData(date, station, component)
    start_trace = obs.UTCDateTime(date)
    end_trace = start_trace + dt.timedelta(days=1) - 1. / cfg.sampling_rate
    
    load_time = (dt.datetime.now() - start_timer).total_seconds()

    # next bit is just to reseparate traces whose gaps were filled with dirty zeros
    nonzero_data = obs.Stream()
    for trace in raw_data.traces:
        if len(nonzero_data) > cfg.max_traces:
            break

        if np.sum(trace.data) == 0.:
            continue

        bounded_iszero = np.hstack([[0], trace.data == 0, [0]])
        iszero_diff = np.diff(bounded_iszero)
        zeros_start, = np.where(iszero_diff > 0)
        zeros_stop, = np.where(iszero_diff < 0)
        seq_len = zeros_stop - zeros_start
        which_seq, = np.where(seq_len > (1. / trace.stats.delta)) # remove zero sequences longer than 1 second

        if which_seq.size > 0:
            mask = np.zeros(trace.data.shape)
            for j in range(which_seq.size):
                start_mask = int(zeros_start[which_seq[j]])
                stop_mask = int(zeros_stop[which_seq[j]])
                mask[start_mask:stop_mask] = 1

            trace.data = np.ma.masked_where(mask, trace.data)
            nonzero_data.extend(trace.split())
        else:
            nonzero_data.extend([trace])

    nonzero_time = (dt.datetime.now() - start_timer).total_seconds() - load_time

    if len(nonzero_data) == 0:
        print('No data, skipping {} {}:{}!'.format(date, station, component))
        return

    if len(nonzero_data) > cfg.max_traces:
        print('More than {} gaps ({} gaps), skipping {} {}:{}!'.format(cfg.max_traces, len(nonzero_data), date, station, component))
        return

    deltas = np.zeros((len(nonzero_data.traces)))
    for i in range(deltas.size):
        deltas[i] = nonzero_data.traces[i].stats.delta

    delta = stats.mode(deltas).mode
    merged = obs.Trace(data=np.zeros(86400 * int(1. / delta),
                                    dtype=np.float32))
    merged.stats.sampling_rate = cfg.sampling_rate
    merged.stats.delta = delta
    merged.stats.starttime = start_trace
    merged.stats.station = station
    merged.stats.channel = channel
    merged.stats.network = header.traces[0].stats['network']
    merged.stats.location = header.traces[0].stats['location']
    merged.mask = np.zeros(merged.stats.npts, dtype=np.int32)
 
    for trace in nonzero_data.traces:
        to_merge = synchronize(trace, start_trace)
        to_merge = to_merge.slice(starttime=start_trace, endtime=end_trace)

        to_merge.mask = np.ones(to_merge.stats.npts,
                                dtype=np.int32)
        if to_merge.stats.delta == delta:
            merged = mergeStrict(merged, to_merge)

    merge_time = (dt.datetime.now() - start_timer).total_seconds() - load_time - nonzero_time
    
    percent_gap = float(np.sum(merged.mask)) / merged.mask.size
    percent_gap = (1 - percent_gap)
    
    if percent_gap < cfg.gap_limit:
        merged = fill_gaps(merged)
        
        to_write = obs.Stream()
        to_write.extend([merged])
        
        if remove_ir:
            inventory = client.get_stations(starttime=obs.UTCDateTime(date), endtime=obs.UTCDateTime(date + dt.timedelta(days=1)),
                network=header_day.traces[0].stats['network'], sta=station, channel=header_day.traces[0].stats['channel'], level="response")
            to_write.attach_response(inventory)
            to_write.remove_response(zero_mean=True, pre_filt=[cfg.pre_min1, cfg.pre_min2, cfg.pre_max1, cfg.pre_max2], taper=True, taper_fraction=(1. / cfg.pre_min2) / 86400)

        sampling_rate = 1. / delta
        if sampling_rate != cfg.sampling_rate:
            if np.fmod(sampling_rate, cfg.sampling_rate) == 0:
                to_write.decimate(int(sampling_rate / cfg.sampling_rate))
            else:
                print('Cannot perform simple decimation for {} {}:{} (now: {}Hz; target: {}Hz), skipping file!'.format(date, station, component, sampling_rate, cfg.sampling_rate))
                return

        preprocess_time = (dt.datetime.now() - start_timer).total_seconds() - load_time - nonzero_time - merge_time

        day_data.components = [component]
        day_data.traces = to_write
        day_data.operational = {station: True}
        day_data.loaded = True
        day_data.write(raw=True)

    else:
        print('Too much of a gap ({}%) for {} {}:{}, skipping file!)'.format(percent_gap, date, station, component))


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

    beginTimeSync = Tref + dt.timedelta(seconds=nbeg * delta)
    dt_toremove = (trace.stats.starttime - beginTimeSync)

    trace.data = timeshift(trace.data, dt_toremove / delta)
    trace.stats.starttime = beginTimeSync

    return trace


def mergeStrict(trace, to_merge):
    if (trace.stats.sampling_rate - to_merge.stats.sampling_rate) > 0.01:
        print('merge: different delta (trace={}, add={})'.format(
            trace.stats.sampling_rate, to_merge.stats.sampling_rate))
        return

    sr = trace.stats.sampling_rate

    ib_cpy = int(np.floor((max(trace.stats.starttime,
                               to_merge.stats.starttime) -
                           trace.stats.starttime) * sr + 0.5))

    if ib_cpy >= trace.stats.npts:
        return
    ib = int(np.floor((max(trace.stats.starttime,
                           to_merge.stats.starttime) -
                       to_merge.stats.starttime) * sr + 0.5))

    npts_cpy = min(to_merge.stats.npts - ib,
                   trace.stats.npts - ib_cpy)

    ie_cpy = ib_cpy + npts_cpy
    ie = ib + npts_cpy

    mask_nt = np.zeros(trace.stats.npts, dtype=np.int32)
    mask_nt[ib_cpy:ie_cpy] = to_merge.mask[ib:ie]

    trace.mask = np.bitwise_or(trace.mask, mask_nt)
    trace.data[np.nonzero(mask_nt)[0]] = np.float32(to_merge.data[np.nonzero(to_merge.mask[ib:ie])[0]])

    return trace


def fill_gaps(trace):
    maskArray = np.ma.array(trace.data, mask=trace.mask)
    slices = np.ma.extras.flatnotmasked_contiguous(maskArray)

    maskArray = np.ma.array(trace.data, mask=np.logical_not(trace.mask))

    if slices is not None:
        demeaned = trace.data - np.mean(trace.data)
        std = np.std(demeaned)

        for slice in slices:
            begin = slice.start-1
            end = slice.stop+1

            if begin < 0:
                yb = trace.data[end]
            else:
                yb = trace.data[begin]

            if end >= trace.stats.npts:
                ye = yb
            else:
                ye = trace.data[end]

            slope = (ye - yb) / (end - begin - 2)
            #noise = np.random.normal(0, 0.05 * std, end - begin - 2)
            #trace.data[slice] = slope * np.arange(0, end-begin-2) + yb + noise
            trace.data[slice] = slope * np.arange(0, end-begin-2) + yb

    trace.mask = np.ones(trace.stats.npts, dtype=np.int32)

    return trace


def prefilter(network, station_idx=None):
    """Preprocess all data files by prefiltering and decimating.
    """

    if not station_idx:
        stations = network.stations
    else:
        stations = network.stations[int(station_idx[0]):int(station_idx[1])]

    # make the necessary directories to avoid process collisions
    # unfiltered
    day_data = DayData(network.start_date,
                            stations,
                            network.components)
    _write_dirs(day_data)

    # filtered downsampled
    day_data = DayData(network.start_date,
                            stations,
                            network.components,
                            'filtered')
    _write_dirs(day_data)

    # generate dates to preprocess
    dates = network.datelist()

    pool = mp.Pool(processes=int(cfg.n_proc))
    pool.map_async(partial(_prefilter, network=network, stations=stations), dates)

    pool.close()
    pool.join()


def _write_dirs(day_data, raw=False):
    path = day_data.where(raw=raw)

    for station in day_data.stations:
        folder = os.path.join(path, station)
        cfg.chk_folder(folder)


def _prefilter(date, network, stations):
    starttime = obs.UTCDateTime(date)
    endtime = obs.UTCDateTime(date + dt.timedelta(days=1)) - 1. / cfg.sampling_rate

    for station in stations:
        idx = network.stations.index(station)
        unfiltered = DayData(date, station, network.components)
        unfiltered.read(raw=True)

        if not unfiltered.operational[station]:
            continue

        for i in range(len(network.components)):
            if 1. / unfiltered.traces[i].stats.delta != cfg.sampling_rate:
                sampling_rate = 1. / unfiltered.traces[i].stats.delta
                if sampling_rate.is_integer():
                    unfiltered.traces[i].decimate(int(sampling_rate / cfg.sampling_rate))
                else:
                    print('Cannot perform simple decimation for {} {}:{} (now: {}Hz; target: {}Hz), skipping file!'.format(date, raw_data.station, raw_data.component, sampling_rate, cfg.sampling_rate))
                    return
            
            unfiltered.traces[i].stats.delta = 1. / cfg.sampling_rate
            unfiltered.traces[i].trim(starttime=starttime,
                                      endtime=endtime)

        unfiltered.traces.detrend('demean')
        unfiltered.traces.detrend('linear')
        unfiltered.write()

        """for f in cfg.freq_bands:

                        # taper and filter
            filtered.traces.taper(
                np.float32((cfg.sampling_rate / np.float32(f[0])) / filtered.traces[0].data.size),
                type='cosine')

            if cfg.filter_type == 'bandpass':
                filtered.traces.filter('bandpass',
                                    freqmin=np.float32(f[0]),
                                    freqmax=np.float32(f[1]))
            elif cfg.filter_type == 'highpass':
                filtered.traces.filter('highpass', freq=np.float32(f[0]))
            elif cfg.filter_type == 'lowpass':
                filtered.traces.filter('lowpass', freq=np.float32(f[0]))
            """
        filtered = DayData(date, station, network.components)
        filtered.traces = unfiltered.traces.copy()
        filtered.operational = unfiltered.operational
        filtered.loaded = True

        filtered.traces.taper(
            np.float32((cfg.sampling_rate / network.bands[idx][0]) / filtered.traces[0].data.size),
            type='cosine')
        filtered.traces.filter('bandpass',
            freqmin=network.bands[idx][0], freqmax=network.bands[idx][1])

        filtered.set_band('filtered')
        filtered.write()


def remove_spikes(day_data):
    """Remove non-physical spikes from data and set them to zero
    """

    for station in day_data.stations:
        for component in day_data.components:
            trace = day_data.get_trace(station, component)
            abs_dev = np.abs(trace - np.median(trace))
            med_dev = np.median(abs_dev)
            norm_abs_dev = abs_dev / med_dev if med_dev else 0

            spike_idx = np.argwhere(norm_abs_dev > cfg.spike)
            spike_idx = spike_idx.flatten()

            if not spike_idx.any():
                continue

            for idx in spike_idx:
                start_idx = np.int32(idx - cfg.spike_buffer)
                end_idx = np.int32(idx + cfg.spike_buffer)

                if start_idx < 0:
                    start_idx = 0
                if end_idx >= trace.size:
                    end_idx = trace.size - 1

                start_value = trace[start_idx]
                end_value = trace[end_idx]

                n_samples = trace[start_idx:end_idx].size
                trace[start_idx:end_idx] = np.linspace(
                    start_value,
                    end_value,
                    n_samples)
                day_data.set_trace(trace, station, component)

    return day_data


def detect_peaks(x, mph=None, mpd=1, threshold=0, edge='rising',
                 kpsh=False, valley=False, show=False, ax=None):

    """Detect peaks in data based on their amplitude and other features.

    Parameters
    ----------
    x : 1D array_like
        data.
    mph : {None, number}, optional (default = None)
        detect peaks that are greater than minimum peak height.
    mpd : positive integer, optional (default = 1)
        detect peaks that are at least separated by minimum peak distance (in
        number of data).
    threshold : positive number, optional (default = 0)
        detect peaks (valleys) that are greater (smaller) than `threshold`
        in relation to their immediate neighbors.
    edge : {None, 'rising', 'falling', 'both'}, optional (default = 'rising')
        for a flat peak, keep only the rising edge ('rising'), only the
        falling edge ('falling'), both edges ('both'), or don't detect a
        flat peak (None).
    kpsh : bool, optional (default = False)
        keep peaks with same height even if they are closer than `mpd`.
    valley : bool, optional (default = False)
        if True (1), detect valleys (local minima) instead of peaks.
    show : bool, optional (default = False)
        if True (1), plot data in matplotlib figure.
    ax : a matplotlib.axes.Axes instance, optional (default = None).

    Returns
    -------
    ind : 1D array_like
        indeces of the peaks in `x`.

    Notes
    -----
    The detection of valleys instead of peaks is performed internally by simply
    negating the data: `ind_valleys = detect_peaks(-x)`

    The function can handle NaN's

    See this IPython Notebook [1]_.

    References
    ----------
    .. [1]:http://nbviewer.ipython.org/github/demotu/BMC/blob/master/
        notebooks/DetectPeaks.ipynb

    Examples
    --------
    >>> from detect_peaks import detect_peaks
    >>> x = np.random.randn(100)
    >>> x[60:81] = np.nan
    >>> # detect all peaks and plot data
    >>> ind = detect_peaks(x, show=True)
    >>> print(ind)

    >>> x = np.sin(2*np.pi*5*np.linspace(0, 1, 200)) + np.random.randn(200)/5
    >>> # set minimum peak height = 0 and minimum peak distance = 20
    >>> detect_peaks(x, mph=0, mpd=20, show=True)

    >>> x = [0, 1, 0, 2, 0, 3, 0, 2, 0, 1, 0]
    >>> # set minimum peak distance = 2
    >>> detect_peaks(x, mpd=2, show=True)

    >>> x = np.sin(2*np.pi*5*np.linspace(0, 1, 200)) + np.random.randn(200)/5
    >>> # detection of valleys instead of peaks
    >>> detect_peaks(x, mph=0, mpd=20, valley=True, show=True)

    >>> x = [0, 1, 1, 0, 1, 1, 0]
    >>> # detect both edges
    >>> detect_peaks(x, edge='both', show=True)

    >>> x = [-2, 1, -2, 2, 1, 1, 3, 0]
    >>> # set threshold = 2
    >>> detect_peaks(x, threshold = 2, show=True)
    """

    x = np.atleast_1d(x).astype('float64')
    if x.size < 3:
        return np.array([], dtype=int)
    if valley:
        x = -x
    # find indices of all peaks
    dx = x[1:] - x[:-1]
    # handle NaN's
    indnan = np.where(np.isnan(x))[0]
    if indnan.size:
        x[indnan] = np.inf
        dx[np.where(np.isnan(dx))[0]] = np.inf
    ine, ire, ife = np.array([[], [], []], dtype=int)
    if not edge:
        ine = np.where((np.hstack((dx, 0)) < 0) & (np.hstack((0, dx)) > 0))[0]
    else:
        if edge.lower() in ['rising', 'both']:
            ire = np.where((np.hstack((dx, 0)) <= 0) &
                           (np.hstack((0, dx)) > 0))[0]
        if edge.lower() in ['falling', 'both']:
            ife = np.where((np.hstack((dx, 0)) < 0) &
                           (np.hstack((0, dx)) >= 0))[0]
    ind = np.unique(np.hstack((ine, ire, ife)))
    # handle NaN's
    if ind.size and indnan.size:
        # NaN's and values close to NaN's cannot be peaks
        ind = ind[np.in1d(ind, np.unique(np.hstack((indnan,
                                                    indnan - 1, indnan + 1))),
                          invert=True)]
    # first and last values of x cannot be peaks
    if ind.size and ind[0] == 0:
        ind = ind[1:]
    if ind.size and ind[-1] == x.size - 1:
        ind = ind[:-1]
    # remove peaks < minimum peak height
    if ind.size and mph is not None:
        ind = ind[x[ind] >= mph]
    # remove peaks - neighbors < threshold
    if ind.size and threshold > 0:
        dx = np.min(np.vstack([x[ind] - x[ind - 1], x[ind] - x[ind + 1]]),
                    axis=0)
        ind = np.delete(ind, np.where(dx < threshold)[0])
    # detect small peaks closer than minimum peak distance
    if ind.size and mpd > 1:
        ind = ind[np.argsort(x[ind])][::-1]  # sort ind by peak height
        idel = np.zeros(ind.size, dtype=bool)
        for i in range(ind.size):
            if not idel[i]:
                # keep peaks with the same height if kpsh is True
                idel = idel | (ind >= ind[i] - mpd) & (ind <= ind[i] + mpd) \
                    & (x[ind[i]] > x[ind] if kpsh else True)
                idel[i] = 0  # Keep current peak
        # remove the small peaks and sort back the indices by their occurrence
        ind = np.sort(ind[~idel])

    if show:
        if indnan.size:
            x[indnan] = np.nan
        if valley:
            x = -x
        _plot_peaks(x, mph, mpd, threshold, edge, valley, ax, ind)

    return ind


def _plot_peaks(x, mph, mpd, threshold, edge, valley, ax, ind):
    """Plot results of the detect_peaks function, see its help."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print('matplotlib is not available.')
    else:
        if ax is None:
            _, ax = plt.subplots(1, 1, figsize=(8, 4))

        ax.plot(x, 'b', lw=1)
        if ind.size:
            label = 'valley' if valley else 'peak'
            label = label + 's' if ind.size > 1 else label
            ax.plot(ind, x[ind], '+', mfc=None, mec='r', mew=2, ms=8,
                    label='%d %s' % (ind.size, label))
            ax.legend(loc='best', framealpha=.5, numpoints=1)
        ax.set_xlim(-.02 * x.size, x.size * 1.02 - 1)
        ymin, ymax = x[np.isfinite(x)].min(), x[np.isfinite(x)].max()
        yrange = ymax - ymin if ymax > ymin else 1
        ax.set_ylim(ymin - 0.1 * yrange, ymax + 0.1 * yrange)
        ax.set_xlabel('Data #', fontsize=14)
        ax.set_ylabel('Amplitude', fontsize=14)
        mode = 'Valley detection' if valley else 'Peak detection'
        ax.set_title("%s (mph=%s, mpd=%d, threshold=%s, edge='%s')"
                     % (mode, str(mph), mpd, str(threshold), edge))
        # plt.grid()
        plt.show()

