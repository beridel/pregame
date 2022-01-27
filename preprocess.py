import numpy as np
import obspy as obs
import os



#for trace in str.traces:
#    trace.file = file
network = 'TO'
station = 'ACAH'
channel = '01'
component = 'HHZ'
path = pwd + station

file_list = create_file_list(path)

def preprocess(network, station, channel, component, files):
    #file_list = create_file_list(path)
    header = reading_headers(files)















def reading_headers(file_list):
    header = obs.Stream() #set header as an empty stream object
    for file in files_list
        str = obs.read(file, headonly=True)
        header.extend(str)
    print('Finished reading {} headers for {} comp {}'.format(len(header),station,component))

def create_file_list(path):
    file_list = os.listdir(path)
    print('Finished creating file list for {} total files'.format(len(file_list)))
