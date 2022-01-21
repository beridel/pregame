import numpy as np
import obspy as obs
import os

def create_file_list(path):
    file_list = os.listdir(path)
    print('Finished creating file list for {} total files'.format(len(file_list)))

def reading_headers(file_list, station, ):
    header = obs.Stream() #set header as an empty stream object
    for file in files_list
        str = obs.read(file, headonly=True)
        for trace in str.traces:
            trace.file = file
        header.extend(str)
    print('Finished reading headers')
