"""
Generates files for fake Jupiter with 40 degrees subtracted from right ascension

Samantha Lomuscio
"""

import math
import pandas as pd
import os
import numpy as np 
from astropy.time import Time
from coordinates import equatorial 
from gt_apps import filter, evtbin, maketime

roi = 20
emin = 100
emax = 10000
num_bins = int(math.log10(emax/emin)*10)
zmax = 90 #deg

#read in data from file
info = pd.read_csv('DataRanges.txt', delimiter=',')

for i in range(len(info)):
    #set number
    setnum = info['SetNumber'][i]

    #make a directory for this set 
    os.mkdir('FakeJupiter-40/'+setnum)

    #make a list of photon files to be used
    start_week = info['StartWeek'][i]
    end_week = info['EndWeek'][i]

    f = open('FakeJupiter-40/' + setnum + '/photon_file.txt','w+')

    #write locations of photon files to a text file
    if start_week != end_week: 
        weeks = list(range(int(start_week),int(end_week)))
        weeks.append(end_week)

        for val in weeks:
            if val < 100:
                f.write('/data/scratch/ysong/fdata/weeklyphoton/lat_photon_weekly_w0' + str(val) + '_p305_v001.fits\n')
            else:
                f.write('/data/scratch/ysong/fdata/weeklyphoton/lat_photon_weekly_w' + str(val) + '_p305_v001.fits\n')

    elif start_week == end_week:
        if start_week < 100:
            f.write('/data/scratch/ysong/fdata/weeklyphoton/lat_photon_weekly_w0' + str(start_week) + '_p305_v001.fits\n')
        else:
	    f.write('/data/scratch/ysong/fdata/weeklyphoton/lat_photon_weekly_w' + str(start_week) + '_p305_v001.fits\n')

    f.close()

    data_file = '@FakeJupiter-40/' + setnum + '/photon_file.txt' 

    time_start = Time(info['StartDate'][i]) 
    time_end = Time(info['EndDate'][i])

    time_start_MET = int(info['StartDateMET'][i])
    time_end_MET = int(info['EndDateMET'][i])

    time_step = 1 #in units of days

    times = np.arange(time_start_MET, time_end_MET, time_step*(24*60*60), dtype=float)
    times = np.append(times, time_end_MET)

    skycoords_eq = equatorial(obj='jupiter', mode='c', tstart=time_start, tend=time_end, tstep=time_step) #used to find ra and dec

    #loop through times to create stacked counts map
    for i in range(1, len(skycoords_eq) + 1):
        ra = skycoords_eq[i-1].ra.degree - 40
        dec = skycoords_eq[i-1].dec.degree
        
        if ra < 0:
	    ra = 360 - abs(ra) 

        tmin = times[i-1]
        tmax = times[i]
     
        #gtselect
        filter['evclass'] = 128
        filter['evtype'] = 3
        filter['rad'] = roi
        filter['emin'] = emin
        filter['emax'] = emax
        filter['zmax'] = zmax
        filter['infile'] = data_file
        filter['ra'] = ra
        filter['dec'] = dec
        filter['tmin'] = tmin #must be in MET (s)
        filter['tmax'] = tmax
        outfile = 'FakeJupiter-40/' + setnum + '/fakejupiter-40_' + setnum + '_' + str(i) + '.fits'
        filter['outfile'] = outfile
        filter.run()
