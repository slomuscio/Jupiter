import math
import matplotlib.pyplot as plt
import numpy as np 
from astropy.io import fits 
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord 
from gt_apps import filter, evtbin, maketime
from coordinates import galactic, equatorial
from bright_sources_filter import bright_source_filter

#files
data_file = '@jupiter_data.txt'
sc_file = 'L1907011408184582BCD965_SC00.fits'

with open('jupiter_data.txt') as j:
     files = j.readlines()
file_start = files[0][:-1]
file_end = files[-1][:-1] 

#get parameters from Fermi photon event file
hdulist_start = fits.open(file_start)
hdu_start = hdulist_start[1]

hdulist_end = fits.open(file_end)
hdu_end = hdulist_end[1]

roi_data = (hdu_start.header['DSVAL3']).split(',')
roi = int(float(((roi_data[2])[:-1]))) #degrees

time_start = Time(hdu_start.header['DATE-OBS'], format='fits') #used to get coords, fits format
time_end = Time(hdu_end.header['DATE-END'], format='fits') #used to get coords, fits format

time_step = 0.5 #in units of days

time_start_MET = hdu_start.header['TSTART'] #used in gtselect MET (s)
time_end_MET = hdu_end.header['TSTOP'] #used in gtselect MET (s)

times = np.arange(time_start_MET, time_end_MET, time_step*(24*60*60), dtype=float)

energyrange = (hdu_start.header['DSVAL5']).split(':')
emin = float(energyrange[0]) #MeV
emax = float(energyrange[1]) #MeV
num_bins = int(math.log10(emax/emin)*10)

zmax = 90 #deg

bin_size = 0.1 #deg/pix
dim = (2*roi)/bin_size #gives dimension of cmap diameter/bin size = num pixels

hdulist_start.close()
hdulist_end.close()


#make array that will hold stacked counts map data
total_cmap_data = np.zeros((int(dim),int(dim)), dtype = float)


#make list of SkyCoord objects with coordinate data
skycoords_eq = equatorial(obj='jupiter', tstart=time_start, tend=time_end, tstep=time_step) #used to find ra and dec
skycoords_gal = galactic(obj='jupiter', tstart=time_start, tend=time_end, tstep=time_step) #used to remove points too close to galactic plane

#removes times/positions when object is near bright sources (sun, moon, others)
data = bright_source_filter(eq=skycoords_eq, gal=skycoords_gal, tstart=time_start,tend=time_end, tstep=time_step, times=times)

#loop through times to create stacked counts map
for i in range(1,len(data)):
     ra = data[i-1][1].ra.degree
     dec = data[i-1][1].dec.degree

     tmin = data[i][2]
     tmax = data[i][2]
     
     if (tmax-tmin) > time_step:
          continue
     
     else:
          #gtselect
          filter['evclass'] = 128
          filter['evtype'] = 3
          filter['rad'] = 'INDEF'
          filter['emin'] = emin
          filter['emax'] = emax
          filter['zmax'] = zmax
          filter['infile'] = data_file
          filter['ra'] = ra
          filter['dec'] = dec
          filter['tmin'] = times[i-1] #must be in MET (s)
          filter['tmax'] = times[i]
          outfile = 'JUPITER_binned_filtered_' + str(i) + '.fits'
          filter['outfile'] = outfile
          filter.run()

          #gtmktime
          maketime['scfile'] = sc_file
          maketime['filter'] = '(DATA_QUAL==1)&&(LAT_CONFIG==1)'
          maketime['roicut'] = 'yes'
          maketime['evfile'] = outfile
          gti_file = 'JUPITER_binned_gti_'+ str(i) + '.fits'
          maketime['outfile'] = gti_file
          maketime.run()

          #counts map
          evtbin['scfile'] = sc_file
          evtbin['nxpix'] = int(dim)
          evtbin['nypix'] = int(dim)
          evtbin['binsz'] = bin_size
          evtbin['axisrot'] = 0
          evtbin['coordsys'] = 'CEL'
          evtbin['proj'] = 'AIT'
          evtbin['algorithm'] = 'CMAP'
          evtbin['evfile'] = gti_file
          cmap_file = 'JUPITER_binned_cmap_'+ str(i) + '.fits'
          evtbin['outfile'] = cmap_file
          evtbin['xref'] = ra
          evtbin['yref'] = dec
          evtbin.run()
     
          #stack each cmap 
          total_cmap_data += fits.getdata(cmap_file)


cmap_outfile = 'stacked_cmap.fits'
hdu = fits.PrimaryHDU(total_cmap_data)
hdu.writeto(cmap_outfile, overwrite=True)

plt.imshow(total_cmap_data, cmap='gray')
plt.colorbar()
plt.show()

