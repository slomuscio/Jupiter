"""
eq must be a list of astropy SkyCoord objects in equitorial coordinates
gal must be a list of astropy SkyCoord objects in galactic functions

both can be obtained using coordinates.py functions equitorial and galactic

Returns a list of coordinates (eq and gal), times, and distances to bright sources only at times where the object in question is far enough away from the bright sources

Samantha Lomuscio
"""

from coordinates import equatorial
from astropy.coordinates import SkyCoord

def bright_source_filter(eq=None,gal=None,tstart=None,tend=None,tstep=None,times=None):
    
    skycoords_sun = equatorial(obj='sun', tstart=tstart, tend=tend, tstep=tstep)
    skycoords_moon = equatorial(obj='moon', tstart=tstart, tend=tend, tstep=tstep)


    #make list of SkyCoord objects for bright sources --- from Fermi LAT 8-Year Point Source Catalog                                                               
    source1 = SkyCoord('21h58m51.4s -30d13m30s')
    source2 = SkyCoord('04h57m02.6s -23d24m54s')
    source3 = SkyCoord('05h38m50.1s -44d05m10s')
    source4 = SkyCoord('07h21m57.2s +71d20m26s')
    source5 = SkyCoord('04h28m41.5s -37d56m25s')
    source6 = SkyCoord('12h56m10.0s -05d47m19s')
    source7 = SkyCoord('11h04m28.5s +38d12m25s')
    source8 = SkyCoord('15h12m51.5s -09d06m23s')
    source9 = SkyCoord('22h53m59.1s +16d09m02s')
    source10 = SkyCoord('18h36m13.3s +59d25m40s')

    skycoords_sources = [source1, source2, source3, source4, source5, source6, source7, source8, source9, source10]


    #find distances to sun, moon, bright sources                                    
    dist_sun = list()
    dist_moon = list()
    dist_sources = list()
    for i in range(0,len(eq)):
        dist_sun.append(eq[i].separation(skycoords_sun[i]))
        dist_moon.append(eq[i].separation(skycoords_moon[i]))
        
        for source in skycoords_sources:
            dist_sources.append(eq[i].separation(source))

    #zip galactic, equitorial, time, distances to sources (sun moon bright) info into one list
    data = zip(gal, eq, times, dist_sun, dist_moon, dist_sources)

    #remove items in list too close to galactic plane, sun, moon, or bright sources
    data2 = [item for item in data if  abs(item[0].b.degree) > 20.0 or item[3].degree > 20 or  item[4].degree > 20 or item[5] > 20]

    return data2
