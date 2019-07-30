"""
eq must be a list of astropy SkyCoord objects in equitorial coordinates

gal must be a list of astropy SkyCoord objects in galactic functions

Both can be obtained and entered as arguments here using coordinates.py functions equitorial and galactic or this function can generate them for you.

Returns a list of coordinates (eq and gal), times, and distances to bright sources only at times where the object in question is far enough away from the bright sources.

Samantha Lomuscio
"""

from coordinates import equatorial, galactic
from astropy.coordinates import SkyCoord

def plane_and_sources_filter(obj=None, frame=None, tstart=None, tend=None, tstep=None, eq=None, gal=None, mode=None): #, times=None):

    if eq is None and gal is None:
        eq = equatorial(obj=obj, frame=frame, tstart=tstart, tend=tend, tstep=tstep, mode=mode)
        gal = galactic(obj=obj, frame=frame, tstart=tstart, tend=tend, tstep=tstep, mode=mode)

    if eq is None and gal is not None:
        eq = equatorial(obj=obj, frame=frame, tstart=tstart, tend=tend, tstep=tstep, mode=mode)

    if eq is not None and gal is None:
        gal = galactic(obj=obj, frame=frame, tstart=tstart, tend=tend, tstep=tstep, mode=mode)

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

    dist_1 = list()
    dist_2 = list()
    dist_3 = list()
    dist_4 = list()
    dist_5 = list()
    dist_6 = list()
    dist_7 = list()
    dist_8 = list()
    dist_9 = list()
    dist_10 = list()

    for i in range(0,len(eq)):
        dist_1.append(eq[i].separation(source1).degree)
        dist_2.append(eq[i].separation(source2).degree)
        dist_3.append(eq[i].separation(source3).degree)
        dist_4.append(eq[i].separation(source4).degree)
        dist_5.append(eq[i].separation(source5).degree)
        dist_6.append(eq[i].separation(source6).degree)
        dist_7.append(eq[i].separation(source7).degree)
        dist_8.append(eq[i].separation(source8).degree)
        dist_9.append(eq[i].separation(source9).degree)
        dist_10.append(eq[i].separation(source10).degree)

    data = zip(eq, gal, dist_1, dist_2, dist_3, dist_4, dist_5, dist_6, dist_7, dist_8, dist_9, dist_10)
    data2 = list()

    for item in data:
        if abs(item[1].b.degree) < 20.0:
            continue
        elif item[2] < 10.0:
            continue
        elif item[3] < 10.0:
            continue
        elif item[4] < 10.0:
            continue
        elif item[5] < 10.0:
            continue
        elif item[6] < 10.0:
            continue
        elif item[7] < 10.0:
            continue
        elif item[8] < 10.0:
            continue
        elif item[9] < 10.0:
            continue
        elif item[10] < 10.0:
            continue
        elif item[11] < 10.0:
            continue
        else:
            data2.append(item)

    return data2




def sun_moon_filter(obj=None, frame=None, tstart=None, tend=None, tstep=None, eq=None, gal=None, mode=None):

    if eq is None and gal is None:
        eq = equatorial(obj=obj, frame=frame, tstart=tstart, tend=tend, tstep=tstep, mode=mode)
        gal = galactic(obj=obj, frame=frame, tstart=tstart, tend=tend, tstep=tstep, mode=mode)

    if eq is None and gal is not None:
        eq = equatorial(obj=obj, frame=frame, tstart=tstart, tend=tend, tstep=tstep, mode=mode)

    if eq is not None and gal is None:
        gal = galactic(obj=obj, frame=frame, tstart=tstart, tend=tend, tstep=tstep, mode=mode)

    skycoords_sun = equatorial(obj='sun', frame=frame, tstart=tstart, tend=tend, tstep=tstep, mode=mode)
    skycoords_moon = equatorial(obj='moon', frame=frame, tstart=tstart, tend=tend, tstep=tstep, mode=mode)

    dist_sun = list()
    dist_moon = list()

    for i in range(0,len(eq)):
        dist_sun.append(eq[i].separation(skycoords_sun[i]).degree)
        dist_moon.append(eq[i].separation(skycoords_moon[i]).degree)

    data = zip(eq, gal, dist_sun, dist_moon)
    data2 = list()

    for item in data:
        if item[2] < 20.0:
            continue
        elif item[3] < 5.0:
            continue
        else:
            data2.append(item)

    return data2
