import numpy as np
import string

import pylab
import matplotlib.pyplot as pl
import matplotlib


__doc__ = """

gpi_utils

Various low level utils used in the GPI python code.
Pretty much a catch-all miscellaneous bin.

"""

# Functions for JD <-> GD conversion,
# courtesy of Ian Crossfield at http://www.astro.ucla.edu/~ianc/python/_modules/date.html
"""
Functions for handling dates.

Contains:
   gd2jd  -- converts gregorian date to julian date
   jd2gd  -- converts julian date to gregorian date

Wish list:
   Function to convert heliocentric julian date!



These functions were taken from Enno Middleberg's site of useful
astronomical python references:
http://www.astro.rub.de/middelberg/python/python.html

"Feel free to download, use, modify and pass on these scripts, but
please do not remove my name from it." --E. Middleberg
"""

# 2009-02-15 13:12 IJC: Converted to importable function


def gd2jd(*date): #, verbose=False):
    """gd2jd.py converts a UT Gregorian date to Julian date.

    Usage: gd2jd.py (2009, 02, 25, 01, 59, 59)

    To get the current Julian date:
        import time
        gd2jd(time.gmtime())

    Hours, minutes and/or seconds can be omitted -- if so, they are
    assumed to be zero.

    Year and month are converted to type INT, but all others can be
    type FLOAT (standard practice would suggest only the final element
    of the date should be float)
    """
    verbose=False
    if verbose: print(date)

    date = list(date)

    if len(date)<3:
        print("You must enter a date of the form (2009, 02, 25)!")
        return -1
    elif len(date)==3:
        for ii in range(3): date.append(0)
    elif len(date)==4:
        for ii in range(2): date.append(0)
    elif len(date)==5:
        date.append(0)

    yyyy = int(date[0])
    mm = int(date[1])
    dd = float(date[2])
    hh = float(date[3])
    min = float(date[4])
    sec = float(date[5])

    UT=hh+min/60+sec/3600


    total_seconds=hh*3600+min*60+sec
    fracday=total_seconds/86400

    if (100*yyyy+mm-190002.5)>0:
        sig=1
    else:
        sig=-1

    JD = 367*yyyy - int(7*(yyyy+int((mm+9)/12))/4) + int(275*mm/9) + dd + 1721013.5 + UT/24 - 0.5*sig +0.5

    months=["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"]


    # Now calculate the fractional year. Do we have a leap year?
    daylist=[31,28,31,30,31,30,31,31,30,31,30,31]
    daylist2=[31,29,31,30,31,30,31,31,30,31,30,31]
    if (yyyy%4 != 0):
        days=daylist2
    elif (yyyy%400 == 0):
        days=daylist2
    elif (yyyy%100 == 0):
        days=daylist
    else:
        days=daylist2

    daysum=0
    for y in range(mm-1):
        daysum=daysum+days[y]
    daysum=daysum+dd-1+UT/24

    if days[1]==29:
        fracyear=yyyy+daysum/366
    else:
        fracyear=yyyy+daysum/365
    if verbose:
        print(yyyy,mm,dd,hh,min,sec)
        print("UT={}".format(UT))
        print("Fractional day: %f" % fracday)
        print("\n"+months[mm-1]+" %i, %i, %i:%i:%i UT = JD %f" % (dd, yyyy, hh, min, sec, JD),)
        print(" = {}\n".format(fracyear))



    return JD


def jd2gd(jd, verbose=False):

    """Task to convert a list of julian dates to gregorian dates
    description at http://mathforum.org/library/drmath/view/51907.html
    Original algorithm in Jean Meeus, "Astronomical Formulae for
    Calculators"

    2009-02-15 13:36 IJC: Converted to importable, callable function
    """
   
    jd=jd+0.5
    Z=int(jd)
    F=jd-Z
    alpha=int((Z-1867216.25)/36524.25)
    A=Z + 1 + alpha - int(alpha/4)

    B = A + 1524
    C = int( (B-122.1)/365.25)
    D = int( 365.25*C )
    E = int( (B-D)/30.6001 )

    dd = B - D - int(30.6001*E) + F

    if E<13.5: mm=E-1

    if E>13.5: mm=E-13

    if mm>2.5: yyyy=C-4716

    if mm<2.5: yyyy=C-4715

    months=["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"]
    daylist=[31,28,31,30,31,30,31,31,30,31,30,31]
    daylist2=[31,29,31,30,31,30,31,31,30,31,30,31]

    h=int((dd-int(dd))*24)
    min=int((((dd-int(dd))*24)-h)*60)
    sec=86400*(dd-int(dd))-h*3600-min*60

    # Now calculate the fractional year. Do we have a leap year?
    if (yyyy%4 != 0): days=daylist2
    elif (yyyy%400 == 0): days=daylist2
    elif (yyyy%100 == 0): days=daylist
    else: days=daylist2

    hh = 24.0*(dd % 1.0)
    min = 60.0*(hh % 1.0)
    sec = 60.0*(min % 1.0)

    dd =  dd-(dd%1.0)
    hh =  hh-(hh%1.0)
    min =  min-(min%1.0)

    if verbose:
        print(str(jd)+" = "+str(months[mm-1])+ ',' + str(dd) +',' +str(yyyy))
        print(string.zfill(h,2)+":"+string.zfill(min,2)+":"+string.zfill(sec,2)+" UTC")

    return (yyyy, mm, dd, hh, min, sec)


def jd2datestr(jd, longform=False):
    """ Convert a JD to a GPI-style datestring YYMMDD
    Note that these are ALWAYS in UTC.

    Parameters
    -----------
    longform : bool
        return 'YYYY-MM-DD' instead of just 'YYMMDD'
    """
    gd = jd2gd(jd)
    if longform:
        datestr = "%04d-%02d-%02d" % ( gd[0], gd[1], gd[2])
    else:
        datestr = "%02d%02d%02d" % ( np.mod(gd[0], 100), gd[1], gd[2])
    return datestr


def expand_path(string):
    import os.path
    return os.path.abspath(os.path.expanduser(os.path.expandvars(string)))



def is_integer(value):
    """ Is something an integer?"""
    try:
        x = int(value)
        return True
    except:
        return False


def imshow_with_mouseover(image, ax=None,  *args, **kwargs):
    """ Wrapper for pyplot.imshow that sets up a custom mouseover display formatter
    so that mouse motions over the image are labeled in the status bar area
    with pixel numerical value as well as X and Y coords.

    Why this behavior isn't the matplotlib default, I have no idea...
    """
    if ax is None: ax = pl.gca()
    myax = ax.imshow(image, *args, **kwargs)
    aximage = ax.images[0].properties()['array']
    # need to account for half pixel offset of array coordinates for mouseover relative to pixel center,
    # so that the whole pixel from e.g. ( 1.5, 1.5) to (2.5, 2.5) is labeled with the coordinates of pixel (2,2)


    # We use the extent and implementation to map back from the data coord to pixel coord
    # There is probably an easier way to do this...
    imext = ax.images[0].get_extent()  # returns [-X, X, -Y, Y]
    imsize = ax.images[0].get_size()   # returns [sY, sX]g
    # map data coords back to pixel coords:
    #pixx = (x - imext[0])/(imext[1]-imext[0])*imsize[1]
    #pixy = (y - imext[2])/(imext[3]-imext[2])*imsize[0]
    # and be sure to clip appropriatedly to avoid array bounds errors
    #report_pixel = lambda x, y :
    def report_pixel(x,y):
        pixvalue = aximage[np.floor( (y - imext[2])/(imext[3]-imext[2])*imsize[0]  ).clip(0,imsize[0]-1),\
                        np.floor( (x - imext[0])/(imext[1]-imext[0])*imsize[1]  ).clip(0,imsize[1]-1)]
        pixstr = "%g"% pixvalue if np.isfinite(pixvalue) else str(pixvalue)
        return "(%6.3f, %6.3f)     %s    " %  (x,y,   pixstr)

        #(x,y, aximage[np.floor(y+0.5),np.floor(x+0.5)])   # this works for regular pixels w/out an explicit extent= call
    ax.format_coord = report_pixel

    return ax


def idl_path_find_file(functionname):
    """ Search the $IDL_PATH to find the source code for a given
    function name. Like 'which' in IDL..."""
    import os
    idlpaths = os.getenv('IDL_PATH').split(':')
    looking_for = functionname.lower() + ".pro"

    for top_path in idlpaths:
        if top_path.startswith('+'):
            top_path = top_path[1:]
            top_path = expand_path(top_path)
            for dirpath, subdirnames, filenames in os.walk(top_path):
                if looking_for in filenames:
                    return os.path.join(dirpath, looking_for)
        else:
            top_path = expand_path(top_path)
            filename = glob.glob(os.path.join(top_path,'*.pro'))
            if looking_for in filenames:
                return os.path.join(dirpath, looking_for)
       
    return None


__all__ = ['gd2jd', 'jd2gd', 'jd2datestr', 'expand_path', 'is_integer', 'imshow_with_mouseover', 'idl_path_find_file']
