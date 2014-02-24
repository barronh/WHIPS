'''
Module for utility functions useful in a variety
of GIS programs

@author: Jacob Oberman
'''

import time
import calendar
from itertools import izip

import numpy
import netCDF4

import filetypes

def wrap_lon_0_360(lon):
    '''Wrap a single longitude to the interval [0, 360)'''
    while lon < 0:
        lon += 360
    while lon >= 360:
        lon -= 360
    return lon

def wrap_lon_neg180_180(lon):
    '''Wrap a single longitude to the interval (-180, 180]'''
    while lon <= -180:
        lon += 360
    while lon > 180:
        lon -= 360
    return lon

def timestr_to_nsecs(timestr, 
                     epoch='00:00:00 01-01-1970', 
                     format='%H:%M:%S %m-%d-%Y'):
    '''
    Convert a given string to the number of seconds
    since a specified epoch.
    
    Defaults to converting into Unix time
    
    Inputs:
        timestr - a string representing the time to convert
        epoch - a string representing desired epoch time
        format - a format string as in time.strptime used for timestr and epoch
        
    Outputs: 
        double - number of seconds since epoch
    '''
    desStructTime = time.strptime(timestr, format)
    epochStructTime = time.strptime(epoch, format)
    desUnixTime = calendar.timegm(desStructTime)
    epochUnixTime = calendar.timegm(epochStructTime)
    return desUnixTime - epochUnixTime

def nsecs_to_timestr(nSecsSinceEpoch,
                     epoch='00:00:00 01-01-1970',
                     format='%H:%M:%S %m-%d-%Y'):
    '''
    Convert the number of seconds since a specified epoch
    to a human-readable time string.
    
    Defaults to converting out of Unix time
    
    Inputs:
        nSecsSinceEpoch - number of seconds since specified epoch
        epoch - a string in the style of format representing epoch time
        format - a format string to convert into as in time.strptime
        
    Outputs:
        string - Human readable time string
    '''
    tStructAtEpoch = time.strptime(epoch, format)
    nSecsAtEpoch = calendar.timegm(tStructAtEpoch)
    nSecsUnix = nSecsSinceEpoch + nSecsAtEpoch
    tStructUnix = time.gmtime(nSecsUnix)
    return time.strftime(format, tStructUnix)
    
def UTCoffset_from_lon(lon):
    '''
    Calculate the approximate offset from UTC based on longitude.
    Offset is returned in seconds and will be positive
    or negative as appropriate
    
    Does not account for daylight savings time.  Uses an 
    algorithm based on simply the longitude, and therefore
    does not account for actual political timezones.
    
    Longitude assumed to be in degrees.
    '''
    hours2secs = 60*60
    return hours2secs*round(wrap_lon_neg180_180(lon)/15.0)

def find_occurences(superArray, subArray):
    '''
    Find the occurrences of a particular subarray within a superarray, 
    searching along the rightmost axis of the array.  Returns a logical array
    with of rank 1 less than the superarray.  Output is true where the subarray
    matched, false otherwise.
    '''
    return numpy.apply_along_axis(numpy.array_equal, -1, superArray, subArray)

def write_grid_to_netcdf(griddef, outFname):
    '''
    Function to create netCDF files that contain
    the lat/lon data needed to plot for a given grid 
    definition that uses rows/columns
    
    Inputs:
        griddef - an instantiated griddef object
        outFname - a path to the outfile (will be clobbered)
    '''
    (minRow, maxRow, minCol, maxCol) = griddef.indLims()
    nRows = maxRow-minRow+1
    nCols = maxCol-minCol+1
    # create the index vectors we need to make grids
    (cols, rows) = numpy.meshgrid(numpy.arange(minCol, maxCol+1), 
                                  numpy.arange(minRow, maxRow+1))
    cols = cols.astype(numpy.float32) # cast as precaution
    rows = rows.astype(numpy.float32)
    # write out netcdf
    fid = netCDF4.Dataset(outFname, 'w', format='NETCDF3_CLASSIC')
    # create dimensions and variables
    fid.createDimension('row', nRows)
    fid.createDimension('col', nCols)
    dims = ('row', 'col')
    # create the 5 grid definitions
    offsets = [(0,0), (1,0), (1,1), (0,1), (.5,.5)] # row,col
    labels = ['ll', 'ul', 'ur', 'lr', 'cent']
    for (lbl, (rowOff, colOff)) in izip(labels,offsets):
        lon = fid.createVariable(lbl+'_lon', 'f', dims)
        setattr(lon, 'Units', 'degrees_east')
        lat = fid.createVariable(lbl+'_lat', 'f', dims)
        setattr(lat, 'Units', 'degrees_north')
        (lat[:], lon[:]) = griddef.griddedToGeo(rows+rowOff, cols+colOff)
    # write grid parameters to file as global attributes
    setattr(fid, 'Projection', griddef.__class__.__name__[:-8])
    for (k,v) in griddef.parms.iteritems():
        setattr(fid, k, v)
    fid.close()

def parse_fromFile_input_file(inFileName, dryRun):
    '''
    Open and read the file and parse it 
    for various parameter values.
    Generate and return a command line call
    that can be used to call whips
    with those parameter values.
    '''
    call = []
    attrs = []
    f = open(inFileName, 'r')
    s = f.readline()
    try:
        while(s != ""):
            if(s == "BEGIN\n"):
                while(s != ""):
                    s = f.readline()
                    if(s == "END\n"):
                        return call + ["--projAttrs"] + attrs
                    if(s[0] == '.'):
                        continue
                    s = s.split('"')
                    if(len(s) > 1 and dryRun):
                        print "Warning: Found line containing "\
                              "stray quotation marks... stripping."
                    words = ("".join(s)).split()
                    '''Determine what the line is supposed to do'''
                    if(words == []):
                        continue
                    elif(words[1] != "="):
                        break
                    elif(words[0] == "DIRECTORY"):
                        call += ["--directory", 
                                 "{0}".format(' '.join(words[2:]))]
                    elif(words[0] == "FILELIST"):
                        call += ["--fileList"] + words[2:]
                    elif(words[0] == "FILETYPE"):
                        call += ["--filetype", 
                                 "{0}".format(' '.join(words[2:]))]
                    elif(words[0] == "GRIDPROJ"):
                        call += ["--gridProj", 
                                 "{0}".format(' '.join(words[2:]))]
                    elif(words[0] == "MAPFUNC"):
                        call += ["--mapFunc", 
                                 "{0}".format(' '.join(words[2:]))]
                    elif(words[0] == "OUTFUNC"):
                        call += ["--outFunc",
                                 "{0}".format(' '.join(words[2:]))]
                    elif(words[0] == "OUTDIRECTORY"):
                        call += ["--outDirectory", 
                                 "{0}".format(' '.join(words[2:]))]
                    elif(words[0] == "OUTFILENAME"):
                        call += ["--outFileName", 
                                 "{0}".format(' '.join(words[2:]))]
                    elif(words[0] == "INCLUDEGRID"):
                        call += ["--includeGrid",
                                 "{0}".format(' '.join(words[2:]))]
                    elif(words[0] == "VERBOSE"):
                        call += ["--verbose", 
                                 "{0}".format(' '.join(words[2:]))]
                    elif(words[0] == "INTERACTIVE"):
                        call += ["--interactive", 
                                 "{0}".format(' '.join(words[2:]))]
                    elif(dryRun):
                        attrs += ['"{0}:{1}"'.format(words[0], \
                                                         ' '.join(words[2:]))]
                    else:
                        attrs += ['{0}:{1}'.format(words[0], \
                                                         ' '.join(words[2:]))]
                break
            s = f.readline()
    except:
        pass
    if(s != ""):
        raise SyntaxError ("Invalid input file.  File must be formatted " \
                      "correctly.\nCheck line '{0}' and try again".format(s))
    raise SyntaxError ("Invalid input file.  Check that the file matches the "\
                      "format described in the documentation and try again")


def parse_filetype(namespace):
    '''
    Open and read the filetype file
    to add associated parameters to namespace 
    '''
    filetype = getattr(filetypes, namespace.filetype + "_filetype")
    printfiletype = namespace.filetype
    wng = ""

    for attr in dir(filetype):
        if attr == "parser":
            setattr(namespace, "filetype", getattr(filetype, attr))
        elif attr == "doutf":
            setattr(namespace, "outFunc", getattr(filetype, attr))
        elif attr[0] == "_":
            pass
        else:
            try:
                wng += "Warning: Value {0} supplied for attribute {1} " \
                      "has been ignored. Filetype {2} does not support "\
                      "custom values for this parameter.\n".format(\
                      getattr(namespace, attr), attr, printfiletype)
            except AttributeError:
                pass
            setattr(namespace, attr, getattr(filetype, attr))
            wng += "Value for {0} given by filetype information\n".format(attr)

    try:
        if (namespace.verbose != 'False'):
            print wng
    except:
        pass
    return namespace
