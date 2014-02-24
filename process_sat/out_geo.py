'''
Framework for creating output from satellite files

All function classes are named after the convention 
func_out_geo.  The function ValidOutfuncs 
automatically maintains a list of implemented
output functions

When INITIATED, output functions may accept a dictionary of
fieldnames and a dictionary of parameters.  If these parameters
are not passed, the constructor is expected to ask the user for
input, but no actual IO should take place in this class (all
io should be piped through the designated IO class)

All output function classes must accept a list of the outputs
from map_geo functions, a GridDef instance, and a base 
output filename when CALLED, in that order.  Note that all 
output functions are not required to use all these
inputs, but they *do* have to accept them without throwing
errors.  Also note that this means all valid outfunc classes
will have a __call__ method with the following syntax.  Many call functions
also implement an additional parameter "verbose" that, when set, prevents the 
function from dumping to the command line.

Output functions have a significant amount of freedom in
what they can do.  They may write files and/or return
data.  At this point, the specifications for returned
data have not yet been decided upon.  returned files
must be saved into the filenames based on the base
filename passed (IE if output0 is passed, a function
may choose to call outputs output0-0, output0-1, etc...
though something more descriptive is preferable)
'''
import sys
from itertools import izip
import datetime
import warnings
import pdb


import numpy
import netCDF4

import utils

def vsnmsg(version): 
    return "This file was generated using WHIPS v{0}".format(version)

def ValidOutfuncs():
    '''Return a list of valid output function names'''
    currentModule = sys.modules[__name__]
    names = dir(currentModule)
    return [el[:-9] for el in names if el.endswith("_out_func")]

class out_func:
    '''Abstract class to for <>_out_geo classes'''
    def __init__(self, parmDict=None):
        self.parmDict = parmDict
    def __call__(self, map_geo, griddef, outfilenames, verbose, version):
        raise NotImplementedError
    @staticmethod
    def parm_list():
        raise NotImplementedError
    @staticmethod
    def required_parms():
        raise NotImplementedError

def _OMNO2e_formula(cloudFrac, fieldOfView):
    eps = 1.5*pow(10,15)*(1+3*cloudFrac)
    capE = pow(10,-16)*(18.5+2.8*pow(10,-4)*pow(abs(fieldOfView-29.7), 3.5))
    return pow((eps*capE), -2)

class invalidPixCeption(Exception):
    pass
    
tai93conv = lambda(timestring):utils.timestr_to_nsecs(timestring, 
                               '00:00:00_01-01-1993', '%H:%M:%S_%m-%d-%Y')

def boolCaster(boolStr):
    if boolStr == 'True':
        return True
    elif boolStr == 'False':
        return False
    else:
        msg = 'Attempt to cast invalid string %s to boolean' % boolStr
        raise TypeError(msg)

# currently borked.  No immediate plans to fix
#class OMNO2e_wght_avg_out_func(out_func):
class OMNO2e_wght_avg_BORKED(out_func): 
    '''
    Weighted avg based on OMNO2e algorithm

    Note: this algorithm doesn't note what to do when the weight for a term is
    nonzero, but the term itself contains a fillValue.  This assumption
    is checked by an assert statement in the code, so it won't be checked
    if optimization is requested
    
    OMNO2e algorithm and theoretical basis laid out
    at <http://disc.sci.gsfc.nasa.gov/Aura/data-holdings/OMI/omno2e_v003.shtml>
    
    Set up to work specifically for OMI instruments.
    
    parameters dict must contain keys:
        toAvg
        customCriteria
        cloudFrac
        pixIndXtrackAxis
        fillVal
    
    writes out a single file of name outfilename
    when called.  That file is an ASCII representation
    weighted averages for each gridcell.  It is a csv 
    file with all numbers in e-format scientific 
    notation.  Cells without valid measurements contain the fillvalue.
    '''  
    @staticmethod
    def parm_list():
        return ['toAvg', 'customCritiera', 'cloudFrac', 
                'pixIndXtrackAxis', 'fillVal']

    @staticmethod
    def required_parms():
        return {'toAvg' : ('The name of the field to be averaged',None),
                'customCriteria' : ('A string that can be evaluated' \
                                    'according to the syntax in the OMNO2d' \
                                    'Description field as described in the' \
                                    'file specification document, except that'\
                                    'data is pre-scaled. For example:\n'\
                                    '\t"SolarZenithAngle=[0:85],CloudFraction=[0:0.3],VcdQualityFlags=~19,XTrackQualityFlags=0,RootMeanSquareErrorOfFit=[0:0.0003],TerrainReflectivity=[0:0.3]"\n\n' \
                                    'would be interpreted as zenith between 0 and 85;'\
                                    'cloud fraction and terrain reflectivity must '\
                                    'be between  0 and 0.3; VcdQualityFlags'\
                                    'exclude 1, 2, 16, 17, 18, and 19; XTrackQualityFlags'\
                                    'only include zeros; and RMS must be between 0 and 0.0003', None),
                'cloudFrac' : ('The name of the field containing the ' \
                               'cloud fractions.\n{ OMI KNMI - CloudFraction' \
                               '\n  OMI NASA - CloudFraction }',None),
                'pixIndXtrackAxis' : ('The dimension order (0 based) of the ' \
                                      'of the "cross-track" dimension (which' \
                                      'ever dimension has size 60).  For all' \
                                      ' currently known cases should be 1 ' \
                                      ' (may change in future versions of ' \
                                      'OMI products).','int'),
                'fillVal' : ('The value to use as a fill value in the output ' \
                             'netCDF file.  This value will replace any missing '\
                             ' or invalid output values','decimal')}
    
    # userKeys not necessary, so 'filler' field used instead
    __userKeys__ = "filetype"

    def __call__(self, maps, griddef, outfilename, verbose, version):

        # function is broken until it can be refactored such that
        # _OMNO2e_func doesn't require totFlag.  Needs to have pixel
        # loop to check each 
        BORKED = 2
        print('THIS FUNCTION IS BORKED.  ABORT! ABORT! ABORT!')
        sys.exit(BORKED)

        # even though IO interface handles casting already
        # a catchblock has been added here for safety
        # in case someone wants to use this class directly
        castDict = {'toAvg':str, 'customCriteria':str,
                    'cloudFrac':str,
                    'pixIndXtrackAxis':int, 'fillVal':float}
        for (k,func) in castDict.items():
            self.parmDict[k] = func(self.parmDict[k])

        '''Write out single weighted-avg file'''
        numpy.seterr(over='raise')
        nRows = griddef.indLims()[1] - griddef.indLims()[0] + 1
        nCols = griddef.indLims()[3] - griddef.indLims()[2] + 1
        sum_weights = numpy.zeros((nRows, nCols))
        sum_weighted_vals = numpy.zeros((nRows, nCols))
        if not isinstance(maps, list):
            maps = [maps] # create list if we didn't get one
        for map in maps:
            with map.pop('parser') as p: # pop out so we can loop
                if verbose:
                    print('Processing {0} for output at {1}.'.format(\
                            p.name, str(datetime.datetime.now())))
                for (k,v) in map.iteritems():
                    # Replace with custom criteria
                    #sumFlag = numpy.array([p.get_cm(self.parmDict['overallQualFlag'], pxind)
                    #                        for (pxind, unused_weight) in v])
                    sumFlag = numpy.mod(sumFlag, 2)
                    cFrac = numpy.array([p.get_cm(self.parmDict['cloudFrac'], pxind)
                                          for (pxind, unused_weight) in v])
                    totFlag = numpy.logical_or(numpy.logical_or(sumFlag, cFracFlag), solZenFlag)
                    fov = numpy.array([pxind[self.parmDict['pixIndXtrackAxis']]
                                        for (pxind, unused_weight) in v])
                    toAvg = numpy.array([p.get_cm(self.parmDict['toAvg'], pxind)
                                          for (pxind, unused_weight) in v])
                    # BORKED
                    weights = 0
                    #   weights = _OMNO2e_formula(cFrac, fov)
                    assert ~any(numpy.logical_and(~numpy.isnan(weights), numpy.isnan(toAvg)))
                    sumWeight = numpy.nansum(weights)
                    sumWeightVals = numpy.nansum(toAvg*weights)
                    # only add if we had some element (otherwise we fill
                    # sum_weights with nans)
                    if ~numpy.isnan(sumWeight) and ~numpy.isnan(sumWeightVals): 
                        sum_weights[k] += sumWeight
                        sum_weighted_vals[k] += sumWeightVals
                map['parser'] = p  # return parser to map
        oldsettings = numpy.seterr(divide='ignore')  # ignore any div by zero errors
        avgs = numpy.where(sum_weights != 0, 
                           numpy.divide(sum_weighted_vals, sum_weights), 
                           self.parmDict['fillVal'])
        numpy.seterr(divide=oldsettings['divide'])  # reset to default
        numpy.savetxt(outfilename, avgs, delimiter=',', fmt='%7e')
        return avgs

class OMNO2e_netCDF_avg_out_func(out_func):
    '''
    Weighted average for a given set of filtered values
    based on OMNO2e algorithm.
    
    Assumptions:
        - this algorithm assumes that fields to be averaged will have,
        at most, 1 extra dimension.  If not, an assertion error is raised.
        - this algorithm is undefined for cases where the weight of a term
        is nonzero, but the term contains a fillValue.  If this condition
        is met, unexpected results may occur.  This assumption is NOT checked
        - The timestamp of the file is assumed to be in the TAI93 standard.
    
    OMNO2e algorithm and theoretical basis laid out
    at <http://disc.sci.gsfc.nasa.gov/Aura/data-holdings/OMI/omno2e_v003.shtml>
    
    Set up to work specifically for OMI instruments.

    parameters dict must contain keys:
        customCriteria:
            Criteria for pixel 
            selection. Each ciriteria is separated
            by a comma, each critieria follows the 
            OMI NO2 field "Description" attribute
            as documented in section 3.4.1 of the 
            OMNO2d File Specification document
            version 1.1 (Jan 10, 2013). For example,            
            { OMI NASA OMNO2d - SolarZenithAngle=[0:85],VcdQualityFlags=~1,CloudFraction=[0:0.3]
              OMI NASA OMNO2d - SolarZenithAngle=[0:85],VcdQualityFlags=~19,XTrackQualityFlags=0,RootMeanSquareErrorOfFit=[0:0.0003],TerrainReflectivity=[0:0.3] }
        cloudFrac:
            field with cloud fraction (0 to 1).  When this
            field is GREATER than the cutoff value the
            pixel is ignored.
        time:
            Field with the timestamps for each pixel. 
            Assumed to be in TAI-93 format.  When
            this field is less than the timeStart 
            parameter or greater than the timeStop
            parameter the pixel is ignored.
        longitude:
            Field with the longitudes at cell centers.
            Used to estimate timezones of the pixels if
            'local' is selected for timeComparison.  Not
            used when timeComparison is 'UTC'
        inFieldNames:
            List of fields to process.  Each of these
            is output as a separate variable in the 
            netcdf output file.
        outFieldNames:
            List of desired output names.  Must be of the
            same length and co-indexed to the list above.
            These will be the actual variable names in
            the netcdf file.
        outUnits:
            List of string labels for the units of each
            output quantity.  Must be of the same length
            and co-indexed to the lists above.
        extraDimLabel: 
            List of the labels for the above extra 
            dimensions.  1 per variable.  Only used if the
            coindexed field has an extra dim. Must be of the
            same length and co-indexed to the lists above.
        extraDimSize:
            List of the sizes of the above extra dimensions.
            1 per variable.  If the coindexed field does not
            have an extra dim, put in 0 or 1.  Must be
            of the same length and co-indexed to the lists 
            above.
        timeComparison:
            Determines how the timeStart and timeStop 
            arguments are interpreted.  If the user selects
            'local', these arguments are interpreted as local
            times.  Only those pixels whose timestamps 
            indicate they lie in the desired span in local
            time will be included.  Daylight savings time
            is not considered and time zone calculations are
            only approximate.  If the users selects 'UTC'
            a straight comparison is done between the pixel
            timestamp and the timeStart and timeStop 
            arguments to see if the pixel is valid.
        timeStart:
            Initial time we want included in file.  
            Times must be in TAI93 standard format.
            *format hh:mm:ss_MM-DD-YYYY will also be converted automatically.
        timeStop:
            Final Time we want in included in file.
            Times must be in TAI93 standard format.
            *format hh:mm:ss_MM-DD-YYYY will also be converted automatically.
        pixIndXtrackAxis:
            The axis (IE which dimension in memory order)
            that specifies the pixels cross-track position.
            This way, the cross-track index number can
            be retrieved safely and used in the weighting
            function.
        fillVal:
            The value we want to use to denote missing data
            in the output file.  This will be documented
            within the output file itself.
        includePixelCount:
            If this parameter is True, WHIPS will include a field
            'ValidPixelCount' in the output file that will include
            the number of valid pixels for each grid cell.

    Outputs a netcdf file with name determined by outFileName
    parameter.  This netcdf file contains as many variables
    as there are inFieldNames passed.  Each variable
    is output as an average over the range of values
    where it was valid acccording to the averaging
    scheme dedfined in the NASA document linked above.
    '''
    @staticmethod
    def parm_list():
        return ['customCriteria', 'cloudFrac',
                'time', 'longitude', 'inFieldNames', 'outFieldNames',
                'outUnits', 'extraDimLabel', 'extraDimSize',
                'timeComparison', 'timeStart', 'timeStop',
                'pixIndXtrackAxis', 'fillVal', 'includePixelCount']
    @staticmethod
    def required_parms():
        return {'customCriteria' : ('A string that can be evaluated' \
                                    'according to the syntax in the OMNO2d' \
                                    'Description field as described in the' \
                                    'file specification document, except that'\
                                    'data is pre-scaled. For example:\n'\
                                    '\t"SolarZenithAngle=[0:85],CloudFraction=[0:0.3],VcdQualityFlags=~19,XTrackQualityFlags=0,RootMeanSquareErrorOfFit=[0:0.0003],TerrainReflectivity=[0:0.3]"\n\n' \
                                    'would be interpreted as zenith between 0 and 85;'\
                                    'cloud fraction and terrain reflectivity must '\
                                    'be between  0 and 0.3; VcdQualityFlags'\
                                    'exclude 1, 2, 16, 17, 18, and 19; XTrackQualityFlags'\
                                    'only include zeros; and RMS must be between 0 and 0.0003', None),
                'cloudFrac' : ('The name of the field containing the ' \
                               'cloud fractions\n{ OMI KNMI - CloudFraction' \
                               '\n  OMI NASA - CloudFraction }',None),
                'time' : ('The name of the field containing the timestamps. '\
                          ' Timestamps are assumed to be the in TAI-93 ' \
                          'format.\n{ OMI KNMI - Time\n  OMI NASA - TIME }', \
                          None),
                'longitude' : ('The name of the field containing longitudes ' \
                               'at cell centers.  Longitudes should be in ' \
                               'degrees east.\n{ OMI KNMI - Longitude\n  ' \
                               'OMI NASA - Longitude }',None),
                'inFieldNames' : ('The names of the fields desired to be ' \
                                  'output.  Input as comma-delimited list ', \
                                  'list'),
                'outFieldNames' : ('The names of the output variables (even ' \
                                   'if they are to be the same as the input ' \
                                   'variables).  Should be a comma-delimited ' \
                                   'list co-indexed to inFieldNames','list'),
                'outUnits' : ('The units of the variables to be written out. ' \
                              'Should be a comma-delimited list co-indexed ' \
                              'to inFieldNames','list'),
                'extraDimLabel' : ('The label for the extra dimension '  \
                                   '(should the variable have an extra ' \
                                   'dimension).  Ignored in the case of a ' \
                                   '2D variable.  Should be a comma-delimited '\
                                   'list co-indexed to inFieldNames','list'),
                'extraDimSize' : ('The size of the extra dimensions (should ' \
                                  'the variable have an extra dimension).  ' \
                                  'For 2D variables, must be set to 0. (zero)' \
                                  '  Should be a comma-delimited list ' \
                                  'co-indexed to inFieldNames.','list'),
                'timeComparison' : ('Must be set to either "local" or "UTC". ' \
                                    ' Determines how the file timestamps are ' \
                                    'compared to the start/stop time.  If set '\
                                    'to "local", then the file timestamps are '\
                                    'converted to local time on a pixel-by-'\
                                    'pixel basis (using longitude to estimate '\
                                    'time zone) before being compared to time '\
                                    'boundaries.  If set to "UTC" the file ' \
                                    'timestamps (which are assumed to be in ' \
                                    'UTC) are compared against the start/stop '\
                                    'time directly.',None),
                'timeStart' : ('The earliest time for which data should be ' \
                               'recorded into the output file.  All times in ' \
                               'input files before this time will be filtered '\
                               'out.  Must be in the format hh:mm:ss_MM-DD-' \
                               'YYYY','time'),
                'timeStop' : ('The latest time for which data should be ' \
                               'recorded into the output files.  All times in '\
                               'input files after this time will be filtered ' \
                               'out.  Must be in the format hh:mm:ss_MM-DD-' \
                               'YYYY','time'),
                'pixIndXtrackAxis' : ('The dimension order (0 based) of the ' \
                                      '"cross-track" dimension (whichever ' \
                                      'dimension has size 60).  For all ' \
                                      'currently known cases set equal to 1 ' \
                                      '(depends on the construction of the ' \
                                      'parser function.  If you rewrite the ' \
                                      'parser, check this).','int'),
                'fillVal' : ('The value to use as a fill value in the output ' \
                             'netCDF file.  This value will replace any missing '\
                             'or invalid output values','decimal'), 
                'includePixelCount' : ('If set to true, the output will include '\
                                       'a field "ValidPixelCount" that contains '\
                                       'the number of valid pixels in each grid '\
                                       'cell.  Only pixels with nonzero weight  '\
                                       'are considered valid.', 'bool')}
    # variable signifying which list is to act as the master list index
    __userKeys__ = "inFieldNames"

    def __call__(self, maps, griddef, outfilename, verbose, version):
        '''Write out a weighted-average file in netcdf format.'''
        #Make sure non-string parameters are in the correct format
        dimsizes = self.parmDict['extraDimSize']
        for i in range(len(dimsizes)):
            try:
                dimsizes[i] = int(dimsizes[i])
            except ValueError:
                print ("Warning: {0} is not a valid extraDimSize value.  " \
                      "Using 0 instead").format(dimsizes[i])
                dimsizes[i] = 0
                continue
        self.parmDict['extraDimSize'] = dimsizes

        # even though IO interface handles casting already,
        # a catchblock has been added here for safety
        # in case someone wants to use this class directly
        castDict = {'customCriteria':str, 'cloudFrac':str,
                    'time':str,
                    'longitude':str, 'inFieldNames':list,
                    'outFieldNames':list, 'outUnits':list,
                    'extraDimLabel':list, 'extraDimSize':list,
                    'timeComparison':str, 'timeStart':tai93conv,
                    'timeStop':tai93conv, 'pixIndXtrackAxis':int,
                    'fillVal':float, 'includePixelCount':boolCaster}
        for (k,func) in castDict.items():
            try:
                self.parmDict[k] = func(self.parmDict[k])
            except TypeError:
                pass

        #Perform some basic sanity checks with parameters
        if self.parmDict['timeStart'] > self.parmDict['timeStop']:
            msg = 'Input start time must come before stop time.'
            raise IOError(msg)
        if (len(self.parmDict['inFieldNames']) !=  \
            len(self.parmDict['outFieldNames']) or  
            len(self.parmDict['inFieldNames']) !=  \
            len(self.parmDict['outUnits']) or      
            len(self.parmDict['inFieldNames']) !=  \
            len(self.parmDict['extraDimLabel'])):
            msg = 'All field/unit inputs ' + \
                'should have the same number of elements.'
            raise IOError(msg)

        # create numpy arrays to hold our data
        (minRow, maxRow, minCol, maxCol) = griddef.indLims()
        nRows = maxRow - minRow + 1
        nCols = maxCol - minCol + 1
        nValidPixels = numpy.zeros((nRows, nCols))
        sumWght = numpy.zeros((nRows, nCols, 1))  # needs extra dim to generalize for 3D vars
        sumVars = dict()
        for field, size in zip(self.parmDict['inFieldNames'], self.parmDict['extraDimSize']):
            if size:
                sumVars[field] = numpy.zeros((nRows, nCols, size))
            else:
                # pad with a singlet dim if it was 2D
                sumVars[field] = numpy.zeros((nRows, nCols, 1))
        
        # loop over maps
        if not isinstance(maps, list):
            maps = [maps] # create list if we only got a single map
        
        for map in maps:
            cldskips = 0
            sumskips = 0
            zenskips = 0
            cstskips = 0
            checked = 0
            # open up context manager
            with map.pop('parser') as parser: # remove parser for looping
                if verbose:
                    print('Processing {0} for output at {1}.'.format(\
                            parser.name, str(datetime.datetime.now())))

                # loop over gridboxes in map and calculate weights
                for (gridCell, pixTup) in map.iteritems():
                    # translate gridCell to account 
                    # for possible non-zero ll corner
                    gridRow = gridCell[0]
                    gridCol = gridCell[1]
                    gridInd = (gridRow - minRow, gridCol - minCol)
                    # get the values needed to calculate weight
                    for (pxInd, unused_weight) in pixTup:
                        checked += 1
                        custskip = False
                        for customCheck in self.parmDict['customCriteria'].split(','):
                            customVar, crit = [x.strip() for x in  customCheck.split('=')]
                            custField = parser.get_cm(customVar, pxInd)
                            if ':' in crit:
                                low, hi = [eval(x.strip()) for x in crit.replace('[', '').replace(']', '').split(':')]
                                if not (custField >= low and custField <= hi):
                                    custskip = True
                            elif '~' in crit:
                                if numpy.bitwise_and(custField, eval(crit.replace('~', '').strip())) > 0:
                                    custskip = True
                            else:
                                if not (custField == eval(crit)):
                                    custskip = True
                        if custskip:
                            cstskips += 1
                            continue                                                
                        cFrac = parser.get_cm(\
                                self.parmDict['cloudFrac'], pxInd)
                        time = parser.get_cm(self.parmDict['time'], pxInd)
                        # calculate and factor in offset 
                        # if the user wanted us to
                        if self.parmDict['timeComparison'] == 'local':
                            pixLon = parser.get_cm(\
                                         self.parmDict['longitude'], pxInd)
                            offset = utils.UTCoffset_from_lon(pixLon)
                            time += offset
                        if time < self.parmDict['timeStart'] \
                               or time > self.parmDict['timeStop']:
                            continue
                        # read in all the data 
                        # abandon ship if data is all NaN
                        rawDataDict = {}
                        try:
                            for field in self.parmDict['inFieldNames']:
                                rawData = parser.get_cm(field, pxInd)
                                if numpy.isnan(rawData).all():
                                    raise invalidPixCeption
                                rawDataDict[field] = rawData
                        except invalidPixCeption:
                            continue
                        # compute the weight
                        fov = pxInd[self.parmDict['pixIndXtrackAxis']]
                        weight = _OMNO2e_formula(cFrac, fov)
                        assert weight != numpy.NaN
                        if weight > 0:
                            nValidPixels[gridInd] += 1

                        # add the weight to the total for this cell 
                        sumWght[gridInd] += weight
                        for field in self.parmDict['inFieldNames']:
                            weightVals = rawDataDict[field] * weight
                            if weightVals.size > 1:
                                sumVars[field][gridInd] = \
                                     numpy.nansum([sumVars[field][gridInd],\
                                                   weightVals], axis=0)        
                            else:
                                sumVars[field][gridInd] = \
                                     numpy.nansum([sumVars[field][gridInd][0],\
                                                   weightVals])
                print 'Checked:', checked
                print 'Excluded Sum,Cld,Zen,Cst', sumskips, cldskips, zenskips, cstskips
                print 'Accepted:', checked - sumskips - cldskips - zenskips - cstskips
                map['parser'] = parser  # return parser to map
                
        # divide out variables by weights to get avgs. 
        oldSettings = numpy.seterr(divide='ignore')
        avgs = dict()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for (field,var) in sumVars.iteritems():
                unfiltAvgs = var/sumWght
                filtAvgs = numpy.where(sumWght != 0, unfiltAvgs, \
                           self.parmDict['fillVal'])
                # strip trailing singlet for 2D arrays
                if filtAvgs.shape[-1] == 1:
                    avgs[field] = filtAvgs.reshape(filtAvgs.shape[0:2])
                else:
                    avgs[field] = filtAvgs
        numpy.seterr(divide=oldSettings['divide'])

        
        # associate coindexed parameters into dicts 
        # so we can loop by field
        outFnames = dict(izip(self.parmDict['inFieldNames'], self.parmDict['outFieldNames']))
        units = dict(izip(self.parmDict['inFieldNames'], self.parmDict['outUnits']))
        extraDim = dict(izip(self.parmDict['inFieldNames'], self.parmDict['extraDimLabel']))
                
        # write out results to a netcdf file
        outFid = netCDF4.Dataset(outfilename, 'w', format='NETCDF3_CLASSIC')
        # create the 2 dimensions all files use
        outFid.createDimension('row', nRows)
        outFid.createDimension('col', nCols)
        # write global attributes
        setattr(outFid, 'Version', vsnmsg(version))
        setattr(outFid, 'File_start_time', utils.nsecs_to_timestr(self.parmDict['timeStart'], '00:00:00 01-01-1993'))
        setattr(outFid, 'File_end_time', utils.nsecs_to_timestr(self.parmDict['timeStop'], '00:00:00 01-01-1993'))
        setattr(outFid, 'FilterDefinition', self.parmDict['customCriteria'])
        setattr(outFid, 'Time_comparison_scheme', self.parmDict['timeComparison'])
        fileListStr = ' '.join([map['parser'].name for map in maps])
        setattr(outFid, 'Input_files', fileListStr)
        setattr(outFid, 'Projection', griddef.__class__.__name__[:-8])
        for (k,v) in griddef.parms.iteritems():
            setattr(outFid, k, v)
        # loop over fields and write variables
        for field in self.parmDict['inFieldNames']:
            # create tuple of dimensions, defining new dim
            # if necessary
            if len(avgs[field].shape) == 2:
                # only row/cols
                varDims = ('row', 'col')
            elif len(avgs[field].shape) == 3:
                # has extra dim
                dimName = extraDim[field]
                dimSize = avgs[field].shape[2]
                if dimName not in outFid.dimensions.keys():
                    outFid.createDimension(dimName, dimSize)
                varDims = ('row', 'col', dimName)
            # create and assign value to variable
            varHandle = outFid.createVariable(outFnames[field], 'd', varDims, fill_value=self.parmDict['fillVal'])
            varHandle[:] = avgs[field]
            # assign variable attributes
            setattr(varHandle, 'Units', units[field])
        # Write out the pixel counts if the user requested them
        if self.parmDict['includePixelCount']:
            varDims = ('row', 'col')
            varHandle = outFid.createVariable('ValidPixelCount', 'i', varDims, 
                                              fill_value=self.parmDict['fillVal'])
            varHandle[:] = nValidPixels
        outFid.close()
        # create a dict with the same data as avgs, but diff names
        outAvg = dict()
        for (k,v) in avgs.iteritems():
            outAvg[outFnames[k]] = v
        if self.parmDict['includePixelCount']:
            outAvg['ValidPixelCount'] = nValidPixels
        return outAvg
    
class wght_avg_netCDF(out_func):
    '''
    Generalized weighted average algorithm
    
    Designed to compute the average of an arbitrary number of desired 
    parameters, with the value weights based on an arbitrary number of input
    parameters.  Note that values may be weighted according to their own value.
    
    The output will be in the form of a netcdf file with name determined by the
    outFileName parameter.  This netCDF file will have dimensions determined by
    the grid_geo file, as well as additional dimensions as required by the 
    input fields.  
    
    Owing to the complexity of the inputs required for this function and the 
    security problems posed by allowing users to input functions to be 
    evaluated, this output function does not support the I/O interface at this
    time.  It is designed to subclassed.

    This function (and therefore subclasses of this function) at present can
    only handle a single input map.  It may be extended to properly handle 
    multiple input maps at some point in the future, but this is difficult
    because the filter function is expected to apply to all pixels in the cell
    (which would require looping over all the maps to find all the pixels)
    but also requires a reference to the parser (which would require those 
    parsers be held open)
    
    parmDict must contain the following keys:  
        time:
            The field associated with the timestamps.  Timestamps may be in any
            format so long as a function is provided to convert them to Unix 
            timestamp values (as this is what the function will use internally)
        longitude:
            Field with the longitudes at cell centers.  Used to estimate
            timezones of the pixels if local is selected for timeComparison.  
            Not used when timeComparison is 'UTC'
        inFieldNames:
            List of strings corresponding to fields for which output is 
            desired.  These must be valid keys for the parser.  Each is output
            as a seperate variable in the netcdf output file.
        outFieldNames:
            List of strings corresponding to the names the output variables
            should have in the final netCDF file.  Must be of the same length
            and co-indexed to the list above.  
        outUnits:
            List of strings corresponding to the labels for the units of each 
            output variable.  These will be attached as the "units" attribute
            for each variable in the output netCDF file.  Must be of the same 
            length and co-indexed to the lists above.
        logNormal:
            Vector indicating whether or not we want to take the
            lognormal mean (as opposed to the simple, arithmetic mean).  If 
            this parameter is set to "True", the mean will be taken as follows:
                logData = log(data)
                logAvg = sum(logData*wghts)/sum(wghts)
                avg = 10^logAvg
            whereas if this parameter is set to "False" the mean will be simply:
                avg = sum(data*wghts)/sum(wghts)
            To allow finer-grained control of the output, logNormal must be 
            set individually for each output field (as it may be appropriate
            to use the log normal distribution only for select fields).  Thus,
            logNormal must be of the same length and co-indexed to the lists
            above.
        dimLabels:
            List of tuple-like strings(delimited by periods with no whitespace),
            each of which contains as many strings as there are
            extra dimensions in the corresponding field.  IE, if a field has 
            dimensions (xtrack, time, layer, quant) and we allocate along
            xtrack and time, then the tuple-like string for that field should be 
            "(layer.quant)" so that the dimensions are named properly in the
            otuput netCDF file.  The list (not the individual tuples) must be
            of the same length and co-indexed to the lists above.  Note that 
            the dimensions looped over and allocated to cells in the map_geo
            function DO NOT need to be represented here.
        dimSizes:
            List of tuple-like strings (delimited by periods with no whitespace),
            each of which contains as many strings (castable to ints) as there
            are extra dimensions in the corresponding field.  IE if the field
            described in dimLabels had dimensions (60, 1300, 9, 4) then the 
            tuple-like string for that field should be "(9.4)".  The list (not the 
            individual tuples) must be of the same length and co-indexed to the
            lists above.  If the field has no extra dimensions, then an empty 
            tuple should be used as a placeholder.Note that the dimensions 
            looped over and allocated to cells in the map_geo function DO NOT
            need to be represented here.
        timeComparison:
            Determines how the timeStart and timeStop parameters are 
            interpreted.  Must be either 'local' or 'UTC'.  If 'local' is 
            selected, only those pixels in the desired timespan in in local 
            time will be included.  Daylight savings time is not considered and
            time zone calculations are approximate based on longitude.  If 
            'UTC' is selected, the raw timestamps are compared directly to the 
            timeStart and timeStop parameters, without attempting to account 
            for timezone.
        timeStart:
            Initial time we want included in the output file.  All measurements
            taken before this time will be discarded, even if they are included
            in the files passed to output function.  Must be a string of the 
            the format hh:mm:ss_MM-DD-YYYY.
        timeStop:
            Final time to be included in the output file.  All measurements
            taken after this time will be discarded.  Must be a string of the
            format hh:mm:ss_MM-DD-YYYY.
        timeConv:
            Function that converts timeStart and timeStop (both of which are 
            strings of the format hh:mm:ss_MM-DD-YYYY) into the format used to
            store time within the parser.  IE if the parser returns time in 
            TAI93, then this function should convert a string
            hh:mm:ss_MM-DD-YYYY to TAI93.
        fillVal:
            The value with which we want to denote missing data in the output
            file.  This value must be castable to the type of all output 
            variables.  Each variable will document the fill value in its
            attributes.
        notes:
            String to be included in the output attributes of the final file.  
            Use this to hold any extra if you you'd like to be packaged with
            the data.
        weightFunction:
            Function that computes the weight of a value.  This can be as 
            simple as assigning a weight of 1 to that value (creating an 
            unweighted average) or using multiple fields to generate a weight
            for each pixel.  The function must be of the form:
                weight = weightFunction(parser, index, prevWght)
            where parser is the parser for the file under consideration, ind
            is the tuple of indices for the pixel under consideration, and
            prevWght is the weight as computed by the mapping function.  Note
            that in authoring these functions it is safe to use both the get()
            and get_cm() functions of the parser argument.  The function may 
            return either 0 or NaN for pixels that should not be considered in
            the average.  Either will be handled appropriately.  The docstring 
            for this function will be included as a global attribute in the 
            output file, so the docstring should be sufficient to describe the
            function in it's entirety.
        filterFunction:
            Function that looks at the entire stack of pixels for a cell and 
            selects any pixels that need to be filtered out.  Note that for 
            efficiency reasons this should not be used to catch filter
            conditions unique to each pixel.  Those sort of filters should
            be applied in the weightFunction.  This function is strictly
            intended for operations that can only be performed once the entire
            stack for a particular cell is available (IE if a value is checked
            for each pixel in the cell and those pixels not matching the 
            majority value are rejected).  The function must be of the form
                flagVec = filterFunction(parser, indStack)
            where indStack is an iterable of index tuples and parser is the 
            parser for the file under consideration.  flagVec (the return
            vector) should be boolean and true for those values that should
            NOT be included in the final average.  It should be the same length
            as indStack.  To reiterate: flagVec should true for those values
            that should be REMOVED, and false for those values to be kept.
            The docstring of this function will be included as a global 
            attribute in the final output file, so the docstring should be
            sufficient to describe the function in it's entirety.  Note that it
            is safe to use both get and get_cm functions within this function -
            it is guaranteed to be called within a context manager.
    '''

    def __init__(self, parmDict=None):
        # call ancestor method
        out_func.__init__(self, parmDict)

        # check that all the lists are the same length
        lists = ['outFieldNames', 'outUnits', 'dimLabels', 'dimSizes', 'logNormal']
        canonicalLength = len(self.parmDict['inFieldNames'])
        isCorrectLength = [len(self.parmDict[list]) == canonicalLength for list in lists]
        if not all(isCorrectLength):
            wrongLength = [name for (name, corr) in izip(lists,isCorrectLength) if not corr]
            msg = "All lists must be the same length.  The following list lengths do " \
                "not match the length of inFieldNames: " + ' '.join(wrongLength)
            raise ValueError(msg)

        # process lists of tuple-like objects into lists of appropriately-typed tuples
        labelTups = self.parmDict['dimLabels']
        sizeTups = self.parmDict['dimSizes']
        # confirm that each tuple is the same size
        tupsMatch = [len(l) == len(s) for (l,s) in izip(labelTups, sizeTups)]
        if not all(tupsMatch):
            misMatch = [l + ' does not match ' + s for (l,s,m) in izip(labelTups,sizeTups,tupsMatch) if not m]
            msg = "All tuple-like strings must correspond to tuples of corresponding size. " \
                "The following sets do not correspond: \n" + '\n'.join(misMatch)
            raise ValueError(msg)
        # convert sizes to integers
        try:
            sizeIntTups = [[int(s) for s in strTup] for strTup in sizeTups]
        except ValueError as err:
            messageWords = err.message.split()
            uncastable = messageWords[-1].strip("'")
            msg = "The value %s in the dimSizes argument was not castable to " \
                "an integer." % uncastable
            raise ValueError(msg)
        # put back into parameter dictionary
        self.parmDict['dimLabels'] = labelTups
        self.parmDict['dimSizes'] = sizeIntTups

        # process logNormal
        try:
            self.parmDict['logNormal'] = [boolCaster(el) for el in self.parmDict['logNormal']]
        except TypeError:
            print('Bad string in logNormal.  Must be either "True" or "False". Exiting.')
            raise

        # convert all the parameters co-indexed to inFieldNames to dictionaries
        # keyed off of inFieldNames
        inFnames = self.parmDict['inFieldNames']  
        for key in lists:
            self.parmDict[key] = dict(zip(inFnames, self.parmDict[key]))
        
    def __call__(self, maps, griddef, outfilename, verbose, version):
        '''Write out a weighted-average file in netcdf format.'''
        
        # create a dictionary of numpy arrays that will hold the data for all 
        # our variables, keyed to inFieldNames
        (minRow, maxRow, minCol, maxCol) = griddef.indLims()
        nRows = maxRow - minRow + 1
        nCols = maxCol - minCol + 1
        outputArrays = dict()
        for field in self.parmDict['inFieldNames']:
            dims = [nRows, nCols] + self.parmDict['dimSizes'][field]
            outputArrays[field] = numpy.zeros(dims)
            
        # prep for computing weights
        wghtDict = dict()  # we only want to compute each weight once
        wghtFunc = self.parmDict['weightFunction']
        filtFunc = self.parmDict['filterFunction']
        
        # convert the times to the proper format
        tConvFunc = self.parmDict['timeConv']
        timeStart = tConvFunc(self.parmDict['timeStart'])
        timeStop = tConvFunc(self.parmDict['timeStop'])
    
        # loop over maps
        if not isinstance(maps, list):
            maps = [maps] # create list if we didn't get one

        # check to see that we were given exactly one map.  If we aren't
        # explain that this function only works for one map and exit
        if len(maps) != 1:
            msg = "Though in general output functions are designed to handle" \
                " multiple input files, this function currently can only " \
                "process individual input files.  Please only provide one " \
                "input file or use a different output function.  This " \
                "limitation may be fixed if there is a convincing reason " \
                "to rewrite the function to accomodate more inputs"
            raise NotImplementedError(msg)
            
        for map in maps:
            with map.pop('parser') as p: 
                if verbose:
                    print('Processing %s for output at %s' %
                          (p.name, str(datetime.datetime.now())))
                    
                # loop over the cells in the map, processing each
                for (cellInd, pixTups) in map.iteritems():
                    
                    # compute the weight only if we haven't already.  In either
                    # case, put the weights in array.
                    wghts = [wghtDict.setdefault(ind, wghtFunc(p, ind, wgt))
                            for (ind, wgt) in pixTups]
                    
                    # create the time array we'll be using to filter
                    tArray = numpy.array([p.get_cm(self.parmDict['time'], ind)
                                          for (ind, wgt) in pixTups]).squeeze()
                    if self.parmDict['timeComparison'] == 'local':
                        offsets = numpy.array([utils.UTCoffset_from_lon(
                                               p.get_cm(self.parmDict['longitude'], ind)) 
                                               for (ind, wgt) in pixTups]).squeeze()
                        tArray += offsets
                    tFlag = numpy.logical_or(tArray < timeStart, tArray > timeStop)
                    
                    # use the filter function on the stack to apply user-defined
                    # filter conditions
                    pixIndStack = [pInd for (pInd, unused_weight) in pixTups]
                    uFlag = numpy.array(filtFunc(p, pixIndStack))

                    # combine time filter and user filter into a single, global flag 
                    gFlag = numpy.logical_or(uFlag, tFlag)

                    # filter the weights so that values that will be rejected don't
                    # have their weights included in the denominator of the final
                    # average.
                    wghts = numpy.where(gFlag, numpy.NaN, wghts)
                    
                    # loop over fields.  For each, compute avg and save
                    for field in self.parmDict['inFieldNames']:
                        
                        # create the array of weighted values    
                        vals = numpy.array([p.get_cm(field, ind)
                                            for (ind, wgt) in pixTups]).squeeze()
                        if self.parmDict['logNormal'][field]:
                            vals = numpy.log(vals) # work with logarithm of data

                        # create a slice object that will allow us to broadcast
                        # weights against the values
                        extraDims = self.parmDict['dimSizes'][field]
                        # handle the special case where we put in 0 for extradims
                        if len(extraDims) == 1 and extraDims == [0]:
                            nExtraDims = 0
                        else:
                            nExtraDims = len(extraDims)
                        wghtSlice = [Ellipsis]+[numpy.newaxis]*nExtraDims

                        # handle special case where there were no pixels in
                        # cell
                        if vals.size == 0:
                            vals = vals.reshape([0]+extraDims)

                        # compute weighted values
                        wghtVals = vals*wghts[wghtSlice]
                        
                        # average the weighted Values
                        wghtValSum = numpy.nansum(wghtVals, axis=0)
                        wghtSum = numpy.nansum(wghts, axis=0)
                        # avoid hassle with div/0 warnings
                        if wghtSum != 0:
                            wghtValAvg = wghtValSum/wghtSum
                        else:
                            wghtValAvg = numpy.NaN
                        
                        # re-exponentiate if we took log average
                        if self.parmDict['logNormal'][field]:
                            wghtValAvg = numpy.exp(wghtValAvg)
                        
                        # mask nan's with fillVal, then slot into output array
                        wghtValAvg = numpy.where(numpy.isnan(wghtValAvg),
                                                 self.parmDict['fillVal'], 
                                                 wghtValAvg)
                        outputArrays[field][cellInd] = wghtValAvg
        
                    # done looping over fields
                # done looping over cells
            # done with context manager on parser
                        
            # return the parser to the map so it can be used elsewhere
            map['parser'] = p
            if verbose:
                print('Done processing %s at %s' %
                      (p.name, str(datetime.datetime.now())))
            
        # done looping over maps
                
        # set up the parts of the netcdf file that AREN'T field specific
        outFid = netCDF4.Dataset(outfilename, 'w', format='NETCDF3_CLASSIC')
        outFid.createDimension('row', nRows)
        outFid.createDimension('col', nCols)

        # set global attributes
        setattr(outFid, 'Version', vsnmsg(version))
        setattr(outFid, 'File_start_time', 
                utils.nsecs_to_timestr(self.parmDict['timeStart'], 
                                       epoch='00:00:00 01-01-1993', 
                                       format='%H:%M:%S %m-%d-%Y'))
        setattr(outFid, 'File_stop_time', 
                utils.nsecs_to_timestr(self.parmDict['timeStop'], 
                                       epoch='00:00:00 01-01-1993', 
                                       format='%H:%M:%S %m-%d-%Y'))
        setattr(outFid, 'Time_comparison_scheme', self.parmDict['timeComparison'])
        flistStr = ' '.join([map['parser'].name for map in maps])
        setattr(outFid, 'Input_files', flistStr)
        setattr(outFid, 'Weighting_function_description', wghtFunc.__doc__)
        setattr(outFid, 'Filter_function_description', filtFunc.__doc__)
        # add in attributes for the projection
        setattr(outFid, 'Projection', griddef.__class__.__name__[:-8])
        setattr(outFid, 'Notes', self.parmDict['notes'])
        for k in griddef.parm_list():
            setattr(outFid, k, griddef.parms[k])

        # loop over fields and write all information for each field
        for field in self.parmDict['inFieldNames']:

            # create the dimensions in the file
            extraDimSizes = self.parmDict['dimSizes'][field]
            extraDimLabels = self.parmDict['dimLabels'][field]
            for (size, label) in zip(extraDimSizes, extraDimLabels):
                if label not in outFid.dimensions.keys():
                    outFid.createDimension(label, size)
                
            # write the variable to file
            vDims = ['row', 'col'] + extraDimLabels
            outFieldName = self.parmDict['outFieldNames'][field]
            varHand = outFid.createVariable(outFieldName, 'd', vDims, fill_value=self.parmDict['fillVal'])
            varHand[:] = outputArrays[field]
            
            # write variable attributes
            setattr(varHand, 'Units', self.parmDict['outUnits'][field])

        # close the output file
        outFid.close()

        # create an output dictionary keyed to output field names

        finalOutArrays = dict()
        for (k,v) in outputArrays.iteritems():
            finalOutArrays[self.parmDict['outFieldNames'][k]] = v
        return finalOutArrays 
            
class unweighted_filtered_MOPITT_avg_netCDF_out_func(wght_avg_netCDF):
    '''
    Output averager designed to work with MOPITT Level 2 data.  Emulates the 
    NASA developed averaging algorithm for level 3 (mostly) faithfully.  

    Following the NASA precedent, data are separated into either daytime or 
    nighttime values according to solar zenith angle.  Unfortunately, since 
    none of the NASA documentation actually specifies what cutoff value was
    used for solar zenith angle, the value is left up to the user with a 
    default of 85.  Also, in contrast to the NASA product, only one time (day
    or night) is actually included.  Which one is left to the user and noted
    in the attributes of the output file.

    Also following NASA precedent, data are filtered based on surface type.
    For cells where one surface type makes up more than 75% of the pixels,
    that surface type is used exclusively.  For cells where no surface type
    reaches the 75% threshold, all pixels are included.

    Again following the NASA algorithm, if pixels containing differing numbers
    of valid levels are present in a single grid cell, only the pixels 
    comprising the majority are retained.  This is tested using the retrieved
    CO mixing ratio profile field.  Note that this applies to BOTH 3D AND 
    NON-3D FIELDS.  That is, a surface measurement will still be filtered
    if the corresponding column measurement does not have the majority number
    of layers present, EVEN IF NO COLUMN MEASUREMENTS ARE REQUESTED.  Furthermore,
    all columns are filtered based on the chosen column measurement, even
    if they have information present in an 'invalid' layer

    The user is given the option of averaging each field assuming either a
    normal or log-normal distribution.  This is left to the user's discretion
    so make sure you know which fields it is appropriate to average and which
    should be log-averaged.

    For further details (don't get your hopes up) on the NASA algorithm, refer
    to 
        Deeter, Merritt N (2009). MOPITT (Measurements of Pollution in the 
          Troposphere) Validated Version 4 Product Users Guide.  Available
          from <http://www.acd.ucar.edu/mopitt/products.shtml>

    The output will be a netCDF file with the name determined in the standard
    way.  Will have the appropriate dimensions for those input fields being
    processed, as well as the fundamental rows/cols determined by grid_geo

    This class subclasses the generic wght_avg_netCDF.  It handles all 
    changes in interface in it's __init__ method and lets any actual calls
    filter up to the parent.

    The parameters dictionary must contain the following keys:
        time:           SEE DOCUMENTATION FOR wght_avg_vals_netCDF_out_func
                            NOTE: must be in TAI93 format
        longitude:      SEE DOCUMENTATION FOR wght_avg_vals_netCDF_out_func
        inFieldNames:   SEE DOCUMENTATION FOR wght_avg_vals_netCDF_out_func
        outFieldNames:  SEE DOCUMENTATION FOR wght_avg_vals_netCDF_out_func
        outUnits:       SEE DOCUMENTATION FOR wght_avg_vals_netCDF_out_func
        logNormal:      SEE DOCUMENTATION FOR wght_avg_vals_netCDF_out_func
        dimLabels:      SEE DOCUMENTATION FOR wght_avg_vals_netCDF_out_func
        dimSizes:       SEE DOCUMENTATION FOR wght_avg_vals_netCDF_out_func
        timeStart:      SEE DOCUMENTATION FOR wght_avg_vals_netCDF_out_func
        timeStop:       SEE DOCUMENTATION FOR wght_avg_vals_netCDF_out_func
        timeComparison: SEE DOCUMENTATION FOR wght_avg_vals_netCDF_out_func
        fillVal:        SEE DOCUMENTATION FOR wght_avg_vals_netCDF_out_func
        solZenAng: The string for the field associated with the solar zenith
            angle.  Must be in degrees
        solZenAngCutoff: The cutoff solar zenith angle value dividing night
            from day.  Pixels with SZA > solZenAngCutoff will be considered
            nighttime pixels.  90 is mathematically correct, values between 
            are typeically used in practice.  In degrees.  If SZA is exactly
            equal to the cutoff, it is included regardless of whether day
            or night was selected.
        dayTime: Boolean variable setting whether the desired output file 
            will be for the daytime or nighttime.  If set to "True", the output
            file will feature daylight retrievals only.  If set to "False", the
            output will feature night retrievals only.  Note that by most 
            estimates, daylight retrievals are higher quality.
        surfTypeField: The string for the field associated with the surface
            type.  This field is assumed to have integers corresponding to 
            different surface types.  No effort is made to distinguish 
            between surface types (only to ensure they are consistent) so 
            the mapping of integers to physical surface types is irrelevant.
        colMeasField: The string for the field associated with a column 
            measurement.  This can be any field with exactly one extra 
            dimension, provided it has NaN's at the same levels as other
            fields where appropriate.  The canonical field to use here is
            the retrieved CO mixing ratio profile.
    '''
    @staticmethod
    def parm_list():
        return ['time', 'longitude', 'inFieldNames', 'outFieldNames',
                'outUnits', 'logNormal', 'dimLabels', 'dimSizes', 'timeStart', 
                'timeStop', 'timeComparison', 'fillVal', 'solZenAngCutoff', 
                'solZenAng', 'dayTime', 'surfTypeField', 'colMeasField']
    @staticmethod
    def required_parms():
        return {'time' : ('The name of the field containing timestamps.  ' \
                          'Timestamps are assumed to be in the TAI-93 format.' \
                          '\n{ MOPITT - TIME }', None),
                'longitude' : ('The name of the field containing longitudes ' \
                               'at cell centers.  Longitudes should be in ' \
                               'degrees east.\n{ MOPITT - Longitude }', None),
                'inFieldNames' : ('The names of the fields desired to be ' \
                                  'output.  Input as comma-delimited list.', \
                                  'list'),
                'outFieldNames': ('The names of the output variables. (even ' \
                                  'if they are to be the same as input ' \
                                  'variables).  Should be a comma-delimited ' \
                                  'list co-indexed to inFieldNames', 'list'),
                'outUnits' : ('The units of the variables to be written out.' \
                              '  Should be a comma-delimited list co-indexed '\
                              'to inFieldNames', 'list'),
                'logNormal' : ('List of boolean strings that specify how to ' \
                               'take the averages of the corresponding fields.'\
                               '  If the string is "True" that field is ' \
                               'averaged assuming a lognormal distribution.  ' \
                               'If the string is "False" that field is ' \
                               'averaged assuming a normal distribution.  ' \
                               'Should be a comma-delimited list co-indexed ' \
                               'to inFieldNames', 'list'),                
                'dimLabels' : ('List of names of the extra dimensions in the ' \
                               'output file.  Must be a semicolon-delimited ' \
                               'list of comma-delimited strings. Fields with no'\
                               'extra dimensions may be left blank.  ' \
                               'For example, if there are four inFields, the ' \
                               'first and third of which have no extra ' \
                               'dimensions, the second of which has one ("foo"),'\
                               ' and the fourth has two ("foo" and "bar"), the '\
                               'dimLabels entry should look like this: '\
                               ';foo;;foo,bar  The outer (semicolon-delimited) '\
                               'list must be co-indexed to inFieldNames', 
                               'listoflists'),
                'dimSizes' : ('List of the sizes of the extra dimensions in the' \
                              ' output file.  Must be a semicolon-delimited list'\
                              ' of comma-delimited lists of integers.  Fields'\
                              'with no extra dimensions may be left blank.  ' \
                              'For example, if there are four inFields, the ' \
                              'first and third of which have no extra ' \
                              'dimensions, the second of which has one (which ' \
                              'has length four), and the fourth has two (which '\
                              'have lengths four and five, respectively), the '\
                              'dimSizes entry should look like this: ;4;;4,5 ' \
                              'The outer (semicolon-delimited) list must be ' \
                              'co-indexed to inFieldNames and all sub-lists ' \
                              'should be the same size as the corresponding ' \
                              'sublist in dimLabels.', 'listoflists'),
                'timeStart' : ('The earliest time for which data should be ' \
                               'recorded into the output file.  All times ' \
                               'before this time in the input file(s) will ' \
                               'be filtered out.  Must be in the format:  hh:' \
                               'mm:ss_MM-DD-YYYY', 'time'),
                'timeStop' : ('The latest time for which data should be ' \
                              'recorded into the output file.  All times after'\
                              ' this time in the input file(s) will be ' \
                              'filtered out.  Must be in the format: ' \
                              'hh:mm:ss_MM-DD-YYYY','time'),
                'timeComparison' : ('Must be set to either "local" or "UTC".  '\
                                    'Determines how the file timestamps are ' \
                                    'compared to the start/stop time.  If set '\
                                    'to "local", then the file timestamps are ' \
                                    'converted to local time on a pixel-by-pixel'\
                                    ' basis (using longitude to estimate time ' \
                                    'zone) before being compared to time ' \
                                    'boundaries.  If set to "UTC" the file ' \
                                    'timestamps (which are assumed to be in UTC)'\
                                    ' are compared against the start/stop time '\
                                    'directly.', None),
                'fillVal' : ('The value to use as a fill value in the output '\
                             'netCDF file.  This value will replace any '\
                             'missing or invalid output values', 'decimal'),
                'solZenAngCutoff' : ('The solar zenith angle that defines the '\
                                     'day to night transition (we use the SZA '\
                                     'to separate day and night pixels, which '\
                                     'should not be averaged together), in ' \
                                     'degrees.  The geometric value here would ' \
                                     'be 90.  Recommended value is 85.', 
                                     'decimal'),
                'solZenAng' : ('The name of the field containing the solar' \
                               ' zenith angle in degrees.  { MOPITT - Solar ' \
                               'Zenith Angle }', None),
                'dayTime' : ('Boolean variable that indicates ' \
			     'whether the output file should contain ' \
			     'values from day or night.  If set to ' \
			     '"True" the output file will have ' \
			     'daylight values.  If set to "False" ' \
			     'the output file will have night ' \
			     'values.', 'bool'),
                'surfTypeField' : ('The name of the field containing the ' \
                                   'surface type index.\n{ MOPITT - Surface ' \
                                   'Index }', None),
                'colMeasField' : ('The name of the field containing the ' \
                                  'column measurement that will be used to ' \
                                  'determine which levels are valid in a ' \
                                  'cell.  Canonically the retrieved CO mixing' \
                                  ' ratio profile field.  It is assumed that ' \
                                  'the field will have a layer dimension first' \
                                  ' and a 2-element second dimension (for ' \
                                  'values and std devs) of which we want the ' \
                                  'first slice.\n{ MOPITT - Retrieved CO Mixing '\
                                  'Ratio Profile }', None)}
    # variable signifying which list is to act as the master list index 
    __userKeys__ = "inFieldNames"
    def __init__(self, pDict):
        '''Convert input to format of parent input'''

        # make a shallow copy to the parameter dict, as we'll be making changes
        # and we don't want to mutate the argument
        parmDict = dict(pDict)

        # even though IO interface handles casting already,
        # a catchblock has been added here for safety
        # in case someone wants to use this class directly
        castDict = {'time':str, 'longitude':str,
                    'inFieldNames':list, 'outFieldNames':list,
                    'outUnits':list, 'logNormal':list,
                    'dimLabels':list, 'dimSizes':list,
                    'timeStart':tai93conv, 'timeStop':tai93conv,
                    'timeComparison':str, 'fillVal':float,
                    'solZenAngCutoff':float, 'solZenAng':str,
                    'dayTime':bool, 'surfTypeField':str,
                    'colMeasField':str}
        for (k,func) in castDict.items():
            try:
                parmDict[k] = func(parmDict[k])
            except TypeError:
                pass
                
        # by this point times are already converted to TAI93 standard
        # no need to convert here
        parmDict['timeConv'] = lambda(x):x
        
        # remove extraneous entries in parmDict.  They will be incorporated in
        # weighting and filtering functions
        SZAcut = parmDict.pop('solZenAngCutoff')
        SZAfield = parmDict.pop('solZenAng')
        dayTime = parmDict.pop('dayTime')
        surfField = parmDict.pop('surfTypeField')
        colMeasField = parmDict.pop('colMeasField')

        dayBool = dayTime
        # note which was chosen
        parmDict['notes'] = 'All values %s with cutoff at %6.2f' % \
            ('daytime' if dayBool else 'nighttime', SZAcut)
        
        # create weighting function
        def wghtFunc(parser, index, prevWght):
            '''
            Values not explicitly weighted.  Values not in desired part of 
            diurnal cycle (as determined by solar zenith angle) are given weight
            of 0 and therefore not included in final average
            '''
            SZA = parser.get_cm(SZAfield, index)
            if dayBool and SZA <= SZAcut:
                # we want day and it's day
                return 1
            elif not dayBool and SZA >= SZAcut:
                # we want night and it's night
                return 1
            else:
                return 0
        parmDict['weightFunction'] = wghtFunc

        # create filtering function
        def filterFunc(parser, indStack):
            '''
            Filter is twofold.  First filter checks if any surface type makes
            up 75% of the pixels in the cell.  If it does, all other surface 
            types are rejected.  Second filter checks if column retrievals have
            different numbers of valid retrievals.  If they do, then the pixels
            in the minority are rejected.  In the case of a tie the pixels with
            more levels present are retained.
            '''
            # short-circuit an empty indstack because it's a common case that's
            # difficult to efficiently code for
            if len(indStack) == 0:
                return numpy.array([])

            # first filter
            sTypes = numpy.array([parser.get_cm(surfField, ind) for ind in indStack]).squeeze()
            uniques = numpy.unique(sTypes)
            uniqueFracs = [float((sTypes == un).sum())/sTypes.size for un in uniques]
            cellType = None
            for (type,frac) in izip(uniques,uniqueFracs):
                # at most one value can meet threshold
                if frac >= .75:
                    cellType = type
            if cellType is None:
                # none met threshold, all are used
                sFlag = numpy.array([False]*len(sTypes))
            else:
                # one met threshold
                sFlag = sTypes != cellType
            
            # second filter
            columns = [parser.get_cm(colMeasField, ind)[:,0] for ind in indStack]
            nValidInCol = numpy.array([col.size - numpy.isnan(col).sum() for col in columns])
            uniqueNvals = set(nValidInCol)
            uniqueCounts = numpy.array([(nValidInCol == val).sum() for val in uniqueNvals])
            maxCount = uniqueCounts.max()
            maxNVals = [nv for (nv,c) in izip(uniqueNvals,uniqueCounts) if c == maxCount]
            # if there are multiples with same count, we want the highest number of valid
            # values, so we take the largest
            maxNVal = max(maxNVals)
            cFlag = numpy.array([nValid != maxNVal for nValid in nValidInCol])

            # combine the filters and return
            return numpy.logical_or(cFlag, sFlag)
        parmDict['filterFunction'] = filterFunc

        # invoke parent's constructor
        wght_avg_netCDF.__init__(self, parmDict)



class MODIS_simp_avg_netCDF_out_func(wght_avg_netCDF):

   @staticmethod
   def parm_list():
      # return ordered list of required parameters
      return ['timeStart', 'timeStop', 'timeComparison', 'time', 'fillVal', 'inFieldNames', 'outUnits', 'outFieldNames', 'extraDimSize', 'extraDimLabels']

   @staticmethod
   def required_parms():
      # return dictionary containing parm:('description',type) entries for each parameter
      return {'timeStart':('The earliest time for which data should be '\
                           'recorded into the output file.  All times '\
                           'before this time in the input file(s) will be '\
                           'filtered out.  Must be in the format: hh:mm:ss_'\
                           'MM-DD-YYYY','time'),
              'timeStop':('The latest time for which data should be recorded '\
                          'into the output file.  All times after this time '\
                          'in the input file(s) will be filtered out.   '\
                          'Must be in the format: hh:mm:ss_MM-DD-YYYY', 'time'),
              'timeComparison':('Must be set to either "local" or "UTC".  '\
                                'Determines how the file timestamps are '\
                                'compared to the start/stop time.  If set '\
                                'to "local", then the file timestamps are '\
                                'converted to local time on a pixel-by-pixel '\
                                'basis (using longitude to estimate time zone)'\
                                ' before being compared to time boundaries.  '\
                                'If set to "UTC" the file timestamps '\
                                '(which are assumed to be in UTC) are compared'\
                                ' against the start/stop time directly.', None),
              'time':('The name of the field containing timestamps.  '\
                      'Timestamps are assumed to be in the '\
                      'TAI-93 format.', None),
              'fillVal':('The value to use as a fill value in the output '\
                         'netCDF file.  This value will replace any '\
                         'missing or invalid output values', 'decimal'),
              'inFieldNames':('The names of the fields desired to be output.  '\
                              'Input as a comma-delimited list.', 'list'),
              'outFieldNames':('The names of the output variables. (even if '\
                               'they are to be the same as input variables).  '\
                               'Should be a comma-delimited list co-indexed '\
                               'to inFieldNames', 'list'),
              'outUnits':('The units of the variables to be written out.' \
                          '  Should be a comma-delimited list co-indexed '\
                          'to inFieldNames', 'list'),
              'extraDimSize':('List of the sizes of the extra dimensions in '\
                              'the output file.  Must be a comma-delimited '\
                              'list of integers.  Fields with no extra '\
                              'dimensions should be left blank.  For example, '\
                              'if there are four inFields, the first and '\
                              'third of which have no extra dimensions, the '\
                              'second of which has one (which has length '\
                              'four), and the fourth has one (which has '\
                              'length five), the extraDimSize entry should '\
                              'look like this: ,4,,5 '\
                              'The list must be co-indexed to inFieldNames', 
                              'listoflists'),
              'extraDimLabels':('List of names of the extra dimensions in the '\
                                'output file.  Must be a comma-delimited list '\
                                'of strings.  Fields with no extra dimensions '\
                                'should be left blank.  '\
                                'For example, if there are four inFields, '\
                                'the first and third of which have no extra '\
                                'dimensions, the second of which has one '\
                                '("foo"), and the fourth has one ("bar"), the '\
                                'extraDimLabels entry should look like this: '\
                                ',foo,,bar '\
                                'The list must be co-indexed to inFieldNames',
                                'listoflists'),
               'includePixelCount':('If set to true, the output will include '\
                                    'a field "ValidPixelCount" that contains '\
                                    'the number of valid pixels in each grid '\
                                    'cell.  Only pixels with nonzero weight  '\
                                    'are considered valid.','bool')}

   # variable signifying which list is to act as the master list index
   __userKeys__ = 'inFieldNames'
   def __init__(self, pDict):
      '''Convert input to format of parent input'''
      # make a shallow copy to the parameter dict, as we'll be making changes
      # and we don't want to mutate the argument
      self.parmDict = dict(pDict)

      # even though IO interface handles casting already
      # a catchblock has been added here for safety
      # in case someone wants to use this class directly
      castDict = {'time':str, 'inFieldNames':list,
                  'outUnits':list, 'outFieldNames':list,
                  'extraDimLabels':list, 'extraDimSize':list, 
                  'timeStart':tai93conv, 'timeStop':tai93conv,
                  'timeComparison':str, 'fillVal':float, 
                  'includePixelCount':boolCaster}
      for (k,func) in castDict.items():
         try:
            self.parmDict[k] = func(self.parmDict[k])
         except TypeError:
            pass
   
      # by this point times are already converted to TAI93 standard
      # no need to convert here
      self.parmDict['timeConv'] = lambda(x):x

      # create Weighting function
      def wghtFunc(parser, index, prevWght):
         '''
         There is no weighting function-- NASA has already weighted pixels
         based on their overall quality.  Only the highest quality pixels
         have been included for land-based data, and medium to high quality
         pixels have been weighted appropriately for ocean-based data.
         '''
         return 1
      self.parmDict['weightFunction'] = wghtFunc
  
      # create filtering function
      def filterFunc(parser, indStack):
         '''
         A filtering function is also not necessary for MODIS data.
         No pixels were removed by this filter.
         '''
         return numpy.zeros(indStack.shape)
 

   def __call__(self, maps, griddef, outfilename, verbose, version):
      '''Write out an averaged file in netcdf format.'''

      # Run sanity checks on values (e.g. logical start/stop times) 
      # raise IOerror on failure
      ###
      if self.parmDict['timeStart'] > self.parmDict['timeStop']:
         msg = 'Input start time must come before stop time.'
         raise IOError(msg)
      if (len(self.parmDict['inFieldNames']) !=  \
              len(self.parmDict['outFieldNames']) or  
              len(self.parmDict['inFieldNames']) !=  \
              len(self.parmDict['outUnits']) or      
              len(self.parmDict['inFieldNames']) !=  \
              len(self.parmDict['extraDimLabels'])):
         msg = 'All field/unit inputs ' + \
               'should have the same number of elements.'
         raise IOError(msg)

      # Create numpy arrays to hold data
      ###
      (minRow, maxRow, minCol, maxCol) = griddef.indLims()
      nRows = maxRow - minRow + 1
      nCols = maxCol - minCol + 1
      sumWght = numpy.zeros((nRows, nCols, 1)) # needs extra dim to generalize for 3D vars
      sumVars = dict()
      for field, size in zip(self.parmDict['inFieldNames'], self.parmDict['extraDimSize']):
         if size:
            sumVars[field] = numpy.zeros((nRows, nCols, int(size[0])))
         else:
            # pad with a singlet dim if it was 2D
            sumVars[field] = numpy.zeros((nRows, nCols, 1))
      nValidPixels = numpy.zeros((nRows, nCols))
      nNaNPixels = numpy.zeros((nRows, nCols))

      # loop over maps
      ###
      if not isinstance(maps, list):
         maps = [maps] # create list if we only got a single map
      for map in maps:
         pixSum = 0
         totofNaN = 0
         # open up context manager
         with map.pop('parser') as parser: # remove parser for looping
            if verbose:
               print('\nProcessing {0} for output at {1}.'.format(\
                      parser.name, str(datetime.datetime.now())))
               sys.stdout.write('\nData added from 0 pixels. ')
               sys.stdout.flush()
            # loop over gridboxes in map and calculate weights
            for (gridCell, pixTup) in map.iteritems():
               # translate gridCell to account for possible non-zero ll corner
               gridRow = gridCell[0]
               gridCol = gridCell[1]
               gridInd = (gridRow - minRow, gridCol - minCol)
               # get the values needed to calculate weight
               for (pxInd, unused_weight) in pixTup:


                  # check quality flags 
                  # and continue if pixel fails to meet requirements
                  # ***To be added as an option in a future update***
                  ###
                  '''
                  QAL_flag = parser.get_cm('Quality_Assurance_Land', pxInd)
                  if QAL_flag % 2 != 1:
                     QAO_flag = parser.get_cm('Quality_Assurance_Ocean', pxInd)
                     if QAO_flag % 2 != 1:
                        continue
                  '''

                  # check time within start/stop range
                  ###
                  time = parser.get_cm(self.parmDict['time'], pxInd)
                  if self.parmDict['timeComparison'] == 'local':
                     pixLon = parser.get_cm(self.parmDict['longitude'], pxInd)
                     offset = utils.UTCoffset_from_lon(pixLon)
                     time += offset
                  if time < self.parmDict['timeStart']\
                                   or time > self.parmDict['timeStop']:
                     continue
                                                      
                  # read in all the data, abandon ship if data is all NaN
                  ###
                  rawDataDict = {}
                  try:
                     for field in self.parmDict['inFieldNames']:
                        rawData = parser.get_cm(field, pxInd)
                        rawDataDict[field] = rawData
                     if all(numpy.isnan(el).all()\
                                    for el in rawDataDict.values()):
                        raise invalidPixCeption
                  except invalidPixCeption:
                     nNaNPixels[gridInd] += 1
                     totofNaN += 1
                     continue              
                  # compute the weight and add to the total for this cell
                  ###
                  pixSum = pixSum + 1
                  if verbose:
                      sys.stdout.write("\rData added from {0} pixels. NaNs found for {1} pixels."\
                                                         .format(pixSum, totofNaN))
                      sys.stdout.flush()
                  sumWght[gridInd] += 1
                  nValidPixels[gridInd] += 1
                  for field in self.parmDict['inFieldNames']:
                     if rawDataDict[field].size > 1:
                        sumVars[field][gridInd] = numpy.nansum([sumVars[field][gridInd], rawDataDict[field]], axis=0)
                     else:
                        sumVars[field][gridInd] = numpy.nansum([sumVars[field][gridInd][0], rawDataDict[field]])
            map['parser'] = parser # return parser to map

      # divide out variables by weights to get avgs
      ###
      oldSettings = numpy.seterr(divide='ignore')
      avgs = dict()
      with warnings.catch_warnings():
         warnings.simplefilter("ignore")
         for (field,var) in sumVars.iteritems():
            unfiltAvgs = var/sumWght
            filtAvgs = numpy.where(sumWght != 0, unfiltAvgs, self.parmDict['fillVal'])
            # strip trailing singlet for 2D arrays
            if filtAvgs.shape[-1] == 1:
               avgs[field] = filtAvgs.reshape(filtAvgs.shape[0:2])
            else:
               avgs[field] = filtAvgs
      numpy.seterr(divide=oldSettings['divide'])
     
      if verbose:
          print "\nWriting final output to file at {0}.".format(\
                                  str(datetime.datetime.now()))
      # associate coindexed parameters into dicts keyed on inFieldnames
      # so we can loop by field
      ###
      outFnames = dict(izip(self.parmDict['inFieldNames'], self.parmDict['outFieldNames']))
      units = dict(izip(self.parmDict['inFieldNames'], self.parmDict['outUnits']))
      extraDim = dict(izip(self.parmDict['inFieldNames'], self.parmDict['extraDimLabels']))

      # write out results to a netcdf file
      ###
      outFid = netCDF4.Dataset(outfilename, 'w', format='NETCDF3_CLASSIC')

      # create the 2 dimensions all files use
      ###
      outFid.createDimension('row', nRows)
      outFid.createDimension('col', nCols)

      # write global attributes
      ###
      setattr(outFid, 'Version', vsnmsg(version))
      # File_start_time, File_end_time, requirements for valid data, 
      # time comparison scheme, etc. written in same manner
      setattr(outFid, 'Input_files', ' '.join([map['parser'].name for map in maps]))
      setattr(outFid, 'Projection', griddef.__class__.__name__[:-8])

      for(k,v) in griddef.parms.iteritems():
         setattr(outFid, k, v)
      # loop over fields and write variables
      ###
      for field in self.parmDict['inFieldNames']:
         # create tuple of dimensions, defining new dim if necessary
         if len(avgs[field].shape) == 2:
            # only row/cols
            varDims = ('row', 'col')
         elif len(avgs[field].shape) == 3:
            # has extra dim
            dimName = extraDim[field][0]
            dimSize = avgs[field].shape[2]
            if dimName not in outFid.dimensions.keys():
               outFid.createDimension(dimName, dimSize)
            varDims = ('row', 'col', dimName)
  
         # create and assign value to variable
         ###
         varHandle = outFid.createVariable(outFnames[field], 'd', varDims, \
                                           fill_value=self.parmDict['fillVal'])
         varHandle[:] = avgs[field]
    
         # assign variable attributes
         ###
         setattr(varHandle, 'Units', units[field])

      if self.parmDict['includePixelCount']:
         varDims = ('row', 'col')
         varHandle = outFid.createVariable('ValidPixelCount', 'i', varDims, 
                                           fill_value=self.parmDict['fillVal'])
         varHandle[:] = nValidPixels
         varHandle = outFid.createVariable('NaNPixelCount', 'i', varDims,
                                           fill_value=self.parmDict['fillVal'])
         varHandle[:] = nNaNPixels

      outFid.close()
      # create a dict with the same data as avgs, but diff names
      ###
      outAvg = dict()
      for (k,v) in avgs.iteritems():
         outAvg[outFnames[k]] = v
      return outAvg

