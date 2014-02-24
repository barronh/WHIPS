WHIPS
=====

Wisconsin Horizontal Interpolation Program for Satellites provides custom gridding of satellite data (i.e., Level 2 to custom Level 3)

PROJECT TITLE: WHIPS
PURPOSE OF PROJECT: Provide a well-documented, easy-to-use general-purpose
                    processing module for processing satellite data
VERSION: 1.2.2 (8/28/13)
AUTHORS: oberman, maki, strom
CONTACT: taholloway@wisc.edu; barronh@ufl.edu


ACKNOWLEDGEMENT POLICY
======================
Whenever you publish research or present results generated with WHIPS,
please include the following acknowledgement or an appropriate
equivalent:

    We wish to thank the University of Wisconsin--Madison for the 
    use and development of the Wisconsin Horizontal Interpolation
    Program for Satellites (WHIPS).  WHIPS was developed by Jacob
    Oberman and Tracey Holloway, and ongoing development by Barron
    Henderson with funding from the NASA Air Quality Applied Science
    Team (AQAST) and the Wisconsin Space Grant Consortium 
    Undergraduate Award.


QUICK START
===========


1. Install the program and run the built-in test module to confirm
that it is working properly.  Installation instructions can be found
in the file INSTALL.txt
 

2. Download whatever data you plan to process.  Currently, the program
is designed to process OMI NO2 DOMINO level 2 data, OMI NO2 NASA level
2 data, MOPITT CO data, and MODIS AOD level 2 data.  See the --filelist
argument documentation for more on aquiring data.


3. Navigate to the folder where whips.py installed or add it to
your path.  It should be the /bin folder corresponding to the /lib 
folder where your python packages live.  Invoke it as:
     
     whips.py --help

4. Follow the on-screen instructions, adding each of the required
parameters.  If you need help with the projection attributes or the
output function attributes, invoke the built-in help as:

     whips.py --AttributeHelp <function_name>

For detailed explanations of all parameters and attributes, see the
"Parameter Details" section below.


5. Invoke whips.py once for each output file you'd like to create.
Note that the software creates output files with only a single
timestep, so you'll need to invoke the command once for each timestep
(IE if you want a month with timesteps every day, you'll probably want
to write a shell script that calls the command once for each day)


6. Additionally, you may create a grid file for your chosen grid by 
including the --includeGrid flag followed by the desired filename.  
A file containing the gridcells used by the projection will be written 
to the same directory as the standard output file.


7. Concatenate your outputs if desired (the authors recommend the NCO
operators at http://nco.sourceforge.net/ if you're using a netCDF
output format) and carry on!  


INVOKING WHIPS - METHOD 1
============================================
Calling directly from the command line
--------------------------------------------

There are two primary methods of running whips.  The first method
is to give whips all the information and settings it needs as flags
in a single command line call.  The calls become very long; typing
them in by hand is therefore inadvisible.  The script can easily
be batch run by building the call in a script.

For clarity and readability, line continuation characters are used to
place each attribute on a separate line.  It is not required to break
up attributes like this, but the command line needs to see the
invocation as a single command, so if you want to break it onto
multple lines you must use line-continuation characters.

1A. Process MOPITT level 2 CO data, (Version 5) 
----------------------------------------------

     - Uses a 36km lambert conic conformal grid centered over North
     America
     - Writes out a 2D, 3D, and 4D parameter from the file

       whips.py --directory /where/you/have/input/files \
         --fileList MOP02T-20050101-L2V10.1.1.prov.hdf \
     --filetype MOPITT_CO_NASA_HDF_V5 \
     --gridProj lcc2par \
     --mapFunc point_in_cell \
         --outDirectory where/you/want/output \
         --outFileName descriptive_name.nc \
     --verbose True \
     --interactive True \
     --projAttrs xOrig:-2916000 yCell:36000 \
     refLon:-97 refLat:40 nCols:162 nRows:126 stdPar2:45 \
     stdPar1:33 xCell:36000 earthRadius:6370000 yOrig:-2268000 \
     "inFieldNames:Time,Retrieved CO Mixing Ratio Profile,Retrieved CO Surface Mixing Ratio" \
     outFieldNames:time,COprof,COsurf logNormal:False,True,True \
     timeStart:00:00:00_01-01-2005 timeStop:23:59:59_01-01-2005 \
     timeComparison:UTC fillVal:-9999.0 \
     solZenAngCutoff:85 dayTime:True

2A. Process OMI level 2 DOMINO NO2 data, (Version 2) 
---------------------------------------------------

     - Uses a 36km lambert conic conformal grid centered over North
     America
     - explicitly specifies 3 input files
     - includes an output grid file (this would be possible in any example, 
       but is only shown in this one)
     - Writes out a 2D, and 3D parameter from the file

      whips.py \
           --directory /where/you/have/input/files \
        --filetype OMI_NO2_KNMI_HDF_v2_0_postFeb2006 \
    --fileList OMI-Aura_L2-OMDOMINO_2011m0901t1502-o37926_v003-2011m1012t121342.he5 \
    OMI-Aura_L2-OMDOMINO_2011m0901t1641-o37927_v003-2011m1012t121458.he5 \
    OMI-Aura_L2-OMDOMINO_2011m0901t1820-o37928_v003-2011m1012t121613.he5 \
        --gridProj lcc2par \
    --mapFunc regional_intersect \
    --outDirectory /where/you/want/output \
        --outFileName some_descriptive_name.nc \
    --verbose True \
    --interactive False \
    --includeGrid /where/you/want/gridfile/output \
    --projAttrs xOrig:-48 yCell:36 refLon:-97 refLat:40 \
    nCols:30 nRows:32 stdPar2:45 stdPar1:33 xCell:36 \
    earthRadius:6370 yOrig:-552 \
    --outFuncAttrs \
        inFieldNames:Time,AveragingKernel,TroposphericVerticalColumn \
    outFieldNames:time,avKern,tropVCD \
    timeStart:00:00:00_09-01-2011 \
    timeStop:23:59:59_09-01-2011 timeComparison:UTC \
    fillVal:-9999\
    customCriteria=CloudFraction=[0:0.3],SolarZenithAngle=[0:0.85] \
    includePixelCount:False

3A. Process OMI level 2 NASA NO2 data (Version 1.2)
--------------------------------------------------

      - uses a 36 km lambert conic conformal grid
      - writes out 2 fields to an output file
      - considers 2 input files
      - looks at all the cornerfiles in the directory

      whips.py \
        --directory /where/you/have/input/files \
    --fileList OMI-Aura_L2-OMNO2_2011m0430t1440-o36120_v003-2011m0501t043317.he5 \
           OMI-Aura_L2-OMNO2_2011m0430t1619-o36121_v003-2011m0501t061955.he5 \
        --filetype OMI_NO2_NASA_HDF_v1_2 \
    --gridProj lcc2par \
    --mapFunc regional_intersect \
        --outDirectory /where/you/want/output \ 
    --outFileName some_file_name.nc \
    --verbose True \
    --interactive False \
    --projAttrs cornerDir:/directery/where/you/keep/cornerfiles \
      cornerFileList: \
      stdPar1:33 stdPar2:45 refLat:40 refLon:-97 xOrig:-2916000 \
      yOrig:-2268000 xCell:36000 yCell:36000 \
      nRows:126 nCols:162 earthRadius:6370000 \
      inFieldNames:ColumnAmountNO2Trop,Time \
          outFieldNames:tropVCD,time \
      timeComparison:local timeStart:00:00:00_04-30-2011 \
      timeStop:23:59:59_04-30-2011 \
      customCriteria=SolarZenithAngle=[0:85],CloudFraction=[0:0.3],VcdQualityFlags=~19,XTrackQualityFlags=0,RootMeanSquareErrorOfFit=[0:0.0003],TerrainReflectivity=[0:0.3] \
      fillVal:-9999.0 includePixelCount:False

4A. Process MODIS level 2 NASA AOD data (Version 1)
--------------------------------------------------

     - uses a 36km Lambert Conical projection centered over the US
     - writes out 2 fields to an output file
     - considers all input files in the directory
     - includes the number of level 2 pixels averaged for each level 3 pixel

     whips.py \
       --directory /where/you/have/input/files \
       --filetype MODIS_AOD_NASA_Collection_5 \
       --gridProj lcc2par \
       --mapFunc point_in_cell \
       --outDirectory /where/you/want/output \
       --outFileName some_descriptive_name.nc \
       --verbose True \
       --interactive True \
       --projAttrs xOrig:-2916000 yOrig:-2268000 \
       xCell:36000 yCell:36000 nRows:126 nCols:162 \
       refLat:40 refLon:-97 stdPar1:33 stdPar2:45 \
       earthRadius:6370000 \
       --outFuncAttrs \
       timeStart:00:00:00_07-03-2008 \
       timeStop:23:59:59_07-03-2008 \
       timeComparison:UTC \
       fillVal -9999.0 \
       inFieldNames:Corrected_Optical_Depth_Land,Effective_Optical_Depth_Average_Ocean \
       outFieldNames:AOD_Land,AOD_Ocean \
       includePixelCount:True

INVOKING WHIPS - METHOD 2
============================================
Using a formatted input file
--------------------------------------------

The second method of running WHIPS is to specify all the options
and information it needs in a specially-formatted input file.  The
command line call then consists only of specifying the location of
the input file WHIPS should read from.  

Below is an explanatory example of the format:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
FAKE INPUT FILE

BEGIN

. This line begins with a period.  That means it is a comment
. Comments will be ignored by WHIPS when processing the file

. All required flags must be specified with the flag name
. and the desired value.  These flags must be in capital case.
. All possible flags (including optional flags) shown below:

DIRECTORY = /where/you/had/your/input/files
FILETYPE = some_file_type
GRIDPROJ = some_grid_projection
MAPFUNC = some_map_function
OUTDIRECTORY = /where/you/want/output/files
OUTFILENAME = name_of_output_file
VERBOSE = True_or_False
INTERACTIVE = True_or_False
OUTFUNC = some_output_function
INCLUDEGRID = /absolute/path/to/output/file/for/grid

. FILELIST is delimited by spaces.  Leave out to use all files
. in DIRECTORY
FILELIST = list of files

. Projection and output function attributes are specified
. the same as flags.  Just as with the command line call, 
. the correct attributes must be present or the program
. will exit.

stdPar1 = 33

inFieldNames = ColumnAmountNO2Trop,Time

. Note that unlike with the command line call, quotations are
. not needed when whitespace is part of an attribute value

outUnits = Molecules cm^-2,seconds since epoch

. You can specify parameters in any order. This includes 
. intermixing attributes and flags

time = Time
OUTFILENAME = whatever_we_want

. WHIPS input files must begin with BEGIN and end with END

END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The examples below are equivalent to their counterparts 
above in the section on method 1.

1B. Process MOPITT level 2 CO data, (Version 5) 
----------------------------------------------

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
BEGIN

DIRECTORY = /where/you/have/input/files
FILELIST =  MOP02T-20050101-L2V10.1.1.prov.hdf 
FILETYPE = MOPITT_CO_NASA_HDF_V5 
GRIDPROJ = lcc2par 
MAPFUNC =  point_in_cell 
VERBOSE = True 
INTERACTIVE = True

. Projection Attributes 
xOrig = -2916000
yCell = 36000 
refLon = -97 
refLat = 40 
nCols = 162 
nRows = 126 
stdPar2 = 45 
stdPar1 = 33 
xCell = 36000 
earthRadius = 6370000 
yOrig = -2268000 

. output function attributes
inFieldNames = Time,Retrieved CO Mixing Ratio Profile,Retrieved CO Surface Mixing Ratio
outFieldNames = time,COprof,COsurf 
logNormal = False,True,True 
timeStart = 00:00:00_01-01-2005 
timeStop = 23:59:59_01-01-2005 
timeComparison = UTC 
fillVal = -9999.0
solZenAngCutoff = 85 
dayTime = True 

OUTDIRECTORY =  /where/you/want/output 
OUTFILENAME = descriptive_name.nc 

END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

whips.py --inFromFile nameOfAboveFile.txt


2B. Process OMI level 2 DOMINO NO2 data, (Version 2) 
---------------------------------------------------

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
BEGIN

DIRECTORY = /where/you/have/input/files 
FILETYPE =  OMI_NO2_KNMI_HDF_v2_0_postFeb2006 
FILELIST =  OMI-Aura_L2-OMDOMINO_2011m0901t1502-o37926_v003-2011m1012t121342.he5 OMI-Aura_L2-OMDOMINO_2011m0901t1641-o37927_v003-2011m1012t121458.he5 OMI-Aura_L2-OMDOMINO_2011m0901t1820-o37928_v003-2011m1012t121613.he5 
GRIDPROJ = lcc2par
MAPFUNC = regional_intersect
VERBOSE = True 
INTERACTIVE = False

OUTDIRECTORY = /where/you/want/output
OUTFILENAME = some_descriptive_name.nc
INCLUDEGRID = /where/you/want/output/grid/full/path.nc

. Proj attrs
xOrig = -48 
yCell = 36 
refLon = -97 
refLat = 40 
nCols = 30 
nRows = 32
stdPar2 = 45 
stdPar1 = 33 
xCell = 36 
earthRadius = 6370 
yOrig = -552 

. Output attrs
inFieldNames = Time,AveragingKernel,TroposphericVerticalColumn 
outFieldNames = time,avKern,tropVCD 
timeStart = 00:00:00_09-01-2011 
timeStop = 23:59:59_09-01-2011 
timeComparison = UTC 
fillVal = -9999 
cloudFractUpperCutoff = 0.3 
solarZenAngUpperCutoff = 85 
includePixelCount = False

END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

whips.py --inFromFile aboveFileName.txt


3B. Process OMI level 2 NASA NO2 data (Version 1.2)
--------------------------------------------------

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
. somefilename.txt
. Input file for whips.py

BEGIN

. Flags
DIRECTORY = /where/you/have/input/files
FILETYPE = OMI_NO2_NASA_HDF_v1_2
FILELIST = OMI-Aura_L2-OMNO2_2011m0430t1440-o36120_v003-2011m0501t043317.he5 OMI-Aura_L2-OMNO2_2011m0430t1619-o36121_v003-2011m0501t061955.he5
GRIDPROJ = lcc2par
MAPFUNC = regional_intersect
OUTDIRECTORY = /where/you/want/output
OUTFILENAME = some_file_name.nc
VERBOSE = True 
INTERACTIVE = False 

. Parser attributes
cornerDir = /directory/where/you/keep/cornerfiles 
cornerFileList = 

. Grid attributes
stdPar1 = 33 
stdPar2 = 45 
refLat = 40 
refLon = -97 
xOrig = -2916000
yOrig = -2268000 
xCell = 36000 
yCell = 36000
nRows = 126 
nCols = 162 
earthRadius = 6370000

. Output function attributes
inFieldNames = ColumnAmountNO2Trop,Time
outFieldNames = tropVCD,time 
timeComparison = local 
timeStart = 00:00:00_04-30-2011 
timeStop = 23:59:59_04-30-2011 
cloudFractUpperCutoff = .3
solarZenAngUpperCutoff = 85 
fillVal = -9999.0 
includePixelCount = False

END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

whips.py --inFromFile aboveFileName.txt

4B. Process MODIS level 2 NASA AOD data (Version 1)
--------------------------------------------------

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
BEGIN

DIRECTORY = /where/you/have/inputfiles
FILETYPE = MODIS_AOD_NASA_Collection_5
GRIDPROJ = lcc2par
MAPFUNC = point_in_cell
VERBOSE = True
INTERACTIVE = True

OUTDIRECTORY = /where/you/want/output
OUTFILENAME = some_descriptive_name.nc

. Proj attrs
xOrig = -2916000
yOrig = -2268000
xCell = 36000
yCell = 36000
nRows = 126
nCols = 162
refLat = 40
refLon = -97
stdPar1 = 33
stdPar2 = 45
earthRadius = 6370000

. Output attrs
inFieldNames = Corrected_Optical_Depth_Land,Effective_Optical_Depth_Average_Ocean
outFieldNames = AOD_Land,AOD_Ocean
timeStart = 00:00:00_07-03-2008
timeStop = 23:59:59_07-03-2008
timeComparison = UTC
fillVal = -9999.0
includePixelCount = True

END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

whips.py --inFromFile aboveFileName.txt



PARAMETER DETAILS
=================

Each input attribute is explained here.  Note that each output
function has a separate set of required inputs and that the
--outFuncAttrs is therefore broken down by output function.  Make sure
you're referencing the details for the attributes relevant to the
function you want to use.  

  --help
    REQUIRED: NO
    DEFAULT: N/A
    - Display the onscreen help message and exit the program.

  --inFromFile /complete/path/to/input_file.txt
      REQUIRED: NO
    DEFAULT: N/A
    - An input file specifying the input options and information
      for WHIPS.  The input file is an alternate interface - 
      everything that can be done with a command line call can
      also be done by specifying the input file and calling the
      input file using this flag.
    - The input file syntax rules are as follows.

        1) The file must have "BEGIN" at the top before specifying
           anything else.  Anything above the BEGIN line will be
           ignored.
        2) Comments are denoted by putting "." at the beginning of
           a line.  Comments are ignored when the input file is 
           parsed.
        3) All required flags must be included.  Values are assigned
           using the following syntax:
               FOO = bar
           All flags must be in caps case.  Optional flags may
           be specified in this fashion as well. 
        4) All attributes for the projection, map and possibly
           parser (right now only MOPITT requires parser attributes)
           must be specified in the same syntax used for flags. 
           Case-sensitive, so the case of the attribute names must
           exactly match that specified below

    See above for both an illustrative example input file as well
    as complete input files for the test cases. 

  --directory /path/to/input/directory
      REQUIRED: NO
    DEFAULT: the current working directory at time of invocation            
      - The input directory that the program will search for whatever
        input files are specified.  If those files are not found in
        this directory, the programs behavior is governed by the
        value of the --interactive flag

  --fileList file1 [file2] [file3] ...
      REQUIRED: NO
    DEFAULT: The list of all files in --directory (non-recursive)
    - The list of files that the program should attempt to
        process for output.  Most output functions (with the 
      exception of the function designed for MOPITT CO) are
        designed to accept an arbitrary number of inputs and only
        use that data which fits the requirements.  
    - filtering out files a priori which cannot contain the
         information required, for instance files known to be on the
      wrong date, saves computational power.  It is advisable to
      use this parameter to include only those files which may
          contain the desired information wherever possible.
    - Locations of data (current as of 09/21/12)
        MOPITT - <http://eosweb.larc.nasa.gov/HPDOCS/datapool/>
        NASA OMI - <http://mirador.gsfc.nasa.gov/>
        KNMI OMI - <http://www.temis.nl/airpollution/no2.html>
        MODIS AOD - <http://ladsweb.nascom.nasa.gov/data/search.html>

  --filetype { HDFknmiomil2_generic, OMI_NO2_KNMI_HDF_v2_0_preFeb2006,
             OMI_NO2_KNMI_v2_0_postFeb2006,
           HDFnasaomil2_generic, OMI_NO2_NASA_HDF_v1_2,
           HDFmopittl2_generic, MOPITT_CO_NASA_HDF_V5,
           HDFmodisl2_generic, MODIS_AOD_NASA_Collection_5 }
      REQUIRED: YES
    DEFAULT: N/A
    - The type of file we're attempting to read in.  Must be one
        of the options listed above.
    - Filetypes have two primary roles.  The first is to specify
      the correct parser so that the file may be read by the
      program.  The second is to provide default values for the
      attributes of the output function.
    - Each filetype has a default output function considered most 
      appropriate for that filetype.  Default attributes specified
      will generally be targeted at this particular output function.
      The default output function can always be overridden with 
      the --outFunc flag.
    - Each parser has an associated filetype that ends in "_generic"
      These filetypes still specify a default output function but do
      not provide default values for any of the output function 
      parameters.
    - Some filetypes require additional "parser parameters".  These 
      can be passed in at the command line under the --projAttrs flag
      and should be formatted just like any projection attribute
      or output function attribute if using an input text file.

          HDFknmiomil2_generic:
              Parser: HDFknmiomil2
          Default output function: OMNO2e_netCDF_avg
          Description: The generic filetype for
            the OMI NO2 data as processed by KNMI (the
            DOMINO retrieval)
          Output attributes supplied: None
          Parser attributes required: None

          OMI_NO2_KNMI_HDF_v2_0_preFeb2006
          Parser: HDFknmiomil2
          Default output function: OMNO2e_netCDF_avg
          Description: The filetype for OMI NO2 as 
            processed by KNMI (version 2.0 of the DOMINO
            retrieval).  All input files should be before
            Feb 2006 (when the vertical grid of the 
            meteorology model used in the retrieval changed).
            Would probably work but has not been tested 
            for previous versions of the product.
          Output attributes supplied: 
               overallQualFlag
             cloudFrac
             solarZenithAngle
             time
             longitude
             pixIndXtrackAxis
             outUnits
             extraDimLabel
             extraDimSize
              Parser attributes required: None

          OMI_NO2_KNMI_HDF_v2_0_postFeb2006
          Parser: HDFknmiomil2
          Default output function: OMNO2e_netCDF_avg
          Description: The filetype for OMI NO2 as 
            processed by KNMI (version 2.0 of the DOMINO
            retrieval).  All input files should be after
            Feb 2006 (when the vertical grid of the 
            meteorology model used in the retrieval changed).
            Would probably work but has not been tested 
            for previous versions of the product.
          Output attributes supplied: 
               overallQualFlag
             cloudFrac
             solarZenithAngle
             time
             longitude
             pixIndXtrackAxis
             outUnits
             extraDimLabel
             extraDimSize
          Parser attributes required: None

          HDFnasaomil2_generic
              Parser: HDFnasaomil2
          Default output function: OMNO2e_netCDF_avg
          Description: The generic filetype for the OMI NO2
            data as processed by NASA (the OMNO2 product).
          Output attributes supplied: None
          Parser attributes required:
               cornerDir      - the directory containing the 
                            auxillary corner files that the
                          parser needs if the selected map
                          function makes use of them .
                 cornerFileList - A list of the corner files in cornerDir
                       that should be searched for corner files
                      that match the input file.  Should be 
                      comma-separated list of arbitrary length.
                      Set to empty string (right hand side of 
                      equals sign or colon blank for input file
                      or command line, respectively) to use
                      all files in the cornerdir.  Order does
                      not matter, files will be matched based
                      on orbit number.

          OMI_NO2_NASA_HDF_v1_2
          Parser: HDFnasaomil2
          Default output function: OMNO2e_netCDF_avg
          Description: The filetype for OMI NO2 as 
            processed by NASA (v1.2 of the OMNO2 product).
            Would probably work but has not been tested 
            for previous versions of the product.
          Output attributes supplied:
               overallQualFlag
             cloudFrac
             solarZenithAngle
             time
             longitude
             pixIndXtrackAxis
             outUnits
             extraDimLabel
             extraDimSize
          Parser attributes required:
               cornerDir      - the directory containing the 
                            auxillary corner files that the
                          parser needs if the selected map
                          function makes use of them .
                 cornerFileList - A list of the corner files in cornerDir
                       that should be searched for corner files
                      that match the input file.  Should be 
                      comma-separated list of arbitrary length

              HDFmopittl2_generic
              Parser: HDFmopittl2
          Default output function: unweighted_filtered_MOPITT_avg_netCDF
          Description: The filetype for MOPITT CO 
            as processed by NASA.
          Output attributes supplied: None
          Parser attributes required: None

          MOPITT_CO_NASA_HDF_V5
              Parser: HDFmopittl2
          Default output function: unweighted_filtered_MOPITT_avg_netCDF
          Description: The filetype for version 5 of the 
            MOPITT CO retrieval as processed by NASA. 
            May work with previous versions of the product
            but has not been tested.
          Output attributes supplied:
               time
             longitude
             solZenAng
             surfTypeField
             colMeasField
             outUnits
             dimLabels
             dimSizes
              Parser attributes required: None

              HDFmodisl2_generic
              Parser: HDFmodisl2
          Default output function: MODIS_simp_avg
          Description: The generic filetype for
            the MODIS AOD data as processed by NASA.
            Support for QA-mean averaging will be
            added in a future update.
          Output attributes supplied: None
          Parser attributes required: None

          MODIS_AOD_NASA_Collection_5
              Parser: HDFmodisl2
          Default output function: MODIS_simp_avg
          Description: The filetype for level 2
            MODIS Aerosol Optical Depth data
            as processed by NASA (Collection 5)
          Output attributes supplied: 
                 time
             extraDimSize
             extraDimLabels
          Parser attributes required: None


  --gridProj {latlon, lcc2par}
      REQUIRED: YES
    DEFAULT: N/A
    - The grid projection used to define the target grid (the grid
      that we wish to regrid our data to).  Must be one of the above
      options.  Further details on these options are as follows:
     
        latlon   - The Plate Caree projection, also known as
                    the "unprojected" projection.  x and y are
                    mapped directly to longitude and latitude,
                    respectively.  
            
            REQUIRED PARAMETERS:
              xCell     - The size of a gridcell in the x
                   (longitude) direction.  In degrees.
              yCell  - The size of a gridcell in the y
                 (latitude) direction.  In degrees.
              xOrig  - The longitude of the lower-left
                   corner of the domain
              yOrig     - The latitude of the lower-left
                     corner of the domain
              nRows     - The number of rows in the grid.
              nCols     - The number of columns in the grid.

        lcc2par     - The Lambert Conic Conformal projection (2
                 parallel construction).  x and y are
                 transformed and scaled, then mapped to
                 latitude and longitude via the projection.
                 The projection parameters must be accurate to
                 get the correct output grid.  All required
                 parameters are available in the GRIDDESC file
                 associated with the MM3 modeling system. A
                 description of the GRIDDESC format can be
                 found at
                 http://www.baronams.com/products/ioapi/GRIDDESC.html

             REQUIRED PARAMETERS:
             stdPar1 - One of the 2 standard parallels
                     used    to define the Lambert Conic
                      Conformal projection.  Must be a
                   valid latitude, in degrees.
              stdPar2 - The second standard parallel used
                    to define the Lambert Conic
                   Conformal projection.  Set this
                    equal to the same value as stdPar1 
                   if the single-parallel form of the
                   projection is being used. Must be a
                   valid latitude in degrees.
             refLat     - The reference latitude upon which
                    the projection is centered. This is
                    the YCENT value in the GRIDDESC
                   file. In degrees
             refLon  - The reference longitude upon which
                    the  projection is centered.  This
                   is BOTH the XCENT and PROJ_GAMMA
                    values in the GRIDDESC file.  If
                    these values are not identical, do
                    not use this function. In degrees.
             xOrig     - The location of the origin in
                    projected x coordinates.  This
                    is the XORIG value in the GRIDDESC
                    file. In same units as earthRadius.
             yOrig     - The location of the origin in
                    projected y coordinates.  This
                    is the YORIG value in the GRIDDESC
                    file. In same units as earthRadius.
             xCell     - The x dimension of a cell, in
                    projected coordinates.  In the same
                    units as earthRadius.  This is the
                    XCELL value in the GRIDDESC file
             yCell     - The y dimension of a cell, in
                    projected coordinates.  In the same
                    units as earthRadius.  This is the
                    YCELL value in the GRIDDESC file
             nRows     - The number of rows in the grid.
             nCols   - The number of columns in the grid.
             earthRadius - The assumed radius of the Earth
                    (assumed spherical).  Must match
                   units used for xCell and yCell.

  --projAttrs name1:value1 name2:value2 ...
    REQUIRED: YES
    DEFAULT: N/A
    - The attributes required for the chosen projection.  Must
      have all the required attributes for that projection.  The
      required parameters for each projection are listed under
          that projection's name above.
    - Case-sensitive.  Attribute names must EXACTLY match those
      laid out above.

  --mapFunc {point_in_cell, regional_intersect, global_intersect}
      REQUIRED: YES
    DEFAULT: N/A
    - The mapping function that will be used to assign pixels to
        cell(s).  Certain datasets have restrictions on which
      mapping functions are usable.
    - Also responsible for computing the "geometric weight" of
        pixels.  That is, this function computes the weight unique
      to a cell/pixel combination.  At present, neither of the
      functions provide this functionality.
            
        point_in_cell - Maps pixels defined by a single
            lat/lon (usually the cell center) to whatever
            grid cell that point lies inside in projected
            space.  This assigns each pixel to a unique
            grid cell.  Grid cells are open on the top and
            right sides and closed on the left and lower
            sides (with directions defined according to
            the projected coordinate system).  This
            function is supported by all currently
            available input filetypes.  NOTE: for 
            filetypes where regional_intersect is 
            available it is strongly recommended to use
            regional_intersect over point_in_cell.

        regional_intersect - Maps pixels (as defined by 
            pairs of geocoordinates that nominally
            correspond to pixel corners) to ALL gridcells
            intersected.  No geometric weights are
            calculated.  This function is currently
            supported by HDFnasaomil2 and HDFknmiomil2
            filetypes only.  

            NOTE: regional intersect should NOT be
            used for grids that wrap all the way 
            around the world longitudinally (IE the
            east and west edges are the same meridian).
            Such grids may give unexpected results, 
            especially if any pixels span the
            discontinuity.

            Pixels with one ore more vertices with
            fill values coordinates are rejected

            Makes several assumptions:
              - Polar discontinuities not encountered
              - Projection is NOT global
              - grid is rectilinear in projected space.
              - Pixels are convex polygons.

            global_intersect - Maps pixels (as defined by
                      n pairs of geocoordinates that nominally
            correspond to the n corners of the 
            pixel) to ALL gridcells intersected.
            No geometric weights are calculated.
            
            This is roughly the same as 
            the regional_intersect approach
            save that it is slightly slower and 
            assumes that the grid wraps around the 
            world longitudinally. 

            Pixels with one ore more vertices with
            fill values coordinates are rejected

            NOTE: using this function with
            grids where the east and west boundary
            are NOT the same meridian is not
            supported and may lead to unexpected
            results.

            Makes several assumptions:
              - polar discontinuities not encountered
              - Projection is global and cyclizes
                along the east and west edges
              - grid is rectilinear in projected space
              - Pixels are convex polygons
        
  --outFunc {OMNO2e_netCDF_avg,unweighted_filtered_MOPITT_avg_netCDF,
         MODIS_simp_avg_netCDF}
      REQUIRED: NO
    DEFAULT: Depends on the value chose for --filetype
    - The function that computes the output and writes the output
        file.  Functions are given significant freedom, but all
      current functions take some kind of average and write it to
        an output file.
    - These functions are frequently designed around a particular
        instrument or input format.  Efforts are made to make them
        as general as possible, but specialized output functions are
      only guaranteed (and really should only be used) for the
      parser types for which they have been designed.  To this 
      end, the associated function will be chosen by default for 
      all filetypes, though it can always be overridden with this 
      command.
    - See the documentaiton for the --filetypes flag to see what
      the default output function for each filetype is.
    - In all cases where a fieldname must be given for a
         parameter it is the short name (the name used to access the
        field through the parser) that must be given.  The 
      fieldnames must correspond to the official field names 
          in the data file.  These can usually be found in the data
      documentation.  Documentation for a few commonly processed
      formats can be found at:
        OMI (KNMI) - <http://www.temis.nl/docs/OMI_NO2_HE5_1.0.2.pdf>
        OMI (NASA) - <http://toms.gsfc.nasa.gov/omi/no2/OMNO2_data_product_specification.pdf>
        MOPITT - <http://www.acd.ucar.edu/mopitt/v5_users_guide_beta.pdf>

              OMNO2e_netCDF_avg - Averaging algorithm based on the
             NASA OMI level 2 to level 3 processing
             algorithm.  Designed for the OMI level 2
             filetypes (HDFnasaomil2 and HDFknmiomil2) and 
            is the default for those filetypes use with 
            other filetypes is of questionable utility.  

            Outputs results to a netCDF file.
            
            Further details available in the
             official NASA documentation located at 
            <http://disc.sci.gsfc.nasa.gov/Aura/data-holdings/OMI/omno2e_v003.shtml>
            
            Assumptions:
              - Data is at most 3 dimensional.
              - Invalid pixels are marked with an overall
               quality flag.
              - Timestamps are in the TAI93 format.

         unweighted_filtered_MOPITT_avg_netCDF - Averaging
             algorithm based on the NASA algorithm for
             processing level 2 MOPITT CO data to level 3
             MOPITT CO data.  Designed for the MOPITT level
             2 filetypes (HDFmopittl2 only at present) and 
            is the default for that filetype. Use with 
            other filetypes is discouraged.
            
            Outputs results to a netCDF file.

            IMPORTANT: Only 1 input file may be used with
            this function.  Conveniently, NASA currently
            provides data in 1-day granules.

            Further details available in the official NASA
            documentation:
            
            Deeter, Merritt N (2009). MOPITT (Measurements
                of Pollution in the Troposphere) Validated
                Version 4 Product Users Guide.  Available from
                <http://www.acd.ucar.edu/mopitt/products.shtml>

            Assumptions/caveats:
              - All fields are filtered based on the
              number of valid layers present in the
              specified column field.  Including 2D
              fields.
              - Timestamps are in the TAI93 format.

         MODIS_simp_avg_netCDF - Averaging
             algorithm based on the NASA algorithm for
             processing level 2 MODIS AOD data to level 3
             MODIS AOD data.  Designed for the MODIS level
             2 filetypes (HDFmodisl2 and Collection 5) and 
            is the default for those filetypes. Use with 
            other filetypes is discouraged.
            
            Outputs results to a netCDF file.

            Further details available in the official NASA
            documentation:

            Hubanks, et al. (2008) MODIS Atmosphere L3 Gridded
                    Product Algorithm Theoretical Basis Document.
                Available from 
                            <http://modis-atmos.gsfc.nasa.gov/_docs/L3_ATBD_2008_12-04.pdf>

            Assumptions/caveats:
              - This averaging method does not factor in
                QA weights.  A more generalized QA-weighted MODIS
                output function will be added in a future release.
              - Timestamps are in the TAI93 format.


  --outFuncAttrs name1:value1 name2:value2
      REQUIRED: YES
    DEFAULT: N/A
    - The attributes required for the chosen output function 
      specified somehow.  Can be either supplied by the filetype
      or specified by the user.  When the user and filetype both
      supply output function attributes, the filetype takes 
      precedence.  The "_generic" filetypes are supplied to
      allow complete control of all output function attributes
      if desired.
    - To see which attributes are supplied by which filetypes,
      see the documentation for the --filetype flag above
    - If one of the "value" elements contains whitespace, enclose
      the entire name:value pair in double quotes.  For example:
          the:full_monty           <- okay
          "the:full monty"          <- okay
          the:full monty         <- not okay
    - In many cases, a comma delimited list is requested.  Make
      sure that elements of the list are not separated by spaces.
      For example:
          pythons:EricIdle,JohnCleese    <- okay
          pythons:EricIdle, JohnCleese    <- not okay
          "pythons:Eric Idle,John Cleese"    <- okay
          "pythons:Eric Idle, John Cleese"    <- not okay
    - Case-sensitive.  Attribute names must EXACTLY match those
      laid out below.
    - Below are the required parameters for each existing output
      function.  { } contain usable/recommended parameters for 
      applicable filetypes.  These are based on the current 
      versions of the data at the time of writing and should be
      double-checked:

        OMNO2e_netCDF_avg -
            customCriteria - Criteria for pixel 
                selection. Each ciriteria is separated
                by a comma, each critieria follows the 
                OMI NO2 field "Description" attribute
                as documented in section 3.4.1 of the 
                OMNO2d File Specification document
                version 1.1 (Jan 10, 2013). For example,            
                { OMI NASA OMNO2d - SolarZenithAngle=[0:85],VcdQualityFlags=~1,CloudFraction=[0:0.3]
                  OMI NASA OMNO2d - SolarZenithAngle=[0:85],VcdQualityFlags=~19,XTrackQualityFlags=0,RootMeanSquareErrorOfFit=[0:0.0003],TerrainReflectivity=[0:0.3] }
			cloudFrac - The name of the field containing
				the cloud fractions.
				{ OMI KNMI - CloudFraction
				  OMI NASA - CloudFraction }
            time - The name of the field containing the
                timestamps.  Timestamps are assumed to
                be in the TAI-93 format.
                { OMI KNMI - Time
                  OMI NASA - Time }
            longitude - The name of the field containing
                the longitudes at cell centers.
                Longitudes should be in degrees east.
                { OMI KNMI - Longitude
                  OMI NASA - Longitude }
            inFieldNames - The names of the fields desired
                to be output.  Input as comma
                delimited list.
            outFieldNames - The names of the output
                variables (even if they are to be the
                same as input variables).  Should be a
                comma-delimited list co-indexed to
                inFieldNames
            outUnits - The units of the variables to be
                written out.  Should be a
                comma-delimited list co-indexed to
                inFieldNames
            extraDimLabel - Label for the extra dimension
                (should the variable have an extra
                dimension).  Ignored in the case of a
                2D variable.  Should be a
                comma-delimited list co-indexed to
                inFieldNames
            extraDimSizes - The size of the extra
                dimensions (should the variable have
                an extra dimension).  For 2D
                variables, must be set to 0. (zero)
                Should be a comma-delimited list
                co-indexed to inFieldNames.
            timeComparison - Must be set to either "local"
                or "UTC".  Determines how the file
                timestamps are compared to the
                start/stop time.  If set to "local",
                then the file timestamps are converted
                to local time on a pixel-by-pixel
                basis (using longitude to estimate
                time zone) before being compared to
                time boundaries.  If set to "UTC" the
                file timestamps (which are assumed to
                be in UTC) are compared against the
                start/stop time directly.
            timeStart - The earliest time for which data
                should be recorded into the output
                file.  All times in input files before
                this time will be filtered out. Must 
                be in the format:
                    hh:mm:ss_MM-DD-YYYY
            timeStop - The latest time for which data
                should be recorded into the output
                files.  All times in input files 
                after this time will be filtered out.  
                Must be in the format:
                    hh:mm:ss_MM-DD-YYYY
            pixIndXtrackAxis - The dimension order (0
                based) of the "cross-track" dimension
                (whichever dimension has size 60).
                For all currently known cases set 
                equal to 1 (depends on the 
                construction of the parser function.  
                If you rewrite the parser, check 
                this).
            fillVal - The value to use as a fill value in
                the output netCDF file.  This value
                will replace any missing or invalid
                output values.
            includePixelCount - If set to "True", the    
                output file will include the field 
                "ValidPixelCount" which will have a 
                count of the valid pixels by cell. 
                If set to "False" this field will
                not be included
            
        unweighted_filtered_MOPITT_avg_netCDF -
            time - The name of the field containing
                timestamps.  Timestamps are assumed to
                be in the TAI-93 format.
                { MOPITT - Time }
            longitude - The name of the field containing
                the longitudes at cell centers.
                Longitudes should be in degrees east.
                { MOPITT - Longitude }
            inFieldNames - The names of the fields desired
                to be output.  Input as comma
                delimited list.
            outFieldNames - The names of the output
                variables (even if they are to be the
                same as input variables).  Should be a
                comma-delimited list co-indexed to
                inFieldNames
            outUnits - The units of the variables to be
                written out.  Should be a
                comma-delimited list co-indexed to
                inFieldNames
            logNormal - List of boolean strings that
                specify how to take the averages of
                the corresponding fields.  If the
                string is "True" that field is
                averaged assuming a lognormal
                distribution.  If the string is
                "False" that field is averaged
                assuming a normal distribution.
                Official documentation (linked above)
                has further information on when
                log-average is appropriate.  Should be
                a comma-delimited list co-indexed to
                inFieldNames
            dimLabels - List of names of the extra
                dimensions in the output file.  Must
                be a forward-slash-delimited list of
                comma-delimited lists of labels.  
                Fields with no extra dimensions may
                be left blank.  For example, if
                there are four inFields, the first
                and third of which have no extra
                dimensions, the second of which has
                one ("foo"), and the fourth has two
                ("foo" and "bar"), the dimLabels
                entry should look like this:
                    /foo//foo,bar
                The outer (slash-delimited) list
                must be    co-indexed to inFieldNames.
            dimSizes - List of the sizes of the extra
                dimensions in the output file.  Must
                be a forward-slash-delimited list of
                comma-delimited lists of integers.
                Fields with no extra dimensions may
                be left blank.  For example, if there
                are four inFields, the first and
                third of which have no extra
                dimensions, the second of which has
                one (which has length 4), and the
                fourth has two (which have lengths 
                four and five, respectively), the
                dimSizes entry should look like this:
                    /4//4,5
                The outer (slash-delimited list
                must be co-indexed to inFieldNames 
                and each inner (comma-delimited) list 
                should be the same size as the 
                corresponding sublist in dimLabels.
            timeStart - The earliest time for which data
                should be recorded into the output
                file.  All times before this time in
                the input file(s) will be filtered 
                out.  Must be in the format:
                    hh:mm:ss_MM-DD-YYYY
            timeStop - The latest time for which data
                should be recorded into the output
                file.  All times after this time in
                the input file(s) will be filtered
                out.  Must be in the format:
                    hh:mm:ss_MM-DD-YYYY
            timeComparison - Must be set to either "local"
                or "UTC".  Determines how the file
                timestamps are compared to the
                start/stop time.  If set to "local",
                then the file timestamps are converted
                to local time on a pixel-by-pixel
                basis (using longitude to estimate
                time zone) before being compared to
                time boundaries.  If set to "UTC" the
                file timestamps (which are assumed to
                be in UTC) are compared against the
                start/stop time directly.
            fillVal - The value to use as a fill value in
                the output netCDF file.  This value
                will replace any missing or invalid
                output values.
            solZenAngCutoff - The solar zenith angle that
                defines the day to night transition
                (we use the SZA to separate day and
                night pixels, which should not be
                averaged together).  The geometric
                value here would be 90.  Recommended 
                value is 85. In degrees.
            solZenAng - The name of the field containing
                the solar zenith angle (in degrees).
                { MOPITT - Solar Zenith Angle }
            dayTime - Boolean variable that indicates
                whether the output file should contain
                values from day or night.  If set to
                "True" the output file will have
                daylight values.  If set to "False"
                the output file will have night
                values.
            surfTypeField - The name of the field
                containing the surface type index.
                { MOPITT - Surface Index }
            colMeasField - The name of the field
                containing the column measurement that
                will be used to determine how many
                valid layers are present in the cell.
                This field must be 4 dimensional, with
                the first extra dimension being the
                level and the first element of the
                second extra dimension containing
                NaN's at the appropriate levels.
                { MOPITT - Retrieved CO Mixing Ratio Profile }
                MODIS_simp_avg_netCDF -
                        timeStart - The earliest time for which data 
                      should be recorded into the output file.  
                  All times before this time in the input 
                  file(s) will be filtered out.  Must be 
                  in the format: hh:mm:ss_MM-DD-YYYY
            timeStop - The latest time for which data should 
                     be recorded into the output file.  All 
                                 times after this time in the input file(s)
                                 will be filtered out.  Must be in the 
                                 format: hh:mm:ss_MM-DD-YYYY
            timeComparison - Must be set to either "local" 
                                       or "UTC".  Determines how the 
                                       file timestamps are compared to the
                                       start/stop time.  If set to "local", 
                                       timestamps are converted to local time 
                                       on a pixel-by-pixel basis (using 
                                       longitude to estimate time zone) 
                                       before being compared to time 
                                       boundaries.  If set to "UTC" the 
                                       file timestamps (which are assumed to 
                       be in UTC) are compared against the 
                       start/stop time directly.
            time - The name of the field containing timestamps.
                             Timestamps are assumed to be in the TAI-93 format.
                 { MODIS - Scan_Start_Time }
            fillVal - The value to use as a fill value in the 
                                output netCDF file.  This value will 
                                replace any missing or invalid output values
            inFieldNames - The names of the fields desired to be
                                     output.  Input as a comma-delimited list.
            outFieldNames - The names of the output variables. 
                          (even if they are to be the same 
                      as input variables).  Should be a
                      comma-delimited list co-indexed 
                      to inFieldNames
            outUnits - The units of the variables to be 
                     written out.  Should be a comma-delimited 
                 list co-indexed to inFieldNames
            extraDimSize - List of the sizes of the extra 
                         dimensions in the output file.  
                     Must be a comma-delimited list of
                     integers.  Fields with no extra 
                     dimensions should be left blank.
                     For example, if there are four
                     inFields, the first and third of which 
                     have no extra dimensions, the second 
                     of which has one (which has length four),
                     and the fourth has one (which has
                     length five), the extraDimSize entry 
                     should look like this: ,4,,5 
                     The list must be co-indexed 
                     to inFieldNames
            extraDimLabels - List of names of the extra dimensions
                           in the output file.  Must be a 
                       comma-delimited list of strings.
                       Fields with no extra dimensions 
                       should be left blank.  For example, 
                       if there are four inFields, the
                       first and third of which have no 
                       extra dimensions, the second of 
                       which has one ("foo"), and the fourth
                       has one ("bar"), the extraDimLabels 
                       entry should look like this: ,foo,,bar 
                       The list must be co-indexed
                       to inFieldNames

   
  --outDirectory /path/to/output/directory
      REQUIRED: YES
    DEFAULT: N/A
    - The output directory to which the output file(s) should be
      written.  Make sure that you have write permissions to this
      directory (the program will complain and quit out if you do
        not).

  --outFileName FileName
      REQUIRED: NO
    DEFAULT: output1
    - The name of the output file itself.  User is responsible for
      adding any file extensions here (IE .nc if it's a netCDF,
      .txt if it's an ASCII).  Output will be written in
      outDirectory under this name.

  --includeGrid GridFileName
        REQUIRED: NO
        DEFAULT: N/A
        - Supply this flag along with the absolute path to a filename
          to which to write out the latitudes and longitudes of the gridcells
          defined by the selected projection.

  --verbose {True,False}
      REQUIRED: NO
    DEFAULT: True
    - Determines how much command-line output the software
      provides while running.  The default behavior (--verbose
      True) provide command line updates for most major
      subprocesses inside the software.  Setting "False" here will
      cause the software to be completely silent while running.

  --interactive {True,False}
      REQUIRED: NO
    DEFAULT: False
    - Determines how the program will handle invalid/nonexistent
      files.  Under the default behavior (--interactive False)
      the software will automatically ignore any file it can't
      process and continue processing any other files in
      --fileList.  If set to True, exectuion will be suspended and
      the user will be given several options when an invalid file
      is encountered.

  --AttributeHelp ProjectionName/OutputFunctionName/FileType [...]
      REQUIRED: NO
    DEFAULT: N/A
    - Prints a message explaining the required attributes for a
      given output function or projection and then exits the
      program.
    - Inputting a filetype outputs the required attributes for
      the default output function associated with that filetype.
