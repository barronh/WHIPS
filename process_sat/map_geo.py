'''
Framework for remapping data.

All functions are named after the convention 
scheme_map_geo.  The function ValidMaps() automatically
maintains a list of implemented mapping schemes.

The interface for schemes is set up to allow
for maximum efficiency of the computations within
the scheme.  They therefore take only a parser and
a griddef (both instances, in that order)
as arguments, and return a dictionary.

The returned dictionary must be formatted as follows:
    -one key "parser" which corresponds to a pointer
    to the parser used to generate the map.
    -one key for each of the gridboxes, which is a tuple
    (row,col) for that gridbox.  These keys correspond
    to lists.  Each list contains tuples, each of which 
    must contain:
        -a tuple of the SWATH indices (this should
        be able to be fed into the get function of 
        the parser as is)
        - the data-independent weight (pass None if
        no weight is computed from this function)
'''
import sys
from itertools import izip
import datetime
import time
import pdb

import map_helpers

from shapely.prepared import prep
import shapely.geometry as geom
import numpy

def ValidMaps():
    '''Return a list of valid map names'''
    currentModule = sys.modules[__name__]
    names = dir(currentModule)
    return [el[:-8] for el in names if el.endswith("_map_geo")]

def global_intersect_map_geo(parser, griddef, verbose=True):
    '''
    For each pixel, find all gridcells that it intersects

    This function does not compute or save the fractional 
    overlap of individual pixels.  

    This function is designed to handle grids that span 
    the entire globe, with the cyclic point (the point 
    where longitude "wraps") ocurring at index [*,0]

    The function assumes that the projection performs the wrapping.  
    IE any pixels east of the far east edge of the domain will be
    on the west side of the projection.  The function will NOT
    automatically wrap pixels that fall off the east/west edges
    of the projection.

    Pixels with NaN for any vertex are rejected

    Assumptions:
        - Straight lines in projected space adequately 
        approximate the edges of pixels/gridcells.
        - polar discontinuities aren't of concern
        - we ARE dealing with a global projection that 
        cyclizes at index (*,0)
        - grid is rectilinear
        - pixels are convex polygons
    '''
    if verbose:
        print('Mapping '+parser.name+'\nat '+str(datetime.datetime.now()))  
    outer_indices = griddef.indLims()
    outer_poly = (geom.Polygon([(outer_indices[0], outer_indices[2]),
                               (outer_indices[0], outer_indices[3]),
                               (outer_indices[1], outer_indices[3]),
                               (outer_indices[1], outer_indices[2])]))
    outer_polyp = prep(outer_poly)              
    # create the dictionary we'll use as a map
    map = map_helpers.init_output_map(outer_indices)
    map['parser'] = parser
    tpreps = []
    tchecks = []
    tskips = []
    # we're going to hold onto both the prepared and unprepared versions
    # of the polys, so we can access the fully method set in the unprep
    # polys, but still do fast comparisons
    gridPolys = map_helpers.rect_grid_polys(outer_indices)
    prepPolys = map_helpers.rect_grid_polys(outer_indices)
    if verbose: print('prepping polys in grid')
    for polyk in prepPolys.keys():
        prepPolys[polyk] = prep(prepPolys[polyk]) # going to get compared a lot
    if verbose: print('done prepping polys in grid')
    cornersStruct = parser.get_geo_corners()
    (row, col) = griddef.geoToGridded(cornersStruct['lat'], \
                                      cornersStruct['lon']) 
    ind = cornersStruct['ind']
    if cornersStruct['lon'].shape[-1] == 4:
        lon = cornersStruct['lon']
        leftlon = lon[..., [0, 3]].max(-1)
        rightlon = lon[..., [1, 2]].min(-1)
        loncheck = (leftlon > rightlon)

    # reshape the matrixes to make looping workable
    row = row.reshape(-1,4)
    col = col.reshape(-1,4)
    loncheck = loncheck.ravel()
    ind = ind.reshape(row.shape[0],-1)
    if verbose: print('Intersecting pixels')
    # create the appropriate pixel(s) depending on whether
    # the pixels span the cyclic point
    minCol = griddef.indLims()[2]
    maxCol = griddef.indLims()[3] + 1
    midCol = (minCol+maxCol)/2.0
    dummy, (colNeg180, colPos180) = griddef.geoToGridded(numpy.array([0, 0]), numpy.array([-180, 180]))
    for (pxrow, pxcol, pxind, pxcheck) in izip(row, col, ind, loncheck):
        if verbose:
            t0 = time.time()
        if (numpy.any(numpy.isnan(pxrow)) or 
            numpy.any(numpy.isnan(pxcol))):
            continue # skip incomplete pixels
        pointsTup = zip(pxrow, pxcol)
        prelimPoly = geom.MultiPoint(pointsTup).convex_hull
        if not outer_polyp.intersects(prelimPoly.envelope):
            if verbose:
                tskips.append(time.time() - t0)
            continue
        (bbBot, bbLeft, bbTop, bbRight) = prelimPoly.bounds
        if bbLeft < midCol and bbRight > midCol:
            pointsLeft = [ (r,c) for (r,c) in pointsTup if c < midCol]
            pointsRight = [ (r,c) for (r,c) in pointsTup if c >= midCol]
            if pxcheck:
                pointsLeft += [ (bbBot, colNeg180), (bbTop, colNeg180) ]
                pointsRight += [ (bbBot, colPos180), (bbTop, colPos180) ]
            else:
                pointsLeft += [ (bbBot, minCol), (bbTop, minCol) ]
                pointsRight += [ (bbBot, maxCol), (bbTop, maxCol) ]
            polyLeft = geom.MultiPoint(pointsLeft).convex_hull
            polyRight = geom.MultiPoint(pointsRight).convex_hull
            splitArea = polyLeft.area + polyRight.area
            spanArea = prelimPoly.area
            if splitArea < spanArea:
                pixPolys = [polyLeft, polyRight]
            else:
                pixPolys = [prelimPoly]
        else:
            pixPolys = [prelimPoly]
        # try intersecting the poly(s) with all the grid polygons
        if verbose:
            t1 = time.time()
        for poly in pixPolys:
            for key in map_helpers.get_possible_cells(outer_indices, poly):
                # Prepped contains is much faster than 
                # than the combination of prepped intersection
                # and unprepped touches
                if prepPolys[key].intersects(poly) and (prepPolys[key].contains(poly) or not gridPolys[key].touches(poly)):
                        map[key].append((tuple(pxind), None))
        if verbose:
            t2 = time.time()
            tprep, tcheck = t1 - t0, t2 - t1
            tpreps.append(tprep)
            tchecks.append(tcheck)
    if verbose:
        print('Checked pixels:', len(tpreps))
        print('Pixel processing:', sum(tpreps))
        print('Pixel checking:', sum(tchecks))
        print('         Total:', sum(tpreps) + sum(tchecks))
        print('Skipped pixels:', len(tskips))
        print('Pixel skip time:',sum(tskips))
        print('Done intersecting.')
    return map

    
def regional_intersect_map_geo(parser, griddef, verbose=True):
    '''
    For each pixel, find all gridcells that it intersects

    This function does not compute or save the fractional
    overlap of individual pixels.  It simply stores the
    pixel indices themselves.
    
    This function is currently not configured to operate 
    on a global scale, or near discontinuities. 
    The current kludge to handle discontinuities is 
    relies on the assumption that any pixel of interest 
    will have at least one corner within the bounds of 
    the grid

    Pixels with NaN for any vertex are rejected
    
    Several assumptions are made:
        - Straight lines in projected space adequately 
        approximate the edges of pixels/gridcells.
        - polar discontinuities aren't of concern
        - we AREN'T dealing with a global projection
        - grid is rectilinear
        - pixels are convex polygons
    '''
    
    if verbose:
        print('Mapping '+parser.name+'\nat '+str(datetime.datetime.now()))
    outer_indices = griddef.indLims()
    map = map_helpers.init_output_map(outer_indices)
    map['parser'] = parser
    bounds = prep(map_helpers.rect_bound_poly(outer_indices))
    # we're going to hold onto both the prepared and unprepared versions
    # of the polys, so we can access the fully method set in the unprep
    # polys, but still do fast comparisons
    gridPolys = map_helpers.rect_grid_polys(outer_indices)
    prepPolys = map_helpers.rect_grid_polys(outer_indices)
    if verbose: print('prepping polys in grid')
    for poly in prepPolys.itervalues():
        poly = prep(poly)  # prepare these, they're going to get compared a lot
    if verbose: print('done prepping polys in grid')
    cornersStruct = parser.get_geo_corners()
    (row, col) = griddef.geoToGridded(cornersStruct['lat'], \
                                      cornersStruct['lon']) 
    ind = cornersStruct['ind']
    # reshape the matrixes to make looping workable
    row = row.reshape(-1,4)
    col = col.reshape(-1,4)
    ind = ind.reshape(row.shape[0],-1)
    if verbose:
        griddedPix = 0 
        print('Intersecting pixels')
        sys.stdout.write("Approximately 0 pixels gridded. ")
        sys.stdout.flush()
        for (pxrow, pxcol, pxind) in izip(row, col, ind):
            if numpy.any(numpy.isnan(pxrow)) or numpy.any(numpy.isnan(pxcol)):
                continue  # if we have only a partial pixel, skip
            elif not any([bounds.contains(geom.asPoint((r,c))) \
                       for (r,c) in izip(pxrow, pxcol)]):
                continue  # if none of the corners are in bounds, skip
            griddedPix += 1
            sys.stdout.write("\rApproximately {0} pixels gridded. ".\
                             format(griddedPix))
            sys.stdout.flush()
            pixPoly = geom.MultiPoint(zip(pxrow, pxcol)).convex_hull
            
            for key in map_helpers.get_possible_cells(outer_indices, pixPoly):
                # Prepped contains is much faster than 
                # than the combination of prepped intersection
                # and unprepped touches
                if prepPolys[key].intersects(pixPoly) and (prepPolys[key].contains(pixPoly) or not gridPolys[key].touches(pixPoly)):
                    map[key].append((tuple(pxind), None))
        print('Done intersecting.')
    else:
        for (pxrow, pxcol, pxind) in izip(row, col, ind):
            if numpy.any(numpy.isnan(pxrow)) or numpy.any(numpy.isnan(pxcol)):
                continue  # if we have only a partial pixel, skip
            elif not any([bounds.contains(geom.asPoint((r,c))) for (r,c) \
                                        in izip(pxrow, pxcol)]):
                continue  # if none of the corners are in bounds, skip
            pixPoly = geom.MultiPoint(zip(pxrow, pxcol)).convex_hull
            for key in map_helpers.get_possible_cells(outer_indices, pixPoly):
                if prepPolys[key].intersects(pixPoly) and not \
                                  gridPolys[key].touches(pixPoly):
                    map[key].append((tuple(pxind), None))
    return map

def point_in_cell_map_geo(parser, griddef, verbose=True):
    '''
    For each object, find the single cell to which it should be assigned.  This 
    cell is determined as the cell into which the representative lat/lon of the 
    pixel would fit in projected space.
    
    Cells are treated as open on upper right and closed on the lower left.  For
    practical purposes, this means that pixels falling on the boundary between 
    two cells will be assigned to the cell up and/or to their right.  Pixels 
    falling on the lower left boundaries of the griddable area will be assigned
    to a cell and those on the upper right boundaries of the griddable area
    will be discarded.
    
    This function computes no weights.  It simply assigns objects on the basis
    of a representative lat-lon given by the parser's get_geo_centers function.
    
    This can operate on any rectilinear rid.  So long as a lat-lon pair can be
    projected to a unique location, it will assign each pixel to one and only 
    one location.
    
    Several assumptions are made:
        - Straight lines in projected space adequately approximate the edges of
        gridcells.
        - Grid is rectilinear.
        - the lat/lon of the col,row origin is the lower left hand corner of 
        the 0,0 gridbox.
    '''
    if verbose:
        print('Mapping '+parser.name+'\nat '+str(datetime.datetime.now()))
    # Create an empty map to be filled
    mapLims = griddef.indLims()
    map = map_helpers.init_output_map(mapLims)
    map['parser'] = parser
    centersStruct = parser.get_geo_centers()
    # get the data to geolocate the pixels and reshape to make looping feasible
    ind = centersStruct['ind']
    (row, col) = griddef.geoToGridded(centersStruct['lat'], \
                                      centersStruct['lon'])
    row = numpy.floor(row.flatten()) # we floor values to get the cell indices
    col = numpy.floor(col.flatten())
    ind = ind.reshape(row.size, -1)
    # loop over pixels and grid as appropriate
    (minRow, maxRow, minCol, maxCol) = mapLims  # unpack this so we can use it
    if verbose:
        nGriddedPix = 0
        print('Assigning pixels to gridboxes.')
        sys.stdout.write("Approximately 0 pixels gridded. ")
        sys.stdout.flush()
        for (pxrow, pxcol, pxind) in izip(row, col, ind):
            if minRow <= pxrow <= maxRow and minCol <= pxcol <= maxCol:
                map[(pxrow, pxcol)].append((tuple(pxind), None))
                nGriddedPix += 1
                sys.stdout.write("\rApproximately {0} pixels gridded. "\
                                 .format(nGriddedPix))
                sys.stdout.flush()
        print('Done intersecting.')
    else:
        for (pxrow, pxcol, pxind) in izip(row, col, ind):
            if minRow <= pxrow <= maxRow and minCol <= pxcol <= maxCol:
                map[(pxrow, pxcol)].append((tuple(pxind), None))
    return map
    
    
    
        
        
