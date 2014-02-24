'''
Miscellaneous helper functions used in
the map_geo module
'''
import math, itertools, pdb

from shapely.geometry import Polygon

def rect_grid_polys((minRow, maxRow, minCol, maxCol)):
    '''
    Create an dictionary of polygons representing
    all the grid cells for a rectilinear grid.
    
    Dictionary is keyed as (row,col) tuples
    '''
    polys = dict()
    for row in range(minRow, maxRow+1):
        for col in range(minCol, maxCol+1):
            ll = (row,col)
            ul = (row+1,col)
            ur = (row+1,col+1)
            lr = (row, col+1)
            polys[(row,col)] = Polygon([ll, ul, ur, lr])
    return polys

def rect_bound_poly((minRow, maxRow, minCol, maxCol)):
    '''
    Create a rectangular polygon representing the 
    bounding box of a rectilinear grid.
    '''
    ll = (minRow, minCol)
    ul = (maxRow+1, minCol)
    ur = (maxRow+1, maxCol+1)
    lr = (minRow, maxCol+1)
    return Polygon([ll, ul, ur, lr])

def get_possible_cells((minRow, maxRow, minCol, maxCol), testPoly):
    '''
    Given the domain and a test polygon, returns a list of index tuples 
    that represent a subset of the domain which contains the polygon.
    Subset may be any size, from a single cell to the full original 
    domain.  Domain is specified as the indLims format

    NOTE: It is assumed that the test_poly at least intersects the bounding
    box of the full domain.  This must be tested BEFORE using this
    function.  
    '''
    (tp_minRow, tp_minCol, tp_maxRow, tp_maxCol) = tuple(
        [int(math.floor(lim)) for lim in testPoly.bounds])
    true_minRow = max(minRow, tp_minRow)
    true_maxRow = min(maxRow, tp_maxRow)
    true_minCol = max(minCol, tp_minCol)
    true_maxCol = min(maxCol, tp_maxCol)
    row_range = range(true_minRow, true_maxRow + 1)
    col_range = range(true_minCol, true_maxCol + 1)
    return itertools.product(row_range, col_range)

def iter_2_of_3(array3d):
    '''Generator iterating over top 2 dimensions'''
    for row in array3d:
        for elem in row:
            yield elem

def iter_all_but_final(arrayNd):
    '''Generator iterating over all but the final dimension'''
    flatterArray = arrayNd.reshape(-1,arrayNd.shape[-1])
    for row in flatterArray:
        yield row
            
def init_output_map((minRow, maxRow, minCol, maxCol)):
    '''
    Initialize and return a dict intended for
    use as an output from a mapping function with
    empty lists for each of the (row,col) keys
    '''
    map = dict()
    for row in range(minRow, maxRow+1):
        for col in range(minCol, maxCol+1):
            map[(row,col)] = []
    return map
