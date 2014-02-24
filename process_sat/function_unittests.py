'''
Run unit tests on the parse, geo, out, and map functions
'''
import unittest
import os
import sys
import tempfile
from itertools import izip, product
import pdb

import numpy
import netCDF4
import tables

import parse_geo
import grid_geo
import map_geo
import out_geo
import utils

class Helpers:

    @staticmethod 
    def genNan(shape):
        zeros = numpy.zeros(shape)
        zeros[:] = numpy.NaN
        return zeros

def does_pytables_close_all(filename):
    '''Returns true if pytables shows all files as closed after one is closed, false otherwise'''
    fid1 = tables.openFile(filename, mode='r')
    fid2 = tables.openFile(filename, mode='r')
    fid1.close()
    retVal = not fid2.isopen
    fid2.close()
    return retVal
    
def skipUnlessSamples():
    '''Skip this test unless our folder has the sample data'''
    dir = os.path.dirname(__file__)
    sDataPath = os.path.join(dir, 'sample_data')
    reason = "Skipped due to lack of sample data.  Download data to use."
    return unittest.skipUnless(os.path.exists(sDataPath), reason)

class fakeParser(parse_geo.GeoFile):
    """Duck-types a parser for testing purposes"""
    def __init__(self, filename, subtype='', extension=None):
        '''Create hash for data'''
        parse_geo.GeoFile.__init__(self, filename, subtype, extension)
        self._next_data = dict()
    def prime_corners(self, lat, lon, ind):
        """
        A properly structured record array constructed from
        these inputs will be returned the next time
        get_geo_corners() is called
        """
        # figure out dims of ind
        indDim = ind.shape[-1] if len(ind.shape)==len(lat.shape) else 1
        dtype = [('lat', lat.dtype, (4,)),
                 ('ind', ind.dtype, (indDim,)),
                 ('lon', lon.dtype, (4,))]
        self._next_corners = numpy.zeros(lat.shape[:-1], dtype=dtype)
        self._next_corners['lat'] = lat
        self._next_corners['lon'] = lon
        self._next_corners['ind'] = ind
    def prime_centers(self, lat, lon, ind):
        """
        A properly structured record array with
        the correct fields will be created from the
        passed data and will be returned when next
        get_geo_centers() is called
        """
        # figure out dims of ind
        if len(ind.shape)==len(lat.shape)+1:
            indDtype = ('ind', ind.dtype, (ind.shape[-1],))
        else:
            indDtype = ('ind', ind.dtype, tuple())
        dtype = [('lat', lat.dtype, tuple()),
                 indDtype,
                 ('lon', lon.dtype, tuple())]
        self._next_centers = numpy.zeros(lat.shape, dtype=dtype)
        self._next_centers['lat'] = lat
        self._next_centers['lon'] = lon
        self._next_centers['ind'] = ind
    def prime_get(self, key, data):
        """
        The data passed will be retrievable via
        both get and get_cm (both of which 
        will retrieve the exact same values).  Index
        tuples passed will be passed on without change,
        and whatever size array is generated 
        will be returned.  Data can only be retrieved
        with key, which can be any hashable
        """
        self._next_data[key] = data
    def get(self, key, indices=None):
        return self._next_data[key][indices]
    def get_cm(self, key, indices=None):
        return self._next_data[key][indices]
    def __enter__(self):
        return self
    def __exit__(self, unused_one, unused_two, unused_three):
        pass
    def get_geo_corners(self):
        return self._next_corners
    def get_geo_centers(self):
        return self._next_centers


class TestParserParentInstantiated(unittest.TestCase):

    def setUp(self):
        dir = os.path.dirname(__file__)
        self.emptyfile = os.path.join(dir, 'sample_data', 'empty')
        self.obj = parse_geo.GeoFile(self.emptyfile)

    def test_get_name(self):
        self.assertEqual(self.obj.name, self.emptyfile)
        
    def test_get_extension(self):
        self.assertEqual(self.obj.ext, '')
        
    def test_get_subtype(self):
        self.assertEqual(self.obj.sub, '')
        
    def test_no_get_geo_corners(self):
        with self.assertRaises(NotImplementedError):
            self.obj.get_geo_corners()
            
    def test_no_get_geo_centers(self):
        with self.assertRaises(NotImplementedError):
            self.obj.get_geo_centers()
        

class TestParserParentInstantiation(unittest.TestCase):

    def setUp(self):
        dir = os.path.dirname(__file__)
        self.emptytxtfile = os.path.join(dir, 'sample_data', 'empty.foo')

    def test_pulls_extension(self):
        obj = parse_geo.GeoFile(self.emptytxtfile)
        self.assertEqual(obj.name, self.emptytxtfile)
        self.assertEqual(obj.ext, 'foo')
        self.assertEqual(obj.sub, '')
        
    def test_overrides_extension(self):
        obj = parse_geo.GeoFile(self.emptytxtfile, extension='bar')
        self.assertEqual(obj.name, self.emptytxtfile)
        self.assertEqual(obj.ext, 'bar')
        self.assertEqual(obj.sub, '')
        
    def test_takes_subtype(self):
        obj = parse_geo.GeoFile(self.emptytxtfile, subtype='baz')
        self.assertEqual(obj.name, self.emptytxtfile)
        self.assertEqual(obj.ext, 'foo')
        self.assertEqual(obj.sub, 'baz')

@unittest.skip("Skipped until it can be rewritten to reflect changes since UI update")    
class TestGetParserFunc(unittest.TestCase):


    def setUp(self):
        dir = os.path.dirname(__file__)
        self.emptyNoExt = os.path.join(dir, 'sample_data', 'empty')
        self.emptyExt = os.path.join(dir, 'sample_data', 'empty.foo')
        self.hdf = os.path.join(dir, 'sample_data', 'omiknmil2sample.hdf')
        self.he5 = os.path.join(dir, 'sample_data', 'omiknmil2sample.he5')
        self.hdfNasa = os.path.join(dir, 'sample_data', 'ominasal2sample.hdf')

    def test_retrieve_parent(self):
        obj = parse_geo.get_parser(self.emptyNoExt)
        self.assertIsInstance(obj, parse_geo.GeoFile)
        self.assertEqual(obj.name, self.emptyNoExt)
        self.assertEqual(obj.ext, '')
        self.assertEqual(obj.sub, '')
        
    def test_retrieve_nonexistent_ext(self):
        obj = parse_geo.get_parser(self.emptyExt)
        self.assertIsInstance(obj, parse_geo.GeoFile)
        self.assertEqual(obj.name, self.emptyExt)
        self.assertEqual(obj.ext, 'foo')
        self.assertEqual(obj.sub, '')
        
    def test_retrieve_generic_hdf(self):
        obj = parse_geo.get_parser(self.hdf)
        self.assertIsInstance(obj, parse_geo.HDF_File)
        self.assertEqual(obj.name, self.hdf)
        self.assertEqual(obj.ext, 'hdf')
        self.assertEqual(obj.sub, '')
        
    def test_retrieve_subtype_knmiomil2(self):
        obj = parse_geo.get_parser(self.hdf, subtype='knmiomil2')
        self.assertIsInstance(obj, parse_geo.HDFknmiomil2_File)
        self.assertEqual(obj.name, self.hdf)
        self.assertEqual(obj.ext, 'hdf')
        self.assertEqual(obj.sub, 'knmiomil2')
        
    def test_retrieve_subtype_nasaomil2(self):
        obj = parse_geo.get_parser(self.hdfNasa, subtype='nasaomil2')
        self.assertIsInstance(obj, parse_geo.HDFnasaomil2_File)
        self.assertEqual(obj.name, self.hdfNasa)
        self.assertEqual(obj.ext, 'hdf')
        self.assertEqual(obj.sub, 'nasaomil2')
        
    def test_retrieve_sub_override_ext(self):
        obj = parse_geo.get_parser(self.he5, extension='hdf', subtype='knmiomil2')
        self.assertIsInstance(obj, parse_geo.HDFknmiomil2_File)
        self.assertEqual(obj.name, self.he5)
        self.assertEqual(obj.ext, 'hdf')
        self.assertEqual(obj.sub, 'knmiomil2')
        
    def test_retrieve_generic_nonexistent_override(self):
        obj = parse_geo.get_parser(self.hdf, subtype='bar')
        self.assertIsInstance(obj, parse_geo.GeoFile)
        self.assertEqual(obj.name, self.hdf)
        self.assertEqual(obj.ext, 'hdf')
        self.assertEqual(obj.sub, 'bar')

@skipUnlessSamples()    
class TestGeneralHDFParser(unittest.TestCase):


    def setUp(self):
        dir = os.path.dirname(__file__)
        self.HDF_Filename = os.path.join(dir, 'sample_data', 'omiknmil2sample.hdf')
        self.nonHDF_Filename = os.path.join(dir, 'sample_data', 'empty')
        
    def test_instantiate(self):
        obj = parse_geo.HDFFile(self.HDF_Filename)
        self.assertIsInstance(obj, parse_geo.HDFFile)
        
    def test_autofill_extension(self):
        obj = parse_geo.HDFFile(self.HDF_Filename)
        self.assertEqual(obj.ext, 'hdf')
        
    def test_override_extension(self):
        obj = parse_geo.HDFFile(self.HDF_Filename, extension='foo')
        self.assertEqual(obj.ext, 'foo')
        
    def test_override_subtype(self):
        obj = parse_geo.HDFFile(self.HDF_Filename, subtype='foo')
        self.assertEqual(obj.sub, 'foo')
        
    def test_reject_nonHDF(self):
        with self.assertRaises(IOError):
            unused_obj = parse_geo.HDFFile(self.nonHDF_Filename)

@skipUnlessSamples()
class TestKnmiOmiL2Parser(unittest.TestCase):


    def setUp(self):
        dir = os.path.dirname(__file__)
        self.fname = os.path.join(dir, 'sample_data', 'omiknmil2sample.hdf')
        self.parser = parse_geo.HDFknmiomil2_File(self.fname, subtype='knmiomil2')
        self.TM4presLev = numpy.array([3.287814, 171.674, 817.214, 2153.9016,
                                       4216.4746, 6889.528, 9949.71, 13104.406,
                                       16025.621, 18367.186, 19833.895,
                                       20333.965, 20134.355, 19552.172,
                                       18687.787, 17590.826, 16314.338,
                                       14915.031, 13452.953, 11988.191,
                                       10568.022, 9219.092, 7953.047, 
                                       6770.6543, 5431.8076, 4011.101,
                                       2825.2905, 1877.6597, 1160.5549,
                                       654.3158, 327.66162, 140.40256, 
                                       48.790634, 10.706806])
        self.checkKeys = ['AirMassFactor',
                          'AveragingKernel',
                          'CloudFraction',
                          'MeasurementQualityFlags',
                          'TM4PressurelevelA',
                          'TerrainHeight',
                          'TroposphericVerticalColumn',
                          'TroposphericColumnFlag',
                          'LatitudeCornerpoints',
                          'Time']
        
    def valid_keys(self):
        return [  'AirMassFactor',
                  'AirMassFactorGeometric',
                  'AirMassFactorTropospheric', 
                  'AssimilatedStratosphericSlantColumn',
                  'AssimilatedStratosphericVerticalColumn',
                  'AveragingKernel',
                  'CloudFraction',
                  'CloudFractionStd',
                  'CloudRadianceFraction',
                  'GhostColumn',
                  'InstrumentConfigurationId',
                  'MeasurementQualityFlags',
                  'SlantColumnAmountNO2',
                  'SlantColumnAmountNO2Std',
                  'SurfaceAlbedo',
                  'TM4PressurelevelA',
                  'TM4PressurelevelB',
                  'TM4SurfacePressure',
                  'TM4TerrainHeight',
                  'TM4TropoPauseLevel',
                  'TerrainHeight',
                  'TotalVerticalColumn',
                  'TotalVerticalColumnError',
                  'TroposphericColumnFlag',
                  'TroposphericVerticalColumn',
                  'TroposphericVerticalColumnError',
                  'TroposphericVerticalColumnModel',
                  'VCDErrorUsingAvKernel',
                  'VCDTropErrorUsingAvKernel',
                  'GroundPixelQualityFlags',
                  'Latitude',
                  'LatitudeCornerpoints',
                  'Longitude',
                  'LongitudeCornerpoints',
                  'SolarAzimuthAngle',
                  'SolarZenithAngle',
                  'Time',
                  'ViewingAzimuthAngle',
                  'ViewingZenithAngle']
        
    def avKern_6_5(self):
        return numpy.array([0, .335, .728, .741, .807, .838, .863, .898, .937,
                          .986, 1.038, 1.077, 1.092, 1.083, 1.073, 1.068,  
                          1.068, 1.066, 1.056, 1.047, 1.039, 1.031, 1.026,
                          1.018, 1.009, 1.005, .993, .980, .957, .933, .875,
                          .824, .819, .925])
        
    def test_instantiate(self):
        self.assertIsInstance(self.parser, parse_geo.HDFknmiomil2_File)
        
    def test_autofill_extension(self):
        self.assertEqual(self.parser.name, self.fname)
        
    def test_override_extension(self):
        obj = parse_geo.HDFknmiomil2_File(self.fname, extension='foo')
        self.assertEqual(obj.ext, 'foo')
        
    def test_override_subtype(self):
        self.assertEqual(self.parser.sub, 'knmiomil2')
        
    def test_reject_nonHDF(self):
        dir = os.path.split(self.fname)[0]
        emptyFile = os.path.join(dir, 'empty')
        with self.assertRaises(IOError):
            parse_geo.HDFknmiomil2_File(emptyFile)
    
    def test_get_retrieve_all_vars(self):
        validKeys = self.valid_keys()
        for key in validKeys:
            var = self.parser.get(key)
            self.assertIsInstance(var, numpy.ndarray)
        
    def test_get_cm_retrieve_all_vars(self):
        validKeys = self.valid_keys()
        with self.parser as p:
            for key in validKeys:
                var = p.get_cm(key)
                self.assertIsInstance(var, numpy.ndarray)
                
    def test_get_cm_closes_when_done(self):
        closes_all = does_pytables_close_all(self.fname)
        fid = tables.openFile(self.fname, mode='r')
        with self.parser as p:
            unused_var = p.get_cm('Time')
        if closes_all:
            self.assertFalse(fid.isopen)
        else:
            self.assertTrue(fid.isopen)
            fid.close()
        
    def test_get_retrieves_right_sizes(self):
        sizes = [1622*60,
                 34*1622*60,
                 1622*60,
                 1622,
                 34,
                 1622*60,
                 1622*60,
                 1622*60,
                 4*1622*60,
                 1622]
        for (key, size) in izip(self.checkKeys, sizes):
            self.assertEqual(self.parser.get(key).size, size)
            
    def test_get_cm_retrieves_right_sizes(self):
        sizes = [1622*60,
                 34*1622*60,
                 1622*60,
                 1622,
                 34,
                 1622*60,
                 1622*60,
                 1622*60,
                 4*1622*60,
                 1622]
        with self.parser as p:
            for (key, size) in izip(self.checkKeys, sizes):
                self.assertEqual(p.get_cm(key).size, size)                   
            
    def test_get_retrieves_right_type(self):
        types = [numpy.float32,
                 numpy.float32,
                 numpy.float32,
                 numpy.uint8,
                 numpy.float32,
                 numpy.int16,
                 numpy.float32,
                 numpy.int8,
                 numpy.float32,
                 numpy.float64]
        self.longMessage = True
        for (key, type) in izip(self.checkKeys, types):
            self.assertEqual(self.parser.get(key).dtype, type, 
                             msg="Type did not match for field %s." % key)
            
    def test_get_cm_retrieves_right_type(self):
        types = [numpy.float32,
                 numpy.float32,
                 numpy.float32,
                 numpy.uint8,
                 numpy.float32,
                 numpy.int16,
                 numpy.float32,
                 numpy.int8,
                 numpy.float32,
                 numpy.float64]
        self.longMessage = True
        with self.parser as p:
            for (key, type) in izip(self.checkKeys, types):
                self.assertEqual(p.get_cm(key).dtype, type, 
                                 msg="Mismatched type in field %s." % key)
                         
    def test_get_retrieves_first_element(self):
        genNan = Helpers.genNan
        vals = [genNan(1),
                numpy.array([-32.767]*34),
                numpy.array([-32.767]),
                numpy.array([2]),
                self.TM4presLev,
                numpy.array([-32768], 'int16'),
                genNan(1),
                numpy.array([-127], 'int8'),
                numpy.array([69.41081, 69.392624, 68.29076, 68.273315]),
                numpy.array([5.782786531883371*pow(10,8)])]
        for (key, val) in izip(self.checkKeys, vals):
            numpy.testing.assert_array_almost_equal(self.parser.get(key, (0,0)), val, decimal=3)
            
    def test_get_cm_retrieves_first_element(self):
        genNan = Helpers.genNan
        vals = [genNan(1),
                numpy.array([-32.767]*34),
                numpy.array([-32.767]),
                numpy.array([2]),
                self.TM4presLev,
                numpy.array([-32768], 'int16'),
                genNan(1),
                numpy.array([-127], 'int8'),
                numpy.array([69.41081, 69.392624, 68.29076, 68.273315]),
                numpy.array([5.782786531883371*pow(10,8)])]
        with self.parser as p:
            for (key, val) in izip(self.checkKeys, vals):
                numpy.testing.assert_array_almost_equal(p.get_cm(key, (0,0)), val, decimal=3, err_msg='Mismatch in %s' % key)
                            
    def test_get_retrieves_nan_at_189_2(self):
        genNan = Helpers.genNan
        vals = [genNan(1),
                numpy.array([-32.767]*34),
                numpy.array([.759]),
                numpy.array([0]),
                self.TM4presLev,
                numpy.array([52]),
                genNan(1),
                numpy.array([-127], 'int8'),
                numpy.array([65.86792, 65.95058, 65.240585, 65.32061]),
                numpy.array([5.78279031260359*pow(10,8)])]
        for (key, val) in izip(self.checkKeys, vals):
            numpy.testing.assert_array_almost_equal(self.parser.get(key, (189,2)), val, decimal=3)
            
    def test_get_cm_retrieves_nan_at_189_2(self):
        genNan = Helpers.genNan
        vals = [genNan(1),
                numpy.array([-32.767]*34),
                numpy.array([.759]),
                numpy.array([0]),
                self.TM4presLev,
                numpy.array([52]),
                genNan(1),
                numpy.array([-127], 'int8'),
                numpy.array([65.86792, 65.95058, 65.240585, 65.32061]),
                numpy.array([5.78279031260359*pow(10,8)])]
        with self.parser as p:
            for (key, val) in izip(self.checkKeys, vals):
                numpy.testing.assert_array_almost_equal(p.get_cm(key, (189,2)), val, decimal=3)
                    
    def test_get_retrieves_last_element(self):
        genNan = Helpers.genNan
        vals = [genNan(1),
                numpy.array([-32.767]*34),
                numpy.array([-32.767]),
                numpy.array([2]),
                self.TM4presLev,
                numpy.array([-32768], 'int16'),
                genNan(1),
                numpy.array([-127], 'int8'),
                numpy.array([80.10472, 79.99527, 80.49099, 80.378]),
                numpy.array([5.78284540309156*pow(10,8)])]
        for (key, val) in izip(self.checkKeys, vals):
            numpy.testing.assert_array_almost_equal(self.parser.get(key, (1621, 59)), val, decimal=3)
            
    def test_get_cm_retrieves_last_element(self):
        genNan = Helpers.genNan
        vals = [genNan(1),
                numpy.array([-32.767]*34),
                numpy.array([-32.767]),
                numpy.array([2]),
                self.TM4presLev,
                numpy.array([-32768], 'int16'),
                genNan(1),
                numpy.array([-127], 'int8'),
                numpy.array([80.10472, 79.99527, 80.49099, 80.378]),
                numpy.array([5.78284540309156*pow(10,8)])]
        with self.parser as p:
            for (key, val) in izip(self.checkKeys, vals):
                numpy.testing.assert_array_almost_equal(p.get_cm(key, (1621, 59)), val, decimal=3)
            
    def test_get_retrieves_valid_at_6_5(self):
        avKern = self.avKern_6_5
        vals = [numpy.array([3.9495022]),
                avKern(),
                numpy.array([1.0]),
                numpy.array([0]),
                self.TM4presLev,
                numpy.array([0]),
                numpy.array([.58010143*9.9999999*pow(10,14)]),
                numpy.array([-1]),
                numpy.array([73.72836, 73.68924, 73.12749, 73.089355]),
                numpy.array([5.782786651906409*pow(10,8)])]
        for (key, val) in izip(self.checkKeys, vals):
            numpy.testing.assert_allclose(self.parser.get(key, (6,5)), val, rtol=1E-5)
            
    def test_get_cm_retrieves_valid_at_6_5(self):
        avKern = self.avKern_6_5
        vals = [numpy.array([3.9495022]),
                avKern(),
                numpy.array([1.0]),
                numpy.array([0]),
                self.TM4presLev,
                numpy.array([0]),
                numpy.array([.58010143*9.9999999*pow(10,14)]),
                numpy.array([-1]),
                numpy.array([73.72836, 73.68924, 73.12749, 73.089355]),
                numpy.array([5.782786651906409*pow(10,8)])]
        with self.parser as p:
            for (key, val) in izip(self.checkKeys, vals):
                numpy.testing.assert_allclose(p.get_cm(key, (6,5)), val, rtol=1E-5)
                
    def test_get_zero_order_array_from_single_point(self):
        self.assertEqual(numpy.rank(self.parser.get('CloudFraction', (6,5))), 0) #non-Nan, applies scale-offset
        self.assertEqual(numpy.rank(self.parser.get('SolarZenithAngle', (6,5))), 0) #non-Nan, no applies scale-offset
        self.assertEqual(numpy.rank(self.parser.get('CloudFraction', (0,0))), 0) #Nan, applies scale-offset
        self.assertEqual(numpy.rank(self.parser.get('SolarZenithAngle', (6,5))), 0) #Nan, no apply scale offset
        
    def test_get_cm_zero_order_array_from_single_point(self):
        with self.parser as p:
            self.assertEqual(numpy.rank(p.get_cm('CloudFraction', (6,5))), 0) #non-Nan, applies scale-offset
            self.assertEqual(numpy.rank(p.get_cm('SolarZenithAngle', (6,5))), 0) #non-Nan, no applies scale-offset
            self.assertEqual(numpy.rank(p.get_cm('CloudFraction', (0,0))), 0) #Nan, applies scale-offset
            self.assertEqual(numpy.rank(p.get_cm('SolarZenithAngle', (6,5))), 0) #Nan, no apply scale offset
    
@skipUnlessSamples()        
class TestKnmiOmiL2GetGeoCorners(TestKnmiOmiL2Parser):


    def setUp(self):
        TestKnmiOmiL2Parser.setUp(self)
        self.geoarray = self.parser.get_geo_corners()
        
    def assert_array_equal(self, x, y, err_msg='', verbose=True):
        numpy.testing.assert_array_equal(x, y, err_msg, verbose)
        
    def test_has_correct_fields(self):
        self.assertItemsEqual(self.geoarray.dtype.names, ('lat', 'lon', 'ind'))

    def test_geocorners_size(self):
        self.assertEqual(self.geoarray.size, 1622*60)
        
    def test_geocorners_types(self):
        self.assertEqual(self.geoarray['lat'].dtype, numpy.float32)
        self.assertEqual(self.geoarray['lon'].dtype, numpy.float32)
    
    def test_geocorners_lat_first_element(self):
        firstElLat = numpy.array([69.41081, 69.392624, 68.29076, 68.273315])
        aFlat = self.geoarray['lat'].reshape((-1,4))
        for i in range(aFlat.shape[0]):
            if numpy.allclose(aFlat[i],firstElLat): break
        numpy.testing.assert_array_equal((self.geoarray['ind'].reshape((-1,2)))[i],(0,0))
        
    def test_geocorners_lon_first_element(self):
        firstElLon = numpy.array([119.52651, 119.691345, 118.76868, 118.924706])
        aFlat = self.geoarray['lon'].reshape((-1,4))
        for i in range(aFlat.shape[0]):
            if numpy.allclose(aFlat[i],firstElLon): break
        numpy.testing.assert_array_equal((self.geoarray['ind'].reshape((-1,2)))[i],(0,0))
        
    def test_geocorners_lat_last_element(self):
        finalElLat = numpy.array([80.10472, 79.99527, 80.49099, 80.378])
        aFlat = self.geoarray['lat'].reshape((-1,4))
        for i in range(aFlat.shape[0]):
            if numpy.allclose(aFlat[i],finalElLat): break
        numpy.testing.assert_array_equal((self.geoarray['ind'].reshape((-1,2)))[i],(1621, 59))
        
    def test_geocorners_lon_last_element(self):
        finalElLon = numpy.array([-167.0503, -167.30943, -173.45343, -173.64778])
        aFlat = self.geoarray['lon'].reshape((-1,4))
        for i in range(aFlat.shape[0]):
            if numpy.allclose(aFlat[i],finalElLon): break
        numpy.testing.assert_array_equal((self.geoarray['ind'].reshape((-1,2)))[i],(1621, 59))
        
    def test_geocorners_lat_random_element(self):
        randElLat = numpy.array([65.86792, 65.95058, 65.240585, 65.32061])
        aFlat = self.geoarray['lat'].reshape((-1,4))
        for i in range(aFlat.shape[0]):
            if numpy.allclose(aFlat[i],randElLat): break
        numpy.testing.assert_array_equal((self.geoarray['ind'].reshape((-1,2)))[i],(189, 2))
        
    def test_geocorners_lon_random_element(self):
        randElLon = numpy.array([56.252617, 56.47187, 57.96543, 58.184425])
        aFlat = self.geoarray['lon'].reshape((-1,4))
        for i in range(aFlat.shape[0]):
            if numpy.allclose(aFlat[i],randElLon): break
        numpy.testing.assert_array_equal((self.geoarray['ind'].reshape((-1,2)))[i],(189, 2))
        
    def test_geocorners_can_feed_ind_to_get(self):
        latAxis = self.geoarray.dtype.names.index('lat')
        lonAxis = self.geoarray.dtype.names.index('lon')
        indAxis = self.geoarray.dtype.names.index('ind')
        flatarray = self.geoarray.flatten()
        for i in range(0,flatarray.size, 5000):
            row = flatarray[i]
            self.assert_array_equal(row[latAxis], self.parser.get('LatitudeCornerpoints', row[indAxis]))
            self.assert_array_equal(row[lonAxis], self.parser.get('LongitudeCornerpoints', row[indAxis]))
            
    def test_geocorners_can_feed_ind_to_get_cm(self):
        latAxis = self.geoarray.dtype.names.index('lat')
        lonAxis = self.geoarray.dtype.names.index('lon')
        indAxis = self.geoarray.dtype.names.index('ind')
        flatarray = self.geoarray.flatten()
        with self.parser as p:
            for i in range(0, flatarray.size, 500):
                row = flatarray[i]
                self.assert_array_equal(row[latAxis], p.get_cm('LatitudeCornerpoints', row[indAxis]))
                self.assert_array_equal(row[lonAxis], p.get_cm('LongitudeCornerpoints', row[indAxis]))
                
@skipUnlessSamples()                
class TestKnmiOmiL2GetGeoCenters(TestKnmiOmiL2Parser):


    def setUp(self):
        TestKnmiOmiL2Parser.setUp(self)
        self.geoarray = self.parser.get_geo_centers()

    def assert_array_equal(self, x, y, err_msg='', verbose=True):
        numpy.testing.assert_array_equal(x, y, err_msg, verbose)
        
    def test_has_correct_fields(self):
        self.assertItemsEqual(self.geoarray.dtype.names, ('lat', 'lon', 'ind'))

    def test_geocenters_size(self):
        self.assertEqual(self.geoarray.size, 1622*60)
        
    def test_geocenters_types(self):
        self.assertEqual(self.geoarray['lat'].dtype, numpy.float32)
        self.assertEqual(self.geoarray['lon'].dtype, numpy.float32)
        
    def test_geocenters_lat_first_element(self):
        firstElInd = (0,0)
        almostFlat = self.geoarray['ind'].reshape((-1,2))
        for i in range(almostFlat.shape[0]):
            if (almostFlat[i] == firstElInd).all(): break
        geoData = (self.geoarray['lat'].reshape(-1))[i]
        self.assertTrue(numpy.isnan(geoData))
        
    def test_geocenters_lon_first_element(self):
        firstElInd = (0,0)
        aFlat = self.geoarray['ind'].reshape((-1,2))
        for i in range(aFlat.shape[0]):
            if (aFlat[i] == firstElInd).all(): break
        geoData = (self.geoarray['lon'].reshape(-1))[i]
        self.assertTrue(numpy.isnan(geoData))
        
    def test_geocenters_lat_last_element(self):
        finalElInd = (1621,59)
        aFlat = self.geoarray['ind'].reshape((-1,2))
        for i in range(aFlat.shape[0]):
            if (aFlat[i] == finalElInd).all(): break
        geoData = (self.geoarray['lat'].reshape(-1))[i]
        self.assertTrue(numpy.isnan(geoData))
        
    def test_geocenters_lon_last_element(self):
        finalElInd = (1621,59)
        aFlat = self.geoarray['ind'].reshape((-1,2))
        for i in range(aFlat.shape[0]):
            if (aFlat[i] == finalElInd).all(): break
        geoData = (self.geoarray['lon'].reshape(-1))[i]
        self.assertTrue(numpy.isnan(geoData))
        
    def test_geocenters_lat_random_element(self):
        randElLat = 66.44143
        mask = abs(self.geoarray['lat']-randElLat) < 0.0000005
        numpy.testing.assert_array_equal(self.geoarray['ind'][mask],[(186,3)])
        
    def test_geocenters_lon_random_element(self):
        randElLon = 57.16488
        mask = abs(self.geoarray['lon']-randElLon) < .0000005 
        numpy.testing.assert_array_equal(self.geoarray['ind'][mask],[(189,2)])
    
    def test_geocenters_can_feed_ind_to_get(self):
        latAxis = self.geoarray.dtype.names.index('lat')
        lonAxis = self.geoarray.dtype.names.index('lon')
        indAxis = self.geoarray.dtype.names.index('ind')
        flatarray = self.geoarray.flatten()
        for i in range(0,flatarray.size, 5000):
            row = flatarray[i]
            self.assert_array_equal(row[latAxis], self.parser.get('Latitude', row[indAxis]))
            self.assert_array_equal(row[lonAxis], self.parser.get('Longitude', row[indAxis]))

    def test_geocenters_can_feed_ind_to_get_cm(self):
        latAxis = self.geoarray.dtype.names.index('lat')
        lonAxis = self.geoarray.dtype.names.index('lon')
        indAxis = self.geoarray.dtype.names.index('ind')
        flatarray = self.geoarray.flatten()
        with self.parser as p:
            for i in range(0, flatarray.size, 500):
                row = flatarray[i]
                self.assert_array_equal(row[latAxis], p.get_cm('Latitude', row[indAxis]))
                self.assert_array_equal(row[lonAxis], p.get_cm('Longitude', row[indAxis]))

class TestParser(unittest.TestCase):
    '''Superclass containing utilities for parser test classes'''

    def __improved_iadd(self, baselist, extension):
        try:
                baselist.__iadd__(extension)
        except TypeError:
                baselist.__iadd__((extension,))
        return baselist

    def bldFlatList(self, argList):
        return reduce(self.__improved_iadd, argList, list())
    
    def bldFlatArray(self, argList):
        return numpy.array(self.bldFlatList(argList))

@skipUnlessSamples()
class TestNASAOmiL2Parser(TestParser):
    

    def setUp(self):
        dir = os.path.dirname(__file__)
        sample_dir = os.path.join(dir, 'sample_data')
        self.fname = os.path.join(sample_dir, 'ominasal2sample.hdf')
        self.parser = parse_geo.HDFnasaomil2_File(self.fname, subtype='knmiomil2',
                                                  cornerDir=sample_dir,
                                                  cornerFileList=['ominasacornersample.hdf'])
        self.checkKeys = ['ColumnAmountNO2Trop', 'PolynomialCoefficients', 
                          'CloudFraction', 'InstrumentConfigurationId',
                          'FitQualityFlags', 'UnpolFldLatBandQualityFlags', 
                          'Time', 'Latitude', 'SolarZenithAngle']
        self.validKeys = ['ColumnAmountNO2', 'ColumnAmountNO2Std', 'ColumnAmountNO2Initial',
                          'ColumnAmountNO2InitialStd', 'ColumnAmountNO2Trop', 'ColumnAmountNO2TropStd',
                          'ColumnAmountNO2BelowCloud', 'ColumnAmountNO2BelowCloudStd',
                          'ColumnAmountNO2Unpolluted', 'ColumnAmountNO2UnpollutedStd',
                          'ColumnAmountNO2Polluted', 'ColumnAmountNO2PollutedStd',
                          'TropFractionUnpolluted', 'TropFractionUnpollutedStd',
                          'SlantColumnAmountNO2', 'SlantColumnAmountNO2Std',
                          'RingCoefficient', 'RingCoefficientStd', 'SlantColumnAmountO3',
                          'SlantColumnAmountO3Std', 'SlantColumnAmountH2O',
                          'SlantColumnAmountH2OStd', 'SlantColumnAmountO2O2', 
                          'SlantColumnAmountO2O2Std', 'PolynomialCoefficients',
                          'PolynomialCoefficientsStd', 'ChiSquaredOfFit', 'RootMeanSquareErrorOfFit',
                          'AMFInitial', 'AMFInitialClear', 'AMFInitialClearStd', 
                          'AMFInitialCloudy', 'AMFInitialCloudyStd', 'AMFUnpolluted',
                          'AMFUnpollutedStd' , 'AMFUnpollutedClear', 'AMFUnpollutedClearStd',
                          'AMFUnpollutedCloudy', 'AMFUnpollutedCloudyStd', 'AMFPolluted',
                          'AMFPollutedStd', 'AMFPollutedClear', 'AMFPollutedClearStd', 
                          'AMFPollutedCloudy' , 'AMFPollutedCloudyStd' , 'AMFPollutedToGround',
                          'AMFPollutedToGroundStd', 'CloudFraction', 'CloudFractionStd', 
                          'CloudRadianceFraction', 'CloudPressure', 'CloudPressureStd',
                          'TerrainReflectivity', 'TerrainPressure', 'TerrainHeight', 
                          'SmallPixelRadiancePointer', 'InstrumentConfigurationId',
                          'MeasurementQualityFlags', 'FitQualityFlags', 'AMFQualityFlags', 
                          'WavelengthRegistrationCheck', 'WavelengthRegistrationCheckStd', 
                          'UnpolFldLatBandQualityFlags', 'vcdQualityFlags',
                          'Time', 'Latitude', 'Longitude', 'SpacecraftLatitude',
                          'SpacecraftLongitude', 'SpacecraftAltitude', 'SolarZenithAngle',
                          'SolarAzimuthAngle', 'ViewingZenithAngle', 'ViewingAzimuthAngle',
                          'GroundPixelQualityFlags']
        self.validSizes = [1644*60, 1644*60*6, 1644*60, 1644, 1644*60, 
                           1, 1644, 1644*60, 1644*60]
        self.validTypes = [numpy.float32, numpy.float32, numpy.float64, numpy.uint8,
                           numpy.uint16, numpy.uint8, numpy.double, numpy.float32, 
                           numpy.float32]
        self.firstElValues = numpy.array(self.bldFlatList([Helpers.genNan(1), Helpers.genNan(6), -32.767,
                              2, 2104, 255, 5.78275952674192*10**8,
                              -79.71093, 94.31383]))
        self.nanElInd = (1545, 39)
        self.nanElValues = numpy.array(self.bldFlatList([Helpers.genNan(1), Helpers.genNan(6), -32.767, 2,
                            2100, 255, 5.782790432626631*10**8, 71.152504, 88.35141]))
        self.lastElInd = (1643, 59)
        self.lastElValues = numpy.array(self.bldFlatList([Helpers.genNan(1), Helpers.genNan(6), -32.767, 
                                         2, 2104, 255, 5.782792393000011*10**8, 61.2573, 102.92529]))
        self.valElInd = (520, 52)
        self.valElValues = numpy.array(self.bldFlatList([4.80269451*10**14, [6.431637*10**20, -9.0338365*10**19,
                                                                             8.2326397*10**16, 7.3165924*10**18,
                                                                             8.2205615*10**18, -5.8147095*10**17],
                                                         .125, 0, 8192, 255, 5.78276992872236*10**8, -15.731218, 46.450844]))

    
    def test_instantiate(self):
        self.assertIsInstance(self.parser, parse_geo.HDFnasaomil2_File)
        
    def test_name(self):
        self.assertEqual(self.parser.name, self.fname)
        
    def test_autofill_extension(self):
        self.assertEqual(self.parser.ext, 'hdf')
        
    def test_override_extension(self):
        dir = os.path.dirname(__file__)
        obj = parse_geo.HDFnasaomil2_File(self.fname, extension='foo',
                                          cornerDir=os.path.join(dir, 'sample_data'),
                                          cornerFileList=['ominasacornersample.hdf'])
        self.assertEqual(obj.ext, 'foo')
        
    def test_override_subtype(self):
        self.assertEqual(self.parser.sub, 'knmiomil2')
        
    def test_reject_nonHDF(self):
        dir = os.path.split(self.fname)[0]
        empty = os.path.join(dir, 'empty')
        with self.assertRaises(IOError):
            parse_geo.HDFknmiomil2_File(empty)
            
    def test_get_retrieve_all_vars(self):
        for key in self.validKeys:
            var = self.parser.get(key)
            self.assertIsInstance(var, numpy.ndarray)
            
    def test_get_cm_retrieves_all_vars(self):
        with self.parser as p:
            for key in self.validKeys:
                var = p.get_cm(key)
                self.assertIsInstance(var, numpy.ndarray)
            
    def test_get_cm_cleans_up_open_file(self):
        closes_all = does_pytables_close_all(self.fname)
        fid = tables.openFile(self.fname, mode='r')
        with self.parser as p:
            unused_var = p.get_cm('Time')
        if closes_all:
            self.assertFalse(fid.isopen)
        else:
            self.assertTrue(fid.isopen)
            fid.close()
        
    def test_get_retrieves_right_sizes_of_full_vars(self):
        chkSizes = [self.parser.get(key).size for key in self.checkKeys]
        self.assertListEqual(chkSizes, self.validSizes)
        
    def test_get_cm_retrieves_right_sizes_of_full_vars(self):
        with self.parser as p:
            chkSizes = [p.get(key).size for key in self.checkKeys]
        self.assertListEqual(chkSizes, self.validSizes)
        
    def test_get_retrieves_correct_type(self):
        chkTypes = [self.parser.get(key).dtype for key in self.checkKeys]
        self.assertListEqual(chkTypes, self.validTypes)
        
    def test_get_cm_retrieves_correct_type(self):
        with self.parser as p:
            chkTypes = [p.get_cm(key).dtype for key in self.checkKeys]
        self.assertListEqual(chkTypes, self.validTypes)
        
    def test_get_retrieves_first_element(self):
        chkVals = numpy.array(self.bldFlatList([self.parser.get(key, (0,0)) for key in self.checkKeys]))
        numpy.testing.assert_array_almost_equal(chkVals.flat, self.firstElValues.flat, decimal=3)
        
    def test_get_cm_retrieves_first_element(self):
        with self.parser as p:
            chkVals = numpy.array(self.bldFlatList([p.get_cm(key, (0,0)) for key in self.checkKeys]))
        numpy.testing.assert_array_almost_equal(chkVals.flat, self.firstElValues.flat, decimal=3)
        
    def test_get_retrieves_nan_element(self):
        chkVals = numpy.array(self.bldFlatList([self.parser.get(key, self.nanElInd) for key in self.checkKeys]))
        numpy.testing.assert_array_almost_equal(chkVals.flat, self.nanElValues.flat, decimal=3)
        
    def test_get_cm_retrieves_nan_element(self):
        with self.parser as p:
            chkVals = numpy.array(self.bldFlatList([p.get_cm(key, self.nanElInd) for key in self.checkKeys]))
        numpy.testing.assert_array_almost_equal(chkVals.flat, self.nanElValues.flat, decimal=3)
        
    def test_get_retrieves_last_element(self):
        chkVals = numpy.array(self.bldFlatList([self.parser.get(key, self.lastElInd) for key in self.checkKeys]))
        numpy.testing.assert_array_almost_equal(chkVals.flat, self.lastElValues.flat, decimal=3)
        
    def test_get_cm_retrieves_last_element(self):
        with self.parser as p:
            chkVals = numpy.array(self.bldFlatList([p.get_cm(key, self.lastElInd) for key in self.checkKeys]))
        numpy.testing.assert_array_almost_equal(chkVals.flat, self.lastElValues.flat, decimal=3)
        
    def test_get_retrives_valid_element(self):
        chkVals = numpy.array(self.bldFlatList([self.parser.get(key, self.valElInd) for key in self.checkKeys]))
        numpy.testing.assert_allclose(chkVals.flat, self.valElValues.flat, rtol=10**-5)
    
    def test_get_cm_retrieves_valid_element(self):
        with self.parser as p:
            chkVals = numpy.array(self.bldFlatList([p.get_cm(key, self.valElInd) for key in self.checkKeys]))
        numpy.testing.assert_allclose(chkVals.flat, self.valElValues.flat, rtol=10**-5)
        
    def test_get_zero_order_array_from_single_point(self):
        # check both nan and non-nan, for a variable that applies scale-offset
        # as well as one that does not
        self.assertEqual(numpy.rank(self.parser.get('CloudFraction', self.valElInd)), 0)
        self.assertEqual(numpy.rank(self.parser.get('SolarZenithAngle', self.valElInd)), 0)
        self.assertEqual(numpy.rank(self.parser.get('CloudFraction', self.nanElInd)), 0)
        self.assertEqual(numpy.rank(self.parser.get('SolarZenithAngle', self.nanElInd)), 0)
        
    def test_get_cm_zero_order_array_from_single_point(self):
        with self.parser as p:
            self.assertEqual(numpy.rank(p.get_cm('CloudFraction', self.valElInd)), 0)
            self.assertEqual(numpy.rank(p.get_cm('SolarZenithAngle', self.valElInd)), 0)
            self.assertEqual(numpy.rank(p.get_cm('CloudFraction', self.nanElInd)), 0)
            self.assertEqual(numpy.rank(p.get_cm('SolarZenithAngle', self.nanElInd)), 0)

@skipUnlessSamples()    
class TestNASAOmiL2GetGeoCorners(unittest.TestCase):
    
    def setUp(self):
        dir = os.path.dirname(__file__)
        self.sample_dir = os.path.join(dir, 'sample_data')
        self.fname = os.path.join(dir, 'sample_data', 'ominasal2sample.hdf')
        self.cornerName = os.path.join(dir, 'sample_data', 'ominasacornersample.hdf')
        parser = parse_geo.HDFnasaomil2_File(self.fname, cornerDir=self.sample_dir,
                                             cornerFileList=['ominasacornersample.hdf'])
        self.geoarray = parser.get_geo_corners()
        
    def test_corners_have_correct_fields(self):
        self.assertItemsEqual(self.geoarray.dtype.names, ('lat', 'lon', 'ind'))
        
    def test_corners_size(self):
        self.assertEqual(self.geoarray.size, 1644*60)
        
    def test_corners_types(self):
        self.assertEqual(self.geoarray['lat'].dtype, numpy.float32)
        self.assertEqual(self.geoarray['lon'].dtype, numpy.float32)
        
    def test_corners_lat_first_element(self):
        knownLat = [-79.5652, -80.0409, -79.8142, -79.3238]
        latInd = utils.find_occurences(self.geoarray['ind'], [0,0])    
        fileLat = self.geoarray['lat'][latInd].squeeze()
        numpy.testing.assert_allclose(knownLat, fileLat)
        
    def test_corners_lon_first_element(self):
        knownLon = [172.91501, -179.07101, -178.61101, 173.598]
        lonInd = utils.find_occurences(self.geoarray['ind'], [0,0])
        fileLon = self.geoarray['lon'][lonInd].squeeze()
        numpy.testing.assert_allclose(knownLon, fileLon)
        
    def test_corners_lat_last_element(self):
        knownLat = [61.3584, 61.3896, 61.126198, 61.1233]
        latInd = utils.find_occurences(self.geoarray['ind'], [1643, 59])
        fileLat = self.geoarray['lat'][latInd].squeeze()
        numpy.testing.assert_allclose(knownLat, fileLat)
        
    def test_corners_lon_last_element(self):
        knownLon = [.96595997, -1.84995, -1.81073, 0.97933]
        lonInd = utils.find_occurences(self.geoarray['ind'], [1643, 59])
        fileLon = self.geoarray['lon'][lonInd].squeeze()
        numpy.testing.assert_allclose(knownLon, fileLon)
        
    def test_corners_lat_rand_element(self):
        knownLat = [20.7784, 20.8055, 20.934, 20.906399]
        latInd = utils.find_occurences(self.geoarray['ind'], [827, 36])
        fileLat = self.geoarray['lat'][latInd].squeeze()
        numpy.testing.assert_allclose(knownLat, fileLat)
        
    def test_corners_lon_rand_element(self):
        knownLon = [-165.727, -165.49, -165.51901, -165.75601]
        lonInd = utils.find_occurences(self.geoarray['ind'], [827, 36])
        fileLon = self.geoarray['lon'][lonInd].squeeze()
        numpy.testing.assert_allclose(knownLon, fileLon)

    def test_raises_IO_if_no_corner_file_valid_dir(self):
        # we should be able to instantiate the parser with an 
        # empty corner file, but not be able to get_geo_corners
        dir = os.path.dirname(__file__)
        newParser = parse_geo.HDFnasaomil2_File(self.fname, cornerDir=dir,
                                                cornerFileList=[''])
        self.assertRaises(IOError, newParser.get_geo_corners)

    def test_raises_IO_if_no_corner_file_invalid_dir(self):
        dir = ''
        newParser = parse_geo.HDFnasaomil2_File(self.fname, cornerDir=dir,
                                                cornerFileList=[''])
        
    def test_closes_cornerfile_when_done(self):
        # The behavior of pytables is consistent within each
        # version but varies from version to version.  Therefore
        # we test above ane ensure we are getting the desired
        # behavior
        closes_all = does_pytables_close_all(self.fname)
        fid = tables.openFile(self.cornerName, 'r')
        parser = parse_geo.HDFnasaomil2_File(self.fname, cornerDir=self.sample_dir,
                                             cornerFileList=['ominasacornersample.hdf'])
        unused_geoarray = parser.get_geo_corners()
        if closes_all:
            self.assertFalse(fid.isopen)
        else:
            self.assertTrue(fid.isopen)
            fid.close()
        
        
@skipUnlessSamples()
class TestNASAOmiL2GetGeoCenters(TestNASAOmiL2Parser):
    
    
    def setUp(self):
        TestNASAOmiL2Parser.setUp(self)
        self.geoarray = self.parser.get_geo_centers()
        
    def test_centers_have_correct_fields(self):
        self.assertItemsEqual(self.geoarray.dtype.names, ('lat', 'lon', 'ind'))
        
    def test_centers_size(self):
        self.assertEqual(self.geoarray.size, 1644*60)
        
    def test_centers_types(self):
        self.assertEqual(self.geoarray['lat'].dtype, numpy.float32)
        self.assertEqual(self.geoarray['lon'].dtype, numpy.float32)
        
    def test_centers_lat_first_element(self):
        knownLat = -79.71093
        latInd = utils.find_occurences(self.geoarray['ind'], [0, 0])
        fileLat = self.geoarray['lat'][latInd].squeeze()
        numpy.testing.assert_allclose(knownLat, fileLat)
        
    def test_centers_lon_first_element(self):
        knownLon = 177.11975
        lonInd = utils.find_occurences(self.geoarray['ind'], [0, 0])
        fileLon = self.geoarray['lon'][lonInd].squeeze()
        numpy.testing.assert_allclose(knownLon, fileLon)
        
    def test_centers_lat_last_element(self):
        knownLat = 61.2573
        latInd = utils.find_occurences(self.geoarray['ind'], [1643, 59])
        fileLat = self.geoarray['lat'][latInd].squeeze()
        numpy.testing.assert_allclose(knownLat, fileLat)
        
    def test_centers_lon_last_element(self):
        knownLon = -.42845777
        lonInd = utils.find_occurences(self.geoarray['ind'], [1643, 59])
        fileLon = self.geoarray['lon'][lonInd].squeeze()
        numpy.testing.assert_allclose(knownLon, fileLon)
        
    def test_centers_lat_rand_element(self):
        knownLat = 2.7522528
        latInd = utils.find_occurences(self.geoarray['ind'], [678, 33])
        fileLat = self.geoarray['lat'][latInd].squeeze()
        numpy.testing.assert_allclose(knownLat, fileLat)
        
    def test_centers_lon_rand_element(self):
        knownLon = -162.4073
        lonInd = utils.find_occurences(self.geoarray['ind'], [678, 33])
        fileLon = self.geoarray['lon'][lonInd].squeeze()
        numpy.testing.assert_allclose(knownLon, fileLon)
        
    def test_geocenters_can_feed_ind_to_get(self):
        latAxis = self.geoarray.dtype.names.index('lat')
        lonAxis = self.geoarray.dtype.names.index('lon')
        indAxis = self.geoarray.dtype.names.index('ind')
        flatarray = self.geoarray.flatten()
        for i in range(0,flatarray.size, 5000):
            row = flatarray[i]
            numpy.testing.assert_array_equal(row[latAxis], self.parser.get('Latitude', row[indAxis]))
            numpy.testing.assert_array_equal(row[lonAxis], self.parser.get('Longitude', row[indAxis]))
            
    def test_geocenters_can_feed_ind_to_get_cm(self): 
        latAxis = self.geoarray.dtype.names.index('lat')
        lonAxis = self.geoarray.dtype.names.index('lon')
        indAxis = self.geoarray.dtype.names.index('ind')
        flatarray = self.geoarray.flatten()
        with self.parser as p:
            for i in range(0, flatarray.size, 500):
                row = flatarray[i]
                numpy.testing.assert_array_equal(row[latAxis], p.get_cm('Latitude', row[indAxis]))
                numpy.testing.assert_array_equal(row[lonAxis], p.get_cm('Longitude', row[indAxis]))  

@skipUnlessSamples()                
class TestMOPPITl2Parser(TestParser):    
    
    def setUp(self):
        dir = os.path.dirname(__file__)
        self.fname = os.path.join(dir, 'sample_data', 'moppitl2sample.hdf')
        self.parser= parse_geo.HDFmopittl2_File(self.fname, subtype='HDFmoppittl2')
        self.checkKeys = ['Time', 'Pressure Grid', 
                          'Retrieved CO Surface Mixing Ratio',
                          'A Priori CO Mixing Ratio Profile',
                          'Surface Index']
        self.validKeys = ['Time', 'Latitude', 'Longitude', 'Seconds in Day',
                          'Pressure Grid', 'Solar Zenith Angle',
                          'Satellite Zenith Angle', 'A Priori Surface Emissivity',
                          'Surface Pressure', 'Retrieved Surface Temperature',
                          'Retrieved Surface Emissivity', 'Surface Index',
                          'Cloud Description', 'Retrieved CO Total Column',
                          'Retrieved CO Mixing Ratio Profile',
                          'Retrieved CO Surface Mixing Ratio',
                          'A Priori Surface Temperature', 
                          'A Priori CO Mixing Ratio Profile',
                          'A Priori CO Surface Mixing Ratio',
                          'Retrieved CO Total Column Diagnostics',
                          'Retrieval Averaging Kernel Matrix', 
                          'Degrees of Freedom for Signal',
                          'Water Vapor Climatology Content', 'Retrieval Iterations',
                          'Retrieval Error Covariance Matrix', 
                          'Level 1 Radiances and Errors', 'DEM Altitude',
                          'MODIS Cloud Diagnostics', 'Signal Chi2', 'Swath Index']
        self.validSizes = [200032, 9, 200032*2, 200032*2*9, 200032]
        self.validTypes = [numpy.float64, numpy.float32, numpy.float32,
                           numpy.float32, numpy.int32]
        presGrid = [900.0, 800.0, 700.0, 600.0, 500.0, 400.0, 300.0, 200.0, 100.0]
        fElProf = self.bldFlatList([[82.91398, 80.92176, 75.33437, 73.38161, 74.29397,
                                     75.11533, 75.61497, 70.074974, 55.287525], 
                                    [25.30945, 24.276491, 22.600277, 22.014448, 22.288157,
                                     22.534565, 22.684454, 21.022459, 16.586231]])
        fElVals = [3.95020813831*10**8, presGrid, [84.239365, 29.985426], fElProf, 0]
        self.firstElVals = self.bldFlatArray(fElVals)
        self.firstElInd = (0,)
        lElProf = self.bldFlatList([[58.33968, 58.10826, 58.364857, 58.7423, 59.134647,
                                     59.697765, 60.06516, 59.69639, 56.674343],
                                    [17.808157, 17.432451, 17.50943, 17.622662, 17.740366,
                                     17.909302, 18.01952, 17.908888, 17.002277]])
        lElVals = [3.95107203890999*10**8, presGrid, [58.138416, 20.491892], lElProf, 0]
        self.lastElVals = self.bldFlatArray(lElVals);
        self.lastElInd = (200031,)
        nElProf = self.bldFlatList([[numpy.nan, numpy.nan, numpy.nan, 61.714626, 62.36158,
                                     62.36822, 60.415375, 51.314304, 26.068226],
                                    [numpy.nan, numpy.nan, numpy.nan, 17.72849, 18.708445,
                                     18.710438, 18.124584, 15.394268, 7.8204556]])
        nElVals = [3.9510279793499964*10**8, presGrid, [62.32344, 20.00358], nElProf, 1]
        self.nanElVals = self.bldFlatArray(nElVals)
        self.nanElInd = (191000,)
        vElProf = self.bldFlatList([[79.70813, 78.30978, 74.31502, 73.186646, 74.27929,
                                     75.443886, 76.37308, 70.705154, 55.171642],
                                    [24.330866, 23.492897, 22.29447, 21.955961, 22.283752,
                                     22.63313, 22.911888, 21.211514, 16.551466]])
        vElVals = [3.950208313689997*10**8, presGrid, [72.061165, 25.169048], vElProf, 0]
        self.valElVals = self.bldFlatArray(vElVals)
        self.valElInd = (50,)
        self.retKeys = ['Time', 'Pressure Grid', 
                        'Retrieved CO Surface Mixing Ratio',
                        'A Priori CO Mixing Ratio Profile',
                        'A Priori CO Mixing Ratio Profile',
                        'Surface Index'] # need to get 2 slices from profile
        self.retInd = [slice(None), slice(None), slice(None), (slice(None),0),
                       (slice(None),1), slice(None)] # how to cut up keys
        
    def retrieve_values_get(self, fileInd, parser):
        '''This function retrieves a flat array from a file parser with get'''
        chkList = [parser.get(k, fileInd)[i] for (k,i) in izip(self.retKeys, self.retInd)]
        return self.bldFlatArray(chkList)
    
    def retrieve_values_get_cm(self, fileInd, parser):
        '''This function retrieves a flat array from a cm file parser with get_cm'''
        chkList = [parser.get_cm(k, fileInd)[i] for (k,i) in izip(self.retKeys, self.retInd)]
        return self.bldFlatArray(chkList)

    def test_instantiate(self):
        self.assertIsInstance(self.parser, parse_geo.HDFmopittl2_File)
        
    def test_name(self):
        self.assertEqual(self.parser.name, self.fname)
        
    def test_autofill_extension(self):
        self.assertEqual(self.parser.ext, 'hdf')
        
    def test_override_extension(self):
        obj = parse_geo.HDFmopittl2_File(self.fname, extension='bananas')
        self.assertEqual(obj.ext, 'bananas')
        
    def test_override_subtype(self):
        self.assertEqual(self.parser.sub, 'HDFmoppittl2')
        
    def test_reject_nonHDF(self):
        dir = os.path.split(self.fname)[0]
        emptyFile = os.path.join(dir, 'empty')
        with self.assertRaises(IOError):
            parse_geo.HDFmopittl2_File(emptyFile)
            
    def test_get_retrieves_all_vars(self):
        for key in self.validKeys:
            var = self.parser.get(key)
            self.assertIsInstance(var, numpy.ndarray)
            
    def test_get_cm_retrieves_all_vars(self):
        with self.parser as p:
            for key in self.validKeys:
                var = p.get_cm(key)
                self.assertIsInstance(var, numpy.ndarray)
            
    @unittest.skip("We should check this but we don't have the technology")    
    def test_get_cm_cleans_up_open_file(self):
        pass
        
    def test_get_retrieves_correct_size(self):
        sizes = [self.parser.get(key).size for key in self.checkKeys]
        self.assertListEqual(sizes, self.validSizes)
        
    def test_get_cm_retrieves_correct_size(self):
        with self.parser as p:
            sizes = [p.get_cm(key).size for key in self.checkKeys]
        self.assertListEqual(sizes, self.validSizes)

    @unittest.expectedFailure
    # we're upcasting everything.  Not ideal, but it works (as 
    # our floats are still floats and our ints still ints)
    def test_get_retrieves_correct_type(self):
        types = [self.parser.get(key).dtype for key in self.checkKeys]
        self.assertListEqual(types, self.validTypes)
        
    @unittest.expectedFailure
    # we're upcasting everything.  Not ideal, but it works (as 
    # our floats are still floats and our ints still ints)
    def test_get_cm_retrieves_correct_type(self):
        with self.parser as p:
            types = [p.get_cm(key).dtype for key in self.checkKeys]
        self.assertListEqual(types, self.validTypes)
        
    def test_get_retrieves_first_element(self):
        chkVals = self.retrieve_values_get(self.firstElInd, self.parser)
        numpy.testing.assert_array_almost_equal(chkVals, self.firstElVals, decimal=5)
        
    def test_get_cm_retrieves_first_element(self):
        with self.parser as p:
            chkVals = self.retrieve_values_get_cm(self.firstElInd, p)
        numpy.testing.assert_array_almost_equal(chkVals, self.firstElVals, decimal=5)
        
    def test_get_retrieves_last_element(self):
        chkVals = self.retrieve_values_get(self.lastElInd, self.parser)
        numpy.testing.assert_array_almost_equal(chkVals, self.lastElVals, decimal=5)
        
    def test_get_cm_retrieves_last_element(self):
        with self.parser as p:
            chkVals = self.retrieve_values_get_cm(self.lastElInd, p)
        numpy.testing.assert_array_almost_equal(chkVals, self.lastElVals, decimal=5)
        
    def test_get_retrieves_nan_element(self):
        chkVals = self.retrieve_values_get(self.nanElInd, self.parser)
        numpy.testing.assert_array_almost_equal(chkVals, self.nanElVals, decimal=5)
        
    def test_get_cm_retrieves_nan_element(self):
        with self.parser as p:
            chkVals = self.retrieve_values_get_cm(self.nanElInd, p)
        numpy.testing.assert_array_almost_equal(chkVals, self.nanElVals, decimal=5)
        
    def test_get_retrieves_valid_element(self):
        chkVals = self.retrieve_values_get(self.valElInd, self.parser)
        numpy.testing.assert_array_almost_equal(chkVals, self.valElVals, decimal=5)
        
    def test_get_cm_retrieves_valid_element(self):
        with self.parser as p:
            chkVals = self.retrieve_values_get_cm(self.valElInd, p)
        numpy.testing.assert_array_almost_equal(chkVals, self.valElVals, decimal=5)
        
@skipUnlessSamples()
class TestMOPPITL2GetGeoCenters(unittest.TestCase):
    
    def setUp(self):
        dir = os.path.dirname(__file__)
        fname = os.path.join(dir, 'sample_data', 'moppitl2sample.hdf')
        self.parser= parse_geo.HDFmopittl2_File(fname, subtype='HDFmoppittl2')
        self.geoarray = self.parser.get_geo_centers()
        
    def test_centers_have_correct_fields(self):
        self.assertItemsEqual(self.geoarray.dtype.names, ('lat', 'lon', 'ind'))
        
    def test_centers_size(self):
        self.assertEqual(self.geoarray.size, 200032)
        
    @unittest.expectedFailure
    # we're upcasting everything.  Not ideal, but it works (as 
    # our floats are still floats and our ints still ints)
    def test_centers_types(self):
        self.assertEqual(self.geoarray['lat'].dtype, numpy.float32)
        self.assertEqual(self.geoarray['lon'].dtype, numpy.float32)
        
    def test_centers_lat_first_element(self):
        knownLat = -12.894982
        latInd = 0
        fileLat = self.geoarray['lat'][latInd].squeeze()
        numpy.testing.assert_allclose(knownLat, fileLat)
        
    def test_centers_lon_first_element(self):
        knownLon = -19.320307
        lonInd = 0
        fileLon = self.geoarray['lon'][lonInd].squeeze()
        numpy.testing.assert_allclose(knownLon, fileLon)
        
    def test_centers_lat_last_element(self):
        knownLat = -8.817091
        latInd = 200031
        fileLat = self.geoarray['lat'][latInd].squeeze()
        numpy.testing.assert_allclose(knownLat, fileLat)
        
    def test_centers_lon_last_element(self):
        knownLon = 158.08238
        lonInd = 200031
        fileLon = self.geoarray['lon'][lonInd].squeeze()
        numpy.testing.assert_allclose(knownLon, fileLon)
        
    def test_centers_lat_rand_element(self):
        knownLat = -77.91937
        latInd = 40000
        fileLat = self.geoarray['lat'][latInd].squeeze()
        numpy.testing.assert_allclose(knownLat, fileLat)
        
    def test_centers_lon_rand_element(self):
        knownLon = 34.588852
        lonInd = 40000
        fileLon = self.geoarray['lon'][lonInd].squeeze()
        numpy.testing.assert_allclose(knownLon, fileLon)
        
    def test_geocenters_can_feed_ind_to_get(self):
        latAxis = self.geoarray.dtype.names.index('lat')
        lonAxis = self.geoarray.dtype.names.index('lon')
        indAxis = self.geoarray.dtype.names.index('ind')
        tupleArray = self.geoarray.flatten()
        for i in range(0, tupleArray.size, 50000):
            row = tupleArray[i]
            numpy.testing.assert_array_equal(row[latAxis], self.parser.get('Latitude', row[indAxis]))
            numpy.testing.assert_array_equal(row[lonAxis], self.parser.get('Longitude', row[indAxis]))
            
    def test_geocenters_can_feed_ind_to_get_cm(self):
        latAxis = self.geoarray.dtype.names.index('lat')
        lonAxis = self.geoarray.dtype.names.index('lon')
        indAxis = self.geoarray.dtype.names.index('ind')
        tupleArray = self.geoarray.flatten()
        with self.parser as p:
            for i in range(0, tupleArray.size, 5000):
                row = tupleArray[i]
                numpy.testing.assert_array_equal(row[latAxis], p.get_cm('Latitude', row[indAxis]))
                numpy.testing.assert_array_equal(row[lonAxis], p.get_cm('Longitude', row[indAxis]))
            
        
class TestOverallGridGeo(unittest.TestCase):
    
    
    def setUp(self):
        # fakeparms should have every possible keyword and
        # work with every implemented class
        self.fakeParms = {'stdPar1'     : 45,
                          'stdPar2'     : 50,
                          'refLat'      : 47.5,
                          'refLon'      : -90,
                          'earthRadius' : 6370000,
                          'xOrig'       : -10000,
                          'yOrig'       : -50000,
                          'xCell'       : 1000,
                          'yCell'       : 2000,
                          'nRows'       : 150,
                          'nCols'       : 250}
        self.projNames = grid_geo.ValidProjections()
        self.origLat = self.fakeParms['refLat']
        self.origLon = self.fakeParms['refLon']
        
    def test_existence_of_ValidProjections(self):
        for name in self.projNames:
            self.assertTrue(hasattr(grid_geo, name+'_GridDef'))
    
    def test_inheritance_of_ValidProjections(self):
        allClasses = [getattr(grid_geo, name+'_GridDef') for name in self.projNames]
        for cls in allClasses:
            self.assertIsInstance(cls(self.fakeParms), grid_geo.GridDef)
    
    def test_all_classes_implement_requiredParms(self):
        allClasses = [getattr(grid_geo, name+'_GridDef') for name in self.projNames]
        for cls in allClasses:
            for parmName in cls.requiredParms():
                self.assertTrue(parmName in self.fakeParms.keys())

    def test_all_classes_implement_indLims(self):
        instClasses = [getattr(grid_geo, name+'_GridDef')(self.fakeParms) for name in self.projNames]
        for inst in instClasses:
            unused_limsTup = inst.indLims()
            
    def test_all_indLims_right_size(self):
        instClasses = [getattr(grid_geo, name+'_GridDef')(self.fakeParms) for name in self.projNames]
        for inst in instClasses:
            self.assertEqual(len(inst.indLims()), 4)
            
    def test_all_indLims_minRow_lt_maxRow(self):
        instClasses = [getattr(grid_geo, name+'_GridDef')(self.fakeParms) for name in self.projNames]
        for inst in instClasses:     
            limsTup = inst.indLims()   
            self.assertLessEqual(limsTup[0], limsTup[1])
            
    def test_all_indLims_minCol_lt_maxCol(self):
        instClasses = [getattr(grid_geo, name+'_GridDef')(self.fakeParms) for name in self.projNames]
        for inst in instClasses: 
            limsTup = inst.indLims()
            self.assertLessEqual(limsTup[2], limsTup[3])
            
    def test_all_classes_implement_geoToProjected(self):
        instClasses = [getattr(grid_geo, name+'_GridDef')(self.fakeParms) for name in self.projNames]
        for inst in instClasses:
            if not hasattr(inst, "geoToProjected"):
                raise AssertionError

    def test_all_classes_implement_geoToGridded(self):
        instClasses = [getattr(grid_geo, name+'_GridDef')(self.fakeParms) for name in self.projNames]
        for inst in instClasses:
            if not hasattr(inst, "geoToGridded"):
                raise AssertionError

    def test_all_classes_implement_projectedToGeo(self):
        instClasses = [getattr(grid_geo, name+'_GridDef')(self.fakeParms) for name in self.projNames]
        for inst in instClasses:
            if not hasattr(inst, "projectedToGeo"):
                raise AssertionError

    def test_all_classes_implement_griddedToGeo(self):
        instClasses = [getattr(grid_geo, name+'_GridDef')(self.fakeParms) for name in self.projNames]
        for inst in instClasses:
            if not hasattr(inst, "griddedToGeo"):
                raise AssertionError

    def test_to_from_projected_are_inverses(self):
        instClasses = [getattr(grid_geo, name+'_GridDef')(self.fakeParms) for name in self.projNames]        
        for inst in instClasses:
            (y, x) = inst.geoToProjected(self.origLat, self.origLon)
            (lat, lon) = inst.projectedToGeo(y, x)
            self.assertEqual((lat, lon), (self.origLat, self.origLon))

    def test_to_from_gridded_are_inverses(self):
        instClasses = [getattr(grid_geo, name+'_GridDef')(self.fakeParms) for name in self.projNames]        
        for inst in instClasses:
            (row, col) = inst.geoToGridded(self.origLat, self.origLon)
            (lat, lon) = inst.griddedToGeo(row, col)
            self.assertEqual((lat, lon), (self.origLat, self.origLon))

    @unittest.expectedFailure
    # this will fail for lat/lon grid because of how it's defined
    def test_origin_is_0_0(self):
        instClasses = [getattr(grid_geo, name+'_GridDef')(self.fakeParms) for name in self.projNames]                
        for inst in instClasses:
            (y, x) = inst.geoToProjected(self.origLat, self.origLon)
            self.assertEqual((y,x), (0,0))

class Testlcc2par_GridDef(TestOverallGridGeo):
    
    
    def setUp(self):
        TestOverallGridGeo.setUp(self)
        self.instance = grid_geo.lcc2par_GridDef(self.fakeParms)
        self.knownGeoCoords = [ (18,      -50),
                                (37.24,  78.9),
                                (47.5,  -87.9),
                                (10,       10),
                                (18,      310),
                                (47.5,    -90) ]
        self.knownProjCoords = [ (4551772.35960682, -2215457.75967858),
                                 (5743841.64228215,  9786960.03533386),
                                 ( 157562.52186212,     2129.68262648),
                                 (9868815.52687009,  2953339.97045807),
                                 (4551772.35960682, -2215457.75967858),
                                 (      0.00000000,        0.00000000) ]
        knownGridCoords = [ (-1081.72887984,  4562.77235961),
                                 ( 4919.48001767,  5754.84164228),
                                 (   27.06484131,   168.56252186),
                                 ( 1502.66998523,  9879.81552687),
                                 (-1081.72887984,  4562.77235961), 
                                 (   26.00000000,    11.00000000) ]
        # subtract one from all grid coords to account for matlab's
        # one based indexing
        self.knownGridCoords = []
        for item in knownGridCoords:
            self.knownGridCoords.append((item[0]-1,item[1]-1))

    def test_accurate_indLims(self):
        (minRow, maxRow, minCol, maxCol) = self.instance.indLims()
        self.assertEqual(maxRow-minRow+1, self.fakeParms['nRows'])
        self.assertEqual(maxCol-minCol+1, self.fakeParms['nCols'])

    def test_geo2Proj_vs_matlab(self):
        for ((lat, lon), (x, y)) in izip(self.knownGeoCoords, self.knownProjCoords):
            (yCalc, xCalc) = self.instance.geoToProjected(lat, lon)
            # use some logic to make sure we don't have way to tight of 
            # tolerances on identically zero points
            if y != 0:
                self.assertAlmostEqual(y, yCalc, delta=abs(y*1E-6))
            else:
                self.assertAlmostEqual(y, yCalc, places=6)
            if x != 0:
                self.assertAlmostEqual(x, xCalc, delta=abs(x*1E-6))
            else:
                self.assertAlmostEqual(x, xCalc, places=6)

    def test_geo2Grid_vs_matlab(self):
        for ((lat, lon), (row, col)) in izip(self.knownGeoCoords, self.knownGridCoords):
            (rowCalc, colCalc) = self.instance.geoToGridded(lat, lon)
            self.assertAlmostEqual(row, rowCalc)
            self.assertAlmostEqual(col, colCalc)

    def test_proj2Geo_vs_matlab(self):
        for ((x, y), (lat, lon)) in izip(self.knownProjCoords, self.knownGeoCoords):
            (latCalc, lonCalc) = self.instance.projectedToGeo(y, x)
            # wrap to 0-360 standard for comparison
            lonCalc = utils.wrap_lon_0_360(lonCalc)
            lon = utils.wrap_lon_0_360(lon)
            self.assertAlmostEqual(lat, latCalc, delta=abs(lat*1E-6))
            self.assertAlmostEqual(lon, lonCalc, delta=abs(lon*1E-6))

    def test_griddedToGeo_vs_matlab(self):
        for ((row, col), (lat, lon)) in izip(self.knownGridCoords, self.knownGeoCoords):
            (latCalc, lonCalc) = self.instance.griddedToGeo(row, col)
            # wrap to 0-360 standard for comparison
            lonCalc = utils.wrap_lon_0_360(lonCalc)
            lon = utils.wrap_lon_0_360(lon)
            self.assertAlmostEqual(lat, latCalc, delta=abs(lat*1E-5))
            self.assertAlmostEqual(lon, lonCalc, delta=abs(lon*1E-5))

    def test_geoToProj_vs_JPsnyder_example(self):
        JPparms = { 'stdPar1' : 33,
                    'stdPar2' : 45,
                    'refLat'  : 23,
                    'refLon'  : -96,
                    'earthRadius' : 1.0,
                    'xOrig'   : 0,
                    'yOrig'   : 0,
                    'xCell'   : 1,
                    'yCell'   : 1,
                    'nRows'   : 1,
                    'nCols'   : 1 }
        instance = grid_geo.lcc2par_GridDef(JPparms)
        (lat, lon) = (35, -75)
        (y, x) = instance.geoToProjected(lat, lon)
        self.assertAlmostEqual(y, .2462112, places=7)
        self.assertAlmostEqual(x, .2966785, places=7)
        
    def test_projToGeo_vs_JPsnyder_example(self):
        JPparms = { 'stdPar1' : 33,
                    'stdPar2' : 45,
                    'refLat'  : 23,
                    'refLon'  : -96,
                    'earthRadius' : 1.0,
                    'xOrig'   : 0,
                    'yOrig'   : 0,
                    'xCell'   : 1,
                    'yCell'   : 1,
                    'nRows'   : 1,
                    'nCols'   : 1 }
        instance = grid_geo.lcc2par_GridDef(JPparms)
        (y, x) = (.2462112, .2966785)
        (lat, lon) = instance.projectedToGeo(y, x)
        self.assertAlmostEqual(lat, 34.9999974, places=5)
        self.assertAlmostEqual(lon, -74.9999981, places=5)


class Testlatlon_GridDef(unittest.TestCase):
    
    
    def setUp(self):
        self.parms = {'xOrig' : -10,
                 'yOrig' : 10,
                 'xCell' : 2.5,
                 'yCell' : 2,
                 'nRows' : 30,
                 'nCols' : 80}
        self.instance = grid_geo.latlon_GridDef(self.parms)
        self.knownGeoCoords = [ (18., 0.),
                                (20., 185.), 
                                (21., -17.5),
                                (10., 1.25),
                                (-20., -30.) ]
        self.knownProjCoords = self.knownGeoCoords
        self.knownGridCoords = [ (4., 4.),
                                 (5., 78.),
                                 (5.5, -3.),
                                 (0., 4.5),
                                 (-15., -8.) ]

    def test_accurate_indLims(self):
        (minRow, maxRow, minCol, maxCol) = self.instance.indLims()
        self.assertEqual(maxRow-minRow+1, self.parms['nRows'])
        self.assertEqual(maxCol-minCol+1, self.parms['nCols'])

    def test_geo2Proj(self):
        for ((lat, lon), (y, x)) in izip(self.knownGeoCoords, self.knownProjCoords):
            (yCalc, xCalc) = self.instance.geoToProjected(lat, lon)
            self.assertEqual((y, x), (yCalc, xCalc))

    def test_geo2Grid(self):
        for ((lat, lon), (row, col)) in izip(self.knownGeoCoords, self.knownGridCoords):
            (rowCalc, colCalc) = self.instance.geoToGridded(lat, lon)
            self.assertEqual((row, col), (rowCalc, colCalc))

    def test_proj2Geo(self):
        for ((y,x), (lat,lon)) in izip(self.knownProjCoords, self.knownGeoCoords):
            (latCalc, lonCalc) = self.instance.projectedToGeo(y, x)
            self.assertEqual((lat, lon), (latCalc, lonCalc))

    def test_grid2Geo(self):
        for ((row, col), (lat, lon)) in izip(self.knownGridCoords, self.knownGeoCoords):
            (latCalc, lonCalc) = self.instance.griddedToGeo(row, col)
            self.assertEqual((lat, lon), (latCalc, lonCalc))


class TestMapGeo(unittest.TestCase):
    

    def setUp(self):
        '''This setup is meant to for child classes, not this class'''
        parms = { 'xOrig' : 0,
                  'yOrig' : 0,
                  'xCell' : 1,
                  'yCell' : 1,
                  'nRows' : 5,
                  'nCols' : 10 }
        self.grid = grid_geo.latlon_GridDef(parms)
        self.parser = fakeParser('foo.dat')
        # create a dictionary with appropriate keys,
        # initialized to empty lists save for parser
        self.testMapDict = {'parser' : self.parser}
        for indTup in product(range(parms['nRows']), range(parms['nCols'])):
            self.testMapDict[indTup] = []

    def test_ValidMaps_are_valid(self):
        maps = map_geo.ValidMaps()
        for name in maps:
            if not hasattr(map_geo, name+'_map_geo'):
                raise AssertionError
            
    def assertMapsEqual(self, mapper):
        '''
        Function to compare the output of the the mapper function to the test
        dictionary.  Assumes that the parser and test dictionary have been
        prepped and are stored in the variables laid out in setUp()
        '''
        mapDict = mapper(self.parser, self.grid)
        self.assertEqual(mapDict, self.testMapDict)
    

class Test_regional_intersect(TestMapGeo):
    

    def setUp(self):
        TestMapGeo.setUp(self)
        self.mapper = getattr(map_geo, 'regional_intersect_map_geo')

    @unittest.skip("Impossible to have empty record array")
    def test_map_empty_parser(self):
        lat = numpy.array([])
        lon = numpy.array([])
        ind = numpy.arange(lat.shape[0])
        self.parser.prime_corners(lat, lon, ind)
        mapDict = self.mapper(self.parser, self.grid)
        self.assertEqual(mapDict, self.testMapDict)

    def test_poly_olap_interior_cell_edges(self):
        lat = numpy.array([[2.5, 3.5, 3.5, 2.5]])
        lon = numpy.array([[1.5, 1.5, 2.5, 2.5]])
        ind = numpy.arange(lat.shape[0])
        self.parser.prime_corners(lat, lon, ind)
        mapDict = self.mapper(self.parser, self.grid)
        listEntry = [((0,),None)]
        keys = [(2,1), (3,1), (3,2), (2,2)]
        for k in keys:
            self.testMapDict[k] = listEntry
        self.assertEqual(mapDict, self.testMapDict)

    def test_mult_poly_olap_interior_cell_edges(self):
        lat = numpy.array([[2.5, 3.5, 3.5, 2.5],
                           [2.75, 3.25, 3.25, 2.75],
                           [1.5, 2, 1.9, 1.5]])
        lon = numpy.array([[1.5, 1.5, 2.5, 2.5],
                           [1.25, 1.25, 1.75, 1.75],
                           [0.8, 0.8, 1.1, 1.2]])
        ind = numpy.array([[i] for i in range(lat.shape[0])])
        self.parser.prime_corners(lat, lon, ind)
        mapDict = self.mapper(self.parser, self.grid)
        self.testMapDict[(1,0)] = [((2,),None)]
        self.testMapDict[(1,1)] = [((2,),None)]
        self.testMapDict[(2,1)] = [((0,),None),((1,),None)]
        self.testMapDict[(3,1)] = [((0,),None),((1,),None)]
        self.testMapDict[(3,2)] = [((0,),None)]
        self.testMapDict[(2,2)] = [((0,),None)]
        self.assertEqual(mapDict, self.testMapDict)

    def test_poly_inside_cell(self):
        lat = numpy.array([[1.25, 1.75, 1.75, 1.25]])
        lon = numpy.array([[2.25, 2.25, 2.75, 2.75]])
        ind = numpy.array([[i] for i in range(lat.shape[0])])
        self.parser.prime_corners(lat, lon, ind)
        mapDict = self.mapper(self.parser, self.grid)
        self.testMapDict[(1,2)] = [((0,),None)]
        self.assertEqual(mapDict, self.testMapDict)

    def test_mult_poly_inside_cell(self):
        lat = numpy.array([[1.25, 1.75, 1.75, 1.25],
                           [1.20, 1.80, 1.80, 1.20]])
        lon = numpy.array([[2.25, 2.25, 2.75, 2.75],
                           [2.20, 2.20, 2.80, 2.80]])
        ind = numpy.array([[i] for i in range(lat.shape[0])])
        self.parser.prime_corners(lat, lon, ind)
        mapDict = self.mapper(self.parser, self.grid)
        self.testMapDict[(1,2)] = [((0,),None), ((1,),None)]
        self.assertEqual(mapDict, self.testMapDict)

    def test_poly_surrounds_cell(self):
        lat = numpy.array([[1.5, 3.5, 3.5, 1.5]])
        lon = numpy.array([[0.5, 0.5, 2.5, 2.5]])
        ind = numpy.array([[i] for i in range(lat.shape[0])])
        self.parser.prime_corners(lat, lon, ind)
        mapDict = self.mapper(self.parser, self.grid)
        listEntry = [((0,),None)]
        keys = [(1,0), (1,1), (1,2), (2,0), (2,1), (2,2), (3,0), (3,1), (3,2)]
        for k in keys:
            self.testMapDict[k] = listEntry
        self.assertEqual(mapDict, self.testMapDict)

    def test_poly_olap_left_edge(self):
        lat = numpy.array([[2.5, 3.5, 3.5, 2.5]])
        lon = numpy.array([[-.5, -.5, .5, .5]])
        ind = numpy.array([[i] for i in range(lat.shape[0])])
        self.parser.prime_corners(lat, lon, ind)
        mapDict = self.mapper(self.parser, self.grid)
        self.testMapDict[(2,0)] = [((0,),None)]
        self.testMapDict[(3,0)] = [((0,),None)]
        self.assertEqual(mapDict, self.testMapDict)

    def test_poly_olap_top_edge(self):
        lat = numpy.array([[4.5, 5.5, 5.5, 4.5]])
        lon = numpy.array([[.25, .25, .75, .75]])
        ind = numpy.array([[i] for i in range(lat.shape[0])])
        self.parser.prime_corners(lat, lon, ind)
        mapDict = self.mapper(self.parser, self.grid)
        self.testMapDict[(4,0)] = [((0,),None)]
        self.assertEqual(mapDict, self.testMapDict)

    def test_poly_olap_right_edge(self):
        lat = numpy.array([[1.25, 1.75, 1.75, 1.25]])
        lon = numpy.array([[9.5, 9.5, 10.5, 10.5]])
        ind = numpy.array([[i] for i in range(lat.shape[0])])
        self.parser.prime_corners(lat, lon, ind)
        mapDict = self.mapper(self.parser, self.grid)
        self.testMapDict[(1,9)] = [((0,),None)]
        self.assertEqual(mapDict, self.testMapDict)

    def test_poly_olap_bot_edge(self):
        lat = numpy.array([[-0.5, 1.5, 1.5, -0.5]])
        lon = numpy.array([[0.5, 0.5, 1.5, 1.5]])
        ind = numpy.array([[i] for i in range(lat.shape[0])])
        self.parser.prime_corners(lat, lon, ind)
        mapDict = self.mapper(self.parser, self.grid)
        self.testMapDict[(0,0)] = [((0,),None)]
        self.testMapDict[(0,1)] = [((0,),None)]
        self.testMapDict[(1,0)] = [((0,),None)]
        self.testMapDict[(1,1)] = [((0,),None)]
        self.assertEqual(mapDict, self.testMapDict)

    def test_poly_on_corner_case_vertex_inside(self):
        lat = numpy.array([[4.5, 5.5, 5.5, 4.5]])
        lon = numpy.array([[-0.5, -0.5, 0.5, 0.5]])
        ind = numpy.array([[i] for i in range(lat.shape[0])])
        self.parser.prime_corners(lat, lon, ind)
        mapDict = self.mapper(self.parser, self.grid)
        self.testMapDict[(4,0)] = [((0,),None)]
        self.assertEqual(mapDict, self.testMapDict)

    def test_poly_olaps_cell_without_vertex(self):
        # bascially, tests if the algorithm catches a polygon
        # that overlaps a cell but does not have a vertex in it
        lat = numpy.array([[2.75, 3.5, 3.25, 2.5]])
        lon = numpy.array([[0.5, 1.25, 1.5, 0.75]])
        ind = numpy.array([[i] for i in range(lat.shape[0])])
        self.parser.prime_corners(lat, lon, ind)
        mapDict = self.mapper(self.parser, self.grid)
        listEntry = [((0,),None)]
        keys = [(2,0),(3,0),(3,1), (2,1)]
        for k in keys:
            self.testMapDict[k] = listEntry
        self.assertEqual(mapDict, self.testMapDict)

    @unittest.expectedFailure
    def test_poly_on_corner_case_vertices_outside(self):
        # we expect this to fail because this particular mapper
        # tests that at least one vertex is inside to screen
        # out spurious pixels from wrapping.
        lat = numpy.array([[.25, .5, -.25, -.5]])
        lon = numpy.array([[-.5, -.25, .5, .25]])
        ind = numpy.array([[i] for i in range(lat.shape[0])])
        self.parser.prime_corners(lat, lon, ind)
        mapDict = self.mapper(self.parser, self.grid)
        self.testMapDict[(0,0)] = [((0,),None)]
        self.assertEqual(mapDict, self.testMapDict)

    def test_poly_outside_left(self):
        lat = numpy.array([[1.5, 2.5, 2.5, 1.5]])
        lon = numpy.array([[-1.5, -1.5, -0.5, -0.5]])
        ind = numpy.array([[i] for i in range(lat.shape[0])])
        self.parser.prime_corners(lat, lon, ind)
        mapDict = self.mapper(self.parser, self.grid)
        self.assertEqual(mapDict, self.testMapDict)

    def test_poly_outside_top(self):
        lat = numpy.array([[5.5, 7.5, 75, 5.5]])
        lon = numpy.array([[0.5, 0.5, 2.5, 2.5]])
        ind = numpy.array([[i] for i in range(lat.shape[0])])
        self.parser.prime_corners(lat, lon, ind)
        mapDict = self.mapper(self.parser, self.grid)
        self.assertEqual(mapDict, self.testMapDict)

    def test_poly_outside_right(self):
        lat = numpy.array([[3.5, 4.5, 4.5, 3.5]])
        lon = numpy.array([[10.25, 10.25, 10.75, 10.75]])
        ind = numpy.array([[i] for i in range(lat.shape[0])])
        self.parser.prime_corners(lat, lon, ind)
        mapDict = self.mapper(self.parser, self.grid)
        self.assertEqual(mapDict, self.testMapDict)

    def test_poly_outside_bottom(self):
        lat = numpy.array([[-1.5, -0.75, -0.75, -1.5]])
        lon = numpy.array([[6.25, 6.25, 6.75, 6.75]])
        ind = numpy.array([[i] for i in range(lat.shape[0])])
        self.parser.prime_corners(lat, lon, ind)
        mapDict = self.mapper(self.parser, self.grid)
        self.assertEqual(mapDict, self.testMapDict)

    def test_2D_array_of_polys(self):
        lat = numpy.array([[[0.5, 1.5, 1.5, 0.5], [1.5, 2.5, 2.5, 1.5], [2.5, 3.5, 3.5, 2.5]],
                           [[0.5, 1.5, 1.5, 0.5], [1.5, 2.5, 2.5, 1.5], [2.5, 3.5, 3.5, 2.5]]])
        lon = numpy.array([[[0.5, 0.5, 1.5, 1.5], [0.5, 0.5, 1.5, 1.5], [0.5, 0.5, 1.5, 1.5]],
                           [[1.5, 1.5, 2.5, 2.5], [1.5, 1.5, 2.5, 2.5], [1.5, 1.5, 2.5, 2.5]]])
        ind = numpy.array([[[0,0], [0,1], [0,2]],
                            [[1,0], [1,1], [1,2]]])
        self.parser.prime_corners(lat, lon, ind)
        mapDict = self.mapper(self.parser, self.grid)
        for (x,y) in product(range(2), range(3)):
            thisPoly = ((x,y), None)
            self.testMapDict[(y,x)].append(thisPoly)
            self.testMapDict[(y+1,x)].append(thisPoly)
            self.testMapDict[(y,x+1)].append(thisPoly)
            self.testMapDict[(y+1,x+1)].append(thisPoly)
        self.assertEqual(mapDict, self.testMapDict)

    def test_indices_unrelated_to_polys(self):
        lat = numpy.array([[1.5, 2.5, 2.5, 1.5]])
        lon = numpy.array([[0.5, 0.5, 1.5, 1.5]])
        ind = numpy.array([[1,2,3]])
        self.parser.prime_corners(lat, lon, ind)
        mapDict = self.mapper(self.parser, self.grid)
        listEntry = [((1,2,3),None)]
        keys = [(1,0), (2,0), (2,1), (1,1)]
        for k in keys:
            self.testMapDict[k] = listEntry
        self.assertEqual(mapDict, self.testMapDict)
        
class Test_point_in_cell(TestMapGeo):
    
    
    def setUp(self):
        TestMapGeo.setUp(self)
        self.mapper = getattr(map_geo, 'point_in_cell_map_geo')
        
    def test_single_pix_inside_cell(self):
        lat = numpy.array([1.2])
        lon = numpy.array([2.1])
        ind = numpy.arange(lat.shape[0])
        self.parser.prime_centers(lat, lon, ind)
        self.testMapDict[(1,2)] = [((0,), None)]
        self.assertMapsEqual(self.mapper)
        
    def test_mult_pix_inside_cells(self):
        lat = numpy.array([1.2, 1.9, 1.5])
        lon = numpy.array([2.1, 2.9, 7.5])
        ind = numpy.arange(lat.shape[0])
        self.parser.prime_centers(lat, lon, ind)
        self.testMapDict[(1,2)] = [((0,), None), ((1,), None)]
        self.testMapDict[(1,7)] = [((2,), None)]
        self.assertMapsEqual(self.mapper)
        
    def test_pixel_on_cell_col_edge(self):
        '''As laid out in documentation, should go to right-hand cell'''
        lat = numpy.array([3.1])
        lon = numpy.array([6.0])
        ind = numpy.arange(lat.shape[0])
        self.parser.prime_centers(lat, lon, ind)
        self.testMapDict[(3,6)] = [((0,), None)]
        self.assertMapsEqual(self.mapper)
        
    def test_pixel_on_cell_row_edge(self):
        '''As laid out in documentation, should go to upper cell'''
        lat = numpy.array([2.0])
        lon = numpy.array([4.4])
        ind = numpy.arange(lat.shape[0])
        self.parser.prime_centers(lat, lon, ind)
        self.testMapDict[(2, 4)] = [((0,), None)]
        self.assertMapsEqual(self.mapper)
        
    def test_accepts_pixel_on_lower_edge(self):
        lat = numpy.array([0.0])
        lon = numpy.array([2.4])
        ind = numpy.arange(lat.shape[0])
        self.parser.prime_centers(lat, lon, ind)
        self.testMapDict[(0,2)] = [((0,), None)]
        self.assertMapsEqual(self.mapper)
        
    def test_accepts_pixel_on_left_edge(self):
        lat = numpy.array([1.2])
        lon = numpy.array([0.0])
        ind = numpy.arange(lat.shape[0])
        self.parser.prime_centers(lat, lon, ind)
        self.testMapDict[(1,0)] = [((0,), None)]
        self.assertMapsEqual(self.mapper)
        
    def test_rejects_pixel_on_upper_edge(self):
        lat = numpy.array([5])
        lon = numpy.array([8.3])
        ind = numpy.arange(lat.shape[0])
        self.parser.prime_centers(lat, lon, ind)
        self.assertMapsEqual(self.mapper)
        
    def test_rejects_pixel_on_right_edge(self):
        lat = numpy.array([4.2])
        lon = numpy.array([10])
        ind = numpy.arange(lat.shape[0])
        self.parser.prime_centers(lat, lon, ind)
        self.assertMapsEqual(self.mapper)
        
    def test_rejects_exterior_pixels_sides(self):
        '''Test that pixels "straight out" from the sides are properly rejected'''
        lat = numpy.array([2.5, 10, 2.5, -10])
        lon = numpy.array([-15, 5, 15, 5])
        ind = numpy.arange(lat.shape[0])
        self.parser.prime_centers(lat, lon, ind)
        self.assertMapsEqual(self.mapper)
        
    def test_rejects_exterior_pixels_corners(self):
        '''
        Test that pixels beyond corners (out of bounds in 2 categories)
        are properly rejected.
        '''
        lat = numpy.array([10, 10, -10, -10])
        lon = numpy.array([-15, 15, 15, -15])
        ind = numpy.arange(lat.shape[0])
        self.parser.prime_centers(lat, lon, ind)
        self.assertMapsEqual(self.mapper)
        
    def test_arbitrary_index(self):
        lat = numpy.array([2.3])
        lon = numpy.array([3.1])
        ind = numpy.array([[1010, 2020]])
        self.parser.prime_centers(lat, lon, ind)
        self.testMapDict[(2,3)] = [((1010,2020), None)]
        self.assertMapsEqual(self.mapper)
        
    def test_2D_array_of_points(self):
        lat = numpy.array([[0.5, 1.5, 2.5], [3.5, 4.5, 5.5]])
        lon = numpy.array([[10.5, 8.5, 6.5], [4.5, 2.5, 0.5]])
        ind = numpy.array([[[0,0], [0,1], [0,2]], [[1,0], [1,1], [1,2]]])
        self.parser.prime_centers(lat, lon, ind)
        self.testMapDict[(1,8)] = [((0,1), None)]
        self.testMapDict[(2,6)] = [((0,2), None)]
        self.testMapDict[(3,4)] = [((1,0), None)]
        self.testMapDict[(4,2)] = [((1,1), None)]
        self.assertMapsEqual(self.mapper)

class Test_global_intersect(TestMapGeo):

    def setUp(self):
        TestMapGeo.setUp(self)
        self.mapper = getattr(map_geo, 'global_intersect_map_geo')

    def test_single_pix_inside_cell(self):
        lat = numpy.array([[1.1, 1.1, 1.8, 1.8]])
        lon = numpy.array([[2.1, 2.8, 2.8, 2.1]])
        ind = numpy.array([[ind] for ind in range(lat.shape[0])])
        self.parser.prime_corners(lat, lon, ind)
        self.testMapDict[(1,2)] = [((0,), None)]
        self.assertMapsEqual(self.mapper)

    def test_mult_pix_inside_multiple_cells(self):
        lat = numpy.array([[1.1, 1.1, 1.8, 1.8],
                           [2.1, 2.1, 2.8, 2.8]])
        lon = numpy.array([[2.1, 2.8, 2.8, 2.1],
                           [7.1, 7.8, 7.8, 7.1]])
        ind = numpy.array([[ind] for ind in range(lat.shape[0])])
        self.parser.prime_corners(lat, lon, ind)
        self.testMapDict[(1,2)] = [((0,), None)]
        self.testMapDict[(2,7)] = [((1,), None)]
        self.maxDiff = None
        self.assertMapsEqual(self.mapper)

    def test_mult_pix_inside_single_cell(self):
        lat = numpy.array([[1.1, 1.1, 1.8, 1.8],
                           [1.2, 1.2, 1.7, 1.7]])
        lon = numpy.array([[2.1, 2.8, 2.8, 2.1],
                           [2.2, 2.7, 2.7, 2.2]])
        ind = numpy.array([[ind] for ind in range(lat.shape[0])])
        self.parser.prime_corners(lat, lon, ind)
        self.testMapDict[(1,2)] = [((0,), None), ((1,), None)]
        self.maxDiff = None
        self.assertMapsEqual(self.mapper)

    def test_pixel_in_multiple_cells(self):
        lat = numpy.array([[2.5, 2.5, 3.5, 3.5]])
        lon = numpy.array([[1.5, 2.5, 2.5, 1.5]])
        ind = numpy.array([[ind] for ind in range(lat.shape[0])])
        self.parser.prime_corners(lat, lon, ind)
        self.testMapDict[(2,1)] = [((0,), None)]
        self.testMapDict[(2,2)] = [((0,), None)]
        self.testMapDict[(3,1)] = [((0,), None)]
        self.testMapDict[(3,2)] = [((0,), None)]
        self.assertMapsEqual(self.mapper)

    def test_pixel_on_lower_edge(self):
        lat = numpy.array([[-0.5, -0.5, 0.5, 0.5]])
        lon = numpy.array([[5.1, 5.5, 5.5, 5.1]])
        ind = numpy.array([[ind] for ind in range(lat.shape[0])])
        self.parser.prime_corners(lat, lon, ind)
        self.testMapDict[(0, 5)] = [((0,), None)]
        self.assertMapsEqual(self.mapper)

    def test_pixel_on_upper_edge(self):
        lat = numpy.array([[4.5, 4.5, 5.5, 5.5]])
        lon = numpy.array([[7.1, 7.6, 7.6, 7.1]])
        ind = numpy.array([[ind] for ind in range(lat.shape[0])])
        self.parser.prime_corners(lat, lon, ind)
        self.testMapDict[(4,7)] = [((0,), None)]
        self.assertMapsEqual(self.mapper)

    def test_rejects_exterior_pixels(self):
        lat = numpy.array([[5.5, 5.5, 6.5, 6.5],
                           [2.5, 2.5, 3.5, 3.5],
                           [-1.5, -1.5, -0.5, -0.5],
                           [2.5, 2.5, 3.5, 3.5],
                           [-1.5, -1.5, -0.5, -0.5]])
        lon = numpy.array([[3.5, 4.5, 4.5, 3.5],
                           [11.5, 12.5, 12.5, 11.5],
                           [3.5, 4.5, 4.5, 3.5],
                           [-1.5, -0.5, -0.5, -1.5],
                           [-1.5, -0.5, -0.5, -1.5]])
        ind = numpy.array([[ind] for ind in range(lat.shape[0])])
        self.parser.prime_corners(lat, lon, ind)
        self.assertMapsEqual(self.mapper)

    def test_pixel_on_right_edge(self):
        lat = numpy.array([[3.3, 3.3, 3.6, 3.6]])
        lon = numpy.array([[9.1, 10.0, 10.0, 9.1]])
        ind = numpy.array([[ind] for ind in range(lat.shape[0])])
        self.parser.prime_corners(lat, lon, ind)
        self.testMapDict[(3,9)] = [((0,), None)]
        self.assertMapsEqual(self.mapper)

    def test_pixel_on_left_edge(self):
        lat = numpy.array([[3.3, 3.3, 3.6, 3.6]])
        lon = numpy.array([[0.0, 0.4, 0.4, 0.0]])
        ind = numpy.array([[ind] for ind in range(lat.shape[0])])
        self.parser.prime_corners(lat, lon, ind)
        self.testMapDict[(3,0)] = [((0,), None)]
        self.assertMapsEqual(self.mapper)

    def test_continuous_pixel_spanning_mid(self):
        lat = numpy.array([[1.1, 1.1, 1.7, 1.7]])
        lon = numpy.array([[4.8, 5.3, 5.3, 4.8]])
        ind = numpy.array([[ind] for ind in range(lat.shape[0])])
        self.parser.prime_corners(lat, lon, ind)
        self.testMapDict[(1,4)] = [((0,), None)]
        self.testMapDict[(1,5)] = [((0,), None)]
        self.assertMapsEqual(self.mapper)

    def test_discontinuous_pixel_2_cells(self):
        lat = numpy.array([[0.2, 0.2, 0.8, 0.8]])
        lon = numpy.array([[9.8, 0.2, 0.2, 9.8]])
        ind = numpy.array([[ind] for ind in range(lat.shape[0])])
        self.parser.prime_corners(lat, lon, ind)
        self.testMapDict[(0,9)] = [((0,), None)]
        self.testMapDict[(0,0)] = [((0,), None)]
        self.maxDiff = None
        self.assertMapsEqual(self.mapper)

    def test_discontinuous_pixel_vertical_slant(self):
        lat = numpy.array([[3.1, 3.6, 2.9, 2.4]])
        lon = numpy.array([[9.3, 9.3, 0.7, 0.7]])
        ind = numpy.array([[ind] for ind in range(lat.shape[0])])
        self.parser.prime_corners(lat, lon, ind)
        self.testMapDict[(2,9)] = [((0,), None)]
        self.testMapDict[(3,9)] = [((0,), None)]
        self.testMapDict[(2,0)] = [((0,), None)]
        self.testMapDict[(3,0)] = [((0,), None)]
        self.assertMapsEqual(self.mapper)

    def test_discontinuous_pixel_3_1_split(self):
        lat = numpy.array([[1.1, 1.5, 1.9, 1.5]])
        lon = numpy.array([[0.2, 0.7, 0.2, 9.7]])
        ind = numpy.array([[ind] for ind in range(lat.shape[0])])
        self.parser.prime_corners(lat, lon, ind)
        self.testMapDict[(1,0)] = [((0,), None)]
        self.testMapDict[(1,9)] = [((0,), None)]
        self.assertMapsEqual(self.mapper)

    def test_2D_index(self):
        lat = numpy.array([[2.1, 2.1, 2.4, 2.4]])
        lon = numpy.array([[7.1, 7.4, 7.4, 7.1]])
        ind = numpy.array([[3,2,1]])
        self.parser.prime_corners(lat, lon, ind)
        self.testMapDict[(2,7)] = [((3,2,1), None)]
        self.assertMapsEqual(self.mapper)

    def test_3D_array_of_pixels(self):
        lat = numpy.array([[[1.4, 1.4, 1.6, 1.6], [1.4, 1.4, 1.6, 1.6], [1.4, 1.4, 1.6, 1.6]], 
                           [[0.4, 0.4, 0.6, 0.6], [0.4, 0.4, 0.6, 0.6], [0.4, 0.4, 0.6, 0.6]]])
        lon = numpy.array([[[0.3, 0.7, 0.7, 0.3], [1.3, 1.7, 1.7, 1.3], [2.3, 2.7, 2.7, 2.3]], 
                           [[0.3, 0.7, 0.7, 0.3], [1.3, 1.7, 1.7, 1.3], [2.3, 2.7, 2.7, 2.3]]])
        ind = numpy.array([[[1], [2], [3]],
                           [[4], [5], [6]]])
        self.parser.prime_corners(lat, lon, ind)
        self.testMapDict[(0,0)] = [((4,), None)]
        self.testMapDict[(1,0)] = [((1,), None)]
        self.testMapDict[(0,1)] = [((5,), None)]
        self.testMapDict[(0,2)] = [((6,), None)]
        self.testMapDict[(1,1)] = [((2,), None)]
        self.testMapDict[(1,2)] = [((3,), None)]
        self.assertMapsEqual(self.mapper)

class TestUtils(unittest.TestCase):
    
    
    def test_wrap_lon_0_360_known_inputs(self):
        knownIn = [0, 360, 180, 47.5, -47.5, -80, -180, -365, -725, 410]
        knownOut = [0, 0, 180, 47.5, 312.5, 280, 180, 355, 355, 50]
        compOut = [utils.wrap_lon_0_360(lon) for lon in knownIn]
        self.assertListEqual(knownOut, compOut)
        
    def test_wrap_lon_neg180_180(self):
        knownIn = [0, 360, 180, 47.5, -47.5, -80, -180, -365, -725, 410]
        knownOut = [0, 0, 180, 47.5, -47.5, -80, 180, -5, -5, 50]
        compOut = [utils.wrap_lon_neg180_180(lon) for lon in knownIn]
        self.assertListEqual(knownOut, compOut)
        
    def test_lon_0_360_recoverable(self):
        inputData = numpy.random.rand(20)*360
        intermediateData = [utils.wrap_lon_neg180_180(lon) for lon in inputData]
        finalData = numpy.array([utils.wrap_lon_0_360(lon) for lon in intermediateData])
        numpy.testing.assert_array_almost_equal(inputData, finalData, decimal=5)
        
    def test_lon_neg180_180_recoverable(self):
        inputData = numpy.random.rand(20)*360 - 180
        intermediateData = [utils.wrap_lon_0_360(lon) for lon in inputData]
        finalData = numpy.array([utils.wrap_lon_neg180_180(lon) for lon in intermediateData])
        numpy.testing.assert_array_almost_equal(inputData, finalData, decimal=5)
        
    def test_timestr_to_nsecs_known_input_defaults(self):
        inputTstr = ['00:00:00 01-01-1970', '01:00:00 01-02-1970',
                     '08:30:15 01-01-1970', '00:00:00 12-31-1969',
                     '00:00:00 01-01-1971', '00:00:00 01-01-1973']
        knownOut = [0, 90000,
                    30615, -86400,
                    31536000, 94694400]
        calcOut = [utils.timestr_to_nsecs(tstr) for tstr in inputTstr]
        self.assertListEqual(knownOut, calcOut)
        
    def test_timestr_to_nsecs_accepts_epoch(self):
        inputTstr = '12:30:30 08-15-2011'
        epochs = ['00:00:00 08-15-2011', '00:00:00 08-15-2010', '12:15:15 08-15-2011']
        knownOut = [45030, 31581030, 915]
        calcOut = [utils.timestr_to_nsecs(inputTstr, epoch=ep) for ep in epochs]
        self.assertListEqual(knownOut, calcOut)
        
    def test_timestr_to_nsecs_accepts_format(self):
        formats = ['%H:%M:%S %m/%d/%Y', '%m %d %y %H %M %S ', '%I:%M:%S %p %b %d, %Y']
        inputTstr = ['14:22:18 01/02/1970', '01 02 70 14 22 18 ', '02:22:18 PM Jan 02, 1970']
        epochs = ['00:00:00 01/01/1970', '01 01 70 00 00 00 ', '12:00:00 AM Jan 01, 1970']
        knownOut = [138138, 138138, 138138]
        calcOut = [utils.timestr_to_nsecs(tStr, format=fmt, epoch=ep) for (tStr, fmt, ep) in zip(inputTstr, formats, epochs)]
        self.assertListEqual(knownOut, calcOut)
        
    def test_timestr_to_nsecs_format_affects_epoch(self):
        fmt = '%H:%M:%S %m/%d/%y'
        ep = '00:00:00 07/10/90'
        tStr = '14:08:08 07/10/90'
        knownOut = 50888
        calcOut = utils.timestr_to_nsecs(tStr, epoch=ep, format=fmt)
        self.assertEqual(knownOut, calcOut)
        
    def test_nsecs_to_timestr_known_input_defaults(self):
        inputNsecs = [0, 90000, 30615,
                      -86400, 31536000, 94694400]
        knownOut = ['00:00:00 01-01-1970', '01:00:00 01-02-1970', '08:30:15 01-01-1970',
                    '00:00:00 12-31-1969', '00:00:00 01-01-1971', '00:00:00 01-01-1973']
        calcOut = [utils.nsecs_to_timestr(nSecs) for nSecs in inputNsecs]
        self.assertListEqual(knownOut, calcOut)
        
    def test_nsecs_to_timestr_accepts_epoch(self):
        epochs = ['00:00:00 08-15-2011', '00:00:00 08-15-2010', '12:15:15 08-15-2011']
        inputNsecs = [45030, 31581030, 915]
        knownOut = ['12:30:30 08-15-2011', '12:30:30 08-15-2011', '12:30:30 08-15-2011']
        calcOut = [utils.nsecs_to_timestr(nSec, ep) for (nSec, ep) in zip(inputNsecs, epochs)]
        self.assertListEqual(knownOut, calcOut)
        
    def test_nsecs_to_timestr_accepts_format(self):
        formats = ['%H:%M:%S %m/%d/%Y', '%m %d %y %H %M %S ', '%I:%M:%S %p %b %d, %Y']
        inputNsecs = 138138
        knownOut = ['14:22:18 01/02/1970', '01 02 70 14 22 18 ', '02:22:18 PM Jan 02, 1970']
        epochs = ['00:00:00 01/01/1970', '01 01 70 00 00 00 ', '12:00:00 AM Jan 01, 1970']
        calcOut = [utils.nsecs_to_timestr(inputNsecs,format=fmt, epoch=ep) for (fmt, ep) in zip(formats, epochs)]
        self.assertListEqual(knownOut, calcOut)
        
    def test_nsecs_to_timstr_format_affects_epoch(self):
        fmt = '%H:%M:%S %m/%d/%y'
        ep = '00:00:00 07/10/90'
        nSecs = 50888
        knownOut = '14:08:08 07/10/90'
        calcOut = utils.nsecs_to_timestr(nSecs, epoch=ep, format=fmt)
        self.assertEqual(knownOut, calcOut)
        
    def test_nsecs_and_tstr_reversible(self):
        inputNsecs = numpy.round(numpy.random.rand(25)*200000000-100000000)
        intermediateTstr = [utils.nsecs_to_timestr(nSecs) for nSecs in inputNsecs]
        outputNsecs = numpy.array([utils.timestr_to_nsecs(tStr) for tStr in intermediateTstr])
        numpy.testing.assert_allclose(inputNsecs, outputNsecs, atol=.5)
        
    def test_UTCoffset_from_lon_known_input(self):
        inLons = [0, -4.5, -10.5, 17.5, 95, -110]
        knownOut = [0, 0, -1, 1, 6, -7]  # in hours
        knownOut = [3600*out for out in knownOut] # in seconds
        calcOut = [utils.UTCoffset_from_lon(lon) for lon in inLons]
        self.assertListEqual(knownOut, calcOut)
        
    def test_UTC_offset_from_lon_0_360_input(self):
        inLons = [355.5, 349.5, 250, 360]
        knownOut = [0, -1, -7, 0] # in hours
        knownOut = [3600*out for out in knownOut] # in seconds
        calcOut = [utils.UTCoffset_from_lon(lon) for lon in inLons]
        self.assertListEqual(knownOut, calcOut)

class TestOutGeo(unittest.TestCase):
    
    
    def setUp(self):
        self.parser = fakeParser('foo.dat')
    

    def test_ValidOutfuncs_are_valid(self):
        funcs = out_geo.ValidOutfuncs()
        for func in funcs:
            if not hasattr(out_geo, func+'_out_func'):
                raise AssertionError("Function %s_out_func not in module" % func)

    def toTAI93(self, timeStr):
        epoch = '00:00:00 01-01-1993'
        return utils.timestr_to_nsecs(timeStr, epoch=epoch)

@unittest.skip("Skipped until this function is fixed")
class Test_OMNO2e_wght_avg_out_func(TestOutGeo):


    def setUp(self):
        TestOutGeo.setUp(self)
        gridParms = { 'xOrig' : 0, 'yOrig' : 0, 
                      'xCell' : 0, 'yCell' : 0,
                      'nRows' : 1, 'nCols' : 1 }
        self.one_el_grid = grid_geo.latlon_GridDef(gridParms)
        self.mapDict = { 'parser' : self.parser }
        fieldnames = { 'toAvg' : 'toAvg',
                       'overallQualFlag' : 'qualFlag',
                       'cloudFrac' : 'cFrac',
                       'solarZenithAngle' : 'solZenAng' }
        outfuncParms = { 'cloudFractUpperCutoff' : .25,
                         'solarZenAngUpperCutoff' : 80,
                         'pixIndXtrackAxis' : 1,
                         'fillVal' : -99999}
        outfuncParms.update(fieldnames)
        self.outfunc = out_geo.OMNO2e_wght_avg_out_func(outfuncParms)
        self.outFname = (tempfile.mkstemp())[1]
        # create a bunch of empty arrays to put data in.
        self.cfrac = numpy.zeros((20,60))
        self.qualFlag = numpy.zeros((20,60))
        self.solZenAng = numpy.zeros((20,60))
        self.toAvg = numpy.zeros((20,60))
        # bind shallow copies of these arrays into the parser
        self.parser.prime_get('toAvg', self.toAvg)
        self.parser.prime_get('qualFlag', self.qualFlag)
        self.parser.prime_get('solZenAng', self.solZenAng)
        self.parser.prime_get('cFrac', self.cfrac)


    def tearDown(self):
        os.remove(self.outFname)

    def test_parser_still_in_map_at_end(self):
        self.cfrac[(1,30)] = .2
        self.qualFlag[(1,30)] = 0
        self.solZenAng[(1,30)] = 30
        self.toAvg[(1,30)] = 1337
        self.mapDict[(0,0)] = [((1,30),None)]
        self.outfunc(self.mapDict, self.one_el_grid, self.outFname, verbose=False)
        self.assertIs(self.mapDict['parser'], self.parser)

    def test_works_for_single_value(self):
        self.cfrac[(1,30)] = .2
        self.qualFlag[(1,30)] = 0
        self.solZenAng[(1,30)] = 30
        self.toAvg[(1,30)] = 1337
        self.mapDict[(0,0)] = [((1,30),None)]
        avg = self.outfunc(self.mapDict, self.one_el_grid, self.outFname, verbose=False)
        self.assertAlmostEqual(avg[(0,0)], 1337)

    def test_gives_accurate_weight(self):
        # calculation performed by hand with scientific calculator
        self.cfrac[0,29:31] = [.2, .15]
        self.qualFlag[0,29:31] = [0, 0]
        self.solZenAng[0,29:31] = [20, 25]
        self.toAvg[0,29:31] = [1, 2]
        self.mapDict[(0,0)] = [((0,29),None), ((0,30),None)]
        avg = self.outfunc(self.mapDict, self.one_el_grid, self.outFname, verbose=False)
        self.assertAlmostEqual(avg[(0,0)], 1.5490637021)

    def test_zero_weight_if_sum_flag_is_set(self):
        self.cfrac[0,29:31] = [.4, .1]
        self.qualFlag[0,29:31] = [1, 0]
        self.solZenAng[0,29:31] = [30, 30]
        self.toAvg[0,29:31] = [34, 42]
        self.mapDict[(0,0)] = [((0,29),None), ((0,30),None)]
        avg = self.outfunc(self.mapDict, self.one_el_grid, self.outFname, verbose=False)
        self.assertEqual(avg[(0,0)], 42)

    def test_zero_weight_if_cfrac_gt_threshold(self):
        self.cfrac[0,29:31] = [.28, .24]
        self.qualFlag[0,29:31] = [0,0]
        self.solZenAng[0,29:31] = [30, 24]
        self.toAvg[0,29:31] = [1, 421]
        self.mapDict[(0,0)] = [((0,29),None), ((0,30),None)]
        avg = self.outfunc(self.mapDict, self.one_el_grid, self.outFname, verbose=False)
        self.assertEqual(avg[(0,0)], 421)

    def test_zero_weight_if_SZA_gt_threshold(self):
        self.cfrac[0,29:31] = [.1, .1]
        self.qualFlag[0,29:31] = [0,0]
        self.solZenAng[0,29:31] = [81.2, 79]
        self.toAvg[0,29:31] = [1, 421]
        self.mapDict[(0,0)] = [((0,29),None), ((0,30),None)]
        avg = self.outfunc(self.mapDict, self.one_el_grid, self.outFname, verbose=False)
        self.assertEqual(avg[(0,0)], 421)

    def test_zero_weight_if_mult_problems(self):
        self.cfrac[0,28:33] = [.4, .5, .1, .5, .1]
        self.qualFlag[0,28:33] = [1, 0, 1, 1, 0]
        self.solZenAng[0,28:33] = [40, 88, 88, 88, 40]
        self.toAvg[0,28:33] = [10, 22, 36, 49, 24601]
        self.mapDict[(0,0)] = [((0,28),None), ((0,29),None), 
                               ((0,30),None), ((0,31),None),
                               ((0,32),None)]
        avg = self.outfunc(self.mapDict, self.one_el_grid, self.outFname, verbose=False)
        self.assertEqual(avg[(0,0)], 24601)

    def test_all_zero_weight_yields_fillVal(self):
        self.cfrac[0,29:31] = [.1, .1]
        self.qualFlag[0,29:31] = [1, 1]
        self.solZenAng[0,29:31] = [40, 40]
        self.toAvg[0,29:31] = [1, 2]
        self.mapDict[(0,0)] = [((0,29),None), ((0,30),None)]
        avg = self.outfunc(self.mapDict, self.one_el_grid, self.outFname, verbose=False)
        self.assertEqual(avg[(0,0)], -99999)

    def test_NaN_in_toAvg_yields_fillVal(self):
        self.cfrac[0,29:31] = [.1, .1]
        self.qualFlag[0,29:31] = [1, 1]
        self.solZenAng[0,29:31] = [40, 40]
        self.toAvg[0,29:31] = [numpy.nan, numpy.nan]
        self.mapDict[(0,0)] = [((0,29),None), ((0,30),None)]
        avg = self.outfunc(self.mapDict, self.one_el_grid, self.outFname, verbose=False)
        self.assertEqual(avg[(0,0)], -99999)

    def test_NaN_in_toAvg_ignored_when_weight_is_0(self):
        self.cfrac[0,29:31] = [.1, .1]
        self.qualFlag[0,29:31] = [0, 1]
        self.solZenAng[0,29:31] = [40, 40]
        self.toAvg[0,29:31] = [713, numpy.nan]
        self.mapDict[(0,0)] = [((0,29),None), ((0,30),None)]
        avg = self.outfunc(self.mapDict, self.one_el_grid, self.outFname, verbose=False)
        self.assertEqual(avg[(0,0)], 713)

    def test_no_pixels_in_cell_yields_fillVal(self):
        self.mapDict[(0,0)] = []
        avg = self.outfunc(self.mapDict, self.one_el_grid, self.outFname, verbose=False)
        self.assertEqual(avg[(0,0)], -99999)

    def test_multi_element_pix_array(self):
        self.cfrac[0:2, 28:34] = [[.11, .12, .13, .14, .15, .16],
                                  [.21, .22, 23, .24, .25, .26]]
        self.qualFlag[0:2, 28:34] = [[0, 0, 0, 0, 0, 0],
                                     [0, 0, 0, 1, 1, 1]]
        self.solZenAng[0:2, 28:34] = [[10, 11, 12, 13, 14, 15],
                                      [87, 88, 89, 20, 12, 13]]
        self.toAvg[0:2, 28:34] = [[100, 200, 300, 400, 500, 600],
                                    [10, 20, 30, 40, 50, 60]]
        gridParms = { 'xOrig' : 0, 'yOrig' : 0, 
                      'xCell' : 0, 'yCell' : 0,
                      'nRows' : 2, 'nCols' : 3 }
        six_el_grid = grid_geo.latlon_GridDef(gridParms)
        for (i,j) in product(range(2), range(3)):
            oneDind = i*3+j
            self.mapDict[(i,j)] = [((0,28+oneDind),None), ((1,28+oneDind),None)]
        avg = self.outfunc(self.mapDict, six_el_grid, self.outFname, verbose=False)
        expectedOut = [[100, 200, 300], [400, 500, 600]]
        numpy.testing.assert_array_equal(avg, expectedOut)

    def test_multidict_one_zero_weight(self):
        self.cfrac[0, 29:31] = [.1, .5]
        self.qualFlag[0, 29:31] = [0, -127]
        self.solZenAng[0, 29:31] = [10, 89]
        self.toAvg[0, 29:31] = [117, 343]
        secondMapDict = dict(self.mapDict)
        self.mapDict[(0,0)] = [((0,29),None)]
        secondMapDict[(0,0)] = [((0,30),None)]
        dictList = [self.mapDict, secondMapDict]
        avg = self.outfunc(dictList, self.one_el_grid, self.outFname, verbose=False)
        self.assertAlmostEqual(avg[(0,0)], 117)

    def test_multidict_both_with_nonzero_weight(self):
        # calculation performed by hand with scientific calculator
        self.cfrac[0,0:2] = [.05, .10]
        self.qualFlag[0,0:2] = [0, 0]
        self.solZenAng[0,0:2] = [50, 55]
        self.toAvg[0,0:2] = [25000, 30000]
        secondMapDict = dict(self.mapDict)
        self.mapDict[(0,0)] = [((0,0),None)]
        secondMapDict[(0,0)] = [((0,1),None)]
        dictList = [self.mapDict, secondMapDict]
        avg = self.outfunc(dictList, self.one_el_grid, self.outFname, verbose=False)
        self.assertAlmostEqual(27394.49255639, avg[(0,0)])
        
class Test_OMNO2e_netCDF_avg_out_func(TestOutGeo):
    
    def setUp(self):
        TestOutGeo.setUp(self)
        self.version = 'TEST VERSION'
        self.one_el_gridParms = {'xOrig' : 0, 'yOrig' : 0,
                            'xCell' : 0, 'yCell' : 0,
                            'nRows' : 1, 'nCols' : 1}
        self.one_el_grid = grid_geo.latlon_GridDef(self.one_el_gridParms)
        self.six_el_gridParms = {'xOrig' : 0, 'yOrig' : 0,
                            'xCell' : 0, 'yCell' : 0,
                            'nRows' : 2, 'nCols' : 3}
        self.six_el_grid = grid_geo.latlon_GridDef(self.six_el_gridParms)
        self.mapDict = {'parser' : self.parser}
        self.fnames = {'overallQualFlag' : 'qualFlag',
                  'cloudFrac' : 'cfrac',
                  'solarZenithAngle' : 'solZenAng',
                  'time' : 'time',
                  'longitude' : 'lon'}
        self.startTimeStr = '00:00:00 08-30-2011'
        self.stopTimeStr = '23:59:59 08-30-2011'
        startTimeParm = self.startTimeStr.replace(' ','_')
        stopTimeParm = self.stopTimeStr.replace(' ', '_')
        self.defParms = {'inFieldNames' : ['test2D', 'test3D'],
                    'outFieldNames' : ['outTest2D', 'outTest3D'],
                    'outUnits' : ['Jigawatts', 'fathoms'],
                    'extraDimLabel' : ['Irrelevant', 'layer'],
                    'extraDimSize' : [0, 4],
                    'timeComparison' : 'UTC',
                    'timeStart' : startTimeParm,
                    'timeStop' : stopTimeParm,
                    'cloudFractUpperCutoff' : .25,
                    'solarZenAngUpperCutoff' : 80,
                    'pixIndXtrackAxis' : 1,
                    'fillVal' : -99999.0,
                    'includePixelCount' : False}
        self.defParms.update(self.fnames)
        self.fnames.update(self.defParms)
        self.defOutFunc = out_geo.OMNO2e_netCDF_avg_out_func(self.defParms)
        (outFid, self.outFname) = (tempfile.mkstemp())
        os.close(outFid)
        # create empty arrays to hold data
        self.cfrac = numpy.zeros((20,60))
        self.qualFlag = numpy.zeros((20,60))
        self.solZenAng = numpy.zeros((20,60))
        self.time = numpy.zeros((20,60))
        self.lon = numpy.zeros((20,60))
        self.test2D = numpy.zeros((20,60))
        self.test3D = numpy.zeros((20,60,4))
        # bind empty arrays to fake parser
        self.parser.prime_get('qualFlag', self.qualFlag)
        self.parser.prime_get('solZenAng', self.solZenAng)
        self.parser.prime_get('cfrac', self.cfrac)
        self.parser.prime_get('time', self.time)
        self.parser.prime_get('lon', self.lon)
        self.parser.prime_get('test2D', self.test2D)
        self.parser.prime_get('test3D', self.test3D)
        
    def tearDown(self):
        try:
            self.fid.close()
        except(AttributeError):
            pass
        os.remove(self.outFname)
        
    def test_parser_still_in_empty_map_at_end(self):
        self.mapDict[(0,0)] = []
        unused_result = self.defOutFunc(self.mapDict, self.one_el_grid, 
                                        self.outFname, verbose=False, 
                                        version = self.version)
        self.assertIs(self.mapDict['parser'], self.parser)
    
    def test_parser_still_in_nonempty_map_at_end(self):
        self.mapDict[(0,0)] = [((0,31), None)]
        unused_result = self.defOutFunc(self.mapDict, self.one_el_grid,
                                        self.outFname, verbose=False,
                                        version = self.version)
        self.assertIs(self.mapDict['parser'], self.parser)
    
    def test_result_right_shape_emptyMap_2D(self):
        self.mapDict[(0,0)] = []
        resDict = self.defOutFunc(self.mapDict, self.six_el_grid,
                                  self.outFname, verbose=False,
                                  version = self.version)
        self.assertEqual(resDict['outTest2D'].shape, (2,3))
        
    def test_result_right_shape_emptyMap_3D(self):
        self.mapDict[(0,0)] = []
        resDict = self.defOutFunc(self.mapDict, self.six_el_grid,
                                  self.outFname, verbose=False,
                                  version = self.version)
        self.assertEqual(resDict['outTest3D'].shape, (2,3,4))
        
    def test_result_right_shape_nonempty_map_2D(self):
        self.mapDict[(0,0)] = [((0,0), None), ((10,34), None)]
        self.mapDict[(1,1)] = [((3,0), None), ((19, 5), None), ((0,1), None)] 
        resDict = self.defOutFunc(self.mapDict, self.six_el_grid,
                                  self.outFname, verbose=False,
                                  version = self.version)
        self.assertEqual(resDict['outTest2D'].shape, (2,3))
        
    def test_result_right_shape_nonempty_map_3D(self):
        self.mapDict[(0,0)] = [((0,0), None), ((10,34), None)]
        self.mapDict[(1,1)] = [((3,0), None), ((19, 5), None), ((0,1), None)] 
        resDict = self.defOutFunc(self.mapDict, self.six_el_grid,
                                  self.outFname, verbose=False,
                                  version = self.version)
        self.assertEqual(resDict['outTest3D'].shape, (2,3,4))
    
    def test_single_valid_pixel_2D(self):
        self.cfrac[0,30] = .2
        self.qualFlag[0,30] = 0
        self.solZenAng[0,30] = 30
        self.time[0,30] = self.toTAI93('08:00:00 08-30-2011')
        self.lon[0,30] = 13
        data = numpy.random.rand()
        self.test2D[0,30] = data
        self.mapDict[(0,0)] = [((0,30), None)]
        resDict = self.defOutFunc(self.mapDict, self.one_el_grid,
                                  self.outFname, verbose=False,
                                  version = self.version)
        self.assertAlmostEqual(resDict['outTest2D'][0,0], data)
        
    def test_single_valid_pixel_3D(self):
        self.cfrac[0,30] = .2
        self.qualFlag[0,30] = 0
        self.solZenAng[0,30] = 30
        self.time[0,30] = self.toTAI93('08:00:00 08-30-2011')
        self.lon[0,30] = 13
        data = numpy.random.rand(4)
        self.test3D[0,30,:] = data
        self.mapDict[(0,0)] = [((0,30), None)]
        resDict = self.defOutFunc(self.mapDict, self.one_el_grid,
                                  self.outFname, verbose=False,
                                  version = self.version)
        numpy.testing.assert_array_almost_equal(resDict['outTest3D'][0,0,:], data[:])
        
    def test_correct_pixel_count_single_valid_pixel(self):
        self.defParms['includePixelCount'] = True
        newOutFunc = out_geo.OMNO2e_netCDF_avg_out_func(self.defParms)
        self.cfrac[0,30] = .2
        self.qualFlag[0,30] = 0
        self.solZenAng[0,30] = 30
        self.time[0,30] = self.toTAI93('08:00:00 08-30-2011')
        self.lon[0,30] = 13
        data = numpy.random.rand()
        self.test2D[0,30] = data
        self.mapDict[(0,0)] = [((0,30), None)]
        resDict = newOutFunc(self.mapDict, self.one_el_grid,
                             self.outFname, verbose=False,
                             version = self.version)
        self.assertEqual(resDict['ValidPixelCount'][0,0], 1)

    def test_correct_pixel_count_single_invalid_pixel(self):
        self.defParms['includePixelCount'] = True
        newOutFunc = out_geo.OMNO2e_netCDF_avg_out_func(self.defParms)
        self.cfrac[0,30] = .2
        self.qualFlag[0,30] = -1
        self.solZenAng[0,30] = 30
        self.time[0,30] = self.toTAI93('08:00:00 08-30-2011')
        self.lon[0,30] = 13
        data = numpy.random.rand()
        self.test2D[0,30] = data
        self.mapDict[(0,0)] = [((0,30), None)]
        resDict = newOutFunc(self.mapDict, self.one_el_grid,
                             self.outFname, verbose=False,
                             version = self.version)
        self.assertEqual(resDict['ValidPixelCount'][0,0], 0)

    def test_correct_pixel_count_multiple_valid_pixels(self):
        self.defParms['includePixelCount'] = True
        newOutFunc = out_geo.OMNO2e_netCDF_avg_out_func(self.defParms)        
        self.cfrac[0,29:31] = [.2, .15]
        self.qualFlag[0,29:31] = [0, 0]
        self.solZenAng[0,29:31] = [30, 25]
        self.time[0,29:31] = [self.toTAI93('11:14:38 08-30-2011'),
                              self.toTAI93('16:11:00 08-30-2011')]
        self.lon[0,29:31] = [-48, 67]
        data = numpy.random.rand(2)
        self.test2D[0,29:31] = data
        self.mapDict[(0,0)] = [((0,29), None), ((0,30), None)]
        resDict = newOutFunc(self.mapDict, self.one_el_grid,
                             self.outFname, verbose=False,
                             version = self.version)
        self.assertEqual(resDict['ValidPixelCount'][0,0], 2)

    def test_2_valid_pix_accurate_weight_2D(self):
        # calculation of weight performed by hand
        self.cfrac[0,29:31] = [.2, .15]
        self.qualFlag[0,29:31] = [0, 0]
        self.solZenAng[0,29:31] = [30, 25]
        self.time[0,29:31] = [self.toTAI93('11:14:38 08-30-2011'),
                              self.toTAI93('16:11:00 08-30-2011')]
        self.lon[0,29:31] = [-48, 67]
        data = numpy.random.rand(2)
        self.test2D[0,29:31] = data
        self.mapDict[(0,0)] = [((0,29), None), ((0,30), None)]
        resDict = self.defOutFunc(self.mapDict, self.one_el_grid, 
                                  self.outFname, verbose=False,
                                  version = self.version)
        weight = .8212822960 # weight for first element of data
        expected = (data[0]*weight+data[1])/(weight+1)
        self.assertAlmostEqual(resDict['outTest2D'][0,0], expected, delta=10**-7)
        
    def test_2_valid_pix_accurate_weight_3D(self):
        # calculation fo weight performed by hand
        self.cfrac[0,29:31] = [.2, .15]
        self.qualFlag[0,29:31] = [0, 0]
        self.solZenAng[0,29:31] = [30, 25]
        self.time[0,29:31] = [self.toTAI93('11:14:38 08-30-2011'),
                              self.toTAI93('16:11:00 08-30-2011')]
        self.lon[0,29:31] = [-48, 67]
        data = numpy.random.rand(2,4)
        self.test3D[0,29:31,:] = data
        self.mapDict[(0,0)] = [((0,29), None), ((0,30), None)]
        resDict = self.defOutFunc(self.mapDict, self.one_el_grid, 
                                  self.outFname, verbose=False,
                                  version = self.version)
        weight = .8212822960
        expected = (data[0,:]*weight+data[1,:])/(weight+1)
        numpy.testing.assert_array_almost_equal(resDict['outTest3D'][0,0], expected, decimal=7)

    def test_zero_weight_if_sum_flag_is_set_2D(self):
        self.cfrac[0,29:31] = [.2, .15]
        self.qualFlag[0,29:31] = [-1, 0]
        self.solZenAng[0,29:31] = [10, 10]
        self.time[0,29:31] = [self.toTAI93('11:14:38 08-30-2011'),
                              self.toTAI93('16:11:00 08-30-2011')]
        self.lon[0,29:31] = [-48, 67]
        data = numpy.random.rand(2)
        self.test2D[0,29:31] = data
        self.mapDict[(0,0)] = [((0,29), None), ((0,30), None)]
        resDict = self.defOutFunc(self.mapDict, self.one_el_grid,
                                  self.outFname, verbose=False,
                                  version = self.version)
        self.assertAlmostEqual(resDict['outTest2D'][0,0], data[1])
        
    def test_zero_weight_if_sum_flag_is_set_3D(self):
        self.cfrac[0,29:31] = [.2, .15]
        self.qualFlag[0,29:31] = [-127,0]
        self.solZenAng[0,29:31] = [10, 10]
        self.time[0,29:31] = [self.toTAI93('11:14:38 08-30-2011'),
                              self.toTAI93('16:11:00 08-30-2011')]
        self.lon[0,29:31] = [-48, 67]
        data = numpy.random.rand(2,4)
        self.test3D[0,29:31,:] = data
        self.mapDict[(0,0)] = [((0,29), None), ((0,30), None)]
        resDict = self.defOutFunc(self.mapDict, self.one_el_grid,
                                  self.outFname, verbose=False,
                                  version = self.version)
        numpy.testing.assert_array_almost_equal(resDict['outTest3D'][0,0], data[1,:])
        
    def test_zero_weight_if_cfrac_gt_threshold_2D(self):
        self.cfrac[0,29:31] = [.28, .24]
        self.qualFlag[0,29:31] = [0,0]
        self.solZenAng[0,29:31] = [11.1, 13.4]
        self.time[0,29:31] = [self.toTAI93('11:14:38 08-30-2011'),
                              self.toTAI93('16:11:00 08-30-2011')]
        self.lon[0,29:31] = [-48, 67]
        data = numpy.random.rand(2)
        self.test2D[0,29:31] = data
        self.mapDict[(0,0)] = [((0,29), None), ((0,30), None)]
        resDict = self.defOutFunc(self.mapDict, self.one_el_grid, 
                                  self.outFname, verbose=False,
                                  version = self.version)
        self.assertAlmostEqual(resDict['outTest2D'][0,0], data[1])
        
    def test_zero_weight_if_cfrac_gt_threshold_3D(self):
        self.cfrac[0,29:31] = [.28, .24]
        self.qualFlag[0,29:31] = [0,0]
        self.solZenAng[0,29:31] = [19, 58]
        self.time[0,29:31] = [self.toTAI93('11:14:38 08-30-2011'),
                              self.toTAI93('16:11:00 08-30-2011')]
        self.lon[0,29:31] = [-48, 67]
        data = numpy.random.rand(2,4)
        self.test3D[0,29:31,:] = data
        self.mapDict[(0,0)] = [((0,29), None), ((0,30), None)]
        resDict = self.defOutFunc(self.mapDict, self.one_el_grid,
                                  self.outFname, verbose=False,
                                  version = self.version)
        numpy.testing.assert_array_almost_equal(resDict['outTest3D'][0,0,:], data[1,:])
        
    def test_zero_weight_if_SZA_gt_threshold_2D(self):
        self.cfrac[0,29:31] = [.1, .1]
        self.qualFlag[0,29:31] = [0, 0]
        self.solZenAng[0,29:31] = [80.2, 79.79]
        self.time[0,29:31] = [self.toTAI93('11:14:38 08-30-2011'),
                              self.toTAI93('16:11:00 08-30-2011')]
        self.lon[0,29:31] = [-48, 67]
        data = numpy.random.rand(2)
        self.test2D[0,29:31] = data
        self.mapDict[(0,0)] = [((0,29), None), ((0,30), None)]
        resDict = self.defOutFunc(self.mapDict, self.one_el_grid,
                                  self.outFname, verbose=False,
                                  version = self.version)
        self.assertAlmostEqual(resDict['outTest2D'][0,0], data[1])

    def test_zero_weight_if_SZA_gt_threshold_3D(self):
        self.cfrac[0,29:31] = [.1, .1]
        self.qualFlag[0,29:31] = [0, 0]
        self.solZenAng[0,29:31] = [80.2, 79.79]
        self.time[0,29:31] = [self.toTAI93('11:14:38 08-30-2011'),
                              self.toTAI93('16:11:00 08-30-2011')]
        self.lon[0,29:31] = [-48, 67]
        data = numpy.random.rand(2,4)
        self.test3D[0,29:31,:] = data
        self.mapDict[(0,0)] = [((0,29), None), ((0,30), None)]
        resDict = self.defOutFunc(self.mapDict, self.one_el_grid, 
                                  self.outFname, verbose=False,
                                  version = self.version)
        numpy.testing.assert_array_almost_equal(resDict['outTest3D'][0,0,:], data[1,:])
        
    def test_zero_weight_UTC_time_wrong_day_2D(self):
        self.cfrac[0,29:31] = [.1, .1]
        self.qualFlag[0,29:31] = [0, 0]
        self.solZenAng[0,29:31] = [10, 10]
        self.time[0,29:31] = [self.toTAI93('12:01:10 08-15-2011'),
                              self.toTAI93('12:01:10 08-30-2011')]
        self.lon[0,29:31] = [1, -1]
        data = numpy.random.rand(2)
        self.test2D[0,29:31] = data
        self.mapDict[(0,0)] = [((0,29), None), ((0,30), None)]
        resDict = self.defOutFunc(self.mapDict, self.one_el_grid, 
                                  self.outFname, verbose=False,
                                  version = self.version)
        self.assertAlmostEqual(resDict['outTest2D'][0,0], data[1])
        
    def test_zero_weight_UTC_time_wrong_day_3D(self):
        self.cfrac[0,29:31] = [.1, .1]
        self.qualFlag[0,29:31] = [0, 0]
        self.solZenAng[0,29:31] = [10, 10]
        self.time[0,29:31] = [self.toTAI93('12:01:10 08-15-2011'),
                              self.toTAI93('12:01:10 08-30-2011')]
        self.lon[0,29:31] = [1, -1]
        data = numpy.random.rand(2,4)
        self.test3D[0,29:31,:] = data
        self.mapDict[(0,0)] = [((0,29), None), ((0,30), None)]
        resDict = self.defOutFunc(self.mapDict, self.one_el_grid, 
                                  self.outFname, verbose=False,
                                  version = self.version)
        numpy.testing.assert_array_almost_equal(resDict['outTest3D'][0,0,:], data[1,:])
        
    def test_zero_weight_local_time_wrong_day_2D(self):
        self.cfrac[0,29:31] = [.1, .1]
        self.qualFlag[0,29:31] = [0, 0]
        self.solZenAng[0,29:31] = [10, 10]
        self.time[0,29:31] = [self.toTAI93('12:01:10 08-15-2011'),
                              self.toTAI93('12:01:10 08-30-2011')]
        self.lon[0,29:31] = [1, -1]
        data = numpy.random.rand(2)
        self.test2D[0,29:31] = data
        self.mapDict[(0,0)] = [((0,29), None), ((0,30), None)]
        newParms = self.defParms
        newParms['timeComparison'] = 'local'
        newOutFunc = out_geo.OMNO2e_netCDF_avg_out_func(newParms)
        resDict = newOutFunc(self.mapDict, self.one_el_grid, 
                             self.outFname, verbose=False,
                             version = self.version)
        self.assertAlmostEqual(resDict['outTest2D'][0,0], data[1])
        
    def test_zero_weight_local_time_wrong_day_3D(self):
        self.cfrac[0,29:31] = [.1, .1]
        self.qualFlag[0,29:31] = [0, 0]
        self.solZenAng[0,29:31] = [10, 10]
        self.time[0,29:31] = [self.toTAI93('12:01:10 08-15-2011'),
                              self.toTAI93('12:01:10 08-30-2011')]
        self.lon[0,29:31] = [1, -1]
        data = numpy.random.rand(2,4)
        self.test3D[0,29:31,:] = data
        self.mapDict[(0,0)] = [((0,29), None), ((0,30), None)]
        newParms = self.defParms
        newParms['timeComparison'] = 'local'
        newOutFunc = out_geo.OMNO2e_netCDF_avg_out_func(newParms)
        resDict = newOutFunc(self.mapDict, self.one_el_grid, 
                             self.outFname, verbose=False,
                             version = self.version)
        numpy.testing.assert_array_almost_equal(resDict['outTest3D'][0,0,:], data[1,:])
        
    def test_zero_weight_UTC_time_wrong_year_2D(self):
        self.cfrac[0,29:31] = [.1,.1]
        self.qualFlag[0,29:31] = [0,0]
        self.solZenAng[0,29:31] = [4, 5]
        self.time[0,29:31] = [self.toTAI93('12:01:10 08-30-1911'),
                              self.toTAI93('12:01:10 08-30-2011')]
        self.lon[0,29:31] = [-4.5, -1.2]
        data = numpy.random.rand(2)
        self.test2D[0,29:31] = data
        self.mapDict[(0,0)] = [((0,29), None), ((0,30), None)]
        resDict = self.defOutFunc(self.mapDict, self.one_el_grid, 
                                  self.outFname, verbose=False,
                                  version = self.version)
        self.assertAlmostEqual(resDict['outTest2D'][0,0], data[1])
        
    def test_zero_weight_UTC_time_wrong_year_3D(self):
        self.cfrac[0,29:31] = [.1,.1]
        self.qualFlag[0,29:31] = [0,0]
        self.solZenAng[0,29:31] = [4, 5]
        self.time[0,29:31] = [self.toTAI93('12:01:10 08-30-1911'),
                              self.toTAI93('12:01:10 08-30-2011')]
        self.lon[0,29:31] = [-4.5, -1.2]
        data = numpy.random.rand(2,4)
        self.test3D[0,29:31,:] = data
        self.mapDict[(0,0)] = [((0,29), None), ((0,30), None)]        
        resDict = self.defOutFunc(self.mapDict, self.one_el_grid, 
                                  self.outFname, verbose=False,
                                  version = self.version)
        numpy.testing.assert_array_almost_equal(resDict['outTest3D'][0,0,:], data[1,:])
        
    def test_zero_weight_local_time_wrong_year_2D(self):
        self.cfrac[0,29:31] = [.1,.1]
        self.qualFlag[0,29:31] = [0,0]
        self.solZenAng[0,29:31] = [4, 5]
        self.time[0,29:31] = [self.toTAI93('12:01:10 08-30-1911'),
                              self.toTAI93('12:01:10 08-30-2011')]
        self.lon[0,29:31] = [-4.5, -1.2]
        data = numpy.random.rand(2)
        self.test2D[0,29:31] = data
        self.mapDict[(0,0)] = [((0,29), None), ((0,30), None)]
        newParms = self.defParms
        newParms['timeComparison'] = 'local'
        newOutFunc = out_geo.OMNO2e_netCDF_avg_out_func(newParms)        
        resDict = newOutFunc(self.mapDict, self.one_el_grid, 
                             self.outFname, verbose=False,
                             version = self.version)
        self.assertAlmostEqual(resDict['outTest2D'][0,0], data[1])

    def test_zero_weight_local_time_wrong_year_3D(self):
        self.cfrac[0,29:31] = [.1,.1]
        self.qualFlag[0,29:31] = [0,0]
        self.solZenAng[0,29:31] = [4, 5]
        self.time[0,29:31] = [self.toTAI93('12:01:10 08-30-1911'),
                              self.toTAI93('12:01:10 08-30-2011')]
        self.lon[0,29:31] = [-4.5, -1.2]
        data = numpy.random.rand(2,4)
        self.test3D[0,29:31,:] = data
        self.mapDict[(0,0)] = [((0,29), None), ((0,30), None)]        
        newParms = self.defParms
        newParms['timeComparison'] = 'local'
        newOutFunc = out_geo.OMNO2e_netCDF_avg_out_func(newParms)        
        resDict = newOutFunc(self.mapDict, self.one_el_grid, 
                             self.outFname, verbose=False,
                             version = self.version)
        numpy.testing.assert_array_almost_equal(resDict['outTest3D'][0,0,:], data[1,:])
        
    def test_UTC_time_does_not_incorporate_lon(self):
        self.cfrac[0,29:31] = [.1, .2]
        self.qualFlag[0,29:31] = [0, 0]
        self.solZenAng[0,29:31] = [43, 44]
        self.time[0,29:31] = [self.toTAI93('02:24:01 08-30-2011'),
                              self.toTAI93('22:24:01 08-29-2011')]
        self.lon[0,29:31] = [-150, 150]
        data = numpy.random.rand(2)
        self.test2D[0,29:31] = data
        self.mapDict[(0,0)] = [((0,29), None), ((0,30), None)]
        resDict = self.defOutFunc(self.mapDict, self.one_el_grid, 
                                  self.outFname, verbose=False,
                                  version = self.version)
        self.assertAlmostEqual(resDict['outTest2D'][0,0], data[0])
        
    def test_local_time_does_incorporate_lon(self):
        self.cfrac[0,29:31] = [.1, .2]
        self.qualFlag[0,29:31] = [0, 0]
        self.solZenAng[0,29:31] = [43, 44]
        self.time[0,29:31] = [self.toTAI93('02:24:01 08-30-2011'),
                              self.toTAI93('22:24:01 08-29-2011')]
        self.lon[0,29:31] = [-150, 150]
        data = numpy.random.rand(2)
        self.test2D[0,29:31] = data
        self.mapDict[(0,0)] = [((0,29), None), ((0,30), None)]
        newParms = self.defParms
        newParms['timeComparison'] = 'local'
        newOutFunc = out_geo.OMNO2e_netCDF_avg_out_func(newParms)
        resDict = newOutFunc(self.mapDict, self.one_el_grid, 
                             self.outFname, verbose=False,
                             version = self.version)
        self.assertAlmostEqual(resDict['outTest2D'][0,0], data[1])        
        
    def test_UTC_time_does_not_incorporate_lon_360(self):
        self.cfrac[0,29:31] = [.1, .2]
        self.qualFlag[0,29:31] = [0, 0]
        self.solZenAng[0,29:31] = [43, 44]
        self.time[0,29:31] = [self.toTAI93('02:24:01 08-30-2011'),
                              self.toTAI93('22:24:01 08-29-2011')]
        self.lon[0,29:31] = [210, 150]
        data = numpy.random.rand(2)
        self.test2D[0,29:31] = data
        self.mapDict[(0,0)] = [((0,29), None), ((0,30), None)]
        resDict = self.defOutFunc(self.mapDict, self.one_el_grid, 
                                  self.outFname, verbose=False,
                                  version = self.version)
        self.assertAlmostEqual(resDict['outTest2D'][0,0], data[0])
                
    def test_local_time_does_incorporate_lon_360(self):
        self.cfrac[0,29:31] = [.1, .2]
        self.qualFlag[0,29:31] = [0, 0]
        self.solZenAng[0,29:31] = [43, 44]
        self.time[0,29:31] = [self.toTAI93('02:24:01 08-30-2011'),
                              self.toTAI93('22:24:01 08-29-2011')]
        self.lon[0,29:31] = [210, 150]
        data = numpy.random.rand(2)
        self.test2D[0,29:31] = data
        self.mapDict[(0,0)] = [((0,29), None), ((0,30), None)]
        newParms = self.defParms
        newParms['timeComparison'] = 'local'
        newOutFunc = out_geo.OMNO2e_netCDF_avg_out_func(newParms)
        resDict = newOutFunc(self.mapDict, self.one_el_grid, 
                             self.outFname, verbose=False,
                             version = self.version)
        self.assertAlmostEqual(resDict['outTest2D'][0,0], data[1])
                        
    def test_zero_weight_if_mult_problems_2D(self):
        self.cfrac[0,28:37] = [.4, .5, .1, .5, .4, .5, .1, .5, .1]
        self.qualFlag[0,28:37] = [1, 0, 1, 1, 1, 0, 1, 1, 0]
        self.solZenAng[0,28:37] = [40, 88, 88, 88, 40, 88, 88, 88, 40]
        self.time[0, 28:37] = [self.toTAI93('07:59:59 08-30-2011'),
                               self.toTAI93('07:59:59 08-30-2011'),
                               self.toTAI93('07:59:59 08-30-2011'),
                               self.toTAI93('07:59:59 08-30-2011'),
                               self.toTAI93('09:59:59 09-21-2011'),
                               self.toTAI93('09:59:59 09-21-2011'),
                               self.toTAI93('09:59:59 09-21-2011'),
                               self.toTAI93('09:59:59 09-21-2011'),
                               self.toTAI93('07:59:59 08-30-2011')]
        self.lon[0, 28:37] = [-1, 2, 3, 0, -1, 0, 3, 4, 1]
        data = numpy.random.rand(9)
        self.test2D[0,28:37] = data
        self.mapDict[(0,0)] = [((0,28), None), ((0,29), None), 
                               ((0,30), None), ((0,31), None),
                               ((0,32), None), ((0,33), None),
                               ((0,34), None), ((0,35), None),
                               ((0,36), None)]
        resDict = self.defOutFunc(self.mapDict, self.one_el_grid, 
                                  self.outFname, verbose=False,
                                  version = self.version)
        self.assertAlmostEqual(resDict['outTest2D'][0,0], data[8])
        
    def test_zero_weight_if_mult_problems_3D(self):
        self.cfrac[0,28:37] = [.4, .5, .1, .5, .4, .5, .1, .5, .1]
        self.qualFlag[0,28:37] = [1, 0, 1, 1, 1, 0, 1, 1, 0]
        self.solZenAng[0,28:37] = [40, 88, 88, 88, 40, 88, 88, 88, 40]
        self.time[0, 28:37] = [self.toTAI93('07:59:59 08-30-2011'),
                               self.toTAI93('07:59:59 08-30-2011'),
                               self.toTAI93('07:59:59 08-30-2011'),
                               self.toTAI93('07:59:59 08-30-2011'),
                               self.toTAI93('09:59:59 09-21-2011'),
                               self.toTAI93('09:59:59 09-21-2011'),
                               self.toTAI93('09:59:59 09-21-2011'),
                               self.toTAI93('09:59:59 09-21-2011'),
                               self.toTAI93('07:59:59 08-30-2011')]
        self.lon[0, 28:37] = [-1, 2, 3, 0, -1, 0, 3, 4, 1]
        data = numpy.random.rand(9,4)
        self.test3D[0,28:37,:] = data
        self.mapDict[(0,0)] = [((0,28), None), ((0,29), None), 
                               ((0,30), None), ((0,31), None),
                               ((0,32), None), ((0,33), None),
                               ((0,34), None), ((0,35), None),
                               ((0,36), None)]
        resDict = self.defOutFunc(self.mapDict, self.one_el_grid, 
                                  self.outFname, verbose=False,
                                  version = self.version)
        numpy.testing.assert_array_almost_equal(resDict['outTest3D'][0,0,:], data[8,:])
        
    def test_all_zero_weight_yields_fillVal_2D(self):
        self.cfrac[0,29:31] = [.1, .1]
        self.qualFlag[0,29:31] = [1, 1]
        self.solZenAng[0,29:31] = [23, 14]
        self.time[0,29:31] = [self.toTAI93('15:28:53 08-30-2011'),
                              self.toTAI93('15:28:53 08-30-2011')]
        self.lon[0,29:31] = [-7.5, 7.5]
        data = numpy.random.rand(2)
        self.test2D[0,29:31] = data
        self.mapDict[(0,0)] = [((0,29), None), ((0,30), None)]
        resDict = self.defOutFunc(self.mapDict, self.one_el_grid, 
                                  self.outFname, verbose=False,
                                  version = self.version)
        expected = self.defParms['fillVal']
        self.assertAlmostEqual(resDict['outTest2D'][0,0], expected)
        
    def test_all_zero_weight_yields_fillVal_3D(self):
        self.cfrac[0,29:31] = [.1, .1]
        self.qualFlag[0,29:31] = [1, 1]
        self.solZenAng[0,29:31] = [23, 14]
        self.time[0,29:31] = [self.toTAI93('15:28:53 08-30-2011'),
                              self.toTAI93('15:28:53 08-30-2011')]
        self.lon[0,29:31] = [-7.5, 7.5]
        data = numpy.random.rand(2,4)
        self.test3D[0,29:31,:] = data
        self.mapDict[(0,0)] = [((0,29), None), ((0,30), None)]
        resDict = self.defOutFunc(self.mapDict, self.one_el_grid, 
                                  self.outFname, verbose=False,
                                  version = self.version)
        expected = numpy.array(4*[self.defParms['fillVal']])
        numpy.testing.assert_array_almost_equal(resDict['outTest3D'][0,0,:], expected)
        
    def test_NAN_in_data_and_zero_weight_yields_fillVal_2D(self):
        self.cfrac[0,29:31] = [.1, .1]
        self.qualFlag[0,29:31] = [1, 1]
        self.solZenAng[0,29:31] = [40, 40]
        self.time[0,29:31] = [self.toTAI93('15:28:12 08-30-2011'),
                              self.toTAI93('15:28:12 08-30-2011')] 
        self.lon[0,29:31] = [0, 0.3]
        self.test2D[0,29:31] = numpy.NaN
        self.mapDict[(0,0)] = [((0,29), None), ((0,30), None)]
        resDict = self.defOutFunc(self.mapDict, self.one_el_grid, 
                                  self.outFname, verbose=False,
                                  version = self.version)
        expected = self.defParms['fillVal']
        self.assertAlmostEqual(resDict['outTest2D'][0,0], expected)
        
    def test_NAN_in_data_and_zero_weight_yields_fillVal_3D(self):
        self.cfrac[0,29:31] = [.1, .1]
        self.qualFlag[0,29:31] = [1, 1]
        self.solZenAng[0,29:31] = [40, 40]
        self.time[0,29:31] = [self.toTAI93('15:28:12 08-30-2011'),
                              self.toTAI93('15:28:12 08-30-2011')] 
        self.lon[0,29:31] = [0, 0.3]
        self.test3D[0,29:31,:] = numpy.NaN
        self.mapDict[(0,0)] = [((0,29), None), ((0,30), None)]
        resDict = self.defOutFunc(self.mapDict, self.one_el_grid, 
                                  self.outFname, verbose=False,
                                  version = self.version)
        expected = numpy.array(4*[self.defParms['fillVal']])
        numpy.testing.assert_array_almost_equal(resDict['outTest3D'][0,0,:], expected)
        
    def test_NAN_in_data_with_zero_weight_does_not_affect_avg_2D(self):
        self.cfrac[0,29:31] = [.1, .1]
        self.qualFlag[0,29:31] = [1, 0]
        self.solZenAng[0,29:31] = [40, 40]
        self.time[0,29:31] = [self.toTAI93('15:28:12 08-30-2011'),
                              self.toTAI93('15:28:12 08-30-2011')] 
        self.lon[0,29:31] = [0, 0.3]
        self.test2D[0,29] = numpy.NaN
        data = numpy.random.rand()
        self.test2D[0,30] = data
        self.mapDict[(0,0)] = [((0,29), None), ((0,30), None)]
        resDict = self.defOutFunc(self.mapDict, self.one_el_grid, 
                                  self.outFname, verbose=False,
                                  version = self.version)
        self.assertAlmostEqual(resDict['outTest2D'][0,0], data)
        
    def test_NAN_in_data_with_zero_weight_does_not_affect_avg_3D(self):
        self.cfrac[0,29:31] = [.1, .1]
        self.qualFlag[0,29:31] = [1, 0]
        self.solZenAng[0,29:31] = [40, 40]
        self.time[0,29:31] = [self.toTAI93('16:12:19 08-30-2011'),
                              self.toTAI93('16:12:19 08-30-2011')]
        self.lon[0,29:31] = [-2.3, 6.7]
        self.test3D[0,29,:] = numpy.NaN
        data = numpy.random.rand(4)
        self.test3D[0,30,:] = data
        self.mapDict[(0,0)] = [((0,29), None), ((0,30), None)]
        resDict = self.defOutFunc(self.mapDict, self.one_el_grid, 
                                  self.outFname, verbose=False,
                                  version = self.version)
        numpy.testing.assert_array_almost_equal(resDict['outTest3D'][0,0,:], data)
        
    def test_no_pixels_in_cell_yields_correct_count(self):
        self.defParms['includePixelCount'] = True
        newOutFunc = out_geo.OMNO2e_netCDF_avg_out_func(self.defParms)        
        self.mapDict[(0,0)] = []
        resDict = newOutFunc(self.mapDict, self.one_el_grid, 
                             self.outFname, verbose=False,
                             version = self.version)
        self.assertEqual(resDict['ValidPixelCount'][0,0], 0)

    def test_no_pixels_in_cell_yields_fillVal_2D(self):
        self.mapDict[(0,0)] = []
        resDict = self.defOutFunc(self.mapDict, self.one_el_grid, 
                                  self.outFname, verbose=False,
                                  version = self.version)
        self.assertEqual(resDict['outTest2D'][0,0], self.defParms['fillVal'])
        
    def test_no_pixels_in_cell_yields_fillVal_3D(self):
        self.mapDict[(0,0)] = []
        resDict = self.defOutFunc(self.mapDict, self.one_el_grid, 
                                  self.outFname, verbose=False,
                                  version = self.version)
        expected = numpy.array(4*[self.defParms['fillVal']])
        numpy.testing.assert_array_equal(resDict['outTest3D'][0,0,:], expected)
        
    def test_multi_element_pix_array_correct_counts(self):
        self.defParms['includePixelCount'] = True
        newOutFunc = out_geo.OMNO2e_netCDF_avg_out_func(self.defParms)
        self.cfrac[0:2, 28:34] = [[.11, .12, .13, .14, .15, .16],
                                  [.21, .22, 23, .24, .25, .26]]
        self.qualFlag[0:2, 28:34] = [[0, 0, 0, 0, 0, 0],
                                     [0, 0, 0, 1, 1, 1]]
        self.solZenAng[0:2, 28:34] = [[10, 11, 12, 13, 14, 15],
                                      [87, 88, 89, 20, 12, 13]]
        self.time[0:2, 28:34] = 2*[[self.toTAI93('17:20:00 08-30-2011'),
                                    self.toTAI93('17:20:00 08-30-2011'),
                                    self.toTAI93('17:20:00 08-30-2011'),
                                    self.toTAI93('17:20:00 08-30-2011'),
                                    self.toTAI93('17:20:00 08-30-2011'),
                                    self.toTAI93('17:20:00 08-30-2011')]]
        self.lon[0:2, 28:34] = [[1, 2, 3, 4, 5, 6],
                                [1, 2, 3, 4, 5, 7]]
        data = numpy.random.rand(2,6)
        self.test2D[0:2,28:34] = data
        for (i,j) in product(range(2), range(3)):
            oneDind = i*3+j
            self.mapDict[(i,j)] = [((0,28+oneDind), None), ((1,28+oneDind), None)]
        resDict = newOutFunc(self.mapDict, self.six_el_grid, 
                             self.outFname, verbose=False,
                             version = self.version)
        expected = numpy.ones((2,3))
        numpy.testing.assert_array_equal(resDict['ValidPixelCount'][:], expected)

    def test_multi_element_pix_array_2D(self):
        self.cfrac[0:2, 28:34] = [[.11, .12, .13, .14, .15, .16],
                                  [.21, .22, 23, .24, .25, .26]]
        self.qualFlag[0:2, 28:34] = [[0, 0, 0, 0, 0, 0],
                                     [0, 0, 0, 1, 1, 1]]
        self.solZenAng[0:2, 28:34] = [[10, 11, 12, 13, 14, 15],
                                      [87, 88, 89, 20, 12, 13]]
        self.time[0:2, 28:34] = 2*[[self.toTAI93('17:20:00 08-30-2011'),
                                    self.toTAI93('17:20:00 08-30-2011'),
                                    self.toTAI93('17:20:00 08-30-2011'),
                                    self.toTAI93('17:20:00 08-30-2011'),
                                    self.toTAI93('17:20:00 08-30-2011'),
                                    self.toTAI93('17:20:00 08-30-2011')]]
        self.lon[0:2, 28:34] = [[1, 2, 3, 4, 5, 6],
                                [1, 2, 3, 4, 5, 7]]
        data = numpy.random.rand(2,6)
        self.test2D[0:2,28:34] = data
        for (i,j) in product(range(2), range(3)):
            oneDind = i*3+j
            self.mapDict[(i,j)] = [((0,28+oneDind), None), ((1,28+oneDind), None)]
        resDict = self.defOutFunc(self.mapDict, self.six_el_grid, 
                                  self.outFname, verbose=False,
                                  version = self.version)
        expected = data[0,:].reshape((2,3))
        numpy.testing.assert_array_almost_equal(resDict['outTest2D'][:], expected)
        
    def test_multi_element_pix_array_3D(self):
        self.cfrac[0:2, 28:34] = [[.11, .12, .13, .14, .15, .16],
                                  [.21, .22, 23, .24, .25, .26]]
        self.qualFlag[0:2, 28:34] = [[0, 0, 0, 0, 0, 0],
                                     [0, 0, 0, 1, 1, 1]]
        self.solZenAng[0:2, 28:34] = [[10, 11, 12, 13, 14, 15],
                                      [87, 88, 89, 20, 12, 13]]
        self.time[0:2, 28:34] = 2*[[self.toTAI93('17:20:00 08-30-2011'),
                                    self.toTAI93('17:20:00 08-30-2011'),
                                    self.toTAI93('17:20:00 08-30-2011'),
                                    self.toTAI93('17:20:00 08-30-2011'),
                                    self.toTAI93('17:20:00 08-30-2011'),
                                    self.toTAI93('17:20:00 08-30-2011')]]
        self.lon[0:2, 28:34] = [[1, 2, 3, 4, 5, 6],
                                [1, 2, 3, 4, 5, 7]]
        data = numpy.random.rand(2,6,4)
        self.test3D[0:2,28:34, :] = data
        for (i,j) in product(range(2), range(3)):
            oneDind = i*3+j
            self.mapDict[(i,j)] = [((0,28+oneDind), None), ((1,28+oneDind), None)]        
        resDict = self.defOutFunc(self.mapDict, self.six_el_grid, 
                                  self.outFname, verbose=False,
                                  version = self.version)
        expected = data[0,:,:].reshape((2,3,4))
        numpy.testing.assert_array_almost_equal(resDict['outTest3D'][:], expected)

    def test_multiple_dictionaries_yield_correct_count(self):
        self.defParms['includePixelCount'] = True
        newOutFunc = out_geo.OMNO2e_netCDF_avg_out_func(self.defParms)
        self.cfrac[0, 29:31] = [.1, .2]
        self.qualFlag[0, 29:31] = [0, -127]
        self.solZenAng[0, 29:31] = [0, 2]
        self.time[0, 29:31] = [self.toTAI93('05:15:00 08-30-2011'),
                               self.toTAI93('05:15:00 08-30-2011')]
        self.lon[0, 29:31] = [-6, 3]
        data = numpy.random.rand(2)
        self.test2D[0,29:31] = data
        secondMapDict = dict(self.mapDict)
        self.mapDict[(0,0)] = [((0,29), None)]
        secondMapDict[(0,0)] = [((0,30), None)]
        dictList = [self.mapDict, secondMapDict]
        resDict = newOutFunc(dictList, self.one_el_grid,
                             self.outFname, verbose=False,
                             version = self.version)
        self.assertEqual(resDict['ValidPixelCount'][0,0], 1)

    def test_multidict_one_zero_weight_2D(self):
        self.cfrac[0, 29:31] = [.1, .2]
        self.qualFlag[0, 29:31] = [0, -127]
        self.solZenAng[0, 29:31] = [0, 2]
        self.time[0, 29:31] = [self.toTAI93('05:15:00 08-30-2011'),
                               self.toTAI93('05:15:00 08-30-2011')]
        self.lon[0, 29:31] = [-6, 3]
        data = numpy.random.rand(2)
        self.test2D[0,29:31] = data
        secondMapDict = dict(self.mapDict)
        self.mapDict[(0,0)] = [((0,29), None)]
        secondMapDict[(0,0)] = [((0,30), None)]
        dictList = [self.mapDict, secondMapDict]
        resDict = self.defOutFunc(dictList, self.one_el_grid, 
                                  self.outFname, verbose=False,
                                  version = self.version)
        self.assertAlmostEqual(resDict['outTest2D'][0,0], data[0])
        
    def test_multidict_one_zero_weight_3D(self):
        self.cfrac[0, 29:31] = [.1, .2]
        self.qualFlag[0, 29:31] = [0, -127]
        self.solZenAng[0, 29:31] = [0, 2]
        self.time[0, 29:31] = [self.toTAI93('05:15:00 08-30-2011'),
                               self.toTAI93('05:15:00 08-30-2011')]
        self.lon[0, 29:31] = [-6, 3]
        data = numpy.random.rand(2,4)
        self.test3D[0,29:31,:] = data
        secondMapDict = dict(self.mapDict)
        self.mapDict[(0,0)] = [((0,29), None)]
        secondMapDict[(0,0)] = [((0,30), None)]
        dictList = [self.mapDict, secondMapDict]
        resDict = self.defOutFunc(dictList, self.one_el_grid, 
                                  self.outFname, verbose=False,
                                  version = self.version)
        numpy.testing.assert_array_almost_equal(resDict['outTest3D'][0,0,:], data[0,:])
        
    def test_multidict_both_with_nonzero_weight_2D(self):
        # calculation of weight performed by hand
        self.cfrac[0,0:2] = [.05, .10]
        self.qualFlag[0,0:2] = [0, 0]
        self.solZenAng[0,0:2] = [50, 55]
        self.time[0,0:2] = [self.toTAI93('11:14:38 08-30-2011'),
                              self.toTAI93('16:11:00 08-30-2011')]
        self.lon[0,0:2] = [-48, 67]
        data = numpy.random.rand(2)
        self.test2D[0,0:2] = data
        secondMapDict = dict(self.mapDict)
        self.mapDict[(0,0)] = [((0,0), None)]
        secondMapDict[(0,0)] = [((0,1), None)]
        dictList = [self.mapDict, secondMapDict]
        resDict = self.defOutFunc(dictList, self.one_el_grid, 
                                  self.outFname, verbose=False,
                                  version = self.version)
        weight = 1.088125096 # calculated by hand
        expected = (data[0]*weight+data[1])/(weight+1)
        self.assertAlmostEqual(resDict['outTest2D'][0,0], expected)
        
    def test_multidict_both_with_nonzero_weight_3D(self):
        # calculation of weight performed by hand
        self.cfrac[0,0:2] = [.05, .10]
        self.qualFlag[0,0:2] = [0, 0]
        self.solZenAng[0,0:2] = [50, 55]
        self.time[0,0:2] = [self.toTAI93('11:14:38 08-30-2011'),
                              self.toTAI93('16:11:00 08-30-2011')]
        self.lon[0,0:2] = [-48, 67]
        data = numpy.random.rand(2,4)
        self.test3D[0,0:2] = data
        secondMapDict = dict(self.mapDict)
        self.mapDict[(0,0)] = [((0,0), None)]
        secondMapDict[(0,0)] = [((0,1), None)]
        dictList = [self.mapDict, secondMapDict]
        resDict = self.defOutFunc(dictList, self.one_el_grid, 
                                  self.outFname, verbose=False,
                                  version = self.version)
        weight = 1.088125096 # calculated by hand
        expected = (data[0,:]*weight+data[1,:])/(weight+1)
        numpy.testing.assert_array_almost_equal(resDict['outTest3D'][0,0,:], expected[:])
        
    def test_zero_weight_pix_in_multiple_cells_2D(self):
        self.cfrac[0, 29:32] = [.2, .01, 0]
        self.qualFlag[0, 29:32] = [0, 0, 0]
        self.solZenAng[0, 29:32] = [30, 25, 15]
        self.time[0, 29:32] = [self.toTAI93('08:00:00 08-30-2011'),
                               self.toTAI93('23:45:10 08-29-2011'),
                               self.toTAI93('10:00:00 08-30-2011')]
        self.lon[0, 29:32] = [0, 1, -1]
        data = numpy.random.rand(3)
        self.test2D[0, 29:32] = data
        self.mapDict[(1,1)] = [((0,29), None), ((0,30), None)]
        self.mapDict[(1,2)] = [((0,30), None), ((0,31), None)]
        resDict = self.defOutFunc(self.mapDict, self.six_el_grid, 
                                  self.outFname, verbose=False,
                                  version = self.version)
        fv = self.defParms['fillVal']
        expected = numpy.array([[fv, fv, fv],[fv, data[0], data[2]]])
        numpy.testing.assert_array_almost_equal(resDict['outTest2D'][:], expected)
        
    def test_zero_weight_pix_in_multiple_cells_3D(self):
        self.cfrac[0, 29:32] = [.2, .01, 0]
        self.qualFlag[0, 29:32] = [0, 0, 0]
        self.solZenAng[0, 29:32] = [30, 25, 15]
        self.time[0, 29:32] = [self.toTAI93('08:00:00 08-30-2011'),
                               self.toTAI93('23:45:10 08-29-2011'),
                               self.toTAI93('10:00:00 08-30-2011')]
        self.lon[0, 29:32] = [0, 1, -1]
        data = numpy.random.rand(3,4)
        self.test3D[0,29:32,:] = data
        self.mapDict[(1,1)] = [((0,29), None), ((0,30), None)]
        self.mapDict[(1,2)] = [((0,30), None), ((0,31), None)]
        resDict = self.defOutFunc(self.mapDict, self.six_el_grid, 
                                  self.outFname, verbose=False,
                                  version = self.version)
        expected = numpy.empty((2,3,4))
        expected[:] = self.defParms['fillVal']
        expected[1,1,:] = data[0,:]
        expected[1,2,:] = data[2,:]
        numpy.testing.assert_array_almost_equal(resDict['outTest3D'][:], expected)
        
    def test_nonzero_weight_pix_in_multiple_cells_2D(self):
        self.cfrac[0, 29:32] = [.5, .12, 1.0]
        self.qualFlag[0, 29:32] = [0, 0, 0]
        self.solZenAng[0, 29:32] = [30, 20, 10]
        self.time[0, 29:32] = [self.toTAI93('23:45:10 08-30-2011'),
                               self.toTAI93('23:45:10 08-30-2011'),
                               self.toTAI93('23:45:10 08-30-2011')]
        self.lon[0, 29:32] = [0, 0, 0]
        data = numpy.random.rand(3)
        self.test2D[0, 29:32] = data
        self.mapDict[(1,1)] = [((0,29), None), ((0,30), None)]
        self.mapDict[(1,2)] = [((0,30), None), ((0,31), None)]
        resDict = self.defOutFunc(self.mapDict, self.six_el_grid, 
                                  self.outFname, verbose=False,
                                  version = self.version)
        fv = self.defParms['fillVal']
        expected = numpy.array([[fv, fv, fv], [fv, data[1], data[1]]])
        numpy.testing.assert_array_almost_equal(resDict['outTest2D'][:], expected)
        
    def test_nonzero_weight_pix_in_multiple_cells_3D(self):
        self.cfrac[0, 29:32] = [.5, .12, 1.0]
        self.qualFlag[0, 29:32] = [0, 0, 0]
        self.solZenAng[0, 29:32] = [30, 20, 10]
        self.time[0, 29:32] = [self.toTAI93('23:45:10 08-30-2011'),
                               self.toTAI93('23:45:10 08-30-2011'),
                               self.toTAI93('23:45:10 08-30-2011')]
        self.lon[0, 29:32] = [0, 0, 0]
        data = numpy.random.rand(3,4)
        self.test3D[0, 29:32, :] = data
        self.mapDict[(1,1)] = [((0,29), None), ((0,30), None)]
        self.mapDict[(1,2)] = [((0,30), None), ((0,31), None)]
        resDict = self.defOutFunc(self.mapDict, self.six_el_grid, 
                                  self.outFname, verbose=False,
                                  version = self.version)
        expected = numpy.empty((2,3,4))
        expected[:] = self.defParms['fillVal']
        expected[1,1,:] = data[1,:]
        expected[1,2,:] = data[1,:]
        numpy.testing.assert_array_almost_equal(resDict['outTest3D'][:], expected)
        
    def test_output_file_is_netcdf(self):
        for i,j in product(range(2), range(3)):
            self.mapDict[(i,j)] = []
        unused_result = self.defOutFunc(self.mapDict, self.six_el_grid, 
                                        self.outFname, verbose=False,
                                        version = self.version)
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        # passes if no exception is raised in the above 2 lines

    def test_output_file_opens_if_variables_share_extra_dim(self):
        self.defParms['inFieldNames'].append('test3Dagain')
        self.defParms['outFieldNames'].append('outTest3Dagain')
        self.defParms['outUnits'].append('fortnights')
        self.defParms['extraDimLabel'].append('layer')
        self.defParms['extraDimSize'].append(4)
        newOutFunc = out_geo.OMNO2e_netCDF_avg_out_func(self.defParms)
        test3Dagain = numpy.zeros((20,60,4))
        self.parser.prime_get('test3Dagain', test3Dagain)
        for i,j in product(range(2), range(3)):
            self.mapDict[(i,j)] = []
        unused_result = newOutFunc(self.mapDict, self.six_el_grid, 
                                   self.outFname, verbose=False,
                                   version = self.version)
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        # passes if no exception is raised in the above 2 lines

        
    def test_output_file_contains_start_time(self):
        for i,j in product(range(2), range(3)):
            self.mapDict[(i,j)] = []
        unused_result = self.defOutFunc(self.mapDict, self.six_el_grid, 
                                        self.outFname, verbose=False,
                                        version = self.version)
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        self.assertEqual(self.fid.File_start_time, self.startTimeStr)
        
    def test_output_file_contains_stop_time(self):
        for i,j in product(range(2), range(3)):
            self.mapDict[(i,j)] = []
        unused_result = self.defOutFunc(self.mapDict, self.six_el_grid, 
                                        self.outFname, verbose=False,
                                        version = self.version)
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        self.assertEqual(self.fid.File_end_time, self.stopTimeStr)   
        
    def test_output_file_contains_grid_name(self):    
        for i,j in product(range(2), range(3)):
            self.mapDict[(i,j)] = []
        unused_result = self.defOutFunc(self.mapDict, self.six_el_grid, 
                                        self.outFname, verbose=False,
                                        version = self.version)
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        self.assertEqual(self.fid.Projection, 'latlon')
        
    def test_output_file_contains_cfrac_cutoff(self):
        for i,j in product(range(2), range(3)):
            self.mapDict[(i,j)] = []
        unused_result = self.defOutFunc(self.mapDict, self.six_el_grid, 
                                        self.outFname, verbose=False,
                                        version = self.version)
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        self.assertEqual(self.fid.Max_valid_cloud_fraction, self.defParms['cloudFractUpperCutoff'])        
        
    def test_output_file_contains_sza_cutoff(self):
        for i,j in product(range(2), range(3)):
            self.mapDict[(i,j)] = []
        unused_result = self.defOutFunc(self.mapDict, self.six_el_grid, 
                                        self.outFname, verbose=False,
                                        version = self.version)
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        self.assertEqual(self.fid.Max_valid_solar_zenith_angle, self.defParms['solarZenAngUpperCutoff'])        

        
    def test_output_file_contains_tComp(self):
        for i,j in product(range(2), range(3)):
            self.mapDict[(i,j)] = []
        unused_result = self.defOutFunc(self.mapDict, self.six_el_grid, 
                                        self.outFname, verbose=False,
                                        version = self.version)
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        self.assertEqual(self.fid.Time_comparison_scheme, self.defParms['timeComparison'])        
        
    def test_output_file_contains_input_file_list(self):
        self.mapDict[(0,0)] = []
        firstDict = self.mapDict
        secondParser = fakeParser('bar.dat')
        secondDict = {'parser' : secondParser, (0,0) : []}
        secondParser.prime_get('qualFlag', self.qualFlag)
        secondParser.prime_get('solZenAng', self.solZenAng)
        secondParser.prime_get('cfrac', self.cfrac)
        secondParser.prime_get('time', self.time)
        secondParser.prime_get('lon', self.lon)
        secondParser.prime_get('test2D', self.test2D)
        secondParser.prime_get('test3D', self.test3D)
        dictList = [firstDict, secondDict]
        unused_result = self.defOutFunc(dictList, self.one_el_grid, 
                                        self.outFname, verbose=False,
                                        version = self.version)
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        self.assertEqual(self.fid.Input_files, 'foo.dat bar.dat')
        
    def test_output_file_contains_gridParms(self):
        for i,j in product(range(2), range(3)):
            self.mapDict[(i,j)] = []
        unused_result = self.defOutFunc(self.mapDict, self.six_el_grid, 
                                        self.outFname, verbose=False,
                                        version = self.version)
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        outDict = {'xOrig' : self.fid.xOrig, 'yOrig' : self.fid.yOrig,
                   'xCell' : self.fid.xCell, 'yCell' : self.fid.yCell,
                   'nRows' : self.fid.nRows, 'nCols' : self.fid.nCols}
        self.assertDictEqual(self.six_el_gridParms, outDict)
        
    def test_output_file_dims_right_sizes(self):
        for i,j in product(range(2), range(3)):
            self.mapDict[(i,j)] = []
        unused_result = self.defOutFunc(self.mapDict, self.six_el_grid, 
                                        self.outFname, verbose=False,
                                        version = self.version)
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        expectedDims = {'row' : 2, 'col' : 3, 'layer' : 4}       
        actualDims = {}
        for (name,dim) in self.fid.dimensions.items():
            actualDims[name] = len(dim)
        self.assertDictEqual(expectedDims, actualDims)
        
    def test_output_file_right_variables_no_pixel_count(self):
        for i,j in product(range(2), range(3)):
            self.mapDict[(i,j)] = []
        unused_result = self.defOutFunc(self.mapDict, self.six_el_grid, 
                                        self.outFname, verbose=False,
                                        version = self.version)
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        outVars = self.fid.variables.keys()
        self.assertItemsEqual(self.defParms['outFieldNames'], outVars)

    def test_output_file_right_variables_with_pixel_count(self):
        self.defParms['includePixelCount'] = True
        newOutFunc = out_geo.OMNO2e_netCDF_avg_out_func(self.defParms)
        for i,j in product(range(2), range(3)):
            self.mapDict[(i,j)] = []
        unused_result = newOutFunc(self.mapDict, self.six_el_grid, 
                                   self.outFname, verbose=False,
                                   version = self.version)
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        outVars = self.fid.variables.keys()
        correctFields = self.defParms['outFieldNames'] + ['ValidPixelCount']
        self.assertItemsEqual(correctFields, outVars)
        
    def test_output_file_correct_fillVals(self):
        for i,j in product(range(2), range(3)):
            self.mapDict[(i,j)] = []
        unused_result = self.defOutFunc(self.mapDict, self.six_el_grid, 
                                        self.outFname, verbose=False,
                                        version = self.version)
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        fillVals = [self.fid.variables['outTest2D']._FillValue,
                    self.fid.variables['outTest3D']._FillValue]
        expected = 2*[self.defParms['fillVal']]
        self.assertListEqual(expected, fillVals)
        
    def test_output_file_correct_units(self):
        for i,j in product(range(2), range(3)):
            self.mapDict[(i,j)] = []
        unused_result = self.defOutFunc(self.mapDict, self.six_el_grid, 
                                        self.outFname, verbose=False,
                                        version = self.version)
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        units = [self.fid.variables['outTest2D'].Units,
                 self.fid.variables['outTest3D'].Units]
        self.assertListEqual(units, self.defParms['outUnits'])        
        
    def test_output_file_variables_correct_shape(self):
        for i,j in product(range(2), range(3)):
            self.mapDict[(i,j)] = []
        unused_result = self.defOutFunc(self.mapDict, self.six_el_grid, 
                                        self.outFname, verbose=False,
                                        version = self.version)
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        shapes = [self.fid.variables['outTest2D'].shape,
                  self.fid.variables['outTest3D'].shape]
        expected = [(2,3), (2,3,4)]
        self.assertListEqual(expected, shapes)
        
    def test_output_file_correct_vals_3D_one_el_grid(self):
        # calculation of weight performed by hand
        self.cfrac[0,0:2] = [.05, .10]
        self.qualFlag[0,0:2] = [0, 0]
        self.solZenAng[0,0:2] = [50, 55]
        self.time[0,0:2] = [self.toTAI93('11:14:38 08-30-2011'),
                              self.toTAI93('16:11:00 08-30-2011')]
        self.lon[0,0:2] = [-48, 67]
        data = numpy.random.rand(2,4)
        self.test3D[0,0:2] = data
        secondMapDict = dict(self.mapDict)
        self.mapDict[(0,0)] = [((0,0), None)]
        secondMapDict[(0,0)] = [((0,1), None)]
        dictList = [self.mapDict, secondMapDict]
        unused_result = self.defOutFunc(dictList, self.one_el_grid, 
                                        self.outFname, verbose=False,
                                        version = self.version)
        weight = 1.088125096 # calculated by hand
        expected = (data[0,:]*weight+data[1,:])/(weight+1)
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        out = self.fid.variables['outTest3D'][0,0,:]
        numpy.testing.assert_array_almost_equal(expected, out)
        
    def test_output_file_correct_vals_2D_one_el_grid(self):
        # calculation of weight performed by hand
        self.cfrac[0,0:2] = [.05, .10]
        self.qualFlag[0,0:2] = [0, 0]
        self.solZenAng[0,0:2] = [50, 55]
        self.time[0,0:2] = [self.toTAI93('11:14:38 08-30-2011'),
                              self.toTAI93('16:11:00 08-30-2011')]
        self.lon[0,0:2] = [-48, 67]
        data = numpy.random.rand(2)
        self.test2D[0,0:2] = data
        secondMapDict = dict(self.mapDict)
        self.mapDict[(0,0)] = [((0,0), None)]
        secondMapDict[(0,0)] = [((0,1), None)]
        dictList = [self.mapDict, secondMapDict]
        unused_result = self.defOutFunc(dictList, self.one_el_grid, 
                                        self.outFname, verbose=False,
                                        version = self.version)
        weight = 1.088125096 # calculated by hand
        expected = (data[0]*weight+data[1])/(weight+1)
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        out = self.fid.variables['outTest2D'][0,0]
        numpy.testing.assert_array_almost_equal(expected, out)
        
    def test_ouput_file_correct_vals_3D_six_el_grid(self):
        self.cfrac[0:2, 28:34] = [[.11, .12, .13, .14, .15, .16],
                                  [.21, .22, 23, .24, .25, .26]]
        self.qualFlag[0:2, 28:34] = [[0, 0, 0, 0, 0, 0],
                                     [0, 0, 0, 1, 1, 1]]
        self.solZenAng[0:2, 28:34] = [[10, 11, 12, 13, 14, 15],
                                      [87, 88, 89, 20, 12, 13]]
        self.time[0:2, 28:34] = 2*[[self.toTAI93('17:20:00 08-30-2011'),
                                    self.toTAI93('17:20:00 08-30-2011'),
                                    self.toTAI93('17:20:00 08-30-2011'),
                                    self.toTAI93('17:20:00 08-30-2011'),
                                    self.toTAI93('17:20:00 08-30-2011'),
                                    self.toTAI93('17:20:00 08-30-2011')]]
        self.lon[0:2, 28:34] = [[1, 2, 3, 4, 5, 6],
                                [1, 2, 3, 4, 5, 7]]
        data = numpy.random.rand(2,6,4)
        self.test3D[0:2,28:34, :] = data
        for (i,j) in product(range(2), range(3)):
            oneDind = i*3+j
            self.mapDict[(i,j)] = [((0,28+oneDind), None), ((1,28+oneDind), None)]        
        unused_result = self.defOutFunc(self.mapDict, self.six_el_grid, 
                                        self.outFname, verbose=False,
                                        version = self.version)
        expected = data[0,:,:].reshape((2,3,4))
        self.fid = netCDF4.Dataset(self.outFname, 'r')       
        out = self.fid.variables['outTest3D'][:]
        numpy.testing.assert_array_almost_equal(expected, out)
        
    def test_output_file_correct_vals_2D_six_el_grid(self):
        self.cfrac[0:2, 28:34] = [[.11, .12, .13, .14, .15, .16],
                                  [.21, .22, 23, .24, .25, .26]]
        self.qualFlag[0:2, 28:34] = [[0, 0, 0, 0, 0, 0],
                                     [0, 0, 0, 1, 1, 1]]
        self.solZenAng[0:2, 28:34] = [[10, 11, 12, 13, 14, 15],
                                      [87, 88, 89, 20, 12, 13]]
        self.time[0:2, 28:34] = 2*[[self.toTAI93('17:20:00 08-30-2011'),
                                    self.toTAI93('17:20:00 08-30-2011'),
                                    self.toTAI93('17:20:00 08-30-2011'),
                                    self.toTAI93('17:20:00 08-30-2011'),
                                    self.toTAI93('17:20:00 08-30-2011'),
                                    self.toTAI93('17:20:00 08-30-2011')]]
        self.lon[0:2, 28:34] = [[1, 2, 3, 4, 5, 6],
                                [1, 2, 3, 4, 5, 7]]
        data = numpy.random.rand(2,6)
        self.test2D[0:2,28:34] = data
        for (i,j) in product(range(2), range(3)):
            oneDind = i*3+j
            self.mapDict[(i,j)] = [((0,28+oneDind), None), ((1,28+oneDind), None)]
        unused_result = self.defOutFunc(self.mapDict, self.six_el_grid, 
                                        self.outFname, verbose=False,
                                        version = self.version)
        expected = data[0,:].reshape((2,3))
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        out = self.fid.variables['outTest2D'][:]
        numpy.testing.assert_array_almost_equal(expected, out)


class Test_unweighted_filtered_MOPITT_avg_netCDF_out_func(TestOutGeo):
    
    def setUp(self):
        TestOutGeo.setUp(self)
        self.version = 'TEST VERSION'
        self.six_el_gridParms = {'xOrig' : 0, 'yOrig' : 0, 'xCell' : 1,
                            'yCell' : 1, 'nRows' : 2, 'nCols' : 3}
        self.sixElGr = grid_geo.latlon_GridDef(self.six_el_gridParms)
        # initialize the dictionaries empty for all cells.  Only have to change
        # to add actual reference
        self.mapDict = {'parser' : self.parser, (0,0) : [], (1,0) : [],
                        (0,1) : [], (1,1) : [], (0,2) : [], (1,2) : []}
        self.pDict = { 'time' : 'time',
                       'longitude' : 'lon',
                       'inFieldNames' : ['twoDnorm', 'threeDnorm', 'twoDlog', 'threeDlog', 'fourDcol'],
                       'outFieldNames' : ['twoDnorm', 'threeDnorm', 'twoDlog', 'threeDlog', 'fourDcol'],
                       'outUnits' : ['foo', 'bar', 'baz', 'qux', 'spam'],
                       'logNormal' : ['False', 'False', 'True', 'True', 'False'],
                       'dimLabels' : [[], ['layer'], [], ['layer'], ['layer', 'value']],
                       'dimSizes' : [[], ['3'], [], ['3'], ['3','2']],
                       'timeStart' : '00:00:00_01-04-2012',
                       'timeStop' : '23:59:59_01-04-2012',
                       'timeComparison' : 'UTC',
                       'fillVal' : -9999.0,
                       'solZenAngCutoff' : 85,
                       'solZenAng' : 'SZA',
                       'dayTime' : True,
                       'surfTypeField': 'sType',
                       'colMeasField' : 'fourDcol' }
        self.defaultOutClass = out_geo.unweighted_filtered_MOPITT_avg_netCDF_out_func(self.pDict)
        # define a default output function, as we'll be using the same grid, map and outfile 
        # many times.
        (outFid, self.outFname) = (tempfile.mkstemp())
        os.close(outFid)
        self.defaultOutFunc = lambda: self.defaultOutClass(self.mapDict, self.sixElGr, 
                                                           self.outFname, verbose=False,
                                                           version = self.version)
        # create default arrays.  These arrays are defined such that 
        # the values to be averaged are random, and the values that determine
        # the validity of the data show it all to be valid.  Thus, test
        # methods need only make changes to INVALIDATE selected chunks of
        # the data
        self.time = numpy.array([self.toTAI93('12:00:00 01-04-2012')]*20).reshape(4,5)
        self.lon = numpy.zeros_like(self.time)
        self.SZA = numpy.zeros_like(self.time)
        self.sType = numpy.zeros_like(self.time)
        self.twoDnorm = numpy.random.rand(4,5)
        self.twoDlog = numpy.random.rand(4,5)
        self.threeDnorm = numpy.random.rand(4,5,3)
        self.threeDlog = numpy.random.rand(4,5,3)
        self.fourDcol = numpy.random.rand(4,5,3,2)
        # bind the default arrays into the parser.  they can still be 
        # modified, and these changes will be reflected in what comes
        # back out of the parser.
        fNames = ['time', 'lon', 'SZA', 'sType', 'twoDnorm', 'twoDlog', 'threeDnorm', 'threeDlog', 'fourDcol']
        for fn in fNames:
            self.parser.prime_get(fn, getattr(self, fn))
            
    def tearDown(self):
        try:
            self.fid.close()
        except AttributeError:
            pass
        os.remove(self.outFname)

    def test_parser_still_in_empty_map(self):
        unused_result = self.defaultOutFunc()
        self.assertIs(self.mapDict['parser'], self.parser)

    def test_parser_still_in_nonempty_map(self):
        self.mapDict[(0,0)] = [((2,3), None)]
        unused_result = self.defaultOutFunc()
        self.assertIs(self.mapDict['parser'], self.parser)

    def test_2D_norm_results_right_shape_empty_map(self):
        resDict = self.defaultOutFunc()
        self.assertEqual(resDict['twoDnorm'].shape, (2,3))

    def test_2D_log_results_right_shape_empty_map(self):
        resDict = self.defaultOutFunc()
        self.assertEqual(resDict['twoDlog'].shape, (2,3))

    def test_3D_norm_results_right_shape_empty_map(self):
        resDict = self.defaultOutFunc()
        self.assertEqual(resDict['threeDnorm'].shape, (2,3,3))

    def test_3D_log_results_right_shape_empty_map(self):
        resDict = self.defaultOutFunc()
        self.assertEqual(resDict['threeDlog'].shape, (2,3,3))
    
    def test_4D_results_right_shape_empty_map(self):
        resDict = self.defaultOutFunc()
        self.assertEqual(resDict['fourDcol'].shape, (2,3,3,2))

    def test_2D_norm_results_right_shape_non_empty_map(self):
        self.mapDict[(0,0)] = [((2,3), None)]
        self.mapDict[(1,2)] = [((3,4), None), ((0,0), None)]
        resDict = self.defaultOutFunc()
        self.assertEqual(resDict['twoDnorm'].shape, (2,3))

    def test_2D_log_results_right_shape_non_empty_map(self):
        self.mapDict[(0,0)] = [((2,3), None)]
        self.mapDict[(1,2)] = [((3,4), None), ((0,0), None)]
        resDict = self.defaultOutFunc()
        self.assertEqual(resDict['twoDlog'].shape, (2,3))

    def test_3D_norm_results_right_shape_non_empty_map(self):
        self.mapDict[(0,0)] = [((2,3), None)]
        self.mapDict[(1,2)] = [((3,4), None), ((0,0), None)]
        resDict = self.defaultOutFunc()
        self.assertEqual(resDict['threeDnorm'].shape, (2,3,3))

    def test_3D_log_results_right_shape_non_empty_map(self):
        self.mapDict[(0,0)] = [((2,3), None)]
        self.mapDict[(1,2)] = [((3,4), None), ((0,0), None)]
        resDict = self.defaultOutFunc()
        self.assertEqual(resDict['threeDlog'].shape, (2,3,3))

    def test_4D_results_right_shape_non_empty_map(self):
        resDict = self.defaultOutFunc()
        self.assertEqual(resDict['fourDcol'].shape, (2,3,3,2))

    def test_single_valid_pixel_2D_norm(self):
        self.mapDict[(1,2)] = [((2,3), None)]
        expected = self.twoDnorm[2,3]
        resDict = self.defaultOutFunc()
        self.assertAlmostEqual(resDict['twoDnorm'][1,2], expected)

    def test_single_valid_pixel_2D_log(self):
        self.mapDict[(1,2)] = [((2,3), None)]
        expected = self.twoDlog[2,3]
        resDict = self.defaultOutFunc()
        self.assertAlmostEqual(resDict['twoDlog'][1,2], expected)

    def test_single_valid_pixel_3D_norm(self):
        self.mapDict[(1,2)] = [((2,3), None)]
        expected = self.threeDnorm[2,3,:]
        resDict = self.defaultOutFunc()
        numpy.testing.assert_array_almost_equal(resDict['threeDnorm'][1,2,:], expected)
        
    def test_single_valid_pixel_3D_log(self):
        self.mapDict[(1,2)] = [((2,3), None)]
        expected = self.threeDlog[2,3,:]
        resDict = self.defaultOutFunc()
        numpy.testing.assert_array_almost_equal(resDict['threeDlog'][1,2,:], expected)

    def test_single_valid_pixel_4D_norm(self):
        self.mapDict[(1,2)] = [((2,3), None)]
        expected = self.fourDcol[2,3,...]
        resDict = self.defaultOutFunc()
        numpy.testing.assert_array_almost_equal(resDict['fourDcol'][1,2,...], expected)

    def test_two_valid_pixel_2D_norm(self):
        self.mapDict[(1,2)] = [((2,3), None), ((2,4), None)]
        self.twoDnorm[2,3] = 4.2
        self.twoDnorm[2,4] = 4.9
        expected = 4.55
        resDict = self.defaultOutFunc()
        self.assertAlmostEqual(resDict['twoDnorm'][1,2], expected)

    def test_two_valid_pixel_2D_log(self):
        self.mapDict[(1,2)] = [((2,3), None), ((2,4), None)]
        self.twoDlog[2,3] = 4.2
        self.twoDlog[2,4] = 4.9
        expected = 4.53651848888
        resDict = self.defaultOutFunc()
        self.assertAlmostEqual(resDict['twoDlog'][1,2], expected)

    def test_two_valid_pixel_3D_norm(self):
        self.mapDict[(1,2)] = [((2,3), None), ((2,4), None)]
        self.threeDnorm[2,3,:] = numpy.array([1.9, 2.7, 3.1])
        self.threeDnorm[2,4,:] = numpy.array([1.7, 2.4, 5.1])
        expected = [1.8, 2.55, 4.1]
        resDict = self.defaultOutFunc()
        numpy.testing.assert_array_almost_equal(resDict['threeDnorm'][1,2,:], expected)

    def test_two_valid_pixel_3D_log(self):
        self.mapDict[(1,2)] = [((2,3), None), ((2,4), None)]
        self.threeDlog[2,3,:] = numpy.array([1.9, 2.7, 3.1])
        self.threeDlog[2,4,:] = numpy.array([1.7, 2.4, 5.1])
        expected = [1.7972200755611, 2.545584412271, 
                    3.9761790704142]
        resDict = self.defaultOutFunc()
        numpy.testing.assert_array_almost_equal(resDict['threeDlog'][1,2,:], expected)

    def test_two_valid_pixel_4D_norm(self):
        self.mapDict[(1,2)] = [((2,3), None), ((2,4), None)]
        self.fourDcol[2,3,...] = numpy.array([[2.2,3.2], [4.4,5.4], [7.6,8.6]])
        self.fourDcol[2,4,...] = numpy.array([[2.4,3.4], [4.6,5.6], [7.8,8.8]])
        expected = [[2.3,3.3], [4.5,5.5], [7.7,8.7]]
        resDict = self.defaultOutFunc()
        numpy.testing.assert_array_almost_equal(resDict['fourDcol'][1,2,...], expected)

    def test_two_valid_pixel_4D_log(self):
        self.pDict['logNormal'][-1] = 'True'
        newOutClass = out_geo.unweighted_filtered_MOPITT_avg_netCDF_out_func(self.pDict)
        self.mapDict[(1,2)] = [((2,3), None), ((2,4), None)]
        self.fourDcol[2,3,...] = numpy.array([[2,4], [6,8], [10,12]])
        self.fourDcol[2,4,...] = numpy.array([[.5,1], [1.5,2], [2.5,3]])
        expected = [[1,2], [3,4], [5,6]]
        resDict = newOutClass(self.mapDict, self.sixElGr, 
                              self.outFname, verbose=False,
                              version = self.version)
        numpy.testing.assert_array_almost_equal(resDict['fourDcol'][1,2,...], expected)

    def test_zero_weight_if_SZA_gt_thresh_2D_norm(self):
        self.mapDict[(1,2)] = [((2,3), None), ((2,4), None)]
        self.SZA[2,4] = 100
        expected = self.twoDnorm[2,3]
        resDict = self.defaultOutFunc()
        self.assertAlmostEqual(resDict['twoDnorm'][1,2], expected)
        
    def test_zero_weight_if_SZA_gt_thresh_2D_log(self):
        self.mapDict[(1,2)] = [((2,3), None), ((2,4), None)]
        self.SZA[2,4] = 100
        expected = self.twoDlog[2,3]
        resDict = self.defaultOutFunc()
        self.assertAlmostEqual(resDict['twoDlog'][1,2], expected)

    def test_zero_weight_if_SZA_gt_thresh_3D_norm(self):
        self.mapDict[(1,2)] = [((2,3), None), ((2,4), None)]
        self.SZA[2,4] = 100
        expected = self.threeDnorm[2,3,:]
        resDict = self.defaultOutFunc()
        numpy.testing.assert_array_almost_equal(resDict['threeDnorm'][1,2,:], expected)

    def test_zero_weight_if_SZA_gt_thresh_3D_log(self):
        self.mapDict[(1,2)] = [((2,3), None), ((2,4), None)]
        self.SZA[2,4] = 100
        expected = self.threeDlog[2,3,:]
        resDict = self.defaultOutFunc()
        numpy.testing.assert_array_almost_equal(resDict['threeDlog'][1,2,:], expected)

    def test_zero_weight_if_SZA_gt_thresh_4D(self):
        self.mapDict[(1,2)] = [((2,3), None), ((2,4), None)]
        self.SZA[2,4] = 100
        expected = self.fourDcol[2,3,...]
        resDict = self.defaultOutFunc()
        numpy.testing.assert_array_almost_equal(resDict['fourDcol'][1,2,...], expected)

    def test_properly_flips_SZA_when_night_specified_2D(self):
        self.pDict['dayTime'] = False
        newOutClass = out_geo.unweighted_filtered_MOPITT_avg_netCDF_out_func(self.pDict)
        self.mapDict[(1,2)] = [((2,3), None), ((2,4), None)]
        self.SZA[2,4] = 100
        expected = self.twoDnorm[2,4]
        resDict = newOutClass(self.mapDict, self.sixElGr, 
                              self.outFname, verbose=False,
                              version = self.version)
        self.assertAlmostEqual(resDict['twoDnorm'][1,2], expected)

    def test_properly_flips_SZA_when_night_specified_3D(self):
        self.pDict['dayTime'] = False
        newOutClass = out_geo.unweighted_filtered_MOPITT_avg_netCDF_out_func(self.pDict)
        self.mapDict[(1,2)] = [((2,3), None), ((2,4), None)]
        self.SZA[2,4] = 100
        expected = self.threeDnorm[2,4,:]
        resDict = newOutClass(self.mapDict, self.sixElGr, 
                              self.outFname, verbose=False,
                              version = self.version)
        numpy.testing.assert_array_almost_equal(resDict['threeDnorm'][1,2,:], expected)

    def test_screen_out_bad_sTypes_2D(self):
        self.mapDict[(1,2)] = [((0,0), None), ((0,1), None), ((0,2), None), 
                               ((1,0), None), ((1,1), None), ((1,2), None),
                               ((2,0), None), ((2,1), None), ((2,2), None)]
        self.sType[0,0] = 1
        self.sType[0,1] = 2
        expected = numpy.random.rand()
        self.twoDnorm[0,2] = expected
        for i in [0,1,2]:
            self.twoDnorm[1,i] = expected-i
            self.twoDnorm[2,i] = expected+i
        resDict = self.defaultOutFunc()
        self.assertAlmostEqual(resDict['twoDnorm'][1,2], expected)

    def test_screen_out_bad_sTypes_3D(self):
        self.mapDict[(1,2)] = [((0,0), None), ((0,1), None), ((0,2), None), 
                               ((1,0), None), ((1,1), None), ((1,2), None),
                               ((2,0), None), ((2,1), None), ((2,2), None)]
        self.sType[0,0] = 1
        self.sType[0,1] = 2
        expected = numpy.random.rand(3)
        self.threeDnorm[0,2,:] = expected
        for i in [0,1,2]:
            self.threeDnorm[1,i,:] = expected-i
            self.threeDnorm[2,i,:] = expected+i
        resDict = self.defaultOutFunc()
        numpy.testing.assert_array_almost_equal(resDict['threeDnorm'][1,2,:], expected)
        
    def test_screen_out_bad_stypes_4D(self):
        self.mapDict[(1,2)] = [((0,0), None), ((0,1), None), ((0,2), None), 
                               ((1,0), None), ((1,1), None), ((1,2), None),
                               ((2,0), None), ((2,1), None), ((2,2), None)]
        self.sType[0,0] = 1
        self.sType[0,1] = 2
        expected = numpy.random.rand(3,2)
        self.fourDcol[0,2,:,:] = expected
        for i in [0,1,2]:
            self.fourDcol[1,i,...] = expected-i
            self.fourDcol[2,i,...] = expected+i
        resDict = self.defaultOutFunc()
        numpy.testing.assert_array_almost_equal(resDict['fourDcol'][1,2,:,:], expected)

    def test_use_all_sTypes_if_no_threshold_met_2D(self):
        self.mapDict[(1,2)] = [((3,1), None), ((3,2), None), ((3,3), None), ((3,4), None)]
        self.sType[3,1] = 10
        self.sType[3,2] = 100
        expected = numpy.random.rand()
        self.twoDnorm[3,1] = expected+1
        self.twoDnorm[3,2] = expected+1
        self.twoDnorm[3,3] = expected+1
        self.twoDnorm[3,4] = expected-3
        resDict = self.defaultOutFunc()
        self.assertAlmostEqual(resDict['twoDnorm'][1,2], expected)
        
    def test_use_all_sTypes_if_no_threshold_met_3D(self):
        self.mapDict[(1,2)] = [((3,1), None), ((3,2), None), ((3,3), None), ((3,4), None)]
        self.sType[3,1] = 10
        self.sType[3,2] = 100
        expected = numpy.random.rand(3)
        self.threeDnorm[3,1,:] = expected+1
        self.threeDnorm[3,2,:] = expected+1
        self.threeDnorm[3,3,:] = expected+1
        self.threeDnorm[3,4,:] = expected-3
        resDict = self.defaultOutFunc()
        numpy.testing.assert_array_almost_equal(resDict['threeDnorm'][1,2,:], expected)

    def test_screen_properly_on_exactly_threshold_2D(self):
        self.mapDict[(1,2)] = [((3,0), None), ((3,1), None), ((3,2), None), ((3,4), None)]
        self.sType[3,4] = -1
        expected = numpy.random.rand()
        self.twoDnorm[3,0] = expected*1.5
        self.twoDnorm[3,1] = expected*.5
        self.twoDnorm[3,2] = expected
        resDict = self.defaultOutFunc()
        self.assertAlmostEqual(resDict['twoDnorm'][1,2], expected)

    def test_screen_properly_on_exactly_threshold_3D(self):
        self.mapDict[(1,2)] = [((3,0), None), ((3,1), None), ((3,2), None), ((3,4), None)]
        self.sType[3,4] = -1
        expected = numpy.random.rand(3)
        self.threeDnorm[3,0,:] = expected*1.5
        self.threeDnorm[3,1,:] = expected*0.5
        self.threeDnorm[3,2,:] = expected
        resDict = self.defaultOutFunc()
        numpy.testing.assert_array_almost_equal(resDict['threeDnorm'][1,2,:], expected)

    def test_screen_out_minority_levels_when_fewer_2D(self):
        self.mapDict[(1,2)] = [((0,0), None), ((1,0), None), ((2,0), None)]
        self.fourDcol[0,0,2,0] = numpy.NaN
        expected = numpy.random.rand()
        self.twoDnorm[1,0] = expected+3.5
        self.twoDnorm[2,0] = expected-3.5
        resDict = self.defaultOutFunc()
        self.assertAlmostEqual(resDict['twoDnorm'][1,2], expected)

    def test_screen_out_minority_levels_when_fewer_4D(self):
        self.mapDict[(1,2)] = [((0,0), None), ((1,0), None), ((2,0), None)]
        self.fourDcol[0,0,2,0] = numpy.NaN
        expected = numpy.random.rand(3,2)
        self.fourDcol[1,0,:,:] = expected+numpy.ones((3,2))
        self.fourDcol[2,0,:,:] = expected-numpy.ones((3,2))
        resDict = self.defaultOutFunc()
        numpy.testing.assert_array_almost_equal(resDict['fourDcol'][1,2,:,:], expected)

    def test_screen_out_minority_levels_when_fewer_3D(self):
        self.mapDict[(1,2)] = [((0,0), None), ((1,0), None), ((2,0), None)]
        self.fourDcol[0,0,2,0] = numpy.NaN
        expected = numpy.random.rand(3)
        self.threeDnorm[1,0,:] = expected+numpy.array([3.5, 4.5, 5.5])
        self.threeDnorm[2,0,:] = expected-numpy.array([3.5, 4.5, 5.5])
        resDict = self.defaultOutFunc()
        numpy.testing.assert_array_almost_equal(resDict['threeDnorm'][1,2,:], expected)

    def test_screen_out_minority_levels_when_greater_2D(self):
        self.mapDict[(1,2)] = [((0,0), None), ((0,1), None), ((0,2), None)]
        self.fourDcol[0,0,2,0] = numpy.NaN
        self.fourDcol[0,1,2,0] = numpy.NaN
        expected = numpy.random.rand()
        self.twoDnorm[0,0] = expected+1
        self.twoDnorm[0,1] = expected-1
        resDict = self.defaultOutFunc()
        self.assertAlmostEqual(resDict['twoDnorm'][1,2], expected)

    def test_screen_out_minority_levels_when_greater_3D(self):
        self.mapDict[(1,2)] = [((0,0), None), ((0,1), None), ((0,2), None)]
        self.fourDcol[0,0,2,0] = numpy.NaN
        self.fourDcol[0,1,2,0] = numpy.NaN
        expected = numpy.random.rand(3)
        expected[-1] = numpy.NaN
        self.threeDnorm[0,0,:] = expected+1
        self.threeDnorm[0,1,:] = expected-1
        resDict = self.defaultOutFunc()
        expected[-1] = self.pDict['fillVal']
        numpy.testing.assert_array_almost_equal(resDict['threeDnorm'][1,2,:], expected)

    def test_screen_out_minority_levels_when_greater_4D(self):
        self.mapDict[(1,2)] = [((0,0), None), ((0,1), None), ((0,2), None)]
        expected = numpy.random.rand(3,2)
        expected[2,0] = numpy.NaN
        self.fourDcol[0,0,:,:] = expected+1
        self.fourDcol[0,1,:,:] = expected-1
        resDict = self.defaultOutFunc()
        expected[2,0] = self.pDict['fillVal']
        numpy.testing.assert_array_almost_equal(resDict['fourDcol'][1,2,:,:], expected)

    def test_screen_properly_in_tie_2D(self):
        self.mapDict[(1,2)] = [((0,0), None), ((0,1), None), ((0,2), None)]
        self.fourDcol[0,1,2,0] = numpy.NaN
        self.fourDcol[0,2,1,0] = numpy.NaN
        self.fourDcol[0,2,2,0] = numpy.NaN
        expected = self.twoDnorm[0,0]
        resDict = self.defaultOutFunc()
        self.assertAlmostEqual(resDict['twoDnorm'][1,2], expected)

    def test_screen_properly_in_tie_3D(self):
        self.mapDict[(1,2)] = [((0,0), None), ((0,1), None), ((0,2), None)]
        self.fourDcol[0,1,2,0] = numpy.NaN
        self.fourDcol[0,2,1,0] = numpy.NaN
        self.fourDcol[0,2,2,0] = numpy.NaN
        expected = self.threeDnorm[0,0,:]
        resDict = self.defaultOutFunc()
        numpy.testing.assert_array_almost_equal(resDict['threeDnorm'][1,2,:], expected)
        
    def test_screen_properly_in_tie_4D(self):
        self.mapDict[(1,2)] = [((0,0), None), ((0,1), None), ((0,2), None)]
        self.fourDcol[0,1,2,0] = numpy.NaN
        self.fourDcol[0,2,1,0] = numpy.NaN
        self.fourDcol[0,2,2,0] = numpy.NaN
        expected = self.fourDcol[0,0,:,:]
        resDict = self.defaultOutFunc()
        numpy.testing.assert_array_almost_equal(resDict['fourDcol'][1,2,:], expected)

    def test_screen_out_UTC_times_earlier_than_start_2D(self):
        self.mapDict[(1,2)] = [((2,2), None), ((2,3), None)]
        self.time[2,3] = self.toTAI93('12:00:00 01-03-2012')
        expected = self.twoDnorm[2,2]
        resDict = self.defaultOutFunc()
        self.assertAlmostEqual(resDict['twoDnorm'][1,2], expected)

    def test_screen_out_UTC_times_earlier_than_start_3D(self):
        self.mapDict[(1,2)] = [((2,2), None), ((2,3), None)]
        self.time[2,3] = self.toTAI93('12:00:00 01-03-2012')
        expected = self.threeDnorm[2,2,:]
        resDict = self.defaultOutFunc()
        numpy.testing.assert_array_almost_equal(resDict['threeDnorm'][1,2,:], expected)

    def test_screen_out_UTC_times_earlier_than_start_4D(self):
        self.mapDict[(1,2)] = [((2,2), None), ((2,3), None)]
        self.time[2,3] = self.toTAI93('12:00:00 01-03-2012')
        expected = self.fourDcol[2,2,...]
        resDict = self.defaultOutFunc()
        numpy.testing.assert_array_almost_equal(resDict['fourDcol'][1,2,...], expected)

    def test_screen_out_UTC_times_later_than_stop_2D(self):
        self.mapDict[(1,2)] = [((2,2), None), ((2,3), None)]
        self.time[2,3] = self.toTAI93('12:00:00 01-05-2012')
        expected = self.twoDnorm[2,2]
        resDict = self.defaultOutFunc()
        self.assertAlmostEqual(resDict['twoDnorm'][1,2], expected)

    def test_screen_out_UTC_times_later_than_stop_3D(self):
        self.mapDict[(1,2)] = [((2,2), None), ((2,3), None)]
        self.time[2,3] = self.toTAI93('12:00:00 01-05-2012')
        expected = self.threeDnorm[2,2,:]
        resDict = self.defaultOutFunc()
        numpy.testing.assert_array_almost_equal(resDict['threeDnorm'][1,2,:], expected)

    def test_screen_out_UTC_times_later_than_stop_4D(self):
        self.mapDict[(1,2)] = [((2,2), None), ((2,3), None)]
        self.time[2,3] = self.toTAI93('12:00:00 01-05-2012')
        expected = self.fourDcol[2,2,...]
        resDict = self.defaultOutFunc()
        numpy.testing.assert_array_almost_equal(resDict['fourDcol'][1,2,...], expected)

    def test_UTC_time_does_not_account_for_lon(self):
        self.mapDict[(1,2)] = [((2,2), None), ((2,3), None)]
        self.time[2,3] = self.toTAI93('02:00:00 01-04-2012')
        self.time[2,3] = self.toTAI93('22:00:00 01-03-2012')
        self.lon[2,3] = -150
        self.lon[2,3] = 150
        expected = self.twoDnorm[2,2]
        resDict = self.defaultOutFunc()
        self.assertAlmostEqual(resDict['twoDnorm'][1,2], expected)

    def test_local_time_does_account_for_lon(self):
        self.pDict['timeComparison'] = 'local'
        newOutClass = out_geo.unweighted_filtered_MOPITT_avg_netCDF_out_func(self.pDict)
        self.mapDict[(1,2)] = [((2,2), None), ((2,3), None)]
        self.time[2,2] = self.toTAI93('02:00:00 01-04-2012')
        self.time[2,3] = self.toTAI93('22:00:00 01-03-2012')
        self.lon[2,2] = -150
        self.lon[2,3] = 150
        expected = self.twoDnorm[2,3]
        resDict = newOutClass(self.mapDict, self.sixElGr,
                              self.outFname, verbose=False,
                              version = self.version)
        self.assertAlmostEqual(resDict['twoDnorm'][1,2], expected)

    def test_time_and_flag_still_zero_out_pixel_2D(self):
        self.mapDict[(1,2)] = [((0,1), None), ((0,2), None), ((0,3), None), ((0,4), None)]
        self.SZA[0,1] = 100
        self.SZA[0,3] = 100
        self.time[0,2] = self.toTAI93('12:00:00 02-04-2012')
        self.time[0,3] = self.toTAI93('12:00:00 02-04-2012')
        expected = self.twoDnorm[0,4]
        resDict = self.defaultOutFunc()
        self.assertAlmostEqual(resDict['twoDnorm'][1,2], expected)

    def test_time_and_flag_still_zero_out_pixel_3D(self):
        self.mapDict[(1,2)] = [((0,1), None), ((0,2), None), ((0,3), None), ((0,4), None)]
        self.SZA[0,1] = 100
        self.SZA[0,3] = 100
        self.time[0,2] = self.toTAI93('12:00:00 02-04-2012')
        self.time[0,3] = self.toTAI93('12:00:00 02-04-2012')
        expected = self.threeDnorm[0,4,:]
        resDict = self.defaultOutFunc()
        numpy.testing.assert_array_almost_equal(resDict['threeDnorm'][1,2,:], expected)

    def test_time_and_flag_still_zero_out_pixel_4D(self):
        self.mapDict[(1,2)] = [((0,1), None), ((0,2), None), ((0,3), None), ((0,4), None)]
        self.SZA[0,1] = 100
        self.SZA[0,3] = 100
        self.time[0,2] = self.toTAI93('12:00:00 02-04-2012')
        self.time[0,3] = self.toTAI93('12:00:00 02-04-2012')
        expected = self.fourDcol[0,4,:,:]
        resDict = self.defaultOutFunc()
        numpy.testing.assert_array_almost_equal(resDict['fourDcol'][1,2,:,:], expected)

    def test_time_and_filter_still_zero_out_pixel_2D(self):
        self.mapDict[(1,2)] = [((0,0), None), ((0,1), None), ((0,2), None), ((0,3), None), ((0,4), None)]
        self.sType[0,4] = 1
        self.time[0,4] = self.toTAI93('12:00:00 01-05-2012')
        expected = numpy.random.rand()
        self.twoDnorm[0,0] = expected+5
        self.twoDnorm[0,1] = expected-2
        self.twoDnorm[0,2] = expected-2
        self.twoDnorm[0,3] = expected-1
        resDict = self.defaultOutFunc()
        self.assertAlmostEqual(resDict['twoDnorm'][1,2], expected)

    def test_time_and_filter_still_zero_out_pixel_3D(self):
        self.mapDict[(1,2)] = [((0,0), None), ((0,1), None), ((0,2), None), ((0,3), None), ((0,4), None)]
        self.sType[0,4] = 1
        self.time[0,4] = self.toTAI93('12:00:00 01-05-2012')
        expected = numpy.random.rand(3)
        self.threeDnorm[0,0,:] = expected+5
        self.threeDnorm[0,1,:] = expected-2
        self.threeDnorm[0,2,:] = expected-2
        self.threeDnorm[0,3,:] = expected-1
        resDict = self.defaultOutFunc()
        numpy.testing.assert_array_almost_equal(resDict['threeDnorm'][1,2,:], expected)

    def test_time_and_filter_still_zero_out_pixel_4D(self):
        self.mapDict[(1,2)] = [((0,0), None), ((0,1), None), ((0,2), None), ((0,3), None), ((0,4), None)]
        self.sType[0,4] = 1
        self.time[0,4] = self.toTAI93('12:00:00 01-05-2012')
        expected = numpy.random.rand(3,2)
        self.fourDcol[0,0,...] = expected+5
        self.fourDcol[0,1,...] = expected-2
        self.fourDcol[0,2,...] = expected-2
        self.fourDcol[0,3,...] = expected-1
        resDict = self.defaultOutFunc()
        numpy.testing.assert_array_almost_equal(resDict['fourDcol'][1,2,:,:], expected)

    def test_both_filter_still_zero_out_pixel_2D(self):
        self.mapDict[(1,2)] = [((0,0), None), ((0,1), None), ((0,2), None), ((0,3), None), ((0,4), None)]
        self.sType[0,4] = 1
        self.fourDcol[0,4,2,0] = numpy.NaN
        expected = numpy.random.rand()
        self.twoDnorm[0,0] = expected
        self.twoDnorm[0,1] = expected-3
        self.twoDnorm[0,2] = expected+1
        self.twoDnorm[0,3] = expected+2
        resDict = self.defaultOutFunc()
        self.assertAlmostEqual(resDict['twoDnorm'][1,2], expected)

    def test_both_filter_still_zero_out_pixel_3D(self):
        self.mapDict[(1,2)] = [((0,0), None), ((0,1), None), ((0,2), None), ((0,3), None), ((0,4), None)]
        self.sType[0,4] = 1
        self.fourDcol[0,4,2,0] = numpy.NaN
        expected = numpy.random.rand(3)
        self.threeDnorm[0,0,:] = expected
        self.threeDnorm[0,1,:] = expected-3
        self.threeDnorm[0,2,:] = expected+2
        self.threeDnorm[0,3,:] = expected+1
        resDict = self.defaultOutFunc()
        numpy.testing.assert_array_almost_equal(resDict['threeDnorm'][1,2,:], expected)

    def test_both_filter_still_zero_out_pixel_4D(self):
        self.mapDict[(1,2)] = [((0,0), None), ((0,1), None), ((0,2), None), ((0,3), None), ((0,4), None)]
        self.sType[0,4] = 1
        self.fourDcol[0,4,2,0] = numpy.NaN
        expected = numpy.random.rand(3,2)
        self.fourDcol[0,0,...] = expected
        self.fourDcol[0,1,...] = expected-3
        self.fourDcol[0,2,...] = expected+2
        self.fourDcol[0,3,...] = expected+1
        resDict = self.defaultOutFunc()
        numpy.testing.assert_array_almost_equal(resDict['fourDcol'][1,2,...], expected)


    def test_all_zero_weight_yields_fillVal_2D(self):
        self.mapDict[(1,2)] = [((1,0), None), ((1,1), None)]
        self.SZA[1,0] = 130
        self.time[1,1] = self.toTAI93('12:00:00 01-05-2012')
        expected = self.pDict['fillVal']
        resDict = self.defaultOutFunc()
        self.assertAlmostEqual(resDict['twoDnorm'][1,2], expected)

    def test_all_zero_weight_yields_fillVal_3D(self):
        self.mapDict[(1,2)] = [((1,0), None), ((1,1), None)]
        self.SZA[1,0] = 130
        self.time[1,1] = self.toTAI93('12:00:00 01-05-2012')
        expected = numpy.array([self.pDict['fillVal']]*3)
        resDict = self.defaultOutFunc()
        numpy.testing.assert_array_almost_equal(resDict['threeDnorm'][1,2,:], expected)
        
    def test_all_zero_weight_yields_fillVal_4D(self):
        self.mapDict[(1,2)] = [((1,0), None), ((1,1), None)]
        self.SZA[1,0] = 130
        self.time[1,1] = self.toTAI93('12:00:00 01-05-2012')
        expected = numpy.array([self.pDict['fillVal']]*6).reshape((3,2))
        resDict = self.defaultOutFunc()
        numpy.testing.assert_array_almost_equal(resDict['fourDcol'][1,2,...], expected)

    def test_all_zero_weight_plus_NaN_in_data_yields_fillVal_2D(self):
        self.mapDict[(1,2)] = [((1,0), None), ((1,1), None)]
        self.SZA[1,0] = 130
        self.time[1,1] = self.toTAI93('12:00:00 01-05-2012')
        self.twoDnorm[1,0] = numpy.NaN
        self.twoDnorm[2,0] = numpy.NaN
        expected = self.pDict['fillVal']
        resDict = self.defaultOutFunc()
        self.assertAlmostEqual(resDict['twoDnorm'][1,2], expected)

    def test_all_zero_weight_plus_NaN_in_data_yields_fillVal_3D(self):
        self.mapDict[(1,2)] = [((1,0), None), ((1,1), None)]
        self.SZA[1,0] = 130
        self.time[1,1] = self.toTAI93('12:00:00 01-05-2012')
        self.threeDnorm[1,0,2] = numpy.NaN
        self.threeDnorm[1,1,2] = numpy.NaN
        expected = numpy.array([self.pDict['fillVal']]*3)
        resDict = self.defaultOutFunc()
        numpy.testing.assert_array_almost_equal(resDict['threeDnorm'][1,2,:], expected)

    def test_NaN_in_data_with_zero_weight_does_not_affect_avg_2D(self):
        self.mapDict[(1,2)] = [((2,2), None), ((2,3), None)]
        self.SZA[2,3] = 100
        self.twoDnorm[2,3] = numpy.NaN
        expected = self.twoDnorm[2,2]
        resDict = self.defaultOutFunc()
        self.assertAlmostEqual(resDict['twoDnorm'][1,2], expected)

    def test_NaN_in_data_with_zero_weight_does_not_affect_avg_3D(self):
        self.mapDict[(1,2)] = [((2,2), None), ((2,3), None)]
        self.SZA[2,3] = 100
        self.threeDnorm[2,3,:] = numpy.array([numpy.NaN]*3)
        expected = self.threeDnorm[2,2,:]
        resDict = self.defaultOutFunc()
        numpy.testing.assert_array_almost_equal(resDict['threeDnorm'][1,2,:], expected)

    def test_NaN_in_data_with_zero_weight_does_not_affect_avg_4D(self):
        self.mapDict[(1,2)] = [((2,2), None), ((2,3), None)]
        self.SZA[2,3] = 100
        self.fourDcol[2,3,:,:] = numpy.array([numpy.NaN]*6).reshape((3,2))
        expected = self.fourDcol[2,2,...]
        resDict = self.defaultOutFunc()
        numpy.testing.assert_array_almost_equal(resDict['fourDcol'][1,2,...], expected)

    def test_no_pixels_in_cell_yields_fillVal_2D(self):
        expected = self.pDict['fillVal']
        resDict = self.defaultOutFunc()
        self.assertEqual(resDict['twoDnorm'][1,2], expected)
        
    def test_no_pixels_in_cell_yields_fillVal_3D(self):
        expected = numpy.array([self.pDict['fillVal']]*3)
        resDict = self.defaultOutFunc()
        numpy.testing.assert_array_equal(resDict['threeDnorm'][1,2,:], expected)

    def test_correctly_map_to_multiple_cells_simultaneously_2D(self):
        self.mapDict[(1,0)] = [((2,0), None), ((2,1), None)]
        self.mapDict[(1,1)] = [((2,1), None), ((2,2), None)]
        self.mapDict[(1,2)] = [((2,2), None), ((2,3), None)]
        self.SZA[2,1] = 100
        self.time[2,3] = self.toTAI93('12:00:00 01-06-2012')
        expected = numpy.array([self.twoDnorm[2,0], self.twoDnorm[2,2], self.twoDnorm[2,2]])
        resDict= self.defaultOutFunc()
        output = numpy.array([resDict['twoDnorm'][1,0], resDict['twoDnorm'][1,1], resDict['twoDnorm'][1,2]])
        numpy.testing.assert_array_almost_equal(output, expected)

    def test_correctly_map_to_multiple_cells_simultaneously_3D(self):
        self.mapDict[(1,0)] = [((2,0), None), ((2,1), None)]
        self.mapDict[(1,1)] = [((2,1), None), ((2,2), None)]
        self.mapDict[(1,2)] = [((2,2), None), ((2,3), None)]
        self.SZA[2,1] = 100
        self.time[2,3] = self.toTAI93('12:00:00 01-06-2012')
        expected = numpy.array([self.threeDnorm[2,0,:], self.threeDnorm[2,2,:], self.threeDnorm[2,2,:]])
        resDict = self.defaultOutFunc()
        output = numpy.array([resDict['threeDnorm'][1,0,:], resDict['threeDnorm'][1,1,:], resDict['threeDnorm'][1,2,:]])
        numpy.testing.assert_array_almost_equal(output, expected)
        
    @unittest.skip("Skipped until function can be rewritten to accomodate multiple inputs")
    def test_multidict_raises_appropriate_exception(self):
        # This function is currently busted and doesn't work with more
        # than one dictionary.  It raises an error if >1 dictionary is 
        # provided.  This test is in place until there is a convincing
        # reason to rewrite it to accept more than one dictionary.
        mapDict2 = dict(self.mapDict)
        dictList = [self.mapDict, mapDict2]
        self.assertRaises(NotImplementedError, self.defaultOutClass,
                          dictList, self.sixElGr, self.outFname, verbose=False)

    @unittest.skip("Skipped until function can be rewritten to accomodate multiple inputs")
    def test_multidict_one_zero_weight_2D_norm(self):
        mapDict2 = dict(self.mapDict)
        dictList = [self.mapDict, mapDict2]
        self.mapDict[(1,2)] = [((2,2), None)]
        mapDict2[(1,2)] = [((2,3), None)]
        self.SZA[2,3] = 100
        expected = self.twoDnorm[2,2]
        resDict = self.defaultOutClass(dictList, self.sixElGr, self.outFname, verbose=False)
        self.assertAlmostEqual(resDict['twoDnorm'][1,2], expected)
        
    @unittest.skip("Skipped until function can be rewritten to accomodate multiple inputs")
    def test_multi_dict_one_zero_weight_2D_log(self):
        mapDict2 = dict(self.mapDict)
        dictList = [self.mapDict, mapDict2]
        self.mapDict[(1,2)] = [((2,2), None)]
        mapDict2[(1,2)] = [((2,3), None)]
        self.SZA[2,3] = 100
        expected = self.twoDlog[2,2]
        resDict = self.defaultOutClass(dictList, self.sixElGr, self.outFname, verbose=False)
        self.assertAlmostEqual(resDict['twoDlog'][1,2], expected)
        
    @unittest.skip("Skipped until function can be rewritten to accomodate multiple inputs")
    def test_multi_dict_one_zero_weight_3D_norm(self):
        mapDict2 = dict(self.mapDict)
        dictList = [self.mapDict, mapDict2]
        self.mapDict[(1,2)] = [((2,2), None)]
        mapDict2[(1,2)] = [((2,3), None)]
        self.SZA[2,3] = 100
        expected = self.threeDnorm[2,2,:]
        resDict = self.defaultOutClass(dictList, self.sixElGr, self.outFname, verbose=False)
        numpy.testing.assert_array_almost_equal(resDict['threeDnorm'][1,2,:], expected)
        
    @unittest.skip("Skipped until function can be rewritten to accomodate multiple inputs")
    def test_multi_dict_one_zero_weight_3D_log(self):
        mapDict2 = dict(self.mapDict)
        dictList = [self.mapDict, mapDict2]
        self.mapDict[(1,2)] = [((2,2), None)]
        mapDict2[(1,2)] = [((2,3), None)]
        self.SZA[2,3] = 100
        expected = self.threeDlog[2,2,:]
        resDict = self.defaultOutClass(dictList, self.sixElGr, self.outFname, verbose=False)
        numpy.testing.assert_array_almost_equal(resDict['threeDlog'][1,2,:], expected)
        
    @unittest.skip("Skipped until function can be rewritten to accomodate multiple inputs")
    def test_multidict_two_valid_pix_2D_norm(self):
        mapDict2 = dict(self.mapDict)
        dictList = [self.mapDict, mapDict2]
        self.mapDict[(1,2)] = [((2,2), None)]
        mapDict2[(1,2)] = [((2,3), None)]
        self.twoDnorm[2,2] = 1.1
        self.twoDnorm[2,3] = 1.2
        expected = 1.15
        resDict = self.defaultOutClass(dictList, self.sixElGr, self.outFname, verbose=False)
        self.assertAlmostEqual(resDict['twoDnorm'][1,2], expected)
        
    @unittest.skip("Skipped until function can be rewritten to accomodate multiple inputs")
    def test_multidict_two_valid_pix_2D_log(self):
        mapDict2 = dict(self.mapDict)
        dictList = [self.mapDict, mapDict2]
        self.mapDict[(1,2)] = [((2,2), None)]
        mapDict2[(1,2)] = [((2,3), None)]
        self.twoDlog[2,2] = 2
        self.twoDlog[2,3] = 50
        expected = 10
        resDict = self.defaultOutClass(dictList, self.sixElGr, self.outFname, verbose=False)
        self.assertAlmostEqual(resDict['twoDlog'][1,2], expected)
        
    @unittest.skip("Skipped until function can be rewritten to accomodate multiple inputs")
    def test_multidict_two_valid_pix_3D_norm(self):
        mapDict2 = dict(self.mapDict)
        dictList = [self.mapDict, mapDict2]
        self.mapDict[(1,2)] = [((2,2), None)]
        mapDict2[(1,2)] = [((2,3), None)]
        expected = numpy.random.rand(3)
        self.threeDnorm[2,2,:] = expected+2
        self.threeDnorm[2,3,:] = expected-2
        resDict = self.defaultOutClass(dictList, self.sixElGr, self.outFname, verbose=False)
        numpy.testing.assert_array_almost_equal(resDict['threeDnorm'][1,2,:], expected)
        
    @unittest.skip("Skipped until function can be rewritten to accomodate multiple inputs")
    def test_multidict_two_valid_pix_3D_log(self):
        mapDict2 = dict(self.mapDict)
        dictList = [self.mapDict, mapDict2]
        self.mapDict[(1,2)] = [((2,2), None)]
        mapDict2[(1,2)] = [((2,3), None)]
        expected = numpy.random.rand(3)
        self.threeDlog[2,2,:] = expected*2
        self.threeDlog[2,3,:] = expected/2
        resDict = self.defaultOutClass(dictList, self.sixElGr, self.outFname, verbose=False)
        numpy.testing.assert_array_almost_equal(resDict['threeDlog'][1,2,:], expected)

    def test_output_file_is_netcdf(self):
        unused_result = self.defaultOutFunc()
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        # passes if no exception is raised during this step

    def test_output_file_contains_start_time(self):
        unused_result = self.defaultOutFunc()
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        expected = self.pDict['timeStart'].replace('_', ' ')
        self.assertEqual(self.fid.File_start_time, expected)

    def test_output_file_contains_stop_time(self):
        unused_result = self.defaultOutFunc()
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        expected = self.pDict['timeStop'].replace('_', ' ')
        self.assertEqual(self.fid.File_stop_time, expected)

    def test_output_file_contains_grid_name(self):
        unused_result = self.defaultOutFunc()
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        self.assertEqual(self.fid.Projection, 'latlon')
        
    def test_output_file_contains_tComp(self):
        unused_result = self.defaultOutFunc()
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        self.assertEqual(self.fid.Time_comparison_scheme, self.pDict['timeComparison'])

    def test_output_file_contains_notes(self):
        unused_result = self.defaultOutFunc()
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        self.assertEqual(self.fid.Notes, 'All values daytime with cutoff at  85.00')

    def test_output_file_containts_input_file_list_singledict(self):
        unused_result = self.defaultOutFunc()
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        self.assertEqual(self.fid.Input_files, 'foo.dat')

    @unittest.skip("Skipped until function can be rewritten to accomodate multiple inputs")
    def test_output_file_contains_input_file_list_multidict(self):
        mapDict2 = dict(self.mapDict)
        parser2 = fakeParser('bar.dat')
        mapDict2['parser'] = parser2
        dictList = [self.mapDict, mapDict2]
        unused_result = self.defaultOutFunc()
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        self.assertEqual(self.fid.Input_files, 'foo.dat bar.dat')

    def test_output_file_contains_gridParms(self):
        unused_result = self.defaultOutFunc()
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        keys = self.six_el_gridParms.keys()
        expected = dict()
        for key in keys:
            expected[key] = getattr(self.fid, key)
        self.assertDictEqual(self.six_el_gridParms, expected)

    def test_output_file_has_correct_dims_and_dimsizes(self):
        unused_result = self.defaultOutFunc()
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        expectedDims = {'row' : 2, 'col' : 3, 'layer' : 3, 'value' : 2}
        actualDims = {}
        for (name,dim) in self.fid.dimensions.items():
            actualDims[name] = len(dim)
        self.assertDictEqual(actualDims, expectedDims)

    def test_output_file_has_correct_variables(self):
        expected = ['foo', 'bar', 'baz', 'qux', 'spam']
        self.pDict['outFieldNames'] = expected
        newOutClass = out_geo.unweighted_filtered_MOPITT_avg_netCDF_out_func(self.pDict)
        unused_result = newOutClass(self.mapDict, self.sixElGr, 
                                    self.outFname, verbose=False,
                                    version = self.version)
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        self.assertItemsEqual(self.fid.variables.keys(), expected)

    def test_output_file_has_correctly_documented_units(self):
        unused_result = self.defaultOutFunc()
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        varNames = self.pDict['outFieldNames']
        expected = self.pDict['outUnits']
        output = [self.fid.variables[var].Units for var in varNames]
        self.assertListEqual(output, expected)

    def test_output_file_has_correctly_documented_fillVals(self):
        unused_result = self.defaultOutFunc()
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        varNames = self.pDict['outFieldNames']
        expected = [self.pDict['fillVal']]*5
        output = [self.fid.variables[var]._FillValue for var in varNames]
        self.assertListEqual(output, expected)

    def test_output_file_variables_correct_shape(self):
        unused_result = self.defaultOutFunc()
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        varNames = self.pDict['outFieldNames']
        expected = [(2,3), (2,3,3), (2,3), (2,3,3), (2,3,3,2)]
        output = [self.fid.variables[var].shape for var in varNames]
        self.assertListEqual(output, expected)

    def test_output_correctly_writes_valid_values_2D(self):
        self.mapDict[(1,2)] = [((2,2), None), ((2,3), None)]
        expected = numpy.random.rand()
        self.twoDnorm[2,2] = expected+1
        self.twoDnorm[2,3] = expected-1
        unused_result = self.defaultOutFunc()
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        self.assertAlmostEqual(self.fid.variables['twoDnorm'][1,2], expected)

    def test_output_correctly_writes_valid_values_3D(self):
        self.mapDict[(1,2)] = [((2,2), None), ((2,3), None)]
        expected = numpy.random.rand(3)
        self.threeDnorm[2,2,:] = expected+1
        self.threeDnorm[2,3,:] = expected-1
        unused_result = self.defaultOutFunc()
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        numpy.testing.assert_array_almost_equal(self.fid.variables['threeDnorm'][1,2,:], expected)

    def test_output_correctly_writes_valid_values_4D(self):
        self.mapDict[(1,2)] = [((2,2), None), ((2,3), None)]
        expected = numpy.random.rand(3,2)
        self.fourDcol[2,2,...] = expected+1
        self.fourDcol[2,3,...] = expected-1
        unused_result = self.defaultOutFunc()
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        numpy.testing.assert_array_almost_equal(self.fid.variables['fourDcol'][1,2,...], expected)

    def test_output_correctly_writes_fillVal_2D(self):
        self.mapDict[(1,2)] = [((2,3), None)]
        self.SZA[2,3] = 100
        unused_result = self.defaultOutFunc()
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        expected = self.pDict['fillVal']
        self.assertEqual(self.fid.variables['twoDnorm'][1,2].filled(), expected)

    def test_output_correctly_writes_full_fillVal_3D(self):
        self.mapDict[(1,2)] = [((2,3), None)]
        self.SZA[2,3] = 100
        unused_result = self.defaultOutFunc()
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        expected = numpy.array([self.pDict['fillVal']]*3)
        numpy.testing.assert_array_almost_equal(self.fid.variables['threeDnorm'][1,2,:].filled(), expected)

    def test_output_correctly_writes_partial_fillVal_3D(self):
        self.mapDict[(1,2)] = [((2,3), None)]
        self.threeDnorm[2,3,2] = numpy.NaN
        unused_result = self.defaultOutFunc()
        self.fid = netCDF4.Dataset(self.outFname, 'r')
        expected = self.threeDnorm[2,3,:]
        expected[2] = self.pDict['fillVal']
        numpy.testing.assert_array_almost_equal(self.fid.variables['threeDnorm'][1,2,:].filled(), expected)


'''
if __name__ == '__main__':
    foo = '__main__.TestNASAOmiL2GetGeoCorners.test_raises_IO_if_no_corner_file_invalid_dir'
    suite = unittest.defaultTestLoader.loadTestsFromName(foo)
    unittest.TextTestRunner(verbosity=2).run(suite)
'''
if __name__ == '__main__':
    unittest.main(verbosity=2)

