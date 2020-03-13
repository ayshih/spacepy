
import unittest
import spacepy.coordinates as spc
import glob
import os
import datetime
from numpy import array
import numpy as np
from spacepy.time import Ticktock
try:
    import spacepy.irbempy as ib
except ImportError:
    pass #tests will fail, but won't bring down the entire suite
try:
    from astropy.coordinates import SkyCoord
except ImportError:
    pass #tests will fail, but won't bring down the entire suite

__all__ = ['coordsTest']


class coordsTest(unittest.TestCase):
    def setUp(self):
        #super(tFunctionTests, self).setUp()
        try:
            self.cvals = spc.Coords([[1,2,4],[1,2,2]], 'GEO', 'car')
        except ImportError:
            pass #tests will fail, but won't bring down the entire suite

    def tearDown(self):
        #super(tFunctionTests, self).tearDown()
        pass

    def test_coords(self):
        """Coords should create and do simple conversions"""
        np.testing.assert_equal([1,1], self.cvals.x)
        np.testing.assert_equal([2,2], self.cvals.y)
        np.testing.assert_equal([4,2], self.cvals.z)
        self.cvals.ticks = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:00:00'], 'ISO') # add ticktock
        newcoord = self.cvals.convert('GSM', 'sph')

    def test_append(self):
        c2 = spc.Coords([[6,7,8],[9,10,11]], 'GEO', 'car')
        actual = self.cvals.append(c2)
        expected = [[1,2,4],[1,2,2],[6,7,8],[9,10,11]]
        np.testing.assert_equal(expected, actual.data.tolist())

    def test_slice(self):
        expected = spc.Coords([1,2,4], 'GEO', 'car')
        np.testing.assert_equal(expected.data, self.cvals[0].data)

    def test_slice_with_ticks(self):
        self.cvals.ticks = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:00:00'], 'ISO')
        expected = spc.Coords([1,2,4], 'GEO', 'car')
        np.testing.assert_equal(expected.data, self.cvals[0].data)

    def test_to_skycoord_without_ticks(self):
        with self.assertRaises(ValueError) as cm:
            sc = self.cvals.to_skycoord()

    def test_to_skycoord_with_ticks(self):
        self.cvals.ticks = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:00:00'], 'ISO') # add ticktock

        sc = self.cvals.to_skycoord()

        assert isinstance(sc, SkyCoord)
        assert sc.frame.name == 'itrs'

        # Check that the data was loaded correctly
        sc_data = sc.cartesian.xyz.to_value('m').T
        np.testing.assert_allclose(sc_data, self.cvals.data * self.cvals.Re, rtol=1e-10)

        # Check that the time was loaded correctly
        sc_obstime = sc.obstime.utc.mjd
        np.testing.assert_allclose(sc_obstime, self.cvals.ticks.MJD, rtol=0, atol=1e-9)

    def test_to_skycoord_with_ticks_and_conversion(self):
        self.cvals.ticks = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:00:00'], 'ISO') # add ticktock

        # Convert to a frame other than GEO before calling to_skycoord()
        non_geo = self.cvals.convert('GSE', 'sph')

        sc = non_geo.to_skycoord()

        assert isinstance(sc, SkyCoord)
        assert sc.frame.name == 'itrs'

        # Check that the data converts back correctly
        sc_data = sc.cartesian.xyz.to_value('m').T
        np.testing.assert_allclose(sc_data, self.cvals.data * self.cvals.Re, rtol=1e-10)

        # Check that the time was loaded correctly
        sc_obstime = sc.obstime.utc.mjd
        np.testing.assert_allclose(sc_obstime, self.cvals.ticks.MJD, rtol=0, atol=1e-9)

    def test_from_skycoord(self):
        sc = SkyCoord(x=[1, 2], y=[4, 5], z=[7, 8], unit='Mm', frame='itrs',
                      obstime=['2001-02-03T04:00:00', '2005-06-07T08:00:00'])

        coords = spc.Coords.from_skycoord(sc)

        # Check that the data was loaded correctly
        sc_data = sc.cartesian.xyz.to_value('m').T
        np.testing.assert_allclose(coords.data * coords.Re, sc_data, rtol=1e-10)

        # Check that the time was loaded correctly
        sc_obstime = sc.obstime.utc.mjd
        np.testing.assert_allclose(coords.ticks.MJD, sc_obstime, rtol=0, atol=1e-9)

    def test_from_skycoord_with_conversion(self):
        sc = SkyCoord(x=[1, 2], y=[4, 5], z=[7, 8], unit='Mm', frame='itrs',
                      obstime=['2001-02-03T04:00:00', '2005-06-07T08:00:00'])

        # Convert to a frame other than GEO before calling from_skycoord()
        coords = spc.Coords.from_skycoord(sc.gcrs)

        # Check that the data was loaded correctly
        sc_data = sc.cartesian.xyz.to_value('m').T
        np.testing.assert_allclose(coords.data * coords.Re, sc_data, rtol=1e-10)

        # Check that the time was loaded correctly
        sc_obstime = sc.obstime.utc.mjd
        np.testing.assert_allclose(coords.ticks.MJD, sc_obstime, rtol=0, atol=1e-9)


if __name__ == "__main__":
    ## suite = unittest.TestLoader().loadTestsFromTestCase(coordsTest)
    ## unittest.TextTestRunner(verbosity=2).run(suite)

    unittest.main()
