import ctypes
import math
import unittest

from ctypes import c_char, c_float, c_double
from pathlib import Path


PROJECT_ROOT = Path(__file__).parent.resolve() / '..'

PGCATS = {'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6}

# Import compiled c code using ctypes
_disperse = ctypes.CDLL(PROJECT_ROOT / 'build' / 'ctypes' / 'disperse.so')

# Setup get_sigma_y function using ctypes
# double get_sigma_y(char pgcat, double x)
_disperse.get_sigma_y.argtypes = [c_char, c_double]
_disperse.get_sigma_y.restype = c_double

# Setup get_sigma_z function using ctypes
# double get_sigma_z(char pgcat, double x)
_disperse.get_sigma_z.argtypes = [c_char, c_double]
_disperse.get_sigma_z.restype = c_double


def get_sigma_y(pgcat, x):
    return _disperse.get_sigma_y(c_char(PGCATS[pgcat]), c_double(x))


def get_sigma_z(pgcat, x):
    return _disperse.get_sigma_z(c_char(PGCATS[pgcat]), c_double(x))


class TestSigmaY(unittest.TestCase):
    """ Range of testcases for get_sigma_y function. """

    def test_1(self):
        # stability class D, 0.5km downwind, example from:
        # http://faculty.washington.edu/markbenj/CEE357/CEE%20357%20air%20dispersion%20models.pdf
        self.assertAlmostEqual(get_sigma_y('D', 0.5), 36.146193496038)
    
    def test_2(self):
        # stability class A, 0.997km downwind
        self.assertAlmostEqual(get_sigma_y('A', 0.997), 208.157523627706)
    
    def test_3(self):
        # stability class B, 12.345m downwind
        self.assertAlmostEqual(get_sigma_y('B', 0.012345), 2.835970876943)
    
    def test_4(self):
        # stability class C, 27.85km downwind
        self.assertAlmostEqual(get_sigma_y('C', 27.85), 2025.696103458910)

    def test_5(self):
        # stability class D, 5.78m upwind
        self.assertTrue(math.isnan(get_sigma_y('D', -0.00578)))

    def test_6(self):
        # stability class E, 445m downwind
        self.assertAlmostEqual(get_sigma_y('E', 0.445), 24.275915684479)

    def test_7(self):
        # stability class F, 7.5558km downwind
        self.assertAlmostEqual(get_sigma_y('F', 7.5558), 210.931775211803)


class TestSigmaZ(unittest.TestCase):
    """ Range of testcases for get_sigma_z function. """

    def test_1(self):
        # stability class D, 0.5km downwind, example from:
        # http://faculty.washington.edu/markbenj/CEE357/CEE%20357%20air%20dispersion%20models.pdf
        self.assertAlmostEqual(get_sigma_z('D', 0.5), 18.296892641654)
    
    def test_2(self):
        # stability class D, 5.78m upwind
        self.assertTrue(math.isnan(get_sigma_z('D', -0.00578)))
    
    def test_3(self):
        # stability class A, 50m downwind
        self.assertAlmostEqual(get_sigma_z('A', 0.05), 7.246283645973)
        
    def test_4(self):
        # stability class A, 270m downwind
        self.assertAlmostEqual(get_sigma_z('A', 0.27), 41.523682287423)
        
    def test_5(self):
        # stability class A, 2.86km downwind
        self.assertAlmostEqual(get_sigma_z('A', 2.86), 4196.204889704382)
        
    def test_6(self):
        # stability class A, 54km downwind
        self.assertAlmostEqual(get_sigma_z('A', 54.0), 5000.0)
        
    def test_7(self):
        # stability class B, 50m downwind
        self.assertAlmostEqual(get_sigma_z('B', 0.05), 5.558326444834)
        
    def test_8(self):
        # stability class B, 270m downwind
        self.assertAlmostEqual(get_sigma_z('B', 0.27), 27.177523893054)
        
    def test_9(self):
        # stability class B, 2.86km downwind
        self.assertAlmostEqual(get_sigma_z('B', 2.86), 346.177898273921)
        
    def test_10(self):
        # stability class B, 54km downwind
        self.assertAlmostEqual(get_sigma_z('B', 54.0), 5000.0)
        
    def test_11(self):
        # stability class C, 50m downwind
        self.assertAlmostEqual(get_sigma_z('C', 0.05), 3.947711911749)
        
    def test_12(self):
        # stability class C, 270m downwind
        self.assertAlmostEqual(get_sigma_z('C', 0.27), 18.459902569036)
        
    def test_13(self):
        # stability class C, 2.86km downwind
        self.assertAlmostEqual(get_sigma_z('C', 2.86), 159.862915743170)
        
    def test_14(self):
        # stability class C, 54km downwind
        self.assertAlmostEqual(get_sigma_z('C', 54.0), 2348.910612301645)
        
    def test_15(self):
        # stability class D, 50m downwind
        self.assertAlmostEqual(get_sigma_z('D', 0.05), 2.545334368597)
        
    def test_16(self):
        # stability class D, 270m downwind
        self.assertAlmostEqual(get_sigma_z('D', 0.27), 11.034101898944)
        
    def test_17(self):
        # stability class D, 2.86km downwind
        self.assertAlmostEqual(get_sigma_z('D', 2.86), 63.142784897226)
        
    def test_18(self):
        # stability class D, 54km downwind
        self.assertAlmostEqual(get_sigma_z('D', 54.0), 339.310493995667)
        
    def test_19(self):
        # stability class E, 50m downwind
        self.assertAlmostEqual(get_sigma_z('E', 0.05), 1.979015073784)
        
    def test_20(self):
        # stability class E, 270m downwind
        self.assertAlmostEqual(get_sigma_z('E', 0.27), 7.978143439122)
        
    def test_21(self):
        # stability class E, 2.86km downwind
        self.assertAlmostEqual(get_sigma_z('E', 2.86), 41.083717338729)
        
    def test_22(self):
        # stability class E, 54km downwind
        self.assertAlmostEqual(get_sigma_z('E', 54.0), 155.031915174584)
        
    def test_23(self):
        # stability class F, 50m downwind
        self.assertAlmostEqual(get_sigma_z('F', 0.05), 1.321315762922)
        
    def test_24(self):
        # stability class F, 270m downwind
        self.assertAlmostEqual(get_sigma_z('F', 0.27), 5.178781257565)
        
    def test_25(self):
        # stability class F, 2.86km downwind
        self.assertAlmostEqual(get_sigma_z('F', 2.86), 26.282658227590)
        
    def test_26(self):
        # stability class F, 54km downwind
        self.assertAlmostEqual(get_sigma_z('F', 54.0), 80.882017663045)


if __name__ == '__main__':
    unittest.main()
