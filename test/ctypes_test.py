import ctypes
import math
import numpy as np
import unittest

from ctypes import byref, POINTER, Structure
from ctypes import c_char, c_double, c_int, c_ubyte
from pathlib import Path


PROJECT_ROOT = Path(__file__).parent.resolve() / '..'

C_UCHAR_SS = POINTER(POINTER(c_ubyte))
C_DOUBLE_SS = POINTER(POINTER(c_double))

PGCATS = {'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6}
RESOLUTION = {'LOW': 0, 'MEDIUM': 1, 'HIGH': 2, 'EXTREME': 3}
GRIDTYPE = {'PLAN': 0, 'SECTION': 1}
ROUGHNESS = {'urban': 0, 'rural': 1}

class Components(Structure):
    _fields_ = ('x', c_double), ('y', c_double)


class Domain(Structure):
    _fields_ = [
        ('x_min', c_int), ('x_max', c_int),
        ('y_min', c_int), ('y_max', c_int),
        ('z_min', c_int), ('z_max', c_int),

        ('x_spacing', c_int),
        ('y_spacing', c_int),
        ('z_spacing', c_int),

        ('x_points', c_int),
        ('y_points', c_int),
        ('z_points', c_int)
    ]

class Source(Structure):
    _fields_ = [
        ('x', c_double), ('y', c_double),
        ('height', c_double),
        ('diameter', c_double),
        ('velocity', c_double),
        ('temp', c_double),
        ('emission', c_double)
    ]

# Import compiled c code using ctypes
_disperse = ctypes.CDLL(PROJECT_ROOT / 'build' / 'ctypes' / 'disperse.so')

### Setup ctypes functions ###

# double get_sigma_y(char pgcat, double x)
_disperse.get_sigma_y.argtypes = [c_char, c_double]
_disperse.get_sigma_y.restype = c_double

# double get_sigma_z(char pgcat, double x)
_disperse.get_sigma_z.argtypes = [c_char, c_double]
_disperse.get_sigma_z.restype = c_double

# double calc_uz(double uz_ref, double z, double z_ref, char pgcat, char roughness)
_disperse.calc_uz.argtypes = [c_double, c_double, c_double, c_char, c_char]
_disperse.calc_uz.restype = c_double

# void wind_components(double e_r, double n_r, double e_s, double n_s, double sin_phi, double cos_phi)
_disperse.wind_components.argtypes = [POINTER(Components), c_double, c_double, c_double, c_double, c_double, c_double]
_disperse.wind_components.restype = None

# void plume_rise(double* dH, double* Xf, double us, double vs, double ds, double Ts, double Ta, char pgcat)
_disperse.plume_rise.argtypes = [POINTER(c_double), POINTER(c_double), c_double, c_double, c_double, c_double, c_double, c_char]
_disperse.plume_rise.restype = None

# double conc(double x, double y, double z, double u_z, double Q, double H, double s_y, double s_z)
_disperse.conc.argtypes = [c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double]
_disperse.conc.restype = c_double

# void iter_disp(double* rgrid, double* hgrid, Domain* domain, Source* source)
_disperse.iter_disp.argtypes = [C_DOUBLE_SS, C_DOUBLE_SS, POINTER(Domain), POINTER(Source)]
_disperse.iter_disp.restype = None

# void create_image(unsigned char* destgrid, double* grid, Domain* domain, bool vertical)
_disperse.create_image.argtypes = [C_UCHAR_SS, C_DOUBLE_SS, POINTER(Domain), c_int]
_disperse.create_image.restype = None

# Domain new_domain(int resolution)
_disperse.new_domain.argtypes = [c_int]
_disperse.new_domain.restype = Domain

# Source new_source()
_disperse.new_source.argtypes = None
_disperse.new_source.restype = Source

def get_sigma_y(pgcat, x):
    return _disperse.get_sigma_y(c_char(PGCATS[pgcat]), c_double(x))


def get_sigma_z(pgcat, x):
    return _disperse.get_sigma_z(c_char(PGCATS[pgcat]), c_double(x))


def calc_uz(uzref, z, zref, pgcat, roughness):
    return _disperse.calc_uz(c_double(uzref), c_double(z), c_double(zref), c_char(PGCATS[pgcat]), c_char(ROUGHNESS[roughness]))


def wind_components(wind_components, e_r, n_r, e_s, n_s, sin_phi, cos_phi):
    return _disperse.wind_components(wind_components, c_double(e_r), c_double(n_r), c_double(e_s), c_double(n_s), c_double(sin_phi), c_double(cos_phi))


def plume_rise(dH, Xf, us, vs, ds, Ts, Ta, pgcat):
    return _disperse.plume_rise(dH, Xf, c_double(us), c_double(vs), c_double(ds), c_double(Ts), c_double(Ta), c_char(PGCATS[pgcat]))


def conc(x, y, z, u_z, Q, H, s_y, s_z):
    return _disperse.conc(c_double(x), c_double(y), c_double(z), c_double(u_z), c_double(Q), c_double(H), c_double(s_y), c_double(s_z))


def new_domain(resolution):
    return _disperse.new_domain(c_int(RESOLUTION[resolution]))


def new_source():
    return _disperse.new_source()


def iter_disp(r_grid_np, h_grid_np, domain, source):
    return _disperse.iter_disp(r_grid_np.ctypes.data_as(C_DOUBLE_SS), h_grid_np.ctypes.data_as(C_DOUBLE_SS), byref(domain), byref(source))


def create_image(png_grid_np, grid_np, domain, gridtype):
    return _disperse.create_image(png_grid_np.ctypes.data_as(C_UCHAR_SS), grid_np.ctypes.data_as(C_DOUBLE_SS), byref(domain), c_int(GRIDTYPE[gridtype]))


def new_grids(domain):
    r_grid = np.zeros((domain.x_points * domain.y_points), dtype=np.float64)
    h_grid = np.zeros((domain.x_points * domain.z_points), dtype=np.float64)
    return r_grid, h_grid


def new_images(domain):
    r_image = np.zeros((domain.x_points * domain.y_points), dtype=np.uint8)
    h_image = np.zeros((domain.x_points * domain.z_points), dtype=np.uint8)
    return r_image, h_image


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


class TestWindComponents(unittest.TestCase):
    """ Testcase for wind_components function. """

    def test_1(self):
        
        compenents = Components()
        
        source_x = -2.0
        source_y = -3.0
        
        sin_phi = math.sin(math.radians(200.0))
        cos_phi = math.cos(math.radians(200.0))
        
        wind_components(compenents, 10.0, 10.0, source_x, source_y, sin_phi, cos_phi)
        
        self.assertAlmostEqual(compenents.x, 0.016320245790)
        self.assertAlmostEqual(compenents.y, -6.8300495861972)        


class TestCalcUz(unittest.TestCase):
    """ Testcase for calc_uz function. """

    def test_1(self):
        
        uzref = 3.5
        z = 100.0
        zref = 10.0
        pgcat = 'D'
        roughness = 'rural'

        u_adj = calc_uz(uzref, z, zref, pgcat, roughness)

        self.assertAlmostEqual(u_adj, 4.943881406180)
        
        uzref = 10.0
        z = 50.0
        zref = 45.0
        pgcat = 'A'
        roughness = 'urban'

        u_adj = calc_uz(uzref, z, zref, pgcat, roughness)

        self.assertAlmostEqual(u_adj, 10.159296222811)


class TestPlumeRise(unittest.TestCase):
    """ Testcase for plume_rise function. """

    def test_1(self):
        # Example from:
        # https://ceprofs.civil.tamu.edu/qying/cven301_fall2014_arch/lecture7_c.pdf
        vs = 20.0   # m/s
        ds = 5.0    # m
        U = 6.0     # m/s
        Ts = 400.0  # K
        Ta = 280.0  # K
        pgcat = 'D'

        dH = c_double(0.0)
        Xf = c_double(0.0)
        
        plume_rise(byref(dH), byref(Xf), U, vs, ds, Ts, Ta, pgcat)
        
        self.assertAlmostEqual(dH.value, 223.352113600373)
        self.assertAlmostEqual(Xf.value, 1264.034881130080)        


class TestConc(unittest.TestCase):
    """ Testcase for conc function. """

    def test_1(self):
        # Example from:
        # http://faculty.washington.edu/markbenj/CEE357/CEE%20357%20air%20dispersion%20models.pdf

        x = 0.5          # 500 m downwind
        y = 0.0          # along plume centreline
        z = 0.0          # ground level
        u_z = 6.0        # 6 m/s wind speed at height of 50 m
        pgcat = 'D'      # Neutral stability

        # Source centred on (0,0), height 50 m, 10 g/s mass emission rate
        Q = 10.0         # source.emission
        H = 50.0         # source.height

        # Calculate concentration at (x,y,z) == 19.2 ug/m3
        s_y = get_sigma_y(pgcat, x)
        s_z = get_sigma_z(pgcat, x)
        
        test_conc = conc(x, y, z, u_z, Q, H, s_y, s_z)
        
        self.assertAlmostEqual(test_conc, 1.917230120488e-05)


class TestIterDisp(unittest.TestCase):
    """ Testcase for iter_disp function. """

    def test_1(self):
        
        domain = new_domain('MEDIUM')
        source = new_source()

        source.height = 10.0
        source.temp = 100.0
        source.emission = 1.0
        source.velocity = 10.0
        source.diameter = 0.5

        r_grid, h_grid = new_grids(domain)

        iter_disp(r_grid, h_grid, domain, source)
        
        self.assertAlmostEqual(r_grid[23 * 250 + 14], 6.366502967443e-08)
        self.assertAlmostEqual(h_grid[18 * 250 + 181], 4.086979994894e-07)


class TestCreateImage(unittest.TestCase):
    """ Testcase for create_image function. """

    def test_1(self):

        domain = new_domain('MEDIUM')
        source = new_source()

        source.height = 10.0
        source.temp = 100.0
        source.emission = 1.0
        source.velocity = 10.0
        source.diameter = 0.5

        r_grid, h_grid = new_grids(domain)
        
        iter_disp(r_grid, h_grid, domain, source)
        
        r_grid_image, h_grid_image = new_images(domain)
        create_image(r_grid_image, r_grid, domain, 'PLAN')
        create_image(h_grid_image, h_grid, domain, 'SECTION')

        self.assertAlmostEqual(r_grid_image[48 * 250 + 17], 6)
        self.assertAlmostEqual(r_grid_image[87 * 250 + 80], 8)


if __name__ == '__main__':
    unittest.main()
