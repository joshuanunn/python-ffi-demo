import math
import unittest

from ..core import *


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
        methour = new_methour()
        
        source.height = 10.0
        source.temp = 100.0
        source.emission = 1.0
        source.velocity = 10.0
        source.diameter = 0.5

        methour.hours = 1
        methour.wspd = 2.0
        methour.wdir = 130.0 * math.pi / 180.0
        methour.pgcat = PGCATS['A']

        r_grid, h_grid = new_grids(domain)

        iter_disp(r_grid, h_grid, domain, source, methour)
        
        self.assertAlmostEqual(r_grid[23 * 250 + 14], 6.366502967443e-08)
        self.assertAlmostEqual(h_grid[41 * 250 + 181], 3.963714618520e-07)


class TestCreateImage(unittest.TestCase):
    """ Testcase for create_image function. """

    def test_1(self):

        domain = new_domain('MEDIUM')
        source = new_source()
        methour = new_methour()

        source.height = 10.0
        source.temp = 100.0
        source.emission = 1.0
        source.velocity = 10.0
        source.diameter = 0.5

        methour.hours = 1
        methour.wspd = 2.0
        methour.wdir = 130.0 * math.pi / 180.0
        methour.pgcat = PGCATS['A']

        r_grid, h_grid = new_grids(domain)
        
        iter_disp(r_grid, h_grid, domain, source, methour)
        
        r_grid_image, h_grid_image = new_images(domain)
        create_image(r_grid_image, r_grid, domain, 'PLAN')
        create_image(h_grid_image, h_grid, domain, 'SECTION')

        self.assertAlmostEqual(r_grid_image[48 * 250 + 17], 6)
        self.assertAlmostEqual(r_grid_image[87 * 250 + 80], 8)


if __name__ == '__main__':
    unittest.main()
