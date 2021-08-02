import ctypes
import math
import numpy as np
from ctypes import byref, POINTER, Structure
from ctypes import c_char, c_double, c_int, c_ubyte
from pathlib import Path


MODULE_ROOT = Path(__file__).parent.resolve()

C_UCHAR_SS = POINTER(POINTER(c_ubyte))
C_DOUBLE_SS = POINTER(POINTER(c_double))

PGCATS = {'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6}
RESOLUTION = {'LOW': 0, 'MEDIUM': 1, 'HIGH': 2, 'EXTREME': 3}
GRIDTYPE = {'PLAN': 0, 'SECTION': 1}
ROUGHNESS = {'urban': 0, 'rural': 1}


class MetHour(Structure):
    _fields_ = [
        ('hours', c_int),
        ('wspd', c_double),
        ('wdir', c_double),
        #('temp', c_double),
        ('pgcat', c_int),
    ]

class Domain(Structure):
    _fields_ = [
        ('xr_min', c_int), ('xr_max', c_int),
        ('yr_min', c_int), ('yr_max', c_int),
        ('xh_min', c_int), ('xh_max', c_int),
        ('zh_min', c_int), ('zh_max', c_int),

        ('xr_spacing', c_int),
        ('yr_spacing', c_int),
        ('xh_spacing', c_int),
        ('zh_spacing', c_int),

        ('xr_points', c_int),
        ('yr_points', c_int),
        ('xh_points', c_int),
        ('zh_points', c_int)
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
_disperse = ctypes.CDLL(MODULE_ROOT / 'disperse.so')

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

# void plume_rise(double* dH, double* Xf, double us, double vs, double ds, double Ts, double Ta, char pgcat)
_disperse.plume_rise.argtypes = [POINTER(c_double), POINTER(c_double), c_double, c_double, c_double, c_double, c_double, c_char]
_disperse.plume_rise.restype = None

# double conc(double x, double y, double z, double u_z, double Q, double H, double s_y, double s_z)
_disperse.conc.argtypes = [c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double]
_disperse.conc.restype = c_double

# void iter_disp(double* rgrid, double* hgrid, Domain* domain, Source* source, MetHour* met)
_disperse.iter_disp.argtypes = [C_DOUBLE_SS, C_DOUBLE_SS, POINTER(Domain), POINTER(Source), POINTER(MetHour)]
_disperse.iter_disp.restype = None

# void create_image(unsigned char* destgrid, double* grid, Domain* domain, bool vertical)
_disperse.create_image.argtypes = [C_UCHAR_SS, C_DOUBLE_SS, POINTER(Domain), c_int]
_disperse.create_image.restype = None

# MetHour new_methour()
_disperse.new_methour.argtypes = None
_disperse.new_methour.restype = MetHour

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


def new_methour():
    return _disperse.new_methour()


def new_domain(resolution):
    return _disperse.new_domain(c_int(RESOLUTION[resolution]))


def new_source():
    return _disperse.new_source()


def iter_disp(r_grid_np, h_grid_np, domain, source, methour):
    return _disperse.iter_disp(r_grid_np.ctypes.data_as(C_DOUBLE_SS), h_grid_np.ctypes.data_as(C_DOUBLE_SS), byref(domain), byref(source), byref(methour))


def create_image(png_grid_np, grid_np, domain, gridtype):
    return _disperse.create_image(png_grid_np.ctypes.data_as(C_UCHAR_SS), grid_np.ctypes.data_as(C_DOUBLE_SS), byref(domain), c_int(GRIDTYPE[gridtype]))


def new_grids(domain):
    r_grid = np.zeros((domain.xr_points * domain.yr_points), dtype=np.float64)
    h_grid = np.zeros((domain.xh_points * domain.zh_points), dtype=np.float64)
    return r_grid, h_grid


def new_images(domain):
    r_image = np.zeros((domain.xr_points * domain.yr_points), dtype=np.uint8)
    h_image = np.zeros((domain.xh_points * domain.zh_points), dtype=np.uint8)
    return r_image, h_image
