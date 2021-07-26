#include <math.h>
#include <stdint.h>
#include <stdio.h>

const double g = 9.80616; // Gravitational constant

enum pgcat { A, B, C, D, E, F, G };

enum roughness { URBAN, RURAL };

typedef struct {
    double x;
    double y;
} Components;

typedef struct {
    double u;
    double phi;
    //double temp;
    char pgcat;
} MetHour;

double sigma_y(double c, double d, double x) {
    double theta = 0.017453293 * (c - d * log(x));
    return 465.11628 * x * tan(theta);
}

double get_sigma_y(char pgcat, double x) {
    switch (pgcat) {
        case A:
            return sigma_y(24.1670, 2.5334, x);
        case B:
            return sigma_y(18.3330, 1.8096, x);
        case C:
            return sigma_y(12.5000, 1.0857, x);
        case D:
            return sigma_y(8.3330, 0.72382, x);
        case E:
            return sigma_y(6.2500, 0.54287, x);
        case F:
            return sigma_y(4.1667, 0.36191, x);
        default:
            return 0.0;
    }
}

double sigma_z(double a, double b, double x) { return a * pow(x, b); }

double sigma_za(double x) {
    double Sz;
    if (x <= 0.10) Sz = sigma_z(122.800, 0.94470, x);
    else if (x <= 0.15) Sz = sigma_z(158.080, 1.05420, x);
    else if (x <= 0.20) Sz = sigma_z(170.220, 1.09320, x);
    else if (x <= 0.25) Sz = sigma_z(179.520, 1.12620, x);
    else if (x <= 0.30) Sz = sigma_z(217.410, 1.26440, x);
    else if (x <= 0.40) Sz = sigma_z(258.890, 1.40940, x);
    else if (x <= 0.50) Sz = sigma_z(346.750, 1.72830, x);
    else if (x <= 3.11) Sz = sigma_z(453.850, 2.11660, x);
    else Sz = 5000.0;
    
    if (Sz > 5000.0) Sz = 5000.0;

    return Sz;
}

double sigma_zb(double x) {
    double Sz;
    if (x <= 0.20) Sz = sigma_z(90.673, 0.93198, x);
    else if (x <= 0.40) Sz = sigma_z(98.483, 0.98332, x);
    else Sz = sigma_z(109.300, 1.09710, x);

    if (Sz > 5000.0) Sz = 5000.0;

    return Sz;
}

double sigma_zc(double x) {
    double Sz = sigma_z(61.141, 0.91465, x);
    
    if (Sz > 5000.0) Sz = 5000.0;

    return Sz;
}

double sigma_zd(double x) {
    double Sz;
    if (x <= 0.30) Sz = sigma_z(34.459, 0.86974, x);
    else if (x <= 1.00) Sz = sigma_z(32.093, 0.81066, x);
    else if (x <= 3.00) Sz = sigma_z(32.093, 0.64403, x);
    else if (x <= 10.00) Sz = sigma_z(33.504, 0.60486, x);
    else if (x <= 30.00) Sz = sigma_z(36.650, 0.56589, x);
    else Sz = sigma_z(44.053, 0.51179, x);

    return Sz;
}

double sigma_ze(double x) {
    double Sz;
    if (x <= 0.10) Sz = sigma_z(24.260, 0.83660, x);
    else if (x <= 0.30) Sz = sigma_z(23.331, 0.81956, x);
    else if (x <= 1.00) Sz = sigma_z(21.628, 0.75660, x);
    else if (x <= 2.00) Sz = sigma_z(21.628, 0.63077, x);
    else if (x <= 4.00) Sz = sigma_z(22.534, 0.57154, x);
    else if (x <= 10.00) Sz = sigma_z(24.703, 0.50527, x);
    else if (x <= 20.00) Sz = sigma_z(26.970, 0.46713, x);
    else if (x <= 40.00) Sz = sigma_z(35.420, 0.37615, x);
    else Sz = sigma_z(47.618, 0.29592, x);
    
    return Sz;
}

double sigma_zf(double x) {
    double Sz;
    if (x <= 0.20) Sz = sigma_z(15.209, 0.81558, x);
    else if (x <= 0.70) Sz = sigma_z(14.457, 0.78407, x);
    else if (x <= 1.00) Sz = sigma_z(13.953, 0.68465, x);
    else if (x <= 2.00) Sz = sigma_z(13.953, 0.63227, x);
    else if (x <= 3.00) Sz = sigma_z(14.823, 0.54503, x);
    else if (x <= 7.00) Sz = sigma_z(16.187, 0.46490, x);
    else if (x <= 15.00) Sz = sigma_z(17.836, 0.41507, x);
    else if (x <= 30.00) Sz = sigma_z(22.651, 0.32681, x);
    else if (x <= 60.00) Sz = sigma_z(27.074, 0.27436, x);
    else Sz = sigma_z(34.219, 0.21716, x);

    return Sz;
}

double get_sigma_z(char pgcat, double x) {
    switch (pgcat) {
        case A:
            return sigma_za(x);
        case B:
            return sigma_zb(x);
        case C:
            return sigma_zc(x);
        case D:
            return sigma_zd(x);
        case E:
            return sigma_ze(x);
        case F:
            return sigma_zf(x);
        default:
            return 0.0;
    }
}

/*
Calculate effective wind speed, using "power law" method.
arguments:
uz_ref      [m/s]   wind speed of actual measurment
z           [m]     target elevation
z_ref       [m]     elevation of actual measurement
pgcat       []      Pasquill-Gifford stability category
roughness   []      "urban" or "rural"
returns:
Uz			[m/s] 	estimated wind speed at target elevation z
*/
double calc_uz(double uz_ref, double z, double z_ref, char pgcat, char roughness) {
    double p;
    
    if (roughness == URBAN) {
        switch (pgcat) {
            case A:
                p = 0.15;
                break;
            case B:
                p = 0.15;
                break;
            case C:
                p = 0.20;
                break;
            case D:
                p = 0.25;
                break;
            case E:
                p = 0.30;
                break;
            case F:
                p = 0.30;
                break;
            default:
                p = 0.0;
        }
    } else if (roughness == RURAL) {
        switch (pgcat) {
            case A:
                p = 0.07;
                break;
            case B:
                p = 0.07;
                break;
            case C:
                p = 0.10;
                break;
            case D:
                p = 0.15;
                break;
            case E:
                p = 0.35;
                break;
            case F:
                p = 0.55;
                break;
            default:
                p = 0.0;
        }
    }

    return uz_ref * pow(z/z_ref, p);
}

/*
Calculate the downwind (x) and crosswind (y) plume components from rectangular coordinates.
arguments:
components (pointer to struct)
  -> x      [km]    downwind plume receptor distance
  -> y      [m]     crosswind plume receptor distance
e_r         [m]     receptor easting
n_r         [m]     receptor northing
e_s         [m]     source (stack) easting
n_s         [m]     source (stack) northing
sin_phi     []      sine of wind direction in radians
cos_phi     []      cosine of wind direction in radians
*/
void wind_components(Components* components, double e_r, double n_r, double e_s, double n_s, double sin_phi, double cos_phi) {
    components->x = (-1.0*(e_r-e_s)*sin_phi - (n_r-n_s)*cos_phi) / 1000.0;
    components->y = (e_r-e_s)*cos_phi - (n_r-n_s)*sin_phi;
}

/*
Calculates the plume rise (dH) to add to stack height and a downwind plume offset (Xf), using Briggs model.
arguments:
dH    *     [m]     plume rise
Xf    *     [m]     plume rise offset
us          [m/s]   wind velocity at stack tip
vs          [m/s]   stack exit velocity
ds          [m]     stack tip diameter
Ts          [K]     stack tip temperature
Ta          [K]     ambient temperature
pgcat       []      Pasquill-Gifford stability category
*/
void plume_rise(double* dH, double* Xf, double us, double vs, double ds, double Ts, double Ta, char pgcat) {
    // Compute buoyancy flux
    double Fb = g * vs * ds * ds * (Ts - Ta) / (4.0 * Ts);
    // Calculate momentum flux
    double Fm = vs * vs * ds * ds * Ta / (4.0 * Ts);

    // Stable PG categories
    if (pgcat == E || pgcat == F) {
        double eta;
        if (pgcat == E) {
            eta = 0.020;
        } else {
            eta = 0.035;
        }
        double s = g * eta / Ta;
        double dT = 0.019582 * Ts * vs * sqrt(s);
        // Buoyancy dominated
        if ((Ts - Ta) >= dT) {
            *Xf = 2.0715 * us / sqrt(s);
            *dH = 2.6 * pow(Fb/(us*s), 0.333333333333);
            // Momentum dominated
        } else {
            *Xf = 0.0;
            // Calculate unstable/neutral and stable plume rise and take min
            double prUN = 3.0 * ds * vs / us;
            double prS = 1.5 * Fm / pow(us*sqrt(s), 0.333333333333);
            *dH = prUN < prS ? prUN : prS;
        }
    // Unstable or neutral PG categories
    } else {
        // Unstable or neutral
        if (Fb < 55.0) {
            // Check for buoyancy dominated or momentum
            double dT = 0.0297 * Ts * pow(vs, 0.333333333333) / pow(ds, 0.666666666667);
            // Buoyancy dominated
            if ((Ts - Ta) >= dT) {
                *Xf = 49.0 * pow(Fb, 0.625);
                *dH = 21.425 * pow(Fb, 0.75) / us;
                // Momentum dominated
            } else {
                *Xf = 0.0;
                *dH = 3.0 * ds * vs / us;
            }
        } else {
            double dT = 0.00575 * Ts * pow(vs, 0.666666666667) / pow(ds, 0.333333333333);
            if ((Ts - Ta) >= dT) {
                *Xf = 119.0 * pow(Fb, 0.4);
                *dH = 38.71 * pow(Fb, 0.6) / us;
            } else {
                *Xf = 0.0;
                *dH = 3.0 * ds * vs / us;
            }
        }
    }
}

/*
Calculate concentration at distance x along plume, at perpendicular offset y and height z
arguments:
x           [km]    receptor distance downwind along plume centreline
y           [m]     receptor perpendicular offset from plume centreline
z           [m]     receptor height
u_z         [m/s]   wind speed at stack exit
Q           [g/s]   pollutant mass emission rate
H           [m]		effective stack height (includes plume rise)
s_y         [m]     y plume sigma
s_z         [m]     z plume sigma
returns:
r_conc      [g/m3]  calculated receptor concentration
*/
double conc(double x, double y, double z, double u_z, double Q, double H, double s_y, double s_z) {
    // Early return if coordinate upwind, as concentration always zero
    if (x <= 0.0) return 0.0;

    double c1 = Q / (2.0 * M_PI * u_z * s_y * s_z);
    double c2 = pow(M_E, -1.0 * (z - H) * (z - H) / (2.0 * s_z * s_z));
    double c3 = pow(M_E, -1.0 * (z + H) * (z + H) / (2.0 * s_z * s_z));
    double c4 = pow(M_E, -1.0 * y * y / (2.0 * s_y * s_y));

    double r_conc = c1 * (c2 + c3) * c4; // g/m3
    
    if (isnan(r_conc)) return 0.0;
    
    return r_conc;
}

uint32_t cr_to_linear(uint32_t col, uint32_t row, uint32_t x_points, uint32_t y_points) {
    uint32_t row_offset = y_points - 1;
    return x_points * (row_offset - row) + col;
}

void iter_disp(double* rgrid, double* hgrid) {
    double source_height = 10.0;
    double source_temp = 100.0;
    double source_emission = 1.0;

    double source_velocity = 10.0;
    double source_diameter = 0.5;

    int hours = 1;

    int x_min = -2500;
    int x_max = 2500;
    int y_min = -2500;
    int y_max = 2500;
    int z_min = 0;
    int z_max = 1000;
    int x_spacing = 20;
    int y_spacing = 20;
    int z_spacing = 10;
    
    uint32_t x_points = (x_max - x_min) / x_spacing;
    uint32_t y_points = (y_max - y_min) / y_spacing;
    uint32_t z_points = (z_max - z_min) / z_spacing;

    MetHour metline;
    metline.u = 2.0;
    metline.phi = 130.0 * M_PI / 180.0;
    metline.pgcat = A;
    
    double AMBIENT_TEMP = 293.15; // Fixed ambient temperature [K] (20 C)
    char roughness = URBAN;

    for (int _; _ < hours; _++) {
        // Calculate effective wind speed at stack tip (user specified wind speed is for 10 m)
        double Uz = calc_uz(metline.u, source_height, 10.0, metline.pgcat, roughness);
        
        // Calculate plume rise using Briggs equations
        double Ts = source_temp + 273.15;
        double dH, Xf;
        plume_rise(&dH, &Xf, Uz, source_velocity, source_diameter, Ts, AMBIENT_TEMP, metline.pgcat);
        
        double H = source_height + dH;
        double Q = source_emission;

        double sin_phi = sin(metline.phi);
        double cos_phi = cos(metline.phi);

        // Calculate concentrations for plan view grid (fixed grid height of 0 m)
        double Yr = y_min;
        for (uint32_t y = 0; y < y_points; y++) {
            double Xr = x_min;
            for (uint32_t x = 0; x < x_points; x++) {
                if (Uz > 0.5) {
                    double xx = (-1.0 * Xr * sin_phi - Yr * cos_phi - Xf) / 1000.0;
                    double yy = Xr * cos_phi - Yr * sin_phi;
                    
                    double sig_y = get_sigma_y(metline.pgcat, xx);
                    double sig_z = get_sigma_z(metline.pgcat, xx);
                    
                    uint32_t i = cr_to_linear(x, y, x_points, y_points);
                    rgrid[i] += (conc(xx, yy, 0.0, Uz, Q, H, sig_y, sig_z) / (double)hours);
                }
                Xr += x_spacing;
            }
            Yr += y_spacing;
        }

        // Calculate concentrations for 2d slice showing height profile along plume
        double Zr = z_min;
        for (uint32_t z = 0; z < z_points; z++) {
            double Xr = x_min;
            for (uint32_t x = 0; x < x_points; x++) {
                if (Uz > 0.5) {
                    double xx = (Xr - Xf) / 1000.0;
                    
                    double sig_y = get_sigma_y(metline.pgcat, xx);
                    double sig_z = get_sigma_z(metline.pgcat, xx);
                    
                    uint32_t i = cr_to_linear(x, z, x_points, z_points);
                    hgrid[i] += (conc(xx, 0.0, Zr, Uz, Q, H, sig_y, sig_z) / (double)hours);
                }
                Xr += x_spacing;
            }
            Zr += z_spacing;
        }
    }
}
