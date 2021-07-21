#include <math.h>

enum pgcat { A, B, C, D, E, F, G };

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
