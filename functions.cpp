#include <cmath>

#include "functions.h"

#define pi (2 * M_PI)
#define s1 sin(pi * x)
#define c1 cos(pi * x)
#define s2 sin(pi * y)
#define c2 cos(pi * y)
#define v1 (solution_v1(x, y, t))
#define v2 (solution_v2(x, y, t))
#define rho (solution_rho(x, y, t))

static double mu = 0.1;
static double p_coeff = 1;

void set_mu(double new_mu)
{
    mu = new_mu;
}

void set_p_coeff(double new_coeff)
{
    p_coeff = new_coeff;
}

double f_0(double x, double y, double t)
{
    (void) x;
    (void) y;
    (void) t;
    //return 1. + v1 - v2 + pi * c1 * s2 * exp(t) + pi * s1 * c2 * exp(-t);
    return 1 + v1 - v2 + pi * (c1 * s2 * exp(-t) + s1 * c2 * exp(t));
}

double f_1(double x, double y, double t)
{
    (void) x;
    (void) y;
    (void) t;
    //return v1 + v1 * pi * c1 * s2 * exp(t) + v2 * pi * s1 * c2 * exp(t) + 1 * rho - 0.1 / rho * pi * pi * (-4. / 3. * v1 - v1 + 1. / 3. * c1 * c2 * exp(-t));
    return -v1 + p_coeff + pi * v1 * c1 * s2 * exp(-t) + v2 * pi * s1 * c2 * exp(-t) + 0.1 / rho * 7. / 3. * pi * pi * v1 - 0.1 / rho / 3. * pi * pi * c1 * c2 * exp(t);
}

double f_2(double x, double y, double t)
{
    (void) x;
    (void) y;
    (void) t;
    //return -v2 + v1 * pi * c1 * s2 * exp(-t) + v2 * pi * s1 * c2 * exp(-t) - 1 * rho - 0.1 / rho * pi * pi * (-4. / 3. * v2 - v2 + 1. / 3. * c1 * c2 * exp(t));
    return v2 + -p_coeff + v1 * pi * c1 * s2 * exp(t) + v2 * pi * s1 * c2 * exp(t) + 0.1 / rho * 7. / 3. * pi * pi * v2 - 0.1 / rho / 3. * pi * pi * c1 * c2 * exp(-t);
}

double solution_rho(double x, double y, double t)
{
    (void) x;
    (void) y;
    (void) t;
    //return exp(t + x - y);
    return exp(t + x - y);
}

double solution_v1(double x, double y, double t)
{
    (void) x;
    (void) y;
    (void) t;
    return s1 * s2 * exp(-t);
}

double solution_v2(double x, double y, double t)
{
    (void) x;
    (void) y;
    (void) t;
    return s1 * s2 * exp(t);
}

double rho_0(double x, double y)
{
    (void) x;
    (void) y;
    return exp(x - y);
}

double v1_0(double x, double y)
{
    (void) x;
    (void) y;
    return s1 * s2;
}

double v2_0(double x, double y)
{
    (void) x;
    (void) y;
    return s1 * s2;
}

