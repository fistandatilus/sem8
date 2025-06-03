#include <cmath>

#include "functions.h"
#include "data.h"

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

void set_p_coeff(rho_type type)
{
    switch (type)
    {
    case rho_type::lin1:
        p_coeff = 1;
        return;
    case rho_type::lin10:
        p_coeff = 10;
        return;
    case rho_type::lin100:
        p_coeff = 100;
        return;
    case rho_type::gamma:
        p_coeff = 1.4;
    }
}

double f_0(double x, double y, double t)
{
    (void) x;
    (void) y;
    (void) t;
    return 1 + v1 - v2 + pi * (c1 * s2 * exp(-t) + s1 * c2 * exp(t));
    //return 0;
}

double f_1(double x, double y, double t)
{
    (void) x;
    (void) y;
    (void) t;
    return -v1 + p_coeff * (p_coeff > 1.1 && p_coeff < 2 ? exp(0.4 * (t + x - y)) : 1)+ pi * v1 * c1 * s2 * exp(-t) + v2 * pi * s1 * c2 * exp(-t) + mu / rho * 7. / 3. * pi * pi * v1 - mu / rho / 3. * pi * pi * c1 * c2 * exp(t);
    //return 1 / (x+ 1);
}

double f_2(double x, double y, double t)
{
    (void) x;
    (void) y;
    (void) t;
    return v2 + -p_coeff * (p_coeff > 1.1 && p_coeff < 2 ? exp(0.4 * (t + x - y)) : 1) + v1 * pi * c1 * s2 * exp(t) + v2 * pi * s1 * c2 * exp(t) + mu / rho * 7. / 3. * pi * pi * v2 - mu / rho / 3. * pi * pi * c1 * c2 * exp(-t);
    //return 0;
}

double solution_rho(double x, double y, double t)
{
    (void) x;
    (void) y;
    (void) t;
    return exp(t + x - y);
    //return x + 1;
}

double solution_v1(double x, double y, double t)
{
    (void) x;
    (void) y;
    (void) t;
    return s1 * s2 *exp(-t);
    //return 0;
}

double solution_v2(double x, double y, double t)
{
    (void) x;
    (void) y;
    (void) t;
    return s1 * s2 *exp(t);
    //return 0;
}

double rho_0(double x, double y)
{
    (void) x;
    (void) y;
    return exp(x - y);
    //return x + 1;
}

double v1_0(double x, double y)
{
    (void) x;
    (void) y;
    return s1 * s2;
    //return 0;
}

double v2_0(double x, double y)
{
    (void) x;
    (void) y;
    return s1 * s2;
    //return 0;
}

