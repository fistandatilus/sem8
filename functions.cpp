#include <cmath>

#include "functions.h"

double f_0(double x, double y, double t)
{
    (void) x;
    (void) y;
    (void) t;
    return t;
}

double f_1(double x, double y, double t)
{
    (void) x;
    (void) y;
    (void) t;
    return 0;
}

double f_2(double x, double y, double t)
{
    (void) x;
    (void) y;
    (void) t;
    return 0;
}

double solution_rho(double x, double y, double t)
{
    (void) x;
    (void) y;
    (void) t;
    return exp(t);
}

double solution_v1(double x, double y, double t)
{
    (void) x;
    (void) y;
    (void) t;
    return 0;
}

double solution_v2(double x, double y, double t)
{
    (void) x;
    (void) y;
    (void) t;
    return 0;
}

double rho_0(double x, double y)
{
    (void) x;
    (void) y;
    return 1;
}

double v1_0(double x, double y)
{
    (void) x;
    (void) y;
    return 0;
}

double v2_0(double x, double y)
{
    (void) x;
    (void) y;
    return 0;
}

