#ifndef FUNCTIONS_H
#define FUNCTIONS_H

enum class rho_type;

void set_mu(double new_mu);
void set_p_coeff(rho_type new_coeff);

double f_0(double x, double y, double t);
double f_1(double x, double y, double t);
double f_2(double x, double y, double t);

double solution_rho(double x, double y, double t);
double solution_v1(double x, double y, double t);
double solution_v2(double x, double y, double t);

double rho_0(double x, double y);
double v1_0(double x, double y);
double v2_0(double x, double y);

#endif //FUNCTIONS_H
