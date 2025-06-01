#ifndef NORMS_H
#define NORMS_H

class mesh;
class scheme_params;

double C_norm(const mesh &msh1, const mesh &msh2, const double *array1, const double *array2, const unsigned int scale = 1);
double L_norm(const mesh &msh1, const mesh &msh2, const double *array1, const double *array2, const scheme_params &params,  const unsigned int scale = 1);
double W_norm(const mesh &msh1, const mesh &msh2, const double *array1, const double *array2, const scheme_params &params,  const unsigned int scale = 1);

#endif //NORMS_H
