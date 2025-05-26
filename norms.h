#ifndef NORMS_H
#define NORMS_H

#include "data.h"
#include "mesh.h"

double C_norm(const mesh &msh1, const mesh &msh2, const double *array1, const double *array2, const unsigned int scale);
double L_norm(const mesh &msh1, const mesh &msh2, const double *array1, const double *array2, const scheme_params &params,  const unsigned int scale);
double W_norm(const mesh &msh1, const mesh &msh2, const double *array1, const double *array2, const scheme_params &params,  const unsigned int scale);

#endif //NORMS_H
