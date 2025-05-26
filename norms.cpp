#include <cmath>

#include "data.h"
#include "norms.h"

double C_norm(const mesh &msh1, const mesh &msh2, const double *array1, const double *array2, const unsigned int scale)
{
    double norm = 0;
    const unsigned int size = msh1.size;
    for (unsigned int i_array = 0; i_array < size; i_array++)
    {
        unsigned int i, j, i_array2 = i_array;
        if (scale > 1)
        {
            msh1.get_mesh_by_matrix(i_array, i, j);
            i_array2 = msh2.get_matrix_by_mesh(i * scale, j * scale);
        }
        norm += fabs(array1[i_array] - array2[i_array2]);
    }
    return norm;
}

double L_norm(const mesh &msh1, const mesh &msh2, const double *array1, const double *array2, const scheme_params &params,  const unsigned int scale)
{
    double norm = 0;
    const unsigned int size = msh1.size;
    for (unsigned int i_array = 0; i_array < size; i_array++)
    {
        unsigned int i, j, i_array2 = i_array;
        point_type type = msh1.get_type_by_matrix(i_array);
        if (scale > 1)
        {
            msh1.get_mesh_by_matrix(i_array, i, j);
            i_array2 = msh2.get_matrix_by_mesh(i * scale, j * scale);
        }
        double val = array1[i_array] - array2[i_array2];
        if (type == point_type::i)
            norm += val * val;
        else
            norm += 0.5 * val * val;
    }
    return params.h * sqrt(norm);
}

double W_norm(const mesh &msh1, const mesh &msh2, const double *array1, const double *array2, const unsigned int scale)
{
    double norm = 0;
    const unsigned int size = msh1.size;
    for (unsigned int i_array = 0; i_array < size; i_array++)
    {
        unsigned int i, j, i_array2 = i_array;
        msh1.get_mesh_by_matrix(i_array, i, j);
        point_type type = msh1.get_type_by_mesh(i, j);
        if (type == point_type::l || type == point_type::b || type == point_type::lt || type == point_type::lb|| type == point_type::rb)
            continue;
        unsigned int i_arrayL = msh1.get_matrix_by_mesh(i - 1, j);
        unsigned int i_arrayB = msh1.get_matrix_by_mesh(i, j - 1);
        unsigned int i_array2L = i_arrayL;
        unsigned int i_array2B = i_arrayB;
        if (scale > 1)
        {
            i_array2 = msh2.get_matrix_by_mesh(i * scale, j * scale);
            i_array2L = msh2.get_matrix_by_mesh(i * scale - 1, j * scale);
            i_array2B = msh2.get_matrix_by_mesh(i * scale, j * scale - 1);
        }
        double valL = array1[i_array] - array2[i_array2] - array1[i_arrayL] + array2[i_array2L]; 
        double valB = array1[i_array] - array2[i_array2] - array1[i_arrayB] + array2[i_array2B]; 
        norm += valL * valL + valB * valB;
    }
    return sqrt(norm);
}
