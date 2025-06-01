#include <cmath>

#include "data.h"
#include "norms.h"
#include "mesh.h"

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
        double diff = fabs(array1[i_array] - array2[i_array2]);
        norm = diff > norm ? diff : norm;
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

double W_norm(const mesh &msh1, const mesh &msh2, const double *array1, const double *array2, const scheme_params &params,  const unsigned int scale)
{
    double norm = 0;
    const unsigned int size = msh1.size;
    for (unsigned int i_array = 0; i_array < size; i_array++)
    {
        unsigned int i, j, i_array2 = i_array;
        msh1.get_mesh_by_matrix(i_array, i, j);
        point_type type = msh1.get_type_by_mesh(i, j);

        if (scale > 1)
            i_array2 = msh2.get_matrix_by_mesh(i * scale, j * scale);
        double val0 = array1[i_array] - array2[i_array2];
        val0 *= val0;
        if (type != point_type::i)
            val0 *= 0.5;
        norm += val0;

        if (type == point_type::r || type == point_type::t || type == point_type::rt)
            continue;
        unsigned int i_arrayR = msh1.get_matrix_by_mesh(i + 1, j);
        unsigned int i_arrayT = msh1.get_matrix_by_mesh(i, j + 1);
        unsigned int i_array2R = i_arrayR;
        unsigned int i_array2T = i_arrayT;
        if (scale > 1)
        {
            i_array2R = msh2.get_matrix_by_mesh(i * scale + 1, j * scale);
            i_array2T = msh2.get_matrix_by_mesh(i * scale, j * scale + 1);
        }
        double valR = 0;
        double valT = 0;
        if (type != point_type::b)
            valR = array1[i_array] - array2[i_array2] - array1[i_arrayR] + array2[i_array2R]; 
        if (type != point_type::l)
            valT = array1[i_array] - array2[i_array2] - array1[i_arrayT] + array2[i_array2T]; 
        norm += valR * valR + valT * valT;
    }
    return params.h * sqrt(norm);
}
