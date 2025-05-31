#ifndef SOLVE_H
#define SOLVE_H
#include "eigen/Eigen/Eigen"

class data;
class mesh;
class problem_params;
class scheme_params;

void fill_matrix_G(const data &arrays, const mesh &msh, const scheme_params &scheme, Eigen::SparseMatrix<double> &matrix, Eigen::VectorXd &rhs);
void fill_matrix_V1(const data &arrays, const mesh &msh, const scheme_params &scheme, const problem_params &problem, Eigen::SparseMatrix<double> &matrix, Eigen::VectorXd &rhs);
void fill_matrix_V2(const data &arrays, const mesh &msh, const scheme_params &scheme, const problem_params &problem, Eigen::SparseMatrix<double> &matrix, Eigen::VectorXd &rhs);
void fill_G_H_from_matrix(double *g, double *h, const Eigen::VectorXd &x, unsigned int size);
void fill_V_from_matrix(double *v, const Eigen::VectorXd &x, unsigned int size);
void general_loop(problem_params &problem, scheme_params &scheme, data &arrays);

#endif //SOLVE_H
