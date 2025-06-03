#ifndef SOLVE_H
#define SOLVE_H
#include <memory>

#include "eigen/Eigen/Eigen"

class data;
class mesh;
class problem_params;
class scheme_params;

template <class Matrix, class Vector>
void fill_matrix_G(const data &arrays, const mesh &msh, const scheme_params &scheme, const unsigned int n, Matrix &matrix, Vector &rhs, const unsigned int start, const unsigned int stride);

template <class Matrix, class Vector>
void fill_matrix_V1(const data &arrays, const mesh &msh, const scheme_params &scheme, const problem_params &problem, const unsigned int n, Matrix &matrix, Vector &rhs, const unsigned int start, const unsigned int stride);

template <class Matrix, class Vector>
void fill_matrix_V2(const data &arrays, const mesh &msh, const scheme_params &scheme, const problem_params &problem, const unsigned int n, Matrix &matrix, Vector &rhs, const unsigned int start, const unsigned int stride);

template <class Vector>
void fill_G_H_from_matrix(double *g, double *h, const Vector &x, const unsigned int start, const unsigned int stride);

template <class Vector>
void fill_H_from_matrix(double *h, const Vector &x, const unsigned int start, const unsigned int stride);

template <class Vector>
void fill_V_from_matrix(double *v, const Vector &x, const unsigned int start, const unsigned int stride);
void general_loop(const problem_params &problem, const scheme_params &scheme, data &arrays);
void general_loop_my_solve(const problem_params &problem, const scheme_params &scheme, data &arrays, const unsigned int n_threads);
void *inner_loop_my_solve(void *args_ptr);

void set_initial_data(data &arrays, scheme_params &scheme);
void fill_solution(data &arrays, scheme_params &scheme);
void check_norms(const data &arrays1, const data &arrays2, const mesh &msh, const scheme_params &scheme);

#endif //SOLVE_H
