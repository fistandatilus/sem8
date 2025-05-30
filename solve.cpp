#include <cmath>
#include "eigen/Eigen/Eigen"

#include "solve.h"
#include "data.h"
#include "mesh.h"
#include "functions.h"

void fill_matrix_G(const data &arrays, const mesh &msh, const scheme_params &scheme, Eigen::SparseMatrix<double> &matrix, Eigen::VectorXd &rhs)
{
    const double *g = arrays.g.get();
    const double *v1 = arrays.v1.get();
    const double *v2 = arrays.v2.get();

    const unsigned int size = arrays.size;
    const double tau = scheme.tau;
    const double h = scheme.h;
    for (unsigned int i_matrix; i_matrix < size; i_matrix++)
    {
        unsigned int i, j;
        msh.get_mesh_by_matrix(i_matrix, i, j);
        point_type type = msh.get_type_by_mesh(i, j);
        const double th = tau / h;
        switch (type)
        {
        case point_type::i:
        case point_type::ilt:
        case point_type::irb:
            {
                const unsigned int i0 = i_matrix;
                const unsigned int iL = msh.get_matrix_by_mesh(i - 1, j);
                const unsigned int iR = msh.get_matrix_by_mesh(i + 1, j);
                const unsigned int iB = msh.get_matrix_by_mesh(i, j - 1);
                const unsigned int iT = msh.get_matrix_by_mesh(i, j + 1);
                matrix.insert(i0, i0) = 1. + th * (fabs(v1[i0]) + fabs(v2[i0]));
                matrix.insert(i0, iL) = th * (v1[i0] > 0 ? -v1[i0] : 0);
                matrix.insert(i0, iR) = th * (v1[i0] < 0 ? v1[i0] : 0);
                matrix.insert(i0, iB) = th * (v2[i0] > 0 ? -v2[i0] : 0);
                matrix.insert(i0, iT) = th * (v2[i0] < 0 ? v2[i0] : 0);
                rhs[i0] = g[i0] - 0.5 * th * (v1[iR] - v1[iL] + v2[iT] - v2[iB]);
                break;
             }
        case point_type::l:
        case point_type::lt:
        case point_type::lb:
            {
                const unsigned int i0 = i_matrix;
                const unsigned int iR = msh.get_matrix_by_mesh(i + 1, j);
                matrix.insert(i0, i0) = 1.;
                rhs[i0] = g[i0] - th * (v1[iR] - v1[i0]);
                break;
            }
        case point_type::r:
        case point_type::rt:
        case point_type::rb:
            {
                const unsigned int i0 = i_matrix;
                const unsigned int iL = msh.get_matrix_by_mesh(i - 1, j);
                matrix.insert(i0, i0) = 1.;
                rhs[i0] = g[i0] - th * (v1[i0] - v1[iL]);
                break;
            }
        case point_type::b:
            {
                const unsigned int i0 = i_matrix;
                const unsigned int iT = msh.get_matrix_by_mesh(i, j + 1);
                matrix.insert(i0, i0) = 1.;
                rhs[i0] = g[i0] - th * (v2[iT] - v2[i0]);
                break;
            }
        case point_type::t:
            {
                const unsigned int i0 = i_matrix;
                const unsigned int iB = msh.get_matrix_by_mesh(i, j - 1);
                matrix.insert(i0, i0) = 1.;
                rhs[i0] = g[i0] - th * (v2[i0] - v2[iB]);
                break;
            }
        }
    }
}

void fill_matrix_V1(const data &arrays, const mesh &msh, const scheme_params &scheme, const problem_params &problem, Eigen::SparseMatrix<double> &matrix, Eigen::VectorXd &rhs)
{
    const double *H = arrays.h.get();
    const double *v1 = arrays.v1.get();
    const double *v2 = arrays.v2.get();

    auto P = [&] (double x)
    {
        switch (problem.type)
        {
        case rho_type::lin1:
            return x;
        case rho_type::lin10:
            return 10 * x;
        case rho_type::lin100:
            return 100 * x;
        case rho_type::gamma:
            return pow(x, 1.4);
        }
    };

    const unsigned int size = arrays.size;
    const double tau = scheme.tau;
    const double h = scheme.h;
    const double mu = problem.mu;
    for (unsigned int i_matrix; i_matrix < size; i_matrix++)
    {
        unsigned int i, j;
        msh.get_mesh_by_matrix(i_matrix, i, j);
        point_type type = msh.get_type_by_mesh(i, j);
        const double th = tau / h;
        const double muh = mu / h;
        if (type == point_type::i)
        {
            const unsigned int i0 = i_matrix;
            const unsigned int iL = msh.get_matrix_by_mesh(i - 1, j);
            const unsigned int iR = msh.get_matrix_by_mesh(i + 1, j);
            const unsigned int iB = msh.get_matrix_by_mesh(i, j - 1);
            const unsigned int iT = msh.get_matrix_by_mesh(i, j + 1);
            const unsigned int iLB = msh.get_matrix_by_mesh(i - 1, j - 1);
            const unsigned int iLT = msh.get_matrix_by_mesh(i - 1, j + 1);
            const unsigned int iRB = msh.get_matrix_by_mesh(i + 1, j - 1);
            const unsigned int iRT = msh.get_matrix_by_mesh(i + 1, j + 1);
            const double hv1 = v1[i0] * H[i0];
            const double hv2 = v2[i0] * H[i0];
            matrix.insert(i0, i0) = H[i0] + th * (fabs(hv1) + fabs(hv2) + muh * 14. / 3.);
            matrix.insert(i0, iL) = th * ((hv1 > 0 ? -hv1 : 0) - muh * 4. / 3.);
            matrix.insert(i0, iR) = th * ((hv1 < 0 ? hv1 : 0) - muh * 4. / 3.);
            matrix.insert(i0, iB) = th * ((hv2 > 0 ? -hv2 : 0) - muh);
            matrix.insert(i0, iT) = th * ((hv2 < 0 ? hv2 : 0) - muh);
            rhs[i0] = hv1 + th * ( -0.5 * (P(H[iR]) - P(H[iL])) + 1. / 12. * muh * (v2[iRT] - v2[iRB] - v2[iLT] + v2[iLB]) + H[i0] * f_1(i * h, j * h));
        }
        else
        {
            const unsigned int i0 = i_matrix;
            matrix.insert(i0, i0) = 1.;
            rhs[i0] = 0.;
        }
    }
}

void fill_matrix_V2(const data &arrays, const mesh &msh, const scheme_params &scheme, const problem_params &problem, Eigen::SparseMatrix<double> &matrix, Eigen::VectorXd &rhs)
{
    const double *H = arrays.h.get();
    const double *v1 = arrays.v1.get();
    const double *v2 = arrays.v2.get();

    auto P = [&] (double x)
    {
        switch (problem.type)
        {
        case rho_type::lin1:
            return x;
        case rho_type::lin10:
            return 10 * x;
        case rho_type::lin100:
            return 100 * x;
        case rho_type::gamma:
            return pow(x, 1.4);
        }
    };

    const unsigned int size = arrays.size;
    const double tau = scheme.tau;
    const double h = scheme.h;
    const double mu = problem.mu;
    for (unsigned int i_matrix; i_matrix < size; i_matrix++)
    {
        unsigned int i, j;
        msh.get_mesh_by_matrix(i_matrix, i, j);
        point_type type = msh.get_type_by_mesh(i, j);
        const double th = tau / h;
        const double muh = mu / h;
        if (type == point_type::i)
        {
            const unsigned int i0 = i_matrix;
            const unsigned int iL = msh.get_matrix_by_mesh(i - 1, j);
            const unsigned int iR = msh.get_matrix_by_mesh(i + 1, j);
            const unsigned int iB = msh.get_matrix_by_mesh(i, j - 1);
            const unsigned int iT = msh.get_matrix_by_mesh(i, j + 1);
            const unsigned int iLB = msh.get_matrix_by_mesh(i - 1, j - 1);
            const unsigned int iLT = msh.get_matrix_by_mesh(i - 1, j + 1);
            const unsigned int iRB = msh.get_matrix_by_mesh(i + 1, j - 1);
            const unsigned int iRT = msh.get_matrix_by_mesh(i + 1, j + 1);
            const double hv1 = v1[i0] * H[i0];
            const double hv2 = v2[i0] * H[i0];
            matrix.insert(i0, i0) = H[i0] + th * (fabs(hv1) + fabs(hv2) + muh * 14. / 3.);
            matrix.insert(i0, iL) = th * ((hv1 > 0 ? -hv1 : 0) - muh);
            matrix.insert(i0, iR) = th * ((hv1 < 0 ? hv1 : 0) - muh);
            matrix.insert(i0, iB) = th * ((hv2 > 0 ? -hv2 : 0) - muh * 4. / 3.);
            matrix.insert(i0, iT) = th * ((hv2 < 0 ? hv2 : 0) - muh * 4. / 3.);
            rhs[i0] = hv2 + th * ( -0.5 * (P(H[iT]) - P(H[iB])) + 1. / 12. * muh * (v1[iRT] - v1[iRB] - v1[iLT] + v1[iLB]) + H[i0] * f_2(i * h, j * h));
        }
        else
        {
            const unsigned int i0 = i_matrix;
            matrix.insert(i0, i0) = 1.;
            rhs[i0] = 0.;
        }
    }
}

void fill_G_H_from_matrix(double *g, double *h, const Eigen::VectorXd &x, unsigned int size)
{
    for (unsigned int i = 0; i < size; i++)
    {
        g[i] = x[i];
        h[i] = exp(x[i]);
    }
}

void fill_V_from_matrix(double *v, const Eigen::VectorXd &x, unsigned int size)
{
    for (unsigned int i = 0; i < size; i++)
        v[i] = x[i];
}

void general_loop(problem_params &problem, scheme_params &scheme, data &arrays)
{
    const unsigned int N = scheme.N;
    const unsigned int M = scheme.M;

    const mesh msh(M);
    const unsigned int size = msh.size;
    Eigen::VectorXd x(size), x2(size), b(size);
    Eigen::SparseMatrix<double> A(size, size);
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;

    for (unsigned int n = 0; n < N; n++)
    {
        fill_matrix_G(arrays, msh, scheme, A, b);
        solver.compute(A);
        x = solver.solve(b);
        fill_G_H_from_matrix(arrays.g.get(), arrays.h.get(), x, size);
        fill_matrix_V1(arrays, msh, scheme, problem, A, b);
        solver.compute(A);
        x = solver.solve(b);
        fill_matrix_V2(arrays, msh, scheme, problem, A, b);
        solver.compute(A);
        x2 = solver.solve(b);
        fill_V_from_matrix(arrays.v1.get(), x, size);
        fill_V_from_matrix(arrays.v2.get(), x2, size);
    }
}

