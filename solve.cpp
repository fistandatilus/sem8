#include <cmath>

#include "norms.h"
#include "solve.h"
#include "data.h"
#include "mesh.h"
#include "functions.h"

void fill_matrix_G(const data &arrays, const mesh &msh, const scheme_params &scheme, const unsigned int n, Eigen::SparseMatrix<double> &matrix, Eigen::VectorXd &rhs)
{
    const double *g = arrays.g.get();
    const double *v1 = arrays.v1.get();
    const double *v2 = arrays.v2.get();

    const unsigned int size = arrays.size;
    const double tau = scheme.tau;
    const double h = scheme.h;
    for (unsigned int i_matrix = 0; i_matrix < size; i_matrix++)
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
                matrix.coeffRef(i0, i0) = 1. + th * (fabs(v1[i0]) + fabs(v2[i0]));
                matrix.coeffRef(i0, iL) = 0.5 * th * (-v1[i0] - fabs(v1[i0]));
                matrix.coeffRef(i0, iR) = 0.5 * th * (v1[i0] - fabs(v1[i0]));
                matrix.coeffRef(i0, iB) = 0.5 * th * (-v2[i0] - fabs(v2[i0]));
                matrix.coeffRef(i0, iT) = 0.5 * th * (v2[i0] - fabs(v2[i0]));
                rhs[i0] = g[i0] - 0.5 * th * (v1[iR] - v1[iL] + v2[iT] - v2[iB]) + tau * f_0(i * h, j * h, tau * n);
                break;
             }
        case point_type::l:
        case point_type::lt:
        case point_type::lb:
            {
                const unsigned int i0 = i_matrix;
                const unsigned int iR = msh.get_matrix_by_mesh(i + 1, j);
                matrix.coeffRef(i0, i0) = 1.;
                rhs[i0] = g[i0] - th * (v1[iR] - v1[i0]) + tau * f_0(i * h, j * h, tau * n);
                break;
            }
        case point_type::r:
        case point_type::rt:
        case point_type::rb:
            {
                const unsigned int i0 = i_matrix;
                const unsigned int iL = msh.get_matrix_by_mesh(i - 1, j);
                matrix.coeffRef(i0, i0) = 1.;
                rhs[i0] = g[i0] - th * (v1[i0] - v1[iL]) + tau * f_0(i * h, j * h, tau * n);
                break;
            }
        case point_type::b:
            {
                const unsigned int i0 = i_matrix;
                const unsigned int iT = msh.get_matrix_by_mesh(i, j + 1);
                matrix.coeffRef(i0, i0) = 1.;
                rhs[i0] = g[i0] - th * (v2[iT] - v2[i0]) + tau * f_0(i * h, j * h, tau * n);
                break;
            }
        case point_type::t:
            {
                const unsigned int i0 = i_matrix;
                const unsigned int iB = msh.get_matrix_by_mesh(i, j - 1);
                matrix.coeffRef(i0, i0) = 1.;
                rhs[i0] = g[i0] - th * (v2[i0] - v2[iB]) + tau * f_0(i * h, j * h, tau * n);
                break;
            }
        }
    }
}

void fill_matrix_V1(const data &arrays, const mesh &msh, const scheme_params &scheme, const problem_params &problem, const unsigned int n, Eigen::SparseMatrix<double> &matrix, Eigen::VectorXd &rhs)
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
        return 0.;
    };

    const unsigned int size = arrays.size;
    const double tau = scheme.tau;
    const double h = scheme.h;
    const double mu = problem.mu;
    for (unsigned int i_matrix = 0; i_matrix < size; i_matrix++)
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
            matrix.coeffRef(i0, i0) = H[i0] + th * (fabs(hv1) + fabs(hv2) + muh * 14. / 3.);
            matrix.coeffRef(i0, iL) = th * ((hv1 > 0 ? -hv1 : 0) - muh * 4. / 3.);
            matrix.coeffRef(i0, iR) = th * ((hv1 < 0 ? hv1 : 0) - muh * 4. / 3.);
            matrix.coeffRef(i0, iB) = th * ((hv2 > 0 ? -hv2 : 0) - muh);
            matrix.coeffRef(i0, iT) = th * ((hv2 < 0 ? hv2 : 0) - muh);
            rhs[i0] = hv1 + th * ( -0.5 * (P(H[iR]) - P(H[iL])) + 1. / 12. * muh * (v2[iRT] - v2[iRB] - v2[iLT] + v2[iLB])) + tau * H[i0] * f_1(i * h, j * h, n * tau);
        }
        else
        {
            const unsigned int i0 = i_matrix;
            matrix.coeffRef(i0, i0) = 1.;
            rhs[i0] = 0.;
        }
    }
}

void fill_matrix_V2(const data &arrays, const mesh &msh, const scheme_params &scheme, const problem_params &problem, const unsigned int n, Eigen::SparseMatrix<double> &matrix, Eigen::VectorXd &rhs)
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
        return 0.;
    };

    const unsigned int size = arrays.size;
    const double tau = scheme.tau;
    const double h = scheme.h;
    const double mu = problem.mu;
    for (unsigned int i_matrix = 0; i_matrix < size; i_matrix++)
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
            matrix.coeffRef(i0, i0) = H[i0] + th * (fabs(hv1) + fabs(hv2) + muh * 14. / 3.);
            matrix.coeffRef(i0, iL) = th * ((hv1 > 0 ? -hv1 : 0) - muh);
            matrix.coeffRef(i0, iR) = th * ((hv1 < 0 ? hv1 : 0) - muh);
            matrix.coeffRef(i0, iB) = th * ((hv2 > 0 ? -hv2 : 0) - muh * 4. / 3.);
            matrix.coeffRef(i0, iT) = th * ((hv2 < 0 ? hv2 : 0) - muh * 4. / 3.);
            rhs[i0] = hv2 + th * ( -0.5 * (P(H[iT]) - P(H[iB])) + 1. / 12. * muh * (v1[iRT] - v1[iRB] - v1[iLT] + v1[iLB])) + tau * H[i0] * f_2(i * h, j * h, n * tau);
        }
        else
        {
            const unsigned int i0 = i_matrix;
            matrix.coeffRef(i0, i0) = 1.;
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
    solver.setMaxIterations(2000);
    solver.setTolerance(1e-8);
//    data arrays_solution(msh.size);
//    fill_solution(arrays_solution, scheme);
//
//    check_norms(arrays, arrays_solution, msh, scheme);

    for (unsigned int n = 0; n * scheme.tau <= scheme.T; n++)
    {
        fill_matrix_G(arrays, msh, scheme, n, A, b);
        solver.compute(A);
        x = solver.solve(b);
        fill_G_H_from_matrix(arrays.g.get(), arrays.h.get(), x, size);
        fill_matrix_V1(arrays, msh, scheme, problem, n, A, b);
        solver.compute(A);
        x = solver.solve(b);
        fill_matrix_V2(arrays, msh, scheme, problem, n, A, b);
        solver.compute(A);
        x2 = solver.solve(b);
        fill_V_from_matrix(arrays.v1.get(), x, size);
        fill_V_from_matrix(arrays.v2.get(), x2, size);
//        data arrays_solution(msh.size);
//        fill_solution(arrays_solution, scheme);
//
//        check_norms(arrays, arrays_solution, msh, scheme);
    }
}

void set_initial_data(data &arrays, scheme_params &scheme)
{
    double *h = arrays.h.get();
    double *g = arrays.g.get();
    double *v1 = arrays.v1.get();
    double *v2 = arrays.v2.get();
    const unsigned int size = arrays.size;
    const double h_x = scheme.h;
    const mesh msh(scheme.M);

    for (unsigned int i_array = 0; i_array < size; i_array++)
    {
        unsigned int i;
        unsigned int j;
        msh.get_mesh_by_matrix(i_array, i, j);
        double x1 = i * h_x;
        double x2 = j * h_x;

        h[i_array] = rho_0(x1, x2);
        g[i_array] = log(h[i_array]);

        v1[i_array] = v1_0(x1, x2);
        v2[i_array] = v2_0(x1, x2);
    }
}

void fill_solution(data &arrays, scheme_params &scheme)
{
    double *h = arrays.h.get();
    double *g = arrays.g.get();
    double *v1 = arrays.v1.get();
    double *v2 = arrays.v2.get();
    const unsigned int size = arrays.size;
    const double h_x = scheme.h;
    const double T = scheme.T;
    const mesh msh(scheme.M);

    for (unsigned int i_array = 0; i_array < size; i_array++)
    {
        unsigned int i;
        unsigned int j;
        msh.get_mesh_by_matrix(i_array, i, j);
        double x1 = i * h_x;
        double x2 = j * h_x;

        h[i_array] = solution_rho(x1, x2, T);
        g[i_array] = log(h[i_array]);

        v1[i_array] = solution_v1(x1, x2, T);
        v2[i_array] = solution_v2(x1, x2, T);
    }
}

void check_norms(const data &arrays1, const data &arrays2, const mesh &msh, const scheme_params &scheme)
{
    double c = -1, l = -1, w = -1;
    c = C_norm(msh, msh, arrays1.g.get(), arrays2.g.get());
    l = L_norm(msh, msh, arrays1.g.get(), arrays2.g.get(), scheme);
    printf("  g: c = %e, l = %e, w = %e\n", c, l, w);
    c = C_norm(msh, msh, arrays1.h.get(), arrays2.h.get());
    l = L_norm(msh, msh, arrays1.h.get(), arrays2.h.get(), scheme);
    printf("rho: c = %e, l = %e, w = %e\n", c, l, w);
    c = C_norm(msh, msh, arrays1.v1.get(), arrays2.v1.get());
    l = L_norm(msh, msh, arrays1.v1.get(), arrays2.v1.get(), scheme);
    printf(" v1: c = %e, l = %e, w = %e\n", c, l, w);
    c = C_norm(msh, msh, arrays1.v2.get(), arrays2.v2.get());
    l = L_norm(msh, msh, arrays1.v2.get(), arrays2.v2.get(), scheme);
    printf(" v2: c = %e, l = %e, w = %e\n", c, l, w);
}
