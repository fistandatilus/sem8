#include <sys/sysinfo.h>
#include <cmath>

#include "norms.h"
#include "solve.h"
#include "data.h"
#include "mesh.h"
#include "msr.h"
#include "functions.h"

template <class Matrix, class Vector>
void fill_matrix_G(const data &arrays, const mesh &msh, const scheme_params &scheme, const unsigned int n, Matrix &matrix, Vector &rhs, const unsigned int start, const unsigned int stride)
{
    const double *g = arrays.g.get();
    const double *v1 = arrays.v1.get();
    const double *v2 = arrays.v2.get();

    const double tau = scheme.tau;
    const double h = scheme.h;
    for (unsigned int i_matrix = start; i_matrix < start + stride; i_matrix++)
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

template <class Matrix, class Vector>
void fill_matrix_V1(const data &arrays, const mesh &msh, const scheme_params &scheme, const problem_params &problem, const unsigned int n, Matrix &matrix, Vector &rhs, const unsigned int start, const unsigned int stride)
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

    const double tau = scheme.tau;
    const double h = scheme.h;
    const double mu = problem.mu;
    for (unsigned int i_matrix = start; i_matrix < start + stride; i_matrix++)
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
            matrix.coeffRef(i0, iL) = th * (-0.5 * (hv1  + fabs(hv1)) - muh * 4. / 3.);
            matrix.coeffRef(i0, iR) = th * (0.5 * (hv1 - fabs(hv1)) - muh * 4. / 3.);
            matrix.coeffRef(i0, iB) = th * (-0.5 * (hv2 + fabs(hv2)) - muh);
            matrix.coeffRef(i0, iT) = th * (0.5 * (hv2 - fabs(hv2)) - muh);
            rhs[i0] = hv1 + th * ( -0.5 * (P(H[iR]) - P(H[iL])) + 1. / 12. * muh * (v2[iRT] - v2[iRB] - v2[iLT] + v2[iLB])) + tau * H[i0] * f_1(i * h, j * h, n * tau);
            //printf("i0 = %u, iR = %u, iL = %u, (x, y) = (%f, %f), H[i0] = %e, H[iL] = %e, H[iR] = %e, rhs = %e\n", i0, iR, iL, i * h, j * h, H[i0], H[iL], H[iR], rhs[i0]);
        }
        else
        {
            const unsigned int i0 = i_matrix;
            matrix.coeffRef(i0, i0) = 1.;
            rhs[i0] = 0.;

            if (type == point_type::irb || type == point_type::ilt)
            {
                const unsigned int iL = msh.get_matrix_by_mesh(i - 1, j);
                const unsigned int iR = msh.get_matrix_by_mesh(i + 1, j);
                const unsigned int iB = msh.get_matrix_by_mesh(i, j - 1);
                const unsigned int iT = msh.get_matrix_by_mesh(i, j + 1);
                matrix.coeffRef(i0, iL) = 0.;
                matrix.coeffRef(i0, iR) = 0.;
                matrix.coeffRef(i0, iB) = 0.;
                matrix.coeffRef(i0, iT) = 0.;
            }
        }
    }
}

template <class Matrix, class Vector>
void fill_matrix_V2(const data &arrays, const mesh &msh, const scheme_params &scheme, const problem_params &problem, const unsigned int n, Matrix &matrix, Vector &rhs, const unsigned int start, const unsigned int stride)
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
    for (unsigned int i_matrix = start; i_matrix < start + stride; i_matrix++)
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
            matrix.coeffRef(i0, iL) = th * (-0.5 * (hv1 + fabs(hv1)) - muh);
            matrix.coeffRef(i0, iR) = th * (0.5 * (hv1 - fabs(hv1)) - muh);
            matrix.coeffRef(i0, iB) = th * (-0.5 * (hv2 + fabs(hv2)) - muh * 4. / 3.);
            matrix.coeffRef(i0, iT) = th * (0.5 * (hv2 - fabs(hv2)) - muh * 4. / 3.);
            rhs[i0] = hv2 + th * ( -0.5 * (P(H[iT]) - P(H[iB])) + 1. / 12. * muh * (v1[iRT] - v1[iRB] - v1[iLT] + v1[iLB])) + tau * H[i0] * f_2(i * h, j * h, n * tau);
        }
        else
        {
            const unsigned int i0 = i_matrix;
            matrix.coeffRef(i0, i0) = 1.;
            rhs[i0] = 0.;

            if (type == point_type::irb || type == point_type::ilt)
            {
                const unsigned int iL = msh.get_matrix_by_mesh(i - 1, j);
                const unsigned int iR = msh.get_matrix_by_mesh(i + 1, j);
                const unsigned int iB = msh.get_matrix_by_mesh(i, j - 1);
                const unsigned int iT = msh.get_matrix_by_mesh(i, j + 1);
                matrix.coeffRef(i0, iL) = 0.;
                matrix.coeffRef(i0, iR) = 0.;
                matrix.coeffRef(i0, iB) = 0.;
                matrix.coeffRef(i0, iT) = 0.;
            }
        }
    }
}

template <class Vector>
void fill_G_H_from_matrix(double *g, double *h, const Vector &x, const unsigned int start, const unsigned int stride)
{
    for (unsigned int i = start; i < start + stride; i++)
    {
        g[i] = x[i];
        h[i] = exp(x[i]);
    }
}

template <class Vector>
void fill_H_from_matrix(double *h, const Vector &x, const unsigned int start, const unsigned int stride)
{
    for (unsigned int i = start; i < start + stride; i++)
    {
        h[i] = exp(x[i]);
    }
}

template <class Vector>
void fill_V_from_matrix(double *v, const Vector &x, const unsigned int start, const unsigned int stride)
{
    for (unsigned int i = start; i < start + stride; i++)
        v[i] = x[i];
}

void general_loop(const problem_params &problem, const scheme_params &scheme, data &arrays)
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

    const unsigned int start = 0;
    const unsigned int stride = size;


//    data arrays_solution(msh.size);
//    fill_solution(arrays_solution, scheme);
//
//    check_norms(arrays, arrays_solution, msh, scheme);

    for (unsigned int n = 0; n * scheme.tau < scheme.T; n++)
    {
        fill_matrix_G(arrays, msh, scheme, n, A, b, start, stride);
        solver.compute(A);
        x = solver.solve(b);
        fill_G_H_from_matrix(arrays.g.get(), arrays.h.get(), x, start, stride);
        fill_matrix_V1(arrays, msh, scheme, problem, n, A, b, start, stride);
        solver.compute(A);
        x = solver.solve(b);
        fill_matrix_V2(arrays, msh, scheme, problem, n, A, b, start, stride);
        solver.compute(A);
        x2 = solver.solve(b);
        fill_V_from_matrix(arrays.v1.get(), x, start, stride);
        fill_V_from_matrix(arrays.v2.get(), x2, start, stride);
//        data arrays_solution(msh.size);
//        fill_solution(arrays_solution, scheme);
//
//        check_norms(arrays, arrays_solution, msh, scheme);
    }
}

struct loop_args
{
    const problem_params *problem{};
    const scheme_params *scheme{};
    data *arrays{};
    msr *A;
    double *b{};
    double *s{};
    double *r{};
    double *p{};
    double *v1_tmp{};
    unsigned int thr_num;
    unsigned int thread;
};

void general_loop_my_solve(const problem_params &problem, const scheme_params &scheme, data &arrays, const unsigned int n_threads)
{
    const unsigned int N = scheme.N;
    const unsigned int M = scheme.M;

    const mesh msh(M);
    const unsigned int size = msh.size;

    std::unique_ptr<unsigned int []>msr_template = make_template(msh);
    msr A;
    A.set_template(msr_template.get(), size, msr_size(M, size));

    std::unique_ptr<double []> b = std::make_unique<double []>(size);
    std::unique_ptr<double []> s = std::make_unique<double []>(size);
    std::unique_ptr<double []> r = std::make_unique<double []>(size);
    std::unique_ptr<double []> p = std::make_unique<double []>(size);
    std::unique_ptr<double []> v1_tmp = std::make_unique<double []>(size);

    std::unique_ptr<pthread_t []> tid = std::make_unique<pthread_t []>(n_threads);
    std::unique_ptr<loop_args []> args = std::make_unique<loop_args []>(n_threads);

    for (unsigned int i = 0; i < n_threads; i++)
    {
        args[i].thread = i;
        args[i].thr_num = n_threads;
        args[i].problem = &problem;
        args[i].scheme = &scheme;
        args[i].arrays = &arrays;
        args[i].A = &A;
        args[i].b = b.get();
        args[i].s = s.get();
        args[i].r = r.get();
        args[i].p = p.get();
        args[i].v1_tmp = v1_tmp.get();
        if (i)
        {
            pthread_create(tid.get() + i, 0, inner_loop_my_solve, args.get() + i);
        }
    }
    inner_loop_my_solve(args.get() + 0);

    for (unsigned int i = 1; i < n_threads; i++)
    {
        pthread_join(tid[i], 0);
    }


}

void *inner_loop_my_solve(void *args_ptr)
{
    loop_args &args = *(loop_args *)args_ptr;

    cpu_set_t cpu;
    CPU_ZERO(&cpu);
    int nproc = get_nprocs();
    int cpu_id = nproc - 1 - args.thread%nproc;
    CPU_SET(cpu_id, &cpu);
    pthread_t tid = pthread_self();
    pthread_setaffinity_np(tid, sizeof(cpu_set_t), &cpu);

    const problem_params &problem = *args.problem;
    const scheme_params &scheme = *args.scheme;
    data &arrays = *args.arrays;
    msr &A = *args.A;
    double *b = args.b;
    double *s = args.s;
    double *r = args.r;
    double *p = args.p;
    double *v1_tmp = args.v1_tmp;

    const unsigned int N = scheme.N;
    const unsigned int M = scheme.M;

    const mesh msh(M);
    const unsigned int size = msh.size;

    const unsigned int thr_num = args.thr_num;
    const unsigned int thread = args.thread;
    unsigned int start, stride;
    start_and_size(thr_num, thread, size, start, stride);

    const double desired_eps = 1e-8;
    const unsigned int max_it = 2000;
    unsigned int it;
    for (unsigned int n = 0; n * scheme.tau < scheme.T; n++)
    {
        fill_matrix_G(arrays, msh, scheme, n, A, b, start, stride);
        precond_jacobi(A, b, thr_num, thread);
        solve(A, b, arrays.g.get(), s, r, p, desired_eps, thr_num, thread, max_it, it);
        fill_H_from_matrix(arrays.h.get(), arrays.g.get(), start, stride);
        reduce_sum<int>(thr_num);
        fill_matrix_V1(arrays, msh, scheme, problem, n, A, b, start, stride);
        precond_jacobi(A, b, thr_num, thread);
        memcpy(v1_tmp + start, arrays.v1.get() + start, stride * sizeof(double));
        solve(A, b, v1_tmp, s, r, p, desired_eps, thr_num, thread, max_it, it);
        fill_matrix_V2(arrays, msh, scheme, problem, n, A, b, start, stride);
        precond_jacobi(A, b, thr_num, thread);
        solve(A, b, arrays.v2.get(), s, r, p, desired_eps, thr_num, thread, max_it, it);
        memcpy(arrays.v1.get() + start, v1_tmp + start, stride * sizeof(double));
    }
    return nullptr;
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
    w = W_norm(msh, msh, arrays1.g.get(), arrays2.g.get(), scheme);
    printf("  g: c = %e, l = %e, w = %e\n", c, l, w);
    c = C_norm(msh, msh, arrays1.v1.get(), arrays2.v1.get());
    l = L_norm(msh, msh, arrays1.v1.get(), arrays2.v1.get(), scheme);
    w = W_norm(msh, msh, arrays1.v1.get(), arrays2.v1.get(), scheme);
    printf(" v1: c = %e, l = %e, w = %e\n", c, l, w);
    c = C_norm(msh, msh, arrays1.v2.get(), arrays2.v2.get());
    l = L_norm(msh, msh, arrays1.v2.get(), arrays2.v2.get(), scheme);
    w = W_norm(msh, msh, arrays1.v2.get(), arrays2.v2.get(), scheme);
    printf(" v2: c = %e, l = %e, w = %e\n", c, l, w);
}
