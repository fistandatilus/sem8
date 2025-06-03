#include <cstdio>
#include <cmath>

#include "mesh.h"
#include "msr.h"

unsigned int msr_size(const unsigned int M, const unsigned int size)
{
    return (size + 1 + 4 * ((M - 1) * (M - 1) * 4 + (M - 1) * 3 + 2));
}

void msr::set_template(const unsigned int *ind, const unsigned int n, const unsigned int size)
{
    indexes = ind;
    this->n = n;
    this->size = size;
    data = std::make_unique<double []>(size);

    norm = 0;
}

void msr::copy_template(const msr &x)
{
    set_template(x.indexes, x.n, x.size);
}

double &msr::coeffRef(const unsigned int i, const unsigned int j)
{
    if (i == j)
        return data[i];
    const unsigned int len = indexes[i + 1] - indexes[i] ;
    for (unsigned int k = 0; k < len; k++)
    {
        if (indexes[indexes[i] + k] == j)
            return data[indexes[i] + k];
    }
    printf("invalid msr offdiag access\n");
    return data[i];
}

std::unique_ptr<unsigned int []> make_template(const mesh &msh)
{
    const unsigned int M = msh.width;
    const unsigned int size = msh.size;
    std::unique_ptr<unsigned int []> indexes = std::make_unique<unsigned int []>(msr_size(M, size));

    unsigned int offset = size + 1;

    for (unsigned int i_matrix = 0; i_matrix < size; i_matrix++)
    {
        indexes[i_matrix] = offset;
        unsigned int i, j;
        msh.get_mesh_by_matrix(i_matrix, i, j);
        point_type type = msh.get_type_by_mesh(i, j);
        if (type == point_type::i || type == point_type::ilt || type == point_type::irb)
        {
            const unsigned int iL = msh.get_matrix_by_mesh(i - 1, j);
            const unsigned int iR = msh.get_matrix_by_mesh(i + 1, j);
            const unsigned int iB = msh.get_matrix_by_mesh(i, j - 1);
            const unsigned int iT = msh.get_matrix_by_mesh(i, j + 1);
            indexes[offset] = iL;
            indexes[offset + 1] = iR;
            indexes[offset + 2] = iB;
            indexes[offset + 3] = iT;
            offset += 4;
        }
    }

    indexes[size] = offset;

    if(offset != msr_size(M, size))
        printf("AAA wrong offset!!! excpected = %d, got = %d\n", msr_size(M, size), offset);

    return std::move(indexes);
}

#define UNFOLD 8

void start_and_size(unsigned int p, unsigned int thread, unsigned int n, unsigned int &start, unsigned int &size)
{
  size = n/p;
  unsigned int res = n - size*p;
  start = size*thread + (thread < res ? thread : res);
  if (thread < res) size++;

}

void mul_msr_by_vec(const msr &a, const double *x, double *ret, unsigned int start, unsigned int stride)
{
  for (unsigned int i = start; i < start + stride; i++)
  {
    double s = a.data[i] * x[i];
    unsigned int len = a.indexes[i+1] - a.indexes[i];
    for (unsigned int j = 0; j < len; j++)
      s += a.data[a.indexes[i] + j] * x[a.indexes[a.indexes[i] + j]];
    ret[i] = s;
  }
}

void precond_jacobi(msr &a, double *b, const unsigned int thread_num, const unsigned int thread)
{
    unsigned int stride, start, n = a.n;
    start_and_size(thread_num, thread, n, start, stride);
    for (unsigned int i = start; i < start + stride; i++)
    {
        double c = a.data[i];
        a.data[i] /= c;
        b[i] /= c;
        unsigned int len = a.indexes[i+1] - a.indexes[i];
        for (unsigned int j = 0; j < len; j++)
            a.data[a.indexes[i] + j] /= c;
    }
    
}

void solve(const msr &a, const double *b, double *S, double *s, double *r, double *p, double desired_eps, unsigned int thr_num, unsigned int thread, unsigned int max_it, unsigned int &iter)
{
    unsigned int stride, start, n = a.n;
    start_and_size(thr_num, thread, n, start, stride);

    // S given right from outside
    // z virtual, all same 

    reduce_sum<int>(thr_num);
    mul_msr_by_vec(a, S, r, start, stride);
    double r_norm_prev = 0;
    for (unsigned int i = start; i < start + stride; i++)
    {
        r[i] = b[i] - r[i];
        p[i] = r[i];
        r_norm_prev += r[i] * r[i];
    }
    reduce_sum(thr_num, &r_norm_prev, 1);
    double eps = 0;
    for (unsigned int i = start; i < start + stride; i++)
        eps += b[i] * b[i];
    reduce_sum(thr_num, &eps, 1);
    eps *= desired_eps * desired_eps;
    if (r_norm_prev <= eps * eps)
        return;

    for (iter = 1; iter <= max_it; iter++)
    {
        double ajp[2] = {0, 0};
        for (unsigned int i = start; i < start + stride; i++)
        {
            ajp[0] += r[i];

            double c = a.data[i] * p[i];
            unsigned int len = a.indexes[i+1] - a.indexes[i];
            for (unsigned int j = 0; j < len; j++)
                c += a.data[a.indexes[i] + j] * p[a.indexes[a.indexes[i] + j]];
            ajp[1] += c;
        }
        reduce_sum(thr_num, ajp, 2);
        double aj = ajp[0] / ajp[1];

        double s_norm = 0;
        for (unsigned int i = start; i < start + stride; i++)
        {
            double c = a.data[i] * p[i];
            unsigned int len = a.indexes[i+1] - a.indexes[i];
            for (unsigned int j = 0; j < len; j++)
                c += a.data[a.indexes[i] + j] * p[a.indexes[a.indexes[i] + j]];
            s[i] = r[i] - aj * c;
            s_norm += s[i] * s[i];
        }
        reduce_sum(thr_num, &s_norm, 1);
        if (s_norm <= desired_eps * desired_eps)
        {
            for (unsigned int i = start; i < start + stride; i++)
                S[i] += aj * p[i];
            return;
        }
        double wjp[2] = {0, 0};
        for (unsigned int i = start; i < start + stride; i++)
        {
            double c = a.data[i] * s[i];
            unsigned int len = a.indexes[i+1] - a.indexes[i];
            for (unsigned int j = 0; j < len; j++)
                c += a.data[a.indexes[i] + j] * s[a.indexes[a.indexes[i] + j]];
            wjp[0] += s[i] * c;
            wjp[1] += c * c;
        }
        reduce_sum(thr_num, wjp, 2);
        double wj = wjp[0] / wjp[1];

        double r_norm = 0;
        double bjp[3] = {0, 0, 0};
        for (unsigned int i = start; i < start + stride; i++)
        {
            S[i] += aj * p[i] + wj * s[i];

            double c = a.data[i] * s[i];
            unsigned int len = a.indexes[i+1] - a.indexes[i];
            for (unsigned int j = 0; j < len; j++)
                c += a.data[a.indexes[i] + j] * s[a.indexes[a.indexes[i] + j]];
            bjp[1] += r[i];
            r[i] = s[i] - wj * c;
            bjp[0] += r[i];
            r_norm += r[i] * r[i];
        }
        bjp[2] = r_norm;
        reduce_sum(thr_num, bjp, 3);
        r_norm = bjp[2];
        if (r_norm <= eps * eps || r_norm / r_norm_prev > 0.99)
            return;
        r_norm_prev = r_norm;
        double bj = bjp[0] / bjp[1] * aj / wj;

        for (unsigned int i = start; i < start + stride; i++)
        {
            double c = a.data[i] * p[i];
            unsigned int len = a.indexes[i+1] - a.indexes[i];
            for (unsigned int j = 0; j < len; j++)
                c += a.data[a.indexes[i] + j] * p[a.indexes[a.indexes[i] + j]];
            s[i] = c;
        }
        reduce_sum<int>(thr_num);
        for (unsigned int i = start; i < start + stride; i++)
        {
            p[i] = r[i] + bj * (p[i] - wj * s[i]);
        }
        reduce_sum<int>(thr_num);
    }
}
