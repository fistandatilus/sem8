#ifndef DATA_H
#define DATA_H

#include<memory>

enum class rho_type
{
    lin1,
    lin10,
    lin100,
    gamma,
};

class data
{
public:
    std::unique_ptr<double []> g{};
    std::unique_ptr<double []> h{};
    std::unique_ptr<double []> v1{};
    std::unique_ptr<double []> v2{};

    unsigned int size = 0;

    data(const unsigned int size) : size(size)
    {
        g  = std::make_unique<double []>(size);
        h  = std::make_unique<double []>(size);
        v1 = std::make_unique<double []>(size);
        v2 = std::make_unique<double []>(size);
    }
};

struct problem_params
{
    double mu = 0;
    rho_type type = rho_type::lin1;

    problem_params(double mu, rho_type type) : mu(mu), type(type)
    {};
};

struct scheme_params
{
    unsigned int N = 0;
    unsigned int M = 0;
    double tau = 0;
    double h = 0;
    double T = 0;

    scheme_params(unsigned int N, unsigned int M, double T) : N(N), M(M), T(T), tau(1. / N), h(1. / M)
    {};
};

#endif //DATA_H
