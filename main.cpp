#include <cstdio>

#include "data.h"
#include "mesh.h"

int main(int argc, char *argv[])
{
    double T, mu;
    unsigned int N, M, rho_type_int;
    if (!(argc == 6 && sscanf(argv[1], "%le", &T) == 1 && sscanf(argv[2], "%d", &N) && sscanf(argv[3], "%d", &M) && M >= 3 && sscanf(argv[4], "%le", &mu) && sscanf(argv[5], "%d", &rho_type_int) && rho_type_int >= 1 && rho_type_int <= 4))
        printf("Usage: %s T N M mu rho\nT - maximum time\nN - amount of time steps in 1\nM - amount of space steps in 1\n", argv[0]);

    rho_type type;
    switch (rho_type_int)
    {
    case 1:
        type = rho_type::lin1;
        break;
    case 2:
        type = rho_type::lin10;
        break;
    case 3:
        type = rho_type::lin100;
        break;
    case 4:
        type = rho_type::gamma;
        break;
    }

    problem_params problem(mu, type);
    scheme_params scheme(N, M, T);
    mesh msh(M);
    data arrays(msh.size);
    return 0;
}
