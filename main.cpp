#include <cstdio>
#include <cfenv>

#include "data.h"
#include "mesh.h"
#include "solve.h"
#include "norms.h"


int main(int argc, char *argv[])
{
    //feenableexcept(FE_ALL_EXCEPT ^ FE_INEXACT);
    double T, mu;
    unsigned int N, M, rho_type_int;
    if (!(argc == 6 && sscanf(argv[1], "%le", &T) == 1 && sscanf(argv[2], "%u", &N) == 1 && sscanf(argv[3], "%u", &M) == 1 && M >= 3 && sscanf(argv[4], "%le", &mu) == 1 && sscanf(argv[5], "%u", &rho_type_int) == 1 && rho_type_int >= 1 && rho_type_int <= 4))
    {
        printf("Usage: %s T N M mu rho\nT - maximum time\nN - amount of time steps in 1\nM - amount of space steps in 1\n", argv[0]);
        return -1;
    }

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

    double t = clock();

    set_initial_data(arrays, scheme);
    
    double t_set = clock();
    printf("Inited in %.2f\n", (t_set - t)/CLOCKS_PER_SEC);

    general_loop(problem, scheme, arrays);
    
    double t_solve = clock();
    printf("Solved in %.2f\n", (t_solve - t_set)/CLOCKS_PER_SEC);

    data arrays_solution(msh.size);
    fill_solution(arrays_solution, scheme);

    check_norms(arrays, arrays_solution, msh, scheme);
    printf("Counted norms in %.2f\n", (clock() - t_solve)/CLOCKS_PER_SEC);
    
    t = (clock() - t)/CLOCKS_PER_SEC;
    printf("Elapsed %.2f\n", t);

    return 0;
}
