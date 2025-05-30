#include <cstdio>

int main(int argc, char *argv[])
{
    double T, mu;
    unsigned int N, M, rho_type;
    if (!(argc == 6 && sscanf(argv[1], "%le", &T) == 1 && sscanf(argv[2], "%d", &N) && sscanf(argv[3], "%d", &M) && M >= 3 && sscanf(argv[4], "%le", &mu) && sscanf(argv[5], "%d", &rho_type) && rho_type >= 1 && rho_type <= 4))
        printf("Usage: %s T N M mu rho\nT - maximum time\nN - amount of time steps in 1\nM - amount of space steps in 1\n", argv[0]);

    return 0;
}
