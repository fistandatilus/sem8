#include "eigen/Eigen/src/IterativeLinearSolvers/BiCGSTAB.h"

#include "solve.h"
#include "data.h"

void general_loop(problem_params &problem, scheme_params &scheme)
{
    const unsigned int N = scheme.N;

    for (unsigned int n = 0; n < N; n++)
    {
        fill_matrix_G();
        solve();
        fill_G_matrix();
        fill_matrix_V1();
        solve();
        fill_matrix_V2();
        solve();        // NOTE: mulithread?
        fill_V1_matrix();
        fill_V2_matrix();
    }
}

void fill_matrix_G()
{
    
}
