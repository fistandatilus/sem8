#include <pthread.h>
#include <cstdio>
#include <sys/sysinfo.h>
#include <sys/time.h>
#include <cfenv>

#include "data.h"
#include "mesh.h"
#include "solve.h"
#include "norms.h"
#include "functions.h"

#define N_AMOUNT 6
#define M_AMOUNT 4


void *thread_main(void *args_ptr);
int get_next_task(int task_count = 0);

struct task_t
{
    problem_params problem{};
    scheme_params scheme{};
    data arrays{};
    double time{};

    task_t() = default;

    void init(const unsigned int N, const unsigned int M, const double T, const double mu, const rho_type type)
    {
        problem.mu = mu;
        problem.type = type;
        scheme.N = N;
        scheme.M = M;
        scheme.T = T;
        scheme.tau = 1. / N;
        scheme.h = 1. / M;

        mesh msh(M);
        arrays.init(msh.size);
    };
};

struct args_t{
    task_t *tasks = nullptr;
    int thread = 0;
    int size = 0;
};

double get_full_time() {
    timeval buf;
    gettimeofday(&buf, 0);
    return buf.tv_sec + buf.tv_usec / 1e6;
}

#define MULTITHREAD (true)

int main(int argc, char *argv[])
{
if constexpr (!MULTITHREAD)
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

    set_mu(mu);
    set_p_coeff(type);

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
else
{
    double T;
    unsigned int N0, M0, n_threads;
    if (!(argc == 5 && sscanf(argv[1], "%le", &T) == 1 && sscanf(argv[2], "%u", &N0) == 1 && sscanf(argv[3], "%u", &M0) == 1 && M0 >= 3 && sscanf(argv[4], "%u", &n_threads) == 1))
    {
        printf("Usage: %s T N0 M0 n_threads\nT - maximum time\nN - amount of time steps in 1\nM - amount of space steps in 1\n", argv[0]);
        return -1;
    }

    for (unsigned int rho_type_int = 1; rho_type_int <= 4; rho_type_int++)
    { 
        double mu = 0.1;
        for (unsigned int mu_count = 0; mu_count < 3; mu_count++, mu *= 0.1)
        {
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
            set_mu(mu);
            set_p_coeff(type);

            std::unique_ptr<task_t []> tasks = std::make_unique<task_t []>(N_AMOUNT * M_AMOUNT);
            for (unsigned int i = 0, N = N0; i < N_AMOUNT; i++, N <<= 1u)
                for (unsigned int j = 0, M = M0; j < M_AMOUNT; j++, M <<= 1u)
                {
                    tasks[i * M_AMOUNT + j].init(N, M, T, mu, type);
                }

            std::unique_ptr<pthread_t []> tid = std::make_unique<pthread_t []>(n_threads);
            std::unique_ptr<args_t []> args = std::make_unique<args_t []>(n_threads);

            get_next_task(N_AMOUNT * M_AMOUNT);

            for (unsigned int i = 0; i < n_threads; i++)
            {
                args[i].tasks = tasks.get();
                args[i].thread = i;
                if (i) {
                    pthread_create(tid.get() + i, 0, thread_main, args.get() + i);
                }
            }
            thread_main(args.get() + 0);

            for (unsigned int i = 1; i < n_threads; i++) {
                pthread_join(tid[i], 0);
            }

            fprintf(stderr, "Done!\n");

            printf("\\subsubsection{$\\mu = %.3f, mode = %d$}\n", mu, rho_type_int);
            printf("$$\\Vert g - \\ln(\\rho)\\Vert$$\n");
            printf("\\begin{tabular}{*{%d}{|c}|}\n\\hline\n", M_AMOUNT + 1);
            printf("\\diagbox{$\\tau$}{$h$}");
            for (int i = 1; i <= M_AMOUNT; i++)
                printf("&%#.5g", 1. / (N0 << (i - 1)));
            printf("\\\\\n\\hline\n");
            for (int i = 0; i < N_AMOUNT; i++)
            {
                printf("%#.5g", 1. / (M0 << i));
                for (int j = 0; j < M_AMOUNT; j++)
                {
                    task_t &task = tasks[i * M_AMOUNT + j];
                    mesh msh(task.scheme.M);
                    data solution_arrays(msh.size);
                    fill_solution(solution_arrays, task.scheme);
                    double norm = C_norm(msh, msh, task.arrays.g.get(), solution_arrays.g.get());
                    printf("&$%.3e$", norm);
                }
                printf("\\\\\n");
                for (int j = 0; j < M_AMOUNT; j++)
                {
                    task_t &task = tasks[i * M_AMOUNT + j];
                    mesh msh(task.scheme.M);
                    data solution_arrays(msh.size);
                    fill_solution(solution_arrays, task.scheme);
                    double norm = L_norm(msh, msh, task.arrays.g.get(), solution_arrays.g.get(), task.scheme);
                    printf("&$%.3e$", norm);
                }
                printf("\\\\\n");
                for (int j = 0; j < M_AMOUNT; j++)
                {
                    task_t &task = tasks[i * M_AMOUNT + j];
                    mesh msh(task.scheme.M);
                    data solution_arrays(msh.size);
                    fill_solution(solution_arrays, task.scheme);
                    double norm = W_norm(msh, msh, task.arrays.g.get(), solution_arrays.g.get(), task.scheme);
                    printf("&$%.3e$", norm);
                }
                printf("\\\\\n");
                for (int j = 0; j < M_AMOUNT; j++)
                {
                    task_t &task = tasks[i * M_AMOUNT + j];
                    printf("&$%.3f$", task.time);
                }
                printf("\\\\\n\\hline\n");
            }
            printf("\\end{tabular}\n");
            
            printf("$$\\Vert v_1 - u_1 \\Vert$$\n");
            printf("\\begin{tabular}{*{%d}{|c}|}\n\\hline\n", M_AMOUNT + 1);
            printf("\\diagbox{$\\tau$}{$h$}");
            for (int i = 1; i <= M_AMOUNT; i++)
                printf("&%#.5g", 1. / (N0 << (i - 1)));
            printf("\\\\\n\\hline\n");
            for (int i = 0; i < N_AMOUNT; i++)
            {
                printf("%#.5g", 1. / (M0 << i));
                for (int j = 0; j < M_AMOUNT; j++)
                {
                    task_t &task = tasks[i * M_AMOUNT + j];
                    mesh msh(task.scheme.M);
                    data solution_arrays(msh.size);
                    fill_solution(solution_arrays, task.scheme);
                    double norm = C_norm(msh, msh, task.arrays.v1.get(), solution_arrays.v1.get());
                    printf("&$%.3e$", norm);
                }
                printf("\\\\\n");
                for (int j = 0; j < M_AMOUNT; j++)
                {
                    task_t &task = tasks[i * M_AMOUNT + j];
                    mesh msh(task.scheme.M);
                    data solution_arrays(msh.size);
                    fill_solution(solution_arrays, task.scheme);
                    double norm = L_norm(msh, msh, task.arrays.v1.get(), solution_arrays.v1.get(), task.scheme);
                    printf("&$%.3e$", norm);
                }
                printf("\\\\\n");
                for (int j = 0; j < M_AMOUNT; j++)
                {
                    task_t &task = tasks[i * M_AMOUNT + j];
                    mesh msh(task.scheme.M);
                    data solution_arrays(msh.size);
                    fill_solution(solution_arrays, task.scheme);
                    double norm = W_norm(msh, msh, task.arrays.v1.get(), solution_arrays.v1.get(), task.scheme);
                    printf("&$%.3e$", norm);
                }
                printf("\\\\\n");
                for (int j = 0; j < M_AMOUNT; j++)
                {
                    task_t &task = tasks[i * M_AMOUNT + j];
                    printf("&$%.3f$", task.time);
                }
                printf("\\\\\n\\hline\n");
            }
            printf("\\end{tabular}\n");
            
            printf("$$\\Vert v_2 - u_2\\Vert$$\n");
            printf("\\begin{tabular}{*{%d}{|c}|}\n\\hline\n", M_AMOUNT + 1);
            printf("\\diagbox{$\\tau$}{$h$}");
            for (int i = 1; i <= M_AMOUNT; i++)
                printf("&%#.5g", 1. / (N0 << (i - 1)));
            printf("\\\\\n\\hline\n");
            for (int i = 0; i < N_AMOUNT; i++)
            {
                printf("%#.5g", 1. / (M0 << i));
                for (int j = 0; j < M_AMOUNT; j++)
                {
                    task_t &task = tasks[i * M_AMOUNT + j];
                    mesh msh(task.scheme.M);
                    data solution_arrays(msh.size);
                    fill_solution(solution_arrays, task.scheme);
                    double norm = C_norm(msh, msh, task.arrays.v2.get(), solution_arrays.v2.get());
                    printf("&$%.3e$", norm);
                }
                printf("\\\\\n");
                for (int j = 0; j < M_AMOUNT; j++)
                {
                    task_t &task = tasks[i * M_AMOUNT + j];
                    mesh msh(task.scheme.M);
                    data solution_arrays(msh.size);
                    fill_solution(solution_arrays, task.scheme);
                    double norm = L_norm(msh, msh, task.arrays.v2.get(), solution_arrays.v2.get(), task.scheme);
                    printf("&$%.3e$", norm);
                }
                printf("\\\\\n");
                for (int j = 0; j < M_AMOUNT; j++)
                {
                    task_t &task = tasks[i * M_AMOUNT + j];
                    mesh msh(task.scheme.M);
                    data solution_arrays(msh.size);
                    fill_solution(solution_arrays, task.scheme);
                    double norm = W_norm(msh, msh, task.arrays.v2.get(), solution_arrays.v2.get(), task.scheme);
                    printf("&$%.3e$", norm);
                }
                printf("\\\\\n");
                for (int j = 0; j < M_AMOUNT; j++)
                {
                    task_t &task = tasks[i * M_AMOUNT + j];
                    printf("&$%.3f$", task.time);
                }
                printf("\\\\\n\\hline\n");
            }
            printf("\\end{tabular}\n");
            printf("\n\n");
        }
    }
    return 0;
}
}

void *thread_main(void *args_ptr){
    args_t &args = *((args_t *)args_ptr);
    task_t *tasks = args.tasks;
    cpu_set_t cpu;
    CPU_ZERO(&cpu);
    int nproc = get_nprocs();
    int cpu_id = nproc - 1 - args.thread%nproc;
    CPU_SET(cpu_id, &cpu);
    pthread_t tid = pthread_self();
    pthread_setaffinity_np(tid, sizeof(cpu_set_t), &cpu);
    for (int i = get_next_task(); i >= 0; i = get_next_task())
    {
        task_t &task = tasks[i];
        scheme_params &scheme = task.scheme;
        problem_params &problem = task.problem;
        data &arrays = task.arrays;

        double time = get_full_time();
        set_initial_data(arrays, scheme);
        general_loop(problem, scheme, arrays);
        task.time = get_full_time() - time;

        fprintf(stderr, "thread = %u, N = %u, M = %u, time = %.2f\n", args.thread, scheme.N, scheme.M, task.time);
    }
    return nullptr;
}

int get_next_task(int task_count){
    //инициализируется, если task_count != 0
    static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
    static int tasks_remaining = -1;
    int ret = -2;
    
    if (!task_count && tasks_remaining < 0) {
        return -1;
    }
    pthread_mutex_lock(&mutex);
    if (task_count) {
        tasks_remaining = task_count-1;
    }
    else {
        ret = tasks_remaining--;
    }
    pthread_mutex_unlock(&mutex);
    return ret;
}

