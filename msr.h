#ifndef MSR_H
#define MSR_H

#include <memory>

class mesh;

struct msr
{
    std::unique_ptr<double []> data{};
    const unsigned int *indexes{}; //does not own
    double norm{};
    unsigned int n{};
    unsigned int size{};

  void erase()
  {
    data = nullptr;
    indexes = nullptr;
    n = 0;
    size = 0;
    norm = 0;
  }
  ~msr()
  {
    erase();
  }
    double &coeffRef(const unsigned int i, const unsigned int j);
    void set_template(const unsigned int *ind, const unsigned int n, const unsigned int size);
    msr() = default;
    msr(msr &) = delete;
    msr &operator=(msr &) = delete;
    void copy_template(const msr &x);
};

unsigned int msr_size(const unsigned int M, const unsigned int size);
std::unique_ptr<unsigned int []> make_template(const mesh &msh);

void start_and_size(unsigned int p, unsigned int thread, unsigned int n, unsigned int &start, unsigned int &size);
void mul_msr_by_vec(const msr &a, const double *x, double *ret, unsigned int start, unsigned int stride);
void precond_jacobi(msr &a, double *b, const unsigned int thread_num, const unsigned int thread);
void solve(const msr &a, const double *b, double *S, double *s, double *r, double *p, double desired_eps, unsigned int thr_num, unsigned int thread, unsigned int max_it, unsigned int &iter);

template <class T>
void reduce_sum(int p, T *a = nullptr, unsigned int n = 0);
void start_and_size(unsigned int p, unsigned int thread, unsigned int n, unsigned int &start, unsigned int &size);

template <class T>
void reduce_sum(int p, T *a, unsigned int n)
{
  static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
  static pthread_cond_t condvar_in = PTHREAD_COND_INITIALIZER;
  static pthread_cond_t condvar_out = PTHREAD_COND_INITIALIZER;
  static int threads_in = 0, threads_out = 0;
  static T *pres = nullptr;
  pthread_mutex_lock(&mutex);
  if (pres == nullptr)
  {
    pres = a;
  }
  else
  {
    for (unsigned int i = 0; i < n; i++)
    {
      pres[i] += a[i];
    }
  }
  threads_in++;
  if (threads_in >= p)
  {
    threads_out = 0;
    pthread_cond_broadcast(&condvar_in);
  }
  else
  {
    while (threads_in < p)
      pthread_cond_wait(&condvar_in, &mutex);
  }
  if (pres != a)
  {
    for (unsigned int i = 0; i < n; i++)
      a[i] = pres[i];
  }
  threads_out++;
  if (threads_out >= p)
  {
    pres = nullptr;
    threads_in = 0;
    pthread_cond_broadcast(&condvar_out);
  }
  else
    while (threads_out < p)
      pthread_cond_wait(&condvar_out, &mutex);
  pthread_mutex_unlock(&mutex);
}
#endif //MSR_H
