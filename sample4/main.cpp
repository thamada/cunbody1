/* sample.c
 * by Tsuyoshi Hamada
 */
#include <omp.h>
#include <cstdio>
#include <cstdlib>
#include <cassert>

using namespace std;

void check_Mem()
{
  FILE *fp = popen("bash", "w");
  assert(fp != NULL);
  fprintf(stdout, "Memory Info:\n");  fflush(stdout);
  fprintf(fp, "free -m \n");  fflush(fp);
  pclose(fp);
}

void check_CPUs()
{
  FILE *fp = popen("bash", "w");
  assert(fp != NULL);
  fprintf(stdout, "CPU Info:\n");  fflush(stdout);
  fprintf(fp, "grep 'model name' /proc/cpuinfo\n");  fflush(fp);
  pclose(fp);
}

void check_kernel()
{
  fprintf(stdout, "Linux Kernel Version: ");
  fflush(stdout);
  FILE *fp = popen("bash", "w");
  assert(fp != NULL);
  fprintf(fp, "uname -a\n");
  fflush(fp);
  pclose(fp);
}

void check_nvidia_driver()
{
  fprintf(stdout, "GPU Driver Version: ");
  fflush(stdout);
  FILE *fp = popen("bash", "w");
  assert(fp != NULL);
  fprintf(fp, "ls -F /usr/lib64/libcuda* | grep '*'\n");
  fflush(fp);
  pclose(fp);
}

void check_system()
{
  check_CPUs();
  check_Mem();
  check_kernel();
  check_nvidia_driver();
}

#include <sys/time.h>
#include <sys/resource.h>
double get_time()
{
  static struct timeval tv;
  static struct timezone tz;
  gettimeofday(&tv, &tz);
  return ((double)(tv.tv_sec  + tv.tv_usec*1.0e-6)*1000);
}

#define NMAX (131072)
double xj[NMAX][3];
double mj[NMAX];
double xi[NMAX][3];
double ai[NMAX][3];

/* 
 * cunbody1_force is the only subroutine of libcunbody1.a.
 * cunbody1_force calculates gravitational interactions between
 * i-th and j-th particles.
 */

extern "C"
void cunbody1_force(double xj[][3], // position of j-th particles
                    double mj[],    // mass of j-th particles
                    double xi[][3], // position of i-th particles
                    double eps2,    // softening parameter
                    double ai[][3], // force of i-th particles
                    int ni,         // number of i-th particles
                    int nj);        // number of j-th particles

extern "C"
void cunbody1_force_mp(int devid,
		       double xj[][3], // position of j-th particles
		       double mj[],    // mass of j-th particles
		       double xi[][3], // position of i-th particles
		       double eps2,    // softening parameter
		       double ai[][3], // force of i-th particles
		       int ni,         // number of i-th particles
		       int nj);        // number of j-th particles


/*  main routine */

int main(int argc, char** argv)
{
  check_system();

  int i, dim, ni, nj, nforce;
  double t;
  double eps2 = 0.0001;
  int devid = 0;
  int ngpu = omp_get_max_threads();

  puts("---------------------------------------------------");
  printf("# of OMP threads use %d GPUs: %d\n", ngpu, ngpu);
  puts("---------------------------------------------------");
  assert(ngpu>0);
  srand(0x19740526);

  /* SETUP I,J-PARTICLES */
  ni = nj = 16384; //65536;

  for(i=0; i<ni; i++){
    for(dim=0; dim<3; dim++)
      xi[i][dim] = xj[i][dim] = rand()/(RAND_MAX + 1.0);
    mj[i] = 1.0e-5;
  }

  /* CALL GeForce8800GTX */
#pragma omp parallel for
  for(i=0; i<ngpu; i++){
    cunbody1_force_mp(i%ngpu, xj, mj, xi, eps2, ai, ni, nj); 
  }

  nforce = 60;
  t = get_time();

#pragma omp parallel for
  for(i=0; i< nforce; i++){
    cunbody1_force_mp(i%ngpu, xj, mj, xi, eps2, ai, ni, nj); 
  }

  t = get_time() - t;

  {
    double tsec = t * 1.0e-3;
    double gflops = ((double)ni)*((double)nj)*((double)nforce)*38.0*(1.0e-9)/tsec;
    printf ("%g [sec]\n", tsec);
    printf ("%g [interactions/sec]\n", ((double)ni)*((double)nj) / tsec);
    printf ("%g [forces/sec]\n", nforce/tsec);
    printf ("%g [Gflops]\n", gflops);
  }

  return (0);
}
