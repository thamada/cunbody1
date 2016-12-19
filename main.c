//Time-stamp: <2008-03-02 14:42:42 hamada>

//#define TEST 1  // Force Error Check (with Result File)
//#define TEST 2  // Gflops Benchmark
//#define TEST 3  // Energy Error Check
#define TEST 4  // Visualization Demo with 3plot


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define NMAX (1<<17) // 131072
#include <time.h>
#include <unistd.h> // sleep()



//about timer
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h> 
#include <sys/resource.h>

double get_time(void)
{
  static struct timeval tv;
  static struct timezone tz;
  gettimeofday(&tv, &tz);
  return ((double)(tv.tv_sec  + tv.tv_usec*1.0e-6)); 
}




/* ----------------------------------------------------
 * Force Error Check (with Input/Result File)
 * ---------------------------------------------------- */
#if TEST == 1
#define EPS 0.1
static double r[NMAX][3];
static double m[NMAX];
static double eps2 = (EPS*EPS);
static double acc_d[NMAX][3];
static double acc_h[NMAX][3];
static double pot_h[NMAX];
static double v[NMAX][3];

//#include "force_host.h"
int main( int argc, char** argv) 
{
  int i, n;

  if(argc == 3){
    int nn;
    char ifile[256], rfile[256];
    strcpy(ifile,argv[1]);
    strcpy(rfile,argv[2]);
    init_particles(ifile, &n,  m,r,v);
    read_result   (rfile, &nn, acc_h, pot_h);
    if(n != nn){
      fprintf(stderr, "Abort at %s|%d\n"__FILE__, __LINE__);
      exit(-1);
    }
  }else{
    fprintf(stderr, "cmd <init_filename> <result_file\n");
    fprintf(stderr, "i.e. : run.gpu ../Dfile/init.plum.131072 ../Dfile/result.131072\n");
    fprintf(stderr, "exit program\n");
    exit(-1);
  }

 RESTART_FORCE_ERROR_CHECK:
  cunbody1_force(r, m, r, eps2, acc_d, n, n);
  //  force_host     (r, m, r, eps2, acc_d, n, n);

  /*
  {//--------------- generate result file
    double jk[65536][3];
    force_pot_jerk_host(r, v, m, pot_h, acc_h, jk, n);
    printf("%d\n", n);
    for(i=0;i<n;i++){
      //      printf("%1.16e\t%1.16e\t%1.16e\t%1.16e\t%1.16e\t%1.16e\t%1.16e\t");
      printf("%1.16e\t%1.16e\t%1.16e\t%1.16e\t%1.16e\t%1.16e\t%1.16e\t\n",
	     acc_h[i][0], acc_h[i][1], acc_h[i][2], pot_h[i],
	     jk[i][0], jk[i][1], jk[i][2]);
    }
    exit(0);
  }
  */

  for(i=0;i<n;i++){
    int d;
    static double err_max = 0.0;
    double err=0.0;
    for(d=0;d<3;d++){
      double diff = fabs(acc_d[i][d] - acc_h[i][d]);
      err += (diff*diff);
    }
    err = sqrt(err);
    err = err/sqrt(acc_h[i][0]*acc_h[i][0]+acc_h[i][1]*acc_h[i][1]+acc_h[i][2]*acc_h[i][2]);

    //    printf("HOST %.7e\t%.7e\t%.7e\n", acc_h[i][0], acc_h[i][1], acc_h[i][2]);
    //    printf("GPU  %.7e\t%.7e\t%.7e\n", acc_d[i][0], acc_d[i][1], acc_d[i][2]);

    if (err > err_max){
      err_max = err;
      printf("%.7e\n", err);
    }

    /*
    if (err > 1.0e-5){
      printf("HOST[%d] %.7e\t%.7e\t%.7e\n", i, acc_h[i][0], acc_h[i][1], acc_h[i][2]);
      printf("GPU [%d] %.7e\t%.7e\t%.7e\n", i, acc_d[i][0], acc_d[i][1], acc_d[i][2]);
    }
    */
  }

  goto RESTART_FORCE_ERROR_CHECK;
  return (0);
}
#endif



/* ----------------------------------------------------
 * Gflops Benchmark
 * ---------------------------------------------------- */
#if TEST == 2

//#include <force_host.h>

static double r[NMAX][3];
static double m[NMAX];
static double eps2 = 0.0001;
static double a[NMAX][3];
static double v[NMAX][3];

#define FORCE cunbody1_force
//#define FORCE  force_host

int main(int argc,char** argv)
{
  double tt=0.0;
  int n;
  int step, nstep = 0;

  if((argc == 2)||(argc == 3)){
    char ifile[256];
    strcpy(ifile,argv[1]);
    init_particles(ifile, &n, m, r, v);
    if(argc == 3)
      nstep = atoi(argv[2]);
    else
      nstep = 1;

  }else{
    fprintf(stderr, "argc = %d\n",argc);
    fprintf(stderr, "cmd <init file> <nstep>\n");
    exit(-1);
  }

  FORCE(r, m, r, eps2, a, n, n); // wakeup GPU
RESTART_GFLOPS_BENCH:
  tt = get_time();
  for(step = 0; step < nstep; step++) FORCE(r, m, r, eps2, a, n, n);
  tt = get_time() - tt;
  printf("sec: %g\n", tt);

  {
    double time = tt / ((double)nstep);
    double gflops = 38.0 * ((double)n) * ((double)n) * (1.0e-9) / time;
    printf("N,Gflop/s: %d\t%g\n", n, gflops);
  }

  goto RESTART_GFLOPS_BENCH;
  return (0);
}
#endif

/* -------------------
 * Energy Check 
 * ------------------- */
#if TEST == 3
#include <energy.c>

#define TMAX 1500
#define EPS2 0.0001
#define TIMESTEP 0.005

static double r[NMAX][3];
static double m[NMAX];
static double eps2 = EPS2;
static double a[NMAX][3];
static double v[NMAX][3];
static double a_h[NMAX][3];

#define FORCE  cunbody1_force
//#define FORCE  force_host


int main(int argc,char** argv)
{
  int n;
  double dt = TIMESTEP;
  int time;

  if(argc == 2){
    char ifile[256];
    strcpy(ifile,argv[1]);
    init_particles(ifile, &n, m, r, v);
  }else{
    fprintf(stderr, "cmd <init_filename>\n");
    fprintf(stderr, "exit program\n");
    exit(-1);
  }

  double e_init;
  e_init = energy(m, r, v, eps2, n);
  force (r, m, r, eps2, a, n, n);

  for(time=0;time<10000;time++){
    //    fprintf(stderr, "%d times\n", time);
    leapflog_half (dt,v,a,n);
    leapflog      (dt,r,v,n);
    FORCE(r,m,r,eps2,a,n, n);
    leapflog_half (dt,v,a,n);

    if( (time%10) == 0) {
      double e = energy(m, r, v, eps2, n);
      double err = fabs(e-e_init)/e_init;
      printf("%d : E(int) %e, E(comp) %e, Eerr %e\n", time, e_init, e, err);
      fflush(NULL);
    }
  }

  return (0);
}
#endif


/* -------------------
 * Visualization Demo with 3plot
 * [ The logfile will appear at /dev/shm/xxx.log ]
 * ------------------- */
#if TEST == 4

#include "force_host.h"
#define TMAX 1500
#define EPS2 0.0001
#define TIMESTEP 0.005

static double r[NMAX][3];
static double m[NMAX];
static double eps2 = EPS2;
static double a[NMAX][3];
static double v[NMAX][3];

static double acc_h[NMAX][3];


#define FORCE  cunbody1_force
//#define FORCE  force_host



int main(int argc,char** argv)
{
  int n;
  double dt = TIMESTEP;
  double tt=0.0;
  int time;

  /* initialize particles  */
  
  if(argc == 2){
    char ifile[256];
    strcpy(ifile,argv[1]);
    init_particles(ifile,&n,m,r,v);
  }else{
    fprintf(stderr, "cmd <init_filename>\n");
    fprintf(stderr, "exit program\n");
    exit(-1);
  }

  FORCE (r, m, r, eps2, a, n, n);
  //  force (r, m, r, eps2, a, n, n);


  //  for(time=0;time<1000;time++){
  for(time=0;;time++){
    int nstep = 2;
    int step;
    double flops=0.0;
    tt = get_time();

    for(step = 0; step < nstep; step++){
      leapflog_half (dt, v, a, n);
      leapflog      (dt, r, v, n);
      FORCE(r, m, r, eps2, a, n, n);
      leapflog_half (dt, v, a, n);
    }

    tt = get_time() - tt;
    flops = 38.0 * ((double)n) * ((double)n) * ((double)nstep) * 1.0e-9 / tt; // Gflops
    printf("N=%d, Step=%d, %.1f Gflops\n",n, time, flops);
    debug_position_snap(r, flops ,n, time*nstep);
  }

  debug_position_dump(r, n, 0);

  return (0);
}
#endif

