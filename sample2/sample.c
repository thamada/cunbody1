/* sample.c
 * by Tsuyoshi Hamada
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#define NMAX (131072)

#include <sys/time.h>
#include <sys/resource.h>
double get_time()
{
  static struct timeval tv;
  static struct timezone tz;
  gettimeofday(&tv, &tz);
  return ((double)(tv.tv_sec  + tv.tv_usec*1.0e-6)*1000);
}

double xj[NMAX][3];
double mj[NMAX];
double xi[NMAX][3];
double ai[NMAX][3];

/* 
 * cunbody1_force is the only subroutine of libcunbody1.a.
 * cunbody1_force calculates gravitational interactions between
 * i-th and j-th particles.
 */

void cunbody1_force(double xj[][3], // position of j-th particles
                    double mj[],    // mass of j-th particles
                    double xi[][3], // position of i-th particles
                    double eps2,    // softening parameter
                    double ai[][3], // force of i-th particles
                    int ni,         // number of i-th particles
                    int nj);        // number of j-th particles


void cunbody1_force_mp(int devid,
		    double xj[][3], // position of j-th particles
                    double mj[],    // mass of j-th particles
                    double xi[][3], // position of i-th particles
                    double eps2,    // softening parameter
                    double ai[][3], // force of i-th particles
                    int ni,         // number of i-th particles
                    int nj);        // number of j-th particles


#include <math.h>
void force_host(int devid,
		double xj[][3],
		double m[],
		double xi[][3],
		double eps2,
		double a[][3],
		int ni,
		int nj)
{
  int i,j,d;
  double dx[3];
  for(i=0;i<ni;i++) {
    a[i][0] = 0.0;
    a[i][1] = 0.0;
    a[i][2] = 0.0;
  }

  for(i=0;i<ni;i++){
    double rxi = xi[i][0];
    double ryi = xi[i][1];
    double rzi = xi[i][2];
    double ax,ay,az;
    ax = ay = az = 0.0;
    for(j=0;j<nj;j++){
      double r2,r3,r3i;
      double dx = xj[j][0] - rxi;
      double dy = xj[j][1] - ryi;
      double dz = xj[j][2] - rzi;
      double mj = m[j];
      r2 = eps2+(dx*dx)+(dy*dy)+(dz*dz);
      r3 = sqrt(r2)*r2;
      r3i = 1.0/r3;
      ax += mj*dx*r3i;
      ay += mj*dy*r3i;
      az += mj*dz*r3i;
    }


    a[i][0] =  ax;
    a[i][1] =  ay;
    a[i][2] =  az;

  }
}


/*  main routine */

int main(int argc, char** argv)
{
  int i, dim, ni, nj, nforce;
  double t; 
  double eps2 = 0.0001;
  int devid = 0;

  if(argc==2) devid = atoi(argv[1]) - 1;

  puts("-------------------------");
  printf("Using %d-th GPU.\n", devid+1);
  puts("-------------------------");
  assert(devid>-1);

  srand(0x19740526);

  /* SETUP I,J-PARTICLES */
  //  ni = nj = 65536;
  ni = nj = 8192;
  for(i=0; i<ni; i++){
    for(dim=0; dim<3; dim++)
      xi[i][dim] = xj[i][dim] = rand()/(RAND_MAX + 1.0);
    mj[i] = 1.0e-5;
  }

  /* CALL GeForce8800GTX */
  //  cunbody1_force_mp(devid, xj, mj, xi, eps2, ai, ni, nj); 
  force_host(devid, xj, mj, xi, eps2, ai, ni, nj); 

  nforce = 10;
  t = get_time();

#pragma omp parallel for
  for(i=0; i< nforce; i++){
    //    cunbody1_force_mp(devid, xj, mj, xi, eps2, ai, ni, nj); 
    force_host(devid, xj, mj, xi, eps2, ai, ni, nj); 
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
