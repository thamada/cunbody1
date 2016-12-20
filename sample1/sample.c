/* sample.c
 * by Tsuyoshi Hamada
 */
#include <stdio.h>
#include <stdlib.h>
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

void cunbody1_force(double xj[][3], // position of j-th particles
                    double mj[],    // mass of j-th particles
                    double xi[][3], // position of i-th particles
                    double eps2,    // softening parameter
                    double ai[][3], // force of i-th particles
                    int ni,         // number of i-th particles
                    int nj);        // number of j-th particles


/*  main routine */

int main()
{
  int i, dim, ni, nj;
  double eps2 = 0.001;

  srand(0x19740526);

  /* SETUP I,J-PARTICLES */
  ni = nj = 65536;
  for(i=0; i<ni; i++){
    for(dim=0; dim<3; dim++)
      xi[i][dim] = xj[i][dim] = rand()/(RAND_MAX + 1.0);
    mj[i] = 1.0e-5;
  }

  /* CALL GeForce8800GTX */
  cunbody1_force(xj, mj, xi, eps2, ai, ni, nj); 

  /* OUTPUT RESULTS */
  for(i=ni-10; i<ni; i++)
    for(dim=0; dim<3; dim++)
      printf("ai[%d][%d] = %e\n",i, dim, ai[i][dim]);

  return (0);
}
