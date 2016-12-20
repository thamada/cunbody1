/* sample.c
 * by Tsuyoshi Hamada
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
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

void force_host(double xj[][3],
		double m[],
		double xi[][3],
		double eps2,
		double a[][3],
		int ni,
		int nj);

double xj[NMAX][3];
double mj[NMAX];
double xi[NMAX][3];
double vi[NMAX][3];
double ai[NMAX][3];
double ai_host[NMAX][3];

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

void
init_particles(char* fname,
               int* npar,
               double mass[],
               double posi[][3],
               double veloc[][3])
{
  int idx,n;
  char line[2560];
  FILE *fp;
  fp = fopen(fname,"r");
  if(fp == NULL){
    fprintf(stderr, "Abort: there is not such a file :[%s]\n",fname);
    exit(-1);
  }

  idx=0;

  { // get N (number of particles) 
    char* p;
    fgets(line,100,fp);
    p = (char* )strtok(line,"\n");
    sscanf(p,"%d",&n);
    *npar = n;
  }

  // get all data 
  for(idx=0;idx<n;idx++){
    char* p;
    char col[10][400];
    int i=0;
    if(fgets(line,300,fp)==NULL){printf("Error at loading logfiles\n");exit(0);}
    p = (char* )strtok(line,"\t");
    strcpy(col[i],p);
    i++;
    while((p = (char* )strtok(NULL, "\t"))!=NULL){
      strcpy(col[i],p);
      i++;
    }

    mass[idx]=atof(col[0]);
    posi[idx][0]=atof(col[1]);
    posi[idx][1]=atof(col[2]);
    posi[idx][2]=atof(col[3]);
    veloc[idx][0]=atof(col[4]);
    veloc[idx][1]=atof(col[5]);
    veloc[idx][2]=atof(col[6]);
  }
  fclose(fp);
}

/*  main routine */

int main()
{
  int i, dim, ni, nj, nforce;
  double t; 
  double eps2;
  double emax = 0.0;

  //  init_particles("./init.plum.8192", &ni, mj, xi, vi);
  //  init_particles("./init.plum.2048", &ni, mj, xi, vi);
  init_particles("./init.plum.65536", &ni, mj, xi, vi);
  nj = ni;
  /* SETUP I,J-PARTICLES */

  eps2 = 0.01;

  /* CALL GeForce8800GTX */
  cunbody1_force(xi, mj, xi, eps2, ai, ni, nj);
  force_host(xi, mj, xi, eps2, ai_host, ni, nj);


  for(i=0; i<ni; i++){
    int d;
    double err=0.0;
    /*
    for(d=0; d<3; d++){
      double a_gpu = ai[i][d];
      double a_cpu = ai_host[i][d];
      double diff = a_cpu - a_gpu;
      err += diff*diff;
    }
    err = sqrt(err);
    err = err/sqrt(ai_host[i][0]*ai_host[i][0]+ai_host[i][1]*ai_host[i][1]+ai_host[i][2]*ai_host[i][2]);
    */
    {
      double a_gpu = ai[0][d];
      double diff = ai[i][0] - ai_host[i][0];
      err = sqrt((diff*diff)/(ai_host[i][0]*ai_host[i][0]));
    }
    if(emax < err){
      printf("%e\n",err);
      emax = err;
    }
  }
    
  puts("\n\n\n\n");
  puts("----------------------------------------");
  printf("the worst pair-wise error: %e\n", emax);  
  puts("----------------------------------------");

  if(emax < 4.524054e-07){
    printf("success.\n");
  }else{
    printf("failed.\n");
    exit(-1);
  }
  return (0);
}
