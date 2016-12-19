//Time-stamp: <2007-02-19 11:48:31 hamada>
//Copyright(c) 2000-2006 by Tsuyoshi Hamada. All rights reserved.

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

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

void
read_result(char* fname,      // filename 
	    int* npar,        // number of particles
	    double acc[][3],
	    double pot[])
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

    acc[idx][0]=atof(col[0]);
    acc[idx][1]=atof(col[1]);
    acc[idx][2]=atof(col[2]);
    pot[idx]   =atof(col[3]);
  }
  fclose(fp);
}
