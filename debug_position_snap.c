//Time-stamp: <2007-09-11 17:32:00 hamada>
#include <stdio.h>


void debug_position_snap(double x[][3], double Gflops, int n, int nstep)
{
  //    int i;
  static int flag=0;
  static FILE* fp;

  if(flag==0){
    /*
      {
      char fname[256];
      static int i = 0;
      if(i==0) sprintf(fname,"/mnt/ram/xxx.log%d", time(NULL));
      fp = fopen(fname,"w");
      i=1;
      }
    */
    //      fp = fopen("/mnt/ram/xxx.log","w");
    //      fp = fopen("xxx.log","w");
    //    fp = fopen("/mnt/ram/xxx.log","w");
    fp = fopen("/dev/shm/xxx.log","w");
    flag = 1;
  }else{
    rewind(fp);
  }

  fprintf(fp,"%d\n",nstep);
  fprintf(fp,"%e\n",Gflops);
  fprintf(fp,"%d\n",n);

  fwrite(x, n*24,1,fp);

  /*
    for(i=0;i<n;i++){
    fprintf(fp,"%1.8e\t%1.8e\t%1.8e\n",x[i][0],x[i][1],x[i][2]);
    }
  */

  //    fclose(fp);
  fflush(fp);
}

void debug_position_dump(double x[][3], int n, int nstep)
{
  int i;
  static int flag=0;
  static FILE* fp;

  if(flag==0){
    fp = fopen("/mnt/ram/dump.log","w");
    flag = 1;
  }else{
    rewind(fp);
  }

  fprintf(fp,"%d\n",nstep);
  fprintf(fp,"%d\n",n);

  for(i=0;i<n;i++){
    fprintf(fp,"%1.8e\t%1.8e\t%1.8e\n",x[i][0],x[i][1],x[i][2]);
  }

  //    fclose(fp);
  fflush(fp);
}

