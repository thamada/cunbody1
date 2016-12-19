//Time-stamp: <2007-02-19 11:44:53 hamada>
//Copyright(c) 2006 by Tsuyoshi Hamada. All rights reserved.

#include <stdio.h>

void leapflog(double dt,double x[][3],double v[][3],int n)
{
    int i,d;
    for(i=0;i<n;i++)
        for(d=0;d<3;d++)
            x[i][d] = x[i][d]+v[i][d]*dt;
}

void leapflog_half(double dt, double v[][3], double a[][3], int n)
{
    int i,d;
    for(i=0;i<n;i++)
      for(d=0;d<3;d++)
	v[i][d] = v[i][d] + (a[i][d]*dt)/2.0;
}

