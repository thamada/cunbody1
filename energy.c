//Time-stamp: <2007-02-19 14:02:56 hamada>
#include <stdio.h>

double energy(double m[], double x[][3], double v[][3], double eps2 , int n)
{
    int i,j,d;
    double v2,r2;
    double dx[3];
    double Ep;    /* Potential Energy */
    double Ek;    /* Kinetic Energy   */
    double Ea;    /* Total Energy     */

    /* Compute Potential Energy */ 
    Ep=0.0;
    for(i=0;i<n-1;i++){
        for(j=i+1;j<n;j++){
            r2 = eps2;
            for(d=0;d<3;d++){
                dx[d] = x[j][d] - x[i][d];
                r2 += dx[d] * dx[d];
            }
            Ep -= m[i] * m[j] * ( 1.0/sqrt(r2)  );
        }
    }


    /* Compute Kinetic Energy */ 
    Ek=0.0;
    for(i=0;i<n;i++)
        for(d=0;d<3;d++)
            Ek += 0.5 * m[i] * v[i][d] * v[i][d];

    /* Compute Total Energy */ 
    Ea=0.0;
    Ea = Ek + Ep;

    return Ea;
}

double energy_with_pot(double m[], double x[][3], double v[][3], double p[], double eps2 , int n)
{
    int i,j,d;
    double v2,r2;
    double dx[3];
    double Ep;    /* Potential Energy */
    double Ek;    /* Kinetic Energy   */
    double Ea;    /* Total Energy     */

    /* Compute Potential/Kinetic Energy */ 
    Ep=0.0;
    Ek=0.0;
    for(i=0;i<n;i++){
      Ep += 0.5 * m[i] * p[i];
      for(d=0;d<3;d++) Ek += 0.5 * m[i] * v[i][d] * v[i][d];
    }

    /* Compute Total Energy */ 
    Ea=0.0;
    Ea = Ek + Ep;

    return Ea;
}


double get_err_3D(double host[3], double device[3])
{
  double err = 0.0;
  for(int d=0;d<3;d++){
    double diff = fabs(host[d] - device[d]);
    err += (diff*diff);
  }
  err = sqrt(err);
  err = err / sqrt(host[0]*host[0] + host[1]*host[1] + host[2]*host[2]);

  return (err);
}

double get_err_1D(double host, double device)
{
  double err = 0.0;
  double diff = fabs(host - device);

  err = sqrt((diff*diff)) / sqrt(host*host);

  return (err);
}

