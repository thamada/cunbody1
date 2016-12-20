#include <stdio.h>
#include <math.h>

void force_host(double xj[][3],
		double m[],
		double xi[][3],
		double eps2,
		double a[][3],
		int ni,
		int nj)
{
  int i,j,d;
  double dx[3];
  for(i=0;i<ni;i++) for(d=0;d<3;d++) a[i][d] = 0.0;

  for(i=0;i<ni;i++){
    for(j=0;j<nj;j++){
      double r2,r3;
      r2 = eps2;
      for(d=0;d<3;d++){
	dx[d] = xj[j][d] - xi[i][d];
	r2 += dx[d] * dx[d];
      }
      r3 = sqrt(r2)*r2;
      for(d=0;d<3;d++){
	a[i][d] +=  m[j]*dx[d]/r3;
      }
    }
  }
}

void force_pot_host(double xj[][3],
		    double m[],
		    double xi[][3],
		    double eps2,
		    double a[][3],
		    double pot[],
		    int ni,
		    int nj)
{
  int i,j,d;
  double dx[3];
  for(i=0;i<ni;i++){
    for(d=0;d<3;d++) a[i][d] = 0.0;
    pot[i] = 0.0;
  }
  
  for(i=0;i<ni;i++){
    for(j=0;j<nj;j++){
      double r2,r1,r3;
      r2 = eps2;
      for(d=0;d<3;d++){
	dx[d] = xj[j][d] - xi[i][d];
	r2 += dx[d] * dx[d];
      }
      r1 = sqrt(r2); 
      r3 = r1*r2;
      for(d=0;d<3;d++){
	a[i][d] +=  m[j]*dx[d]/r3;
      }
      pot[i] -= m[j]/r1;
    }
  }
}

void
force_pot_jerk_host(double x[][3], double v[][3], double m[], double p[], double a[][3], double jk[][3], int n)
{
  int    i, j;
  //  const double eps2 = (1.0/256.0)*(1.0/256.0);
  const double eps2 = 0.01;
  double pij;
  double sx_0, sx_1, sx_2; 
  double sjk_0, sjk_1, sjk_2;
  double xj_0, xj_1, xj_2;
  double xi_0, xi_1, xi_2;
  double dx2_0, dx2_1, dx2_2;
  double x2y2, r1, r2, r3, r5, r5i;
  double mj;
  double mf, fx_0, fx_1, fx_2;
  double vj_0, vj_1, vj_2;
  double vi_0, vi_1, vi_2;
  double dv_0, dv_1, dv_2;
  double xv_0, xv_1, xv_2;
  double xv1, xv2, xv2x3;
  double mr5i;
  double jk2a;
  double dx_0, dx_1, dx_2;
  double jk0_0, jk0_1, jk0_2;
  double jk1_0, jk1_1, jk1_2;
  double jk2_0, jk2_1, jk2_2;


  for (i = 0; i < n; i++) {
    xi_0 = x[i][0];
    xi_1 = x[i][1];
    xi_2 = x[i][2];
    vi_0 = v[i][0];
    vi_1 = v[i][1];
    vi_2 = v[i][2];

    pij = 0.0;
    sx_0 = 0.0;
    sx_1 = 0.0;
    sx_2 = 0.0;
    sjk_0 = 0.0;
    sjk_1 = 0.0;
    sjk_2 = 0.0;

    for (j = n - 1; j >= 0; j--) {
 
      xj_0 = x[j][0];
      xj_1 = x[j][1];
      xj_2 = x[j][2];
      vj_0 = v[j][0];
      vj_1 = v[j][1];
      vj_2 = v[j][2];
      mj = m[j];

      dx_0 = xj_0 - xi_0;
      dx_1 = xj_1 - xi_1;
      dx_2 = xj_2 - xi_2;
      dx2_0 = dx_0 * dx_0;
      dx2_1 = dx_1 * dx_1;
      dx2_2 = dx_2 * dx_2;

      x2y2 = dx2_0 + dx2_1;
      r2   = dx2_2 + x2y2 + eps2;
      r1   = sqrt(r2);
      r3   = r1 * r2;
      r5   = r2 * r3;
      r5i  = 1.0/r5;
      //      if(r5 == 0.0) r5i = 0.0;

      mr5i = mj / r5;
      mf   = mj / r3;        // mj / r3 
      pij += mj / r1;        // mj / r1  : Potential

      fx_0 = mf * dx_0;
      fx_1 = mf * dx_1;
      fx_2 = mf * dx_2;

      sx_0 += fx_0;          // : Force.x
      sx_1 += fx_1;          // : Force.y
      sx_2 += fx_2;          // : Force.z
      
      dv_0   = vj_0 - vi_0;
      dv_1   = vj_1 - vi_1;
      dv_2   = vj_2 - vi_2;
      jk1_0  = mf * dv_0;    // mj * Vij.x / r3
      jk1_1  = mf * dv_1;    // mj * Vij.y / r3
      jk1_2  = mf * dv_2;    // mj * Vij.z / r3
      xv_0   = dx_0 * dv_0;
      xv_1   = dx_1 * dv_1;
      xv_2   = dx_2 * dv_2;
      xv1    = xv_0 + xv_1;
      xv2    = xv_2 + xv1;  // Vij * Rij := a
      xv2x3  = 3.0 * xv2;   // 3a
      jk2a   = mr5i * xv2x3;
      jk2_0  = jk2a * dx_0; // 3a * Fij
      jk2_1  = jk2a * dx_1; // 3a * Fij
      jk2_2  = jk2a * dx_2; // 3a * Fij
      jk0_0  = jk1_0 - jk2_0; // Jk.x
      jk0_1  = jk1_1 - jk2_1; // Jk.y
      jk0_2  = jk1_2 - jk2_2; // Jk.z
      sjk_0 += jk0_0;
      sjk_1 += jk0_1;
      sjk_2 += jk0_2;
    }
    p[i] = -pij;
    a[i][0] = sx_0;
    a[i][1] = sx_1;
    a[i][2] = sx_2;
    jk[i][0] = sjk_0;
    jk[i][1] = sjk_1;
    jk[i][2] = sjk_2;
  }
}


double energy_pot(double pot[], int ni)
{
  int i;
  double e = 0.0 ;
  for(i=0;i<ni;i++){
    e += pot[i];
  }
  return (e);
}
