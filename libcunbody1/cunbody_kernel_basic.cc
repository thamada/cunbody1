//Time-stamp: <2008-07-07 20:16:07 hamada>
/*
 * Copyright (C) 2007 
 *      Tsuyoshi Hamada <hamada@progrape.jp>
 *      All rights reserved.
 * This code is released under version 2 of the GNU GPL.
 */

namespace cunbody_kernel_basic{

  __device__ 
  float4 inter(float4 xj,
	       float4 xi,
	       float4 apot)
  {
    float dx,dy,dz;
    float r2,r1i,r2i,r3i;
    float mr3i;
    float mj    = xj.w;
    float ieps2 = xi.w;
    dx = xj.x - xi.x;
    dy = xj.y - xi.y;
    dz = xj.z - xi.z;
    r2 = (dx*dx+ieps2)+ (dy * dy) + (dz * dz); // ** Technic ** 
    r1i = 1/sqrt(r2);
    r2i = r1i * r1i;
    r3i = r2i * r1i;
    mr3i = mj * r3i;
    apot.x += dx * mr3i;
    apot.y += dy * mr3i;
    apot.z += dz * mr3i;
    return (apot);
  }


  // Extended Chamomile-Scheme
  __global__ void
  kernel(float4* g_xj,
	 float* g_xi,
	 float* g_fi,
	 int ni,
	 int nj)
  {
    int tid    = threadIdx.x;
    int bid    = blockIdx.x;

    int n_jblock = (nj+NJ_SHMEM-1)/NJ_SHMEM;

    int i = bid*NPIPE+tid;;
    float4 x_i;
    float4 OUTPUT = make_float4(0.0, 0.0, 0.0, 0.0);

    x_i.x = g_xi[i];
    x_i.y = g_xi[i+ni];
    x_i.z = g_xi[i+ni*2];
    x_i.w = g_xi[i+ni*3];

    __shared__ float4 s_xj[NJ_SHMEM];

    for(int jloop = 0; jloop<n_jblock; jloop++){
      int j_start = NJ_SHMEM * jloop;

      __syncthreads();
      for(int j = 0; j< (NJ_SHMEM/NPIPE); j++){
	s_xj[tid+j*NPIPE] = g_xj[j_start+tid];
      }
      __syncthreads();

      for(int j = 0; j<NJ_SHMEM; j++){
	OUTPUT = inter(s_xj[j], x_i, OUTPUT);
      }
    }

    g_fi[i]      = OUTPUT.x;
    g_fi[i+ni]   = OUTPUT.y;
    g_fi[i+ni*2] = OUTPUT.z;
  }


};
