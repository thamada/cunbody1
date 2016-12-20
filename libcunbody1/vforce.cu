// Time-stamp: <2009-01-17 15:21:58 hamada>

/*
 * Copyright (C) 2007 
 *      Tsuyoshi Hamada <hamada@progrape.jp>
 *      All rights reserved.
 * This code is released under version 2 of the GNU GPL.
 */

#define IDIM  (4)
#define JDIM  (4)
#define FDIM  (3)
#include "vforce.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cutil.h>

//#undef SM_MAX_BYTE
//#define SM_MAX_BYTE (16384-32)
//#define NJ_SHMEM 256
#define NJ_SHMEM 128  // ** Technic ** 
#define NSP 8
#define NVSP 16
//#define NVSP 24 // GTX260
#define NPIPE (NSP*NVSP)
#define   KIRIAGE(x,y)     (((x) % (y)) ?  ((x/y)+1) : (x/y))
#define   MAX(x,y)     (((x) > (y)) ?  (x) : (y))
#define   MIN(x,y)     (((x) < (y)) ?  (x) : (y))

//#include "cunbody_dbg.h"


//#include "cunbody_kernel.cc"
#include "cunbody_kernel_00.cc"
#include "cunbody_kernel_01.cc"
#include "cunbody_kernel_basic.cc"


namespace libcunbody{
  using namespace std;

  class cunbody1_bench
  {
  private:
    bool is_open;
    int devid;
    char gpu_name[128];
    float4* h_xj;
    float* h_xi;
    float* h_fo;
    float4* d_xj;
    float* d_xi;
    float* d_fo;

    void dev_check(void){
      int ndev;
      CUDA_SAFE_CALL(cudaSetDevice(devid));
      CUDA_SAFE_CALL(cudaGetDeviceCount(&ndev));
      if(ndev == 0){
	fprintf(stdout, "ndev = %d @ %s|%d\n", ndev, __FILE__, __LINE__);
	exit(-1);
	/*
	  }else if(ndev > 1){
	  fprintf(stdout, "ndev = %d @ %s|%d\n", ndev, __FILE__, __LINE__);
	  fprintf(stdout, "This library doesn't work with multiple GPUs.\n");
	  exit(-1);
	*/
      }else{
	int dev = devid;
	cudaDeviceProp deviceProp;
	CUDA_SAFE_CALL(cudaGetDeviceProperties(&deviceProp, dev));
	if (deviceProp.major == 9999 && deviceProp.minor == 9999){
	  printf("There is no device supporting CUDA.\n");
	}
	sprintf(gpu_name, "%s", deviceProp.name);
	printf("  GPU : %s\n",                       deviceProp.name);
	printf("  Major # :  %d\n",                  deviceProp.major);
	printf("  Minor # :  %d\n",                  deviceProp.minor);
	printf("  core clock rate :  %.2f GHz\n",    deviceProp.clockRate * 1e-6f);
#if  (CUDART_VERSION >= 2000)
	printf("  Number of cores :  %d\n",      8 * deviceProp.multiProcessorCount);
	printf("  Number of multiprocessors : %d\n", deviceProp.multiProcessorCount);
#endif
	printf("  global memory : %u bytes\n",    deviceProp.totalGlobalMem);
	printf("  constant memory : %u bytes\n",         deviceProp.totalConstMem); 
	printf("  shared memory per block : %u bytes\n", deviceProp.sharedMemPerBlock);
	printf("  registers available per block : %d\n", deviceProp.regsPerBlock);
	printf("  Warp size : %d\n", deviceProp.warpSize);
	printf("  Max # of threads per block : %d\n",   deviceProp.maxThreadsPerBlock);
	printf("  Max sizes of each dimension of a block : %d x %d x %d\n",
	       deviceProp.maxThreadsDim[0],
	       deviceProp.maxThreadsDim[1],
	       deviceProp.maxThreadsDim[2]);
	printf("  Maximum sizes of each dimension of a grid : %d x %d x %d\n",
	       deviceProp.maxGridSize[0],
	       deviceProp.maxGridSize[1],
	       deviceProp.maxGridSize[2]);
	printf("  Maximum memory pitch : %u bytes\n", deviceProp.memPitch);
	printf("  Texture alignment    : %u bytes\n", deviceProp.textureAlignment);
#if  (CUDART_VERSION >= 2000)
	printf("  Concurrent copy and execution: %s\n", deviceProp.deviceOverlap ? "Yes" : "No");
#endif
	printf("  cudaSetDevice to %d \n", dev);
      }
    }

  public:

    cunbody1_bench() {
      is_open = false;
      devid = 0;
    }

    ~cunbody1_bench() {
      this->close();
      devid = 0;
      is_open = false;
    }

    /* no need to call if you wan't*/
    void close(void)
    {
      // cleanup memory
      CUDA_SAFE_CALL(cudaFreeHost(h_xj));
      CUDA_SAFE_CALL(cudaFreeHost(h_xi));
      CUDA_SAFE_CALL(cudaFreeHost(h_fo));
      CUDA_SAFE_CALL(cudaFree(d_xi));
      CUDA_SAFE_CALL(cudaFree(d_fo));
      CUDA_SAFE_CALL(cudaFree(d_xj));
    }

    void set_devid(int id)
    {
      devid = id;
    }

    void force(double xj[][3], double mj[], double xi[][3], double eps2, double a[][3], int ni, int nj)
    {
      int nj1 = ((nj+NJ_SHMEM-1)/NJ_SHMEM)*NJ_SHMEM;

      if( ni > (0x1<<17) ){
	printf("ERROR %s|%d\n",__FILE__, __LINE__);
	printf(" ni > 131072 : ni=%d\n",ni);
	exit(-1);
      }

      if( nj > (0x1<<17) ){
	printf("ERROR %s|%d\n",__FILE__, __LINE__);
	printf(" nj > 131072 : nj=%d\n",nj);
	exit(-1);
      }

      unsigned int ip_size = sizeof(float) * ni * IDIM;
      unsigned int jp_size = sizeof(float4) * nj1;
      unsigned int fo_size = sizeof(float) * ni * FDIM;

      if(is_open == false){
	//    CUT_DEVICE_INIT();
	//    CUT_CHECK_DEVICE();
	unsigned int _nmax = 1<<17; // 131072
	//    unsigned int _nmax = 1<<14; // 16384
	unsigned int _ip_size = sizeof( float) * _nmax * IDIM;
	unsigned int _jp_size = sizeof(float4) * _nmax;
	unsigned int _fo_size = sizeof(float) * _nmax * FDIM;

	dev_check();

	CUDA_SAFE_CALL(  cudaMallocHost( (void**)&h_xj, _jp_size)  );
	CUDA_SAFE_CALL(  cudaMallocHost( (void**)&h_xi, _ip_size)  );
	CUDA_SAFE_CALL(  cudaMallocHost( (void**)&h_fo, _fo_size)  );

	CUDA_SAFE_CALL( cudaMalloc( (void**) &d_xj, _jp_size));
	CUDA_SAFE_CALL( cudaMalloc( (void**) &d_xi, _ip_size));
	CUDA_SAFE_CALL( cudaMalloc( (void**) &d_fo, _fo_size));

	for(int i = 0; i < _nmax ; i++) h_xj[i] = make_float4(0.0, 0.0, 0.0, 0.0);
	CUDA_SAFE_CALL( cudaMemcpy( d_xj, h_xj, _jp_size, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy( d_xi, h_xj, _ip_size, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy( d_fo, h_xj, _fo_size, cudaMemcpyHostToDevice) );
	fprintf(stderr, "open %s by CUNBODY-1 library: rev.hamada20080905  (^<_^)/ %d\n", gpu_name, devid);
	is_open = true;
      }
      // ------------------------------------------------------------ xj
      for(int i = 0; i < nj; i++){
	h_xj[i].x = (float) xj[i][0];
	h_xj[i].y = (float) xj[i][1];
	h_xj[i].z = (float) xj[i][2];
	h_xj[i].w = (float) mj[i];
      }

      if(nj < nj1){
	for(int i = nj; i < nj1; i++)  h_xj[i] = make_float4(0.0, 0.0, 0.0, 0.0);
      }

      CUDA_SAFE_CALL( cudaMemcpy( d_xj, h_xj, jp_size, cudaMemcpyHostToDevice) );

      // ------------------------------------------------------------ xi
      for(int i = 0; i < ni; i++){
	h_xi[i     ] = (float) xi[i][0];
	h_xi[i+ni  ] = (float) xi[i][1];
	h_xi[i+ni*2] = (float) xi[i][2];
	h_xi[i+ni*3] = (float) eps2;
      }
      CUDA_SAFE_CALL( cudaMemcpy( d_xi, h_xi, ip_size, cudaMemcpyHostToDevice) );


      dim3 grid((ni+NPIPE-1)/NPIPE); // ** Technic ** 
      dim3 threads  (NPIPE);
      {
	// type00:Id, N, sec, Gflop/s{max, avg, curr}: 202        131072  4.01055 651.291 651.113 651.118
	// type01:Id, N, sec, Gflop/s{max, avg, curr}: 202        131072  4.01369 650.777 650.627 650.608
	// basic :Id, N, sec, Gflop/s{max, avg, curr}: 202        131072  4.78496 545.754 545.738 545.739
	// nvidia:Id, N, sec, Gflop/s{max, avg, curr}: 202        131072  4.5376  575.526 575.504 575.489

	using namespace cunbody_kernel_type00;
	//	using namespace cunbody_kernel_basic;
	kernel<<< grid, threads >>>(d_xj, d_xi, d_fo, ni, nj1);

	//    using namespace cunbody_kernel_nvidia;
	//    kernel<<< grid, threads, (NJ_SHMEM*sizeof(float4)) >>>(d_xj, d_xi, d_fo, ni, nj1);

	CUT_CHECK_ERROR("Kernel execution failed");
      }

      // ------------------------------------------------------------ fo
      CUDA_SAFE_CALL( cudaMemcpy( h_fo, d_fo, fo_size, cudaMemcpyDeviceToHost) );

      for(int i=0;i<ni; i++){
	a[i][0] = (double)h_fo[i];
	a[i][1] = (double)h_fo[i+ni];
	a[i][2] = (double)h_fo[i+ni*2];
      }
    }


  }; // class cunbody1 __END__
}; // namespace libcunbody __END__  ----------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------------
extern "C" void copyright_cunbody1(void);

void copyright_cunbody1(void)
{
  printf("Copyright(C) 2007 by Tsuyoshi Hamada <hamada@progrape.jp>, All rights reserved.\n");
}


static double T_Force = 0.0;
static double T_Write = 0.0;
static double T_Calc  = 0.0;
static double T_Read  = 0.0;
static unsigned ni_sum = 0;
static unsigned nj_sum = 0;
static double n_inter = 0.0;


extern "C" void cunbody_dumptime(void);

void
_cunbody_dumptime(void)
{
  printf("----------------------------------\n");
  printf("Time: %f |\t W(%f)\t C(%f)\t R(%f)\n",T_Force, T_Write, T_Calc, T_Read);
  printf("MB/s: %f \t %f\n",
	 (1.0e-6)*(sizeof(float)*ni_sum*IDIM+sizeof(float)*nj_sum*JDIM)/T_Write,
	 (1.0e-6)*(sizeof(float)*ni_sum*FDIM)/T_Read);

  printf("Gflop/s: %f \n",(1.0e-9)*n_inter*38.0/T_Calc);
  printf("Gflop/s: %f \n",(1.0e-9)*n_inter*38.0/T_Force);
  printf("n inter %e \n", n_inter);
  printf("\n");

  T_Force = T_Write = T_Calc = T_Read = 0.0;
  ni_sum = 0;
  nj_sum = 0;
  n_inter = 0.0;
}

void
cunbody_dumptime(void){ }

//-------------------------------------------------------------------

#define MAX_OMP_THRE (8)
static libcunbody::cunbody1_bench cunObj[MAX_OMP_THRE];


extern "C" void cunbody1_force_mp(int devid, double xj[][3], double mj[], double xi[][3], double eps2, double a[][3], int ni, int nj)
{
  using namespace std;
  using namespace libcunbody;
  cunObj[devid].set_devid(devid);
  cunObj[devid].force(xj, mj, xi, eps2, a, ni, nj);
}

extern "C" void cunbody1_force(double xj[][3], double mj[], double xi[][3], double eps2, double a[][3], int ni, int nj)
{
  using namespace std;
  using namespace libcunbody;
  int devid = 0;
  cunObj[devid].set_devid(devid);
  cunObj[devid].force(xj, mj, xi, eps2, a, ni, nj);
}



