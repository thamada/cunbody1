#Time-stamp: <2009-01-17 03:00:20 hamada>

ifndef SDK_INSTALL_PATH
SDK_INSTALL_PATH = /home/hamada/NVIDIA_CUDA_SDK
endif

ifndef CUDA_INSTALL_PATH
CUDA_INSTALL_PATH = /usr/local/cuda
endif


LIBCUNBODY = ../libcunbody1/libcunbody1.a

all:
	g++ ./main.cpp $(LIBCUNBODY)  -L $(CUDA_INSTALL_PATH)/lib -L $(SDK_INSTALL_PATH)/lib -lcuda -lcudart -lGL -lGLU -lcutil -fopenmp


clean:
	rm -rf *~ .*~ ./a.out


2gpu:
	export OMP_NUM_THREADS=2


