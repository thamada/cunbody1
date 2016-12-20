#Time-stamp: <2009-01-24 03:46:57 hamada>

ifndef SDK_INSTALL_PATH
SDK_INSTALL_PATH = /home/hamada/NVIDIA_CUDA_SDK
endif

ifndef CUDA_INSTALL_PATH
CUDA_INSTALL_PATH = /usr/local/cuda
endif


LIBCUNBODY = ../libcunbody1/libcunbody1.a

all:
	gcc ./sample.c $(LIBCUNBODY)  -L $(CUDA_INSTALL_PATH)/lib -L $(SDK_INSTALL_PATH)/lib -lcuda -lcudart -lGL -lGLU -lcutil -fopenmp -lm


clean:
	rm -rf *~ .*~ ./a.out


test_1stGPU:
	./a.out 1

test_2ndGPU:
	./a.out 2


