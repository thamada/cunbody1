#Time-stamp: <2009-01-16 22:31:09 hamada>

ifndef SDK_INSTALL_PATH
SDK_INSTALL_PATH = /home/hamada/NVIDIA_CUDA_SDK
endif

ifndef CUDA_INSTALL_PATH
CUDA_INSTALL_PATH = /usr/local/cuda
endif


LIBCUNBODY = ../libcunbody1/libcunbody1.a

all:
	gcc ./main.c force_host.c $(LIBCUNBODY)  -L $(CUDA_INSTALL_PATH)/lib -L $(SDK_INSTALL_PATH)/lib -lcuda -lcudart -lGL -lGLU -lcutil -lm


clean:
	rm -rf *~ .*~ ./a.out

