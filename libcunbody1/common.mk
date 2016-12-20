#Time-stamp: <2009-01-17 15:21:14 hamada>

LIBNAME = libcunbody1.a
TARGET = vforce

ifndef SDK_INSTALL_PATH
SDK_INSTALL_PATH = /home/hamada/NVIDIA_CUDA_SDK
endif

ifndef CUDA_INSTALL_PATH
CUDA_INSTALL_PATH = /usr/local/cuda
endif


all : lib.gcc

c : clean

clean : 
	rm -rf *~ .*~
	rm -rf $(TARGET).{cubin,ptx}
	rm -f $(TARGET).o $(LIBNAME)

# gcc
lib.gcc : 
	nvcc --host-compilation 'C++' -DUNIX -O3 \
	  -Xcompiler "$(ARCHFLAG)" \
          -Xcompiler "-fPIC " \
          -Xcompiler "-O3" \
          -Xcompiler "-ffast-math" \
          -I$(SDK_INSTALL_PATH)/common/inc -I$(CUDA_INSTALL_PATH)/include -I. \
	  -c $(TARGET).cu -o ./$(TARGET).o 
	ar ruv $(LIBNAME) $(TARGET).o
	ranlib $(LIBNAME)


cubin:
	nvcc  -cubin $(TARGET).cu  -I. -I$(CUDA_INSTALL_PATH)/include -I$(SDK_INSTALL_PATH)/common/inc
	grep smem $(TARGET).cubin
	grep reg $(TARGET).cubin

ptx:
	nvcc  --ptx $(TARGET).cu  -I. -I$(CUDA_INSTALL_PATH)/include -I$(SDK_INSTALL_PATH)/common/inc

ins: install

install : lib
	cp $(LIBNAME) /usr/local/lib/

export:
	tar cvfps ../lib2.`date +"%Y%m%d-%H%M%S"`.tar ../lib2


b: export
