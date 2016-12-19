# ----------------------------------- for Intel CC
CC = icc
CFLAGS	= -Wall -static -O3
#CFLAGS	= -Wall -static -O3 -xW -Ob2 -tpp6 -rcd 
#CFLAGS	= -Wall -O4
#CFLAGS	= -Wall -O2

# ----------------------------------- for GNU CC
CC = gcc
#CFLAGS	= -Wall -g
CFLAGS	= -Wall -O3 -pipe



OBJS = main.c  init_particles.c  leapflog.c  debug_position_snap.c libcunbody1.a



CUDA_SDK_PATH = /usr/local/src/NVIDIA_CUDA_SDK
OPT_CUDA = -L/usr/local/cuda/lib -L$(CUDA_SDK_PATH)/lib -lcuda -lcudart -lGL -lGLU  -lcutil 



all : $(OBJS) 
	$(CC) $(OBJS) -o run.x  -lm $(OPT_CUDA)

rmobj : 
	rm -f *.o

c : clean

clean: 
	rm -rf *.o *~ .*~ run.x

t : test

test :
	./run.x ~/Dfile/PL/init.plummer.4096

tt:
	./run.x ./init.plum.65536 ./result.65536

ttt:
	./run.x ~/Dfile/PL/init.plum.1111 ~/Dfile/PL/result.1111

ttx:
	./run.x ~/Dfile/PL/init.plum.2048 ~/Dfile/PL/result.2048
	./run.x ~/Dfile/PL/init.plum.4096 ~/Dfile/PL/result.4096 
	./run.x ~/Dfile/PL/init.plum.8192 ~/Dfile/PL/result.8192
	./run.x ./init.plum.131072 ./result.131072

