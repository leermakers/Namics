SHELL = /bin/sh
# CXX= icpc -c 
# We recommend using "make -D DEBUG" and "make -D LOCALFENK" instead of following options
# Compile sfbox with additional debug code in "#if DEBUG" blocks
# DEBUG = 1               
# Use ~/lib and ~/include instead of normal fenklib locations
# LOCALFENK = 1  

#CXX= icpc -c # if this line is enabled the intel c++ compiler is used instead of Gcc.

LAPACK = 1 

NVCC := nvcc 

DEBUG_STR =              
ifdef DEBUG
    DEBUG_STR = -DDEBUG
endif

	INC_DIR = /usr/local/include
	LIB_DIR = /usr/local/lib
#	MATRIXLIB = -llapack -lblas -lf2c -lm -lpthread
	MATRIXLIB = -lm -lpthread


CPPFLAGS = -I$(INC_DIR) $(DEBUG_STR) -Wall
CPPFLAGS+= -O3

	CPPFLAGS = -I$(INC_DIR) $(DEBUG_STR)  -O3 -Wall -DHAS_LAPACK -DNO_BLAS_WRAP 
#	CPPFLAGS = -I$(INC_DIR) $(DEBUG_STR) -g -O3 -Wall -DHAS_LAPACK -DNO_BLAS_WRAP
CFLAGS = $(CPPFLAGS)

	LIBS +=  $(MATRIXLIB)

ifdef CUDA
	LIBS += -L/usr/include -lcuda -lcudart -lcublas
    
    CPPFLAGS += -DCUDA
    NVCCFLAGS += -DCUDA -arch sm_13  
endif
    
%.o: %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@
   
OBJECTS = namics.o tools.o input.o lattice.o segment.o alias.o molecule.o newton.o engine.o output.o variate.o system.o

#namics : namics.o
#	g++  -o namics namics.o -L$(LIB_DIR) $(LIBS)


namics : $(OBJECTS)
	g++ -o namics $(OBJECTS) -L$(LIB_DIR) $(LIBS)
#	g++ -o namics $(OBJECTS) -L$(LIB_DIR) $(LIBS)
clean : 
	rm *.o
