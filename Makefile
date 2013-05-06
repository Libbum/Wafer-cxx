TARGET = mpisolve 
OBJ = intde2.o random.o paramreader.o grid.o outputroutines.o initialconditions.o potential.o mpisolve.o mexHatPotential.o

PROFILE = #-pg
DEBUG = -Wno-deprecated #-W -Wall -g
OPTIMIZATION = -funroll-loops -finline-functions -O2
INCLUDES = -I include 
CXXFLAGS = $(DEBUG) $(PROFILE) $(OPTIMIZATION) $(FLOWTRACE) $(INCLUDES) $(LIBS) 
MPIFLAGS = $(CXXFLAGS)

HOST_NAME := $(shell hostname)
ifeq ($(HOST_NAME),trifid)
LAPACK = /usr/local/lapack/3.4.2/lib/liblapack.so /usr/local/blas/1.0.248/lib/libblas.so -lm  
LIBS = -L/usr/local/lapack/3.4.2/lib
CXX = g++ #icpc
CPPFLAGS += -DTRIFID=1
else
INCLUDES += -I$(MKL)/include
LAPACK = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread  
LIBS = -L$(MKL)/lib/intel64 
CPPFLAGS += -DVAYU=1
endif

CC = mpicxx

all: $(TARGET)
	@echo 'Built target for '$(HOST_NAME)

mpisolve:	$(OBJ)
	$(CC) $(OBJ) -o $(TARGET) $(LAPACK) 

run:	
	mpirun -np 3 mpisolve

run1:   
	mpirun -np 2 mpisolve

run4:   
	mpirun -np 5 mpisolve

run8:   
	mpirun -np 9 mpisolve

clean:
	rm -f *\.o *~

mrproper: clean
	rm -f $(TARGET)

.c.o:
	$(CC) $(MPIFLAGS) -c $*.c

help:
	@echo 'Build targets:'
	@echo '  all          - Builds mpisolve standalone'
	@echo '  solve        - Builds mpisolve standalone'
	@echo 'Cleanup targets:'
	@echo '  clean        - Remove generated and temp files'
	@echo '  mrproper     - Removes generated targets + all aboves'
	@echo 'Exec targets:'
	@echo '  run          - Run mpisolve (with first cleaning up data)'
