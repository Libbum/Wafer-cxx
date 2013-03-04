TARGET = mpisolve 
OBJ = intde2.o random.o paramreader.o grid.o outputroutines.o initialconditions.o potential.o mpisolve.o mexHatPotential.o

PROFILE = #-pg
DEBUG = -Wno-deprecated #-W -Wall -g
INCLUDES = -I include #-I/usr/local/openmpi/1.6.3-intel/include 
OPTIMIZATION = -funroll-loops -finline-functions -O2
LAPACK = /usr/local/lapack/3.4.2/lib/liblapack.so /usr/local/blas/1.0.248/lib/libblas.so -lm  
LIBS = -L/usr/local/lapack/3.4.2/lib
CXXFLAGS = $(DEBUG) $(PROFILE) $(OPTIMIZATION) $(FLOWTRACE) $(INCLUDES) $(LIBS) 
MPIFLAGS = $(CXXFLAGS)
CC = mpicxx
CXX = g++ #icpc

all: $(TARGET) 

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
	./cleandatafiles.sh

mrproper: clean
	rm -f $(TARGET)

.c.o:
	$(CC) $(MPIFLAGS) -c $*.c

help:
	@echo 'Build targets:'
	@echo '  all          - Builds mpisolve standalone'
	@echo '  solve        - Builds mpisolve standalone'
	@echo 'Cleanup targets:'
	@echo '  clean        - Remove generated files + Emacs leftovers'
	@echo '  mrproper     - Removes generated targets + all aboves'
	@echo 'Exec targets:'
	@echo '  run          - Run mpisolve (with first cleaning up data)'
