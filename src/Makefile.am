DEBUG = -Wno-deprecated #-W -Wall -g
OPTIMIZATION = -funroll-loops -finline-functions -O2
INCLUDES = -I include 
AM_CFLAGS = $(DEBUG) $(OPTIMIZATION) $(INCLUDES)

bin_PROGRAMS = wafer
wafer_SOURCES = intde2.c random.c paramreader.c grid.c outputroutines.c initialconditions.c potential.c mpisolve.c mexHatPotential.c
