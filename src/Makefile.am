DEBUG = -Wno-deprecated #-W -Wall -g
OPTIMIZATION = -funroll-loops -finline-functions -O2
AM_CXXFLAGS = $(DEBUG) $(OPTIMIZATION) $(LIB)

bin_PROGRAMS = wafer 
wafer_SOURCES = intde2.cpp random.cpp paramreader.cpp grid.cpp outputroutines.cpp initialconditions.cpp potential.cpp mexHatPotential.cpp mpisolve.cpp

