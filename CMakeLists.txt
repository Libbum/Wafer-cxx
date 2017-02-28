cmake_minimum_required (VERSION 3.1 FATAL_ERROR)

project (Wafer)
set (Wafer_VERSION_MAJOR 0)
set (Wafer_VERSION_MINOR 9)
set (Wafer_VERSION_REVISION 0)

option (USE_LAPACK "Use LAPACK and BLAS" OFF)
option (USE_MKL "Use MKL Library" OFF)

set (CMAKE_MODULE_PATH "${CMAKE_BINARY_DIR}/Modules/" CACHE PATH "CMake module path" FORCE)

find_package (MPI REQUIRED)

if (USE_MKL)
   find_package (MKL)
   if (MKL_FOUND)
      #If we use the intel compilers, the CMAKE_CXX_STANDARD variable isn't read. Also we add a few more optimisations.
      if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
         #Default IPO, drop superflous warnings
         set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHost -parallel -ansi-alias -DMKL_ILP64 -I${MKL_ROOT}/include -wd11074 -wd11076")
      else()
         #Assuming GNU
         set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -ftree-vectorize -msse2 -ftree-vectorizer-verbose=5 -DMKL_ILP64 -m64 -I${MKL_ROOT}/include")
      endif()
      add_definitions (-DUSINGMKL)
   else()
      message ("-- MKL Library Not Found. Falling back to LAPACK.")
      find_package (LAPACK REQUIRED)
      add_definitions (-DUSINGLAPACK)
   endif()
else()
   if (USE_LAPACK)
      find_package (LAPACK REQUIRED)
      add_definitions (-DUSINGLAPACK)
   endif()
endif()

set (CMAKE_BUILD_TYPE Release)

include (GetGitRevisionDescription)
get_git_head_revision (GIT_REFSPEC GIT_SHA1)
configure_file (include/gitsha1.hpp.in ${CMAKE_BINARY_DIR}/generated/gitsha1.hpp @ONLY)
include_directories (${CMAKE_BINARY_DIR}/generated)

include_directories (include "${MPI_INCLUDE_PATH}" ${CMAKE_BINARY_DIR}/generated)

add_executable (wafer intde2.cpp random.cpp paramreader.cpp grid.cpp outputroutines.cpp initialconditions.cpp potential.cpp mexHatPotential.cpp wafer.cpp)
if (MKL_FOUND)
   if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
      target_link_libraries (wafer "${MPI_LIBRARIES}" -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl)
   else()
      target_link_libraries (wafer "${MPI_LIBRARIES}" -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl)
   endif()
else()
   if (LAPACK_FOUND)
      target_link_libraries (wafer "${MPI_LIBRARIES}" lapack blas)
   else()
      target_link_libraries (wafer "${MPI_LIBRARIES}")
   endif()
endif()

set_target_properties (wafer PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")

set (CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})
install (TARGETS wafer DESTINATION bin)
