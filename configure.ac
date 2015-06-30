#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([wafer], [1.0], [https://github.com/Libbum/Wafer/issues])
AC_CONFIG_SRCDIR([src/mpisolve.cpp])
AM_INIT_AUTOMAKE
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])
AM_MAINTAINER_MODE([disable])

# Checks for programs.
#AC_PROG_CXX
#AC_PROG_CC
#AC_PROG_CXXCPP
AX_PROG_CXX
AC_LANG([C++])

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h sys/time.h unistd.h])
AC_CHECK_HEADERS(iostream)

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([pow sqrt])

# Checks for BLAS/LAPACK libraries:

AX_BLAS([AC_MSG_RESULT()], [AC_MSG_WARN([BLAS library not found])])
AX_LAPACK([AC_MSG_RESULT()], [AC_MSG_WARN([LAPACK library not found])])


# Checks for MPI
AX_MPI([AC_MSG_RESULT()], [AC_MSG_ERROR([MPI not found])])

CXX="$MPICXX"
LIBS="$MPILIBS $LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
