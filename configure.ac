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

#Compile options
#AC_ARG_ENABLE(debug, 
#    [AC_HELP_STRING([--debug-level,Compile with debugging: OFF (default); ON or FULL])],
#    [case "${enableval}" in
         
#THis is auctaully todo. We want to pass through this value to the main file
#enable_debug=$enableval, enable_debug=no)
#if test "$enable_debug" = "yes"; then
#    CFLAGS="-g"
#    CXXFLAGS="-g"
#    FFLAGS="-g"
#    AC_DEFINE(DEBUG,1,[define to enable debugging code])
#fi

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

AC_ARG_ENABLE([blas],
    AS_HELP_STRING([--enable-blas], [Include BLAS library (using an automatic search)]))

AS_IF([test "x$enable_blas" = "xyes"], [
    AC_MSG_NOTICE(BLAS support requested. Configuring...)
    AX_BLAS([AC_MSG_RESULT()], [AC_MSG_ERROR([BLAS library not found])])
])

AC_ARG_ENABLE([lapack],
    AS_HELP_STRING([--enable-lapack], [Include LAPACK library (using an automatic search)]))

AS_IF([test "x$enable_lapack" = "xyes"], [
    AC_MSG_NOTICE(LAPACK support requested. Configuring...)
    AX_LAPACK([AC_MSG_RESULT()], [AC_MSG_ERROR([LAPACK library not found])])
])

# Checks for MPI (required)
AC_MSG_NOTICE(Searching for MPI capability...)
AX_MPI([AC_MSG_RESULT()], [AC_MSG_ERROR([MPI not found])])

CXX="$MPICXX"
LIBS="$MPILIBS $LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
