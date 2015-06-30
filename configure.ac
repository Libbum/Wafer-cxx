#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([wafer], [1.0], [https://github.com/Libbum/Wafer/issues])
AC_CONFIG_SRCDIR([src/mpisolve.c])
AM_INIT_AUTOMAKE
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])


# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([gethostname pow sqrt])

# Checks for BLAS/LAPACK libraries:

AX_BLAS([AC_MSG_RESULT()], [AC_MSG_WARN([BLAS library not found])])
AX_LAPACK([AC_MSG_RESULT()], [AC_MSG_WARN([LAPACK library not found])])

LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
