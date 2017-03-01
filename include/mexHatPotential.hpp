/*

   mexHatPotential.h

   Tim DuBois 24/09/12

*/

#ifndef __mexHatPotential_hpp__
#define __mexHatPotential_hpp__

typedef complex<double> dcomp;
extern double ALX,ALY,ALZ;


#if defined(USINGMKL)
#include <mkl.h>
#define LP_INT lapack_int
#else
#define LP_INT int
#endif

#if defined(USINGLAPACK)
extern "C"
{
    void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv,
            double *b, int *ldb, int *info );

    int dgetrf_(int *m, int *n, double *a, int *
            lda, int *ipiv, int *info);

    int dgetri_(int *n, double *a, int *lda, int
            *ipiv, double *work, int *lwork, int *info);

    void dgemm_(char*,char*,int*,int*,int*,double*,double*,int*,
            double*,int*,double*,double*,int*);
}
#endif


void matRightDivide(double *A, double *B, int M, int N);
void matMultiply(double *A, double *B, double *C, int M, int N);
void matInverse(double *A, int N);
double dist(double *pointsi, double *pointsj, int sizeC);
void gammasm(double *za, double *zb, double *r, double *gambfa, double *gamfafb);
double electroneg(double *system, double *species, int sizeC);
dcomp makepot(double *points, double *species, int sizeC);
dcomp mexHatPotential(double dx, double dy, double dz);

#endif /* __mexHatPotential_hpp__ */
