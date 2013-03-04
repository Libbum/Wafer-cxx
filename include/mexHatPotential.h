/*

   mexHatPotential.h

   Tim DuBois 24/09/12

*/

#ifndef __mexHatPotential_h__
#define __mexHatPotential_h__

typedef complex<double> dcomp;
extern double ALX,ALY,ALZ;

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

void matRightDivide(double *A, double *B, ptrdiff_t M, ptrdiff_t N);
void matMultiply(double *A, double *B, double *C, ptrdiff_t M, ptrdiff_t N);
void matInverse(double *A, ptrdiff_t N);
double dist(double *pointsi, double *pointsj);
void gammasm(double *za, double *zb, double *r, double *gambfa, double *gamfafb);
double electroneg(double *system);
dcomp makepot(double *points);
dcomp mexHatPotential(double dx, double dy, double dz);

#endif /* __mexHatPotential_h__ */
