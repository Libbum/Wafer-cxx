/*

   mexHatPotential.hpp

   Tim DuBois 24/09/12

*/

#ifndef __mexHatPotential_hpp__
#define __mexHatPotential_hpp__

void matRightDivide(double *A, double *B, int M, int N);
void matMultiply(double *A, double *B, double *C, int M, int N);
void matInverse(double *A, int N);
double dist(double *pointsi, double *pointsj, int sizeC);
void gammasm(double *za, double *zb, double *r, double *gambfa, double *gamfafb);
double electroneg(double *system, double *species, int sizeC);
dcomp makepot(double *points, double *species, int sizeC);
dcomp mexHatPotential(double dx, double dy, double dz);

#endif /* __mexHatPotential_hpp__ */
