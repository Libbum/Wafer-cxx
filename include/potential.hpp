/*

   potential.hpp

   Copyright (c) Michael Strickland
   Forked at v2.0; Additions by Tim DuBois

   GNU General Public License (GPLv3)
   See detailed text in license directory

*/

#ifndef __potential_hpp__
#define __potential_hpp__

dcomp potential(int sx, int sy, int sz);
dcomp potentialSub(int sx, int sy, int sz);
double alphas(double mu);
double mu(double t);

double distsq(int sx, int sy, int sz);

double phir(double z);
double psi1(double z);
double psi2(double z);

#endif /* __potential_hpp__ */
