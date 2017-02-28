/*

   intde2.h

   Copyright (c) Michael Strickland
   Forked at v2.0; Additions by Tim DuBois

   GNU General Public License (GPLv3)
   See detailed text in license directory

*/

#ifndef __intde2_h__
#define __intde2_h__

#define SMALL   1.0e-12         // used for upperlimit offset in w (0->k) integration
#define INTACC  1.0e-8          // target integration accuracy
#define TINY    1.0e-300        // tiny number used by integrator
#define BIG     1.0e+300        // big number
#define LENAW   8000            // number of quadrature points
#define IRCUT   1.0e-8          // IR cut on some integrations
#define UVCUT   1.0e+8          // UV cut on inf integrations

//extern int     lenaw;
//extern double  *aw;

void intdeini(int lenaw, double tiny, double eps, double *aw);
void intde(double (*f)(double), double a, double b, double *aw, double *i, double *err);

void intdeiini(int lenaw, double tiny, double eps, double *aw);
void intdei(double (*f)(double), double a, double *aw, double *i, double *err);

void intdeoini(int lenaw, double tiny, double eps, double *aw);
void intdeo(double (*f)(double), double a, double omega, double *aw, double *i, double *err);

#endif /* __intde2_h__ */
