/*

   potential.cpp

   Copyright (c) Michael Strickland

   GNU General Public License (GPLv3)
   See detailed text in license directory

*/

#include <cmath>

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <complex>

using namespace std;

#include "mpisolve.h"
#include "grid.h"
#include "potential.h"
#include "intde2.h"
#include "mexHatPotential.h"

// used for numerical integration when necessary
int     lenaw=LENAW;
double  aw[LENAW];

// global variable useful for subroutines
double  dx,dy,dz;
double  r;
double  md;

// determines square of distance to center of simulation volume in lattice units
double distsq(int sx,int sy, int sz) 
{
	double dx,dy,dz,r2;

	// coordinate system is centered in simulation volume 
	dx = sx - ((double)NUMX+1.)/2.;
	dy = sy - ((double)NUMY+1.)/2.;
	dz = sz - ((double)DISTNUMZ+1.)/2. + ( ((double)nodeID) - ((double)numNodes)/2. )*DISTNUMZ;
	r2 = (dx*dx+dy*dy+dz*dz);
	return r2;
}



dcomp potential(int sx,int sy, int sz) 
{
	double temp,iV,rV;
	double res,err;
	
	// coordinate system is centered in simulation volume 
	dx = ((double) sx) - ((double)NUMX+1.)/2.;
	dy = ((double) sy) - ((double)NUMY+1.)/2.;
	dz = ((double) sz) - ((double)DISTNUMZ+1.)/2. + ( ((double)nodeID) - ((double)numNodes)/2. )*DISTNUMZ;
  r = A*sqrt(dx*dx+dy*dy+dz*dz);

	switch(POTENTIAL) {
	  case 0:
		// none
          	return 0.;
		break;
	  case 1:
		// cubic well
		if ( (sx>NUMX/4 && sx<=3*NUMX/4) && (sy>NUMY/4 && sy<=3*NUMY/4) && (sz>NUMZ/4 && sz<=3*NUMZ/4) )
			return -10.0;
		else
			return 0.0;
		break;
	  case 2:
		// quadrilateral-well in center of cube with short side in z direction
		if ( (sx>NUMX/4 && sx<=3*NUMX/4) && (sy>NUMY/4 && sy<=3*NUMY/4) && (sz>3*NUMZ/8 && sz<=5*NUMZ/8) )
			return -10.0;
		else
			return 0.0;
		break;
	  case 3:
		// 3d periodic
		temp  = sin(2*M_PI*(sx-1)/(NUMX-1))*sin(2*M_PI*(sx-1)/(NUMX-1));
		temp *= sin(2*M_PI*(sy-1)/(NUMY-1))*sin(2*M_PI*(sy-1)/(NUMY-1));
		temp *= sin(2*M_PI*(sz-1)/(NUMZ-1))*sin(2*M_PI*(sz-1)/(NUMZ-1));
		return -temp+1;
		break;
	  case 4:
		// coulomb
		if (r < A)
		  return -1./A;
		else
		  return -1./r; // + 1./A;
		break;
	  case 5:
		// elliptical coulomb
		dz *= 2;
		r = A*sqrt(dx*dx+dy*dy+dz*dz);
		if (r < A)
		  return 0.0;
		else
		  return -1./r + 1./A;
		break;
	  case 6:
		// cornell
		// units here are GeV for energy/momentum and GeV^(-1) for distance
		md = mu(T);
		if (r < A)
		  return 4*MASS;
		else
		return -0.5*(4./3.)/r + SIGMA*r + 4*MASS;
		break;
	  case 7:
		// screened cornell
		// units here are GeV for energy/momentum and GeV^(-1) for distance
		md = mu(T);
		if (r < A)
		  return 4*MASS;
		else
		  return -alphas(2*M_PI*T)*(4./3.)*exp(-md*r)/r + SIGMA*(1. - exp(-md*r))/md + 4*MASS;
		break;
	  case 8:
		// screened cornell + spin correction
		// units here are GeV for energy/momentum and GeV^(-1) for distance
		md = mu(T);
		if (r < A)
		  return 4*MASS;
		else
		  return -alphas(2*M_PI*T)*(4./3.)*exp(-md*r)/r + SIGMA*(1. - exp(-md*r))/md  - 0.8*SIGMA/(4*MASS*MASS*r) + 4*MASS;
		break;
	  case 9:
		// anisotropically screened short distance piece + isotropic cornell + spin correction
		// units here are GeV for energy/momentum and GeV^(-1) for distance
		md = mu(T)*(1 + 0.07*pow(XI,0.2)*(1-A*A*dz*dz/(r*r)))*pow(1+XI,-0.29);
		if (r < A)
		  return 4*MASS;
		else
		  return -alphas(2*M_PI*T)*(4./3.)*exp(-md*r)/r + SIGMA*(1. - exp(-mu(T)*r))/mu(T)  - 0.8*SIGMA/(4*MASS*MASS*r) + 4*MASS;
		break;
	  case 10:
		// anisotropically screened cornell + spin correction
		// units here are GeV for energy/momentum and GeV^(-1) for distance
		md = mu(T)*(1 + 0.07*pow(XI,0.2)*(1-A*A*dz*dz/(r*r)))*pow(1+XI,-0.29);
		if (r < A)
		  return 4*MASS;
		else
		  return -alphas(2*M_PI*T)*(4./3.)*exp(-md*r)/r + SIGMA*(1. - exp(-md*r))/mu(T)  - 0.8*SIGMA/(4*MASS*MASS*r) + 4*MASS;
		break;
	  case 11:
		// fully anisotropic screened cornell + spin correction
		// units here are GeV for energy/momentum and GeV^(-1) for distance
		md = mu(T)*(1 + 0.07*pow(XI,0.2)*(1-A*A*dz*dz/(r*r)))*pow(1+XI,-0.29);
		if (r < A)
		  return 4*MASS;
		else
		  return -alphas(2*M_PI*T)*(4./3.)*exp(-md*r)/r + SIGMA*(1. - exp(-md*r))/md  - 0.8*SIGMA/(4*MASS*MASS*r) + 4*MASS;
		break;
	  case 12:
		// fully anisotropic screened cornell using small xi expression for mu + spin correction
		// units here are GeV for energy/momentum and GeV^(-1) for distance
		md = mu(T)*(1 - 0.125*XI*(A*A*dz*dz/(r*r)+1));
		if (r < A)
		  return 4*MASS;
		else
		  return -alphas(2*M_PI*T)*(4./3.)*exp(-md*r)/r + SIGMA*(1. - exp(-md*r))/md  - 0.8*SIGMA/(4*MASS*MASS*r) + 4*MASS;
		break;
	  case 13:
		// modified fully anisotropic screened cornell + spin correction
		// units here are GeV for energy/momentum and GeV^(-1) for distance
		//md = mu(T)*(pow(1+1.85*pow(XI,1.27),-0.20)+(pow(1+0.74*pow(XI,1.20),-0.23)-pow(1+1.85*pow(XI,1.27),-0.20))*(1-A*A*dz*dz/(r*r)));
		md = mu(T)*(pow(1+1.85*pow(XI,1.27),-0.20)+(pow(1+0.74*pow(XI,1.20),-0.23)-pow(1+1.85*pow(XI,1.27),-0.20))*(1-A*A*dz*dz/(r*r)));
		if (r < A)
			return 4*MASS;
		else
			//return -alphas(2*M_PI*T)*(4./3.)*exp(-md*r)/r + SIGMA*(1. - exp(-md*r))/md  - 0.8*SIGMA/(4*MASS*MASS*r) + 4*MASS;
			return -0.385*exp(-md*r)/r*(1.0 + md*r) + 2.* SIGMA*(1. - exp(-md*r))/md - SIGMA*r*exp(-md*r)  - 0.8*SIGMA/(4*MASS*MASS*r) + 4*MASS;
		break;
	  case 14:
		// newpotential add entropy contribution
		// units here are GeV for energy/momentum and GeV^(-1) for distance
		md = mu(3)*(1 - 0.125*XI*(A*A*dz*dz/(r*r)+1))*T/3;
		if (r < A)
			return 4*MASS;
		else
			return -0.385*exp(-md*r)/r*(1.0 + md*r) + 2.* SIGMA*(1. - exp(-md*r))/md - SIGMA*r*exp(-md*r)  - 0.8*SIGMA/(4*MASS*MASS*r) + 4*MASS;
		break;
	  case 15:
		// 3d harmonic oscillator
          	return r*r/2;
		break;			
	  case 16:
		// Mickey Mouse's Head
		double Dx, Dy, Dz, R;
		if (r/A <= NUMZ/4) return -100.0; // head
		Dx = sx - ((double)NUMX+1.)/2. - (1/sqrt(2.)+0.25)*NUMY/4;
		Dy = sy - ((double)NUMY+1.)/2. - (1/sqrt(2.)+0.25)*NUMX/4;
		Dz = ((double) sz) - ((double)DISTNUMZ+1.)/2. + ( ((double)nodeID) - ((double)numNodes)/2. )*DISTNUMZ;
		R = sqrt(4*Dx*Dx+Dy*Dy+Dz*Dz);
		if (R < NUMZ/8) return -105.0; // ear
		Dy = sy - ((double)NUMY+1.)/2. + (1/sqrt(2.)+0.25)*NUMY/4;
		Dz = ((double) sz) - ((double)DISTNUMZ+1.)/2. + ( ((double)nodeID) - ((double)numNodes)/2. )*DISTNUMZ;
		R = sqrt(4*Dx*Dx+Dy*Dy+Dz*Dz);
		if (R < NUMZ/8) return -105.0; // ear
		Dx = sx - ((double)NUMX+1.)/2.;
		Dy = sy - ((double)NUMY+1.)/2.;
		Dz = ((double) sz) - ((double)DISTNUMZ+1.)/2. - ((double)NUMZ/8.)+ ( ((double)nodeID) - ((double)numNodes)/2. )*DISTNUMZ;
		R = sqrt(Dx*Dx+Dy*Dy+Dz*Dz);
		if (R < NUMZ/6) return -100.0; // nose
		return 0;
		break;
	  case 17:
		// Dodecahedron
		double x, y, z;
		x = dx/((NUMX-1)/2);
		y = dy/((NUMY-1)/2);
		z = dz/((NUMZ-1)/2);
		if (12.70820393249937 + 11.210068307552588*x >= 14.674169922690343*z && 11.210068307552588*x <= 12.70820393249937 + 14.674169922690343*z && 5.605034153776295*(3.23606797749979*x - 1.2360679774997896*z) <= 6.*(4.23606797749979 + 5.23606797749979*y) && 18.1382715378281*x + 3.464101615137755*z <= 12.70820393249937 && 9.06913576891405*x + 15.70820393249937*y <= 12.70820393249937 + 3.464101615137755*z && 9.70820393249937*y <= 12.70820393249937 + 5.605034153776294*x + 14.674169922690343*z && 12.70820393249937 + 5.605034153776294*x + 9.70820393249937*y + 14.674169922690343*z >= 0. && 15.70820393249937*y + 3.464101615137755*z <= 12.70820393249937 + 9.06913576891405*x && 5.605034153776295*(-6.47213595499958*x - 1.2360679774997896*z) <= 25.41640786499874 && 3.464101615137755*z <= 9.06913576891405*x + 3.*(4.23606797749979 + 5.23606797749979*y) && 1.7320508075688772*(3.23606797749979*x + 8.47213595499958*z) <= 3.*(4.23606797749979 + 3.23606797749979*y) && 5.605034153776294*x + 9.70820393249937*y + 14.674169922690343*z <= 12.70820393249937) return -100.0;
		else
		  return 0.;
	  case 18:
		// Complex 3d harmonic oscillator
		return dcomp(1.,1.)*dcomp(r*r/2,0.);
		break;
	  case 19:			
			// complex coulomb
			if (r < A)
				return dcomp(0.,1.)*dcomp(-1./A,0.);
			else
				return dcomp(0.,1.)*dcomp(-1./r,0.);
			break;			
	  case 20:
		// Complex Dodecahedron
		// double x, y, z;
		x = dx/((NUMX-1)/2);
		y = dy/((NUMY-1)/2);
		z = dz/((NUMZ-1)/2);
		if (12.70820393249937 + 11.210068307552588*x >= 14.674169922690343*z && 11.210068307552588*x <= 12.70820393249937 + 14.674169922690343*z && 5.605034153776295*(3.23606797749979*x - 1.2360679774997896*z) <= 6.*(4.23606797749979 + 5.23606797749979*y) && 18.1382715378281*x + 3.464101615137755*z <= 12.70820393249937 && 9.06913576891405*x + 15.70820393249937*y <= 12.70820393249937 + 3.464101615137755*z && 9.70820393249937*y <= 12.70820393249937 + 5.605034153776294*x + 14.674169922690343*z && 12.70820393249937 + 5.605034153776294*x + 9.70820393249937*y + 14.674169922690343*z >= 0. && 15.70820393249937*y + 3.464101615137755*z <= 12.70820393249937 + 9.06913576891405*x && 5.605034153776295*(-6.47213595499958*x - 1.2360679774997896*z) <= 25.41640786499874 && 3.464101615137755*z <= 9.06913576891405*x + 3.*(4.23606797749979 + 5.23606797749979*y) && 1.7320508075688772*(3.23606797749979*x + 8.47213595499958*z) <= 3.*(4.23606797749979 + 3.23606797749979*y) && 5.605034153776294*x + 9.70820393249937*y + 14.674169922690343*z <= 12.70820393249937) return dcomp(-100.,-100.);
		else
		  return 0.;
	  case 21:
		// Anisotropically Screened Quarkonium Potential with Imaginary Part
		if (r < A)
		  return dcomp(4*MASS,0.);
		// real part
		md = mu(3)*(1 - 0.125*XI*(A*A*dz*dz/(r*r)+1))*T/3;
		rV = -0.385*exp(-md*r)/r*(1.0 + md*r) + 2.* SIGMA*(1. - exp(-md*r))/md - SIGMA*r*exp(-md*r)  - 0.8*SIGMA/(4*MASS*MASS*r) + 4*MASS;
		// imaginary part
		md = mu(3); // do not include angular modification of md for imaginary part!
		/* 
		cout << "r = " << r << ";" << endl;
		cout << "md = " << md << ";" << endl;
		cout << "dx = " << A*dx << ";" << endl;
		cout << "dy = " << A*dy << ";" << endl;
		cout << "dz = " << A*dz << ";" << endl;
		*/
		intdeoini(lenaw,TINY,INTACC,aw);
		// phir contribution
		intdeo(phir, TINY, md*r,  aw, &res, &err);
		//cout << "DEBUG: phi = " << res << endl;
		iV = -res;
		// psi1 contribution	
		intdeo(psi1, TINY, md*r, aw, &res, &err);
		//cout << "DEBUG: psi1 = " << res << endl;
		iV += XI*res;
		// psi2 contribution
		intdeo(psi2, TINY, md*r, aw, &res, &err);
		//cout << "DEBUG: psi2 = " << res << endl;
		iV += XI*res;		
		iV *= 0.385*T*TC;
		return dcomp(rV,iV);
		break;
    case 22: //Mexican Hat
    {
      double gr; //gr is now a dynamic version of GR - no longer needed in params.txt
      gr = NUMX*A/2-(A/2);
      double tx = -((double) gr + A) + ((double) sx)*((2.0*((double) gr))/((double) NUMX - 1));
      gr = NUMY*A/2-(A/2);
      double ty = -((double) gr + A) + ((double) sy)*((2.0*((double) gr))/((double) NUMY - 1));
      gr = NUMZ*A/2-(A/2);
      double tz = -((double) gr + A) + ((double) sz+(nodeID-1)*DISTNUMZ)*((2.0*((double) gr))/((double) NUMZ - 1));
      //if (( nodeID == 1) && ( sy == 0 ) && ( sz == 0)) {
      //if ((sx == 0 ) && ( sy == 0 ) && ( sz == 0)) {
        //cout << "tx " << tx << ", ty " << ty << ", tz " << tz << ", sx " << sx << ", numNodes " << numNodes << endl;
      //}

      return mexHatPotential(tx,ty,tz);
    }
    break;
    default:
		return 0.;
		break;
	}
}

// returns value of potential which should be subtracted when computing binding energies
dcomp potentialSub(int sx, int sy, int sz) 
{
	double iV,rV;

	// coordinate system is centered in simulation volume 
	dx = ((double) sx) - ((double)NUMX+1.)/2.;
	dy = ((double) sy) - ((double)NUMY+1.)/2.;
	dz = ((double) sz) - ((double)DISTNUMZ+1.)/2. + ( ((double)nodeID) - ((double)numNodes)/2. )*DISTNUMZ;
	r = A*sqrt(dx*dx+dy*dy+dz*dz);

	switch(POTENTIAL) {
	  case 0:
	  case 1:
	  case 2:
	  case 3:
		return 0.;
		break;
	  case 4:
	  case 5:
		return 1./A;
		break;
	  case 6:
		return 4*MASS;
		break;
	  case 7:
	  case 8:
	  case 9:
	  case 10:
		return SIGMA/mu(T) + 4*MASS;
		break;
	  case 11:
		md = mu(T)*(1 + 0.07*pow(XI,0.2)*(1-A*A*dz*dz/(r*r)))*pow(1+XI,-0.29);
		return SIGMA/md + 4*MASS;
		break;
	  case 12:
		md = mu(T)*(1 - 0.125*XI*(A*A*dz*dz/(r*r)+1));
		//cout << SIGMA <<  ", " << md << ", " << MASS << endl;
		return SIGMA/md + 4*MASS;
		break;
	  case 13:
		md = mu(T)*(pow(1+1.85*pow(XI,1.27),-0.20)+(pow(1+0.74*pow(XI,1.20),-0.23)-pow(1+1.85*pow(XI,1.27),-0.20))*(1-A*A*dz*dz/(r*r)));
		return SIGMA/md + 4*MASS;
		break;
	  case 14:
		md = mu(3)*(1 - 0.125*XI*(A*A*dz*dz/(r*r)+1))*T/3;
		return 2*SIGMA/md + 4*MASS;
		break;
	  case 15:
		return 0.;
		break;			
	  case 16:
		return 0.;
		break;			
	  case 17:
		return 0.;
		break;			
	  case 18:
		return 0.;
		break;	
	  case 19:
		return 0.;
		break;	
	  case 20:
		return 0.;
		break;		
	  case 21:
		md = mu(3)*(1 - 0.125*XI*(A*A*dz*dz/(r*r)+1))*T/3;
		rV = 2*SIGMA/md + 4*MASS;
		iV = -1 + XI/6;
		iV *= 0.385*T*TC;
		return dcomp(rV,iV);
		break;				
    case 22:
    return 0.;
    break;
    default:
		return 0.;
		break;
	}
}

// phir integrand
double phir(double z)
{
	return 2*z*(1-sin(z*r*md)/(z*r*md))/(z*z+1)/(z*z+1);
}

// utility function for ps1 and ps2
double psig(double r, double z)
{	// note r here is rhat
	return (r*z*cos(r*z)-sin(r*z))/(r*r*r*z*z*z);
}

// psi1 integrand
double psi1(double z)
{
	double cos2theta = A*A*dz*dz/(r*r);
	double sin2theta = 1. - cos2theta;
	return z*(1 - 1.5*(sin2theta*sin(z*md*r)/(z*md*r) + (1-3*cos2theta)*psig(r*md, z)))/(z*z+1)/(z*z+1);
}

// psi2 integrand
double psi2(double z)
{
	double cos2theta = A*A*dz*dz/(r*r);
//	double sin2theta = 1. - cos2theta;
	return -4*z*(1 - 3*((2./3. - cos2theta)*sin(z*md*r)/(z*md*r) + (1-3*cos2theta)*psig(r*md, z)))/(z*z+1)/(z*z+1)/(z*z+1)/3;

}


// running coupling
double alphas(double mu)
{
        double b0,b1,b2,L,R;

        b0 = 11 - 2*((double)NF)/3;
        b1 = 51 - 19*((double)NF)/3;
        b2 = 2857 - 5033*((double)NF)/9 + 325*((double)NF)*((double)NF)/27;

        R = 2.3; // scale adjusted to match lattice data from hep-lat/0503017v2

        L = 2*log(mu/R);

        return 4*M_PI*(1 - 2*b1*log(L)/(b0*b0*L) + 4*b1*b1*((log(L)-0.5)*(log(L)-0.5)+b2*b0/(8*b1*b1)-5.0/4.0)/(b0*b0*b0*b0*L*L))/(b0*L);
}

// debye screening mass
double mu(double t)
{
	return 1.4*sqrt( (1+((double)NF)/6)*4*M_PI*alphas(2*M_PI*t) )*t*TC;
}
