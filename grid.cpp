/*

   grid.cpp

   Copyright (c) Michael Strickland

   GNU General Public License (GPLv3)
   See detailed text in license directory

*/

#include <cmath>

#include "mpisolve.h"
#include "grid.h"
#include "outputroutines.h"
#include "initialconditions.h"
#include "potential.h"

#include <complex>

// FDTD 3d Schroedinger Eq Solver

using namespace std;

// this holds the current values of the wavefunction
dcomp ***w;

// this holds the updated values of the wavefunction
dcomp ***W;

// another variable for phi2
dcomp ***W2;

// this holds the snapshots of the wavefunction
dcomp ****wstore;

// this holds the potential array
dcomp ***v;

// these hold the alpha and beta arrays which are used during updates
dcomp ***a,***b;

// temp for point swap
dcomp ***tmp;

//nsnaps = STEPS/SNAPUPDATE;
int nsnaps = 0;

// allocate memory arrays
void allocateMemory() {

    	// cout << "==> Allocating memory\n";

        // indices x,y,z
        
	w = new dcomp**[NUM+2]; 
	for (int sx=0;sx<NUM+2;sx++) w[sx] = new dcomp*[NUM+2];
	for (int sx=0;sx<NUM+2;sx++) for (int sy=0;sy<NUM+2;sy++) w[sx][sy] = new dcomp[NUMZ+2];

	W = new dcomp**[NUM+2]; 
	for (int sx=0;sx<NUM+2;sx++) W[sx] = new dcomp*[NUM+2];
	for (int sx=0;sx<NUM+2;sx++) for (int sy=0;sy<NUM+2;sy++) W[sx][sy] = new dcomp[NUMZ+2];

    W2 = new dcomp**[NUM+2];
    for (int sx=0;sx<NUM+2;sx++) W2[sx] = new dcomp*[NUM+2];
    for (int sx=0;sx<NUM+2;sx++) for (int sy=0;sy<NUM+2;sy++) W2[sx][sy] = new dcomp[NUMZ+2];

    v = new dcomp**[NUM+2]; 
	for (int sx=0;sx<NUM+2;sx++) v[sx] = new dcomp*[NUM+2];
	for (int sx=0;sx<NUM+2;sx++) for (int sy=0;sy<NUM+2;sy++) v[sx][sy] = new dcomp[NUMZ+2];

	a = new dcomp**[NUM+2]; 
	for (int sx=0;sx<NUM+2;sx++) a[sx] = new dcomp*[NUM+2];
	for (int sx=0;sx<NUM+2;sx++) for (int sy=0;sy<NUM+2;sy++) a[sx][sy] = new dcomp[NUMZ+2];

	b = new dcomp**[NUM+2]; 
	for (int sx=0;sx<NUM+2;sx++) b[sx] = new dcomp*[NUM+2];
	for (int sx=0;sx<NUM+2;sx++) for (int sy=0;sy<NUM+2;sy++) b[sx][sy] = new dcomp[NUMZ+2];

	int snaps = 2;
	wstore = new dcomp***[snaps]; 
	for (int n=0;n<snaps;n++) wstore[n] = new dcomp**[NUM+2];
	for (int n=0;n<snaps;n++) for (int sx=0;sx<NUM+2;sx++) wstore[n][sx] = new dcomp*[NUM+2];
	for (int n=0;n<snaps;n++) for (int sx=0;sx<NUM+2;sx++) for (int sy=0;sy<NUM+2;sy++) wstore[n][sx][sy] = new dcomp[NUMZ+2];
		
	return;
}

// "copies" updated arrays using pointer swap
void copyDown() {
	tmp = w;
	w = W;
	W = tmp;
}

// initializes the potentials from potential.cpp 
void loadPotentialArrays()
{
        int sx,sy,sz;

        for (sx=0;sx<=NUM+1;sx++)
        for (sy=0;sy<=NUM+1;sy++)
        for (sz=0; sz<=NUMZ+1;sz++) {
          v[sx][sy][sz] = potential(sx,sy,sz);
          b[sx][sy][sz] = 1./(1.+EPS*v[sx][sy][sz]/((dcomp) 2.));
          a[sx][sy][sz] = (1.-EPS*v[sx][sy][sz]/((dcomp) 2.))*b[sx][sy][sz];
        }  
}

// compute energy of a wave function
dcomp wfncEnergy(dcomp*** wfnc) {
	dcomp res=0;
    	for (int sx=1;sx<=NUM;sx++) 
      	  for (int sy=1;sy<=NUM;sy++) 
       	    for (int sz=1;sz<=NUMZ;sz++) {
		    res += v[sx][sy][sz]*conj(wfnc[sx][sy][sz])*wfnc[sx][sy][sz] - 
			   conj(wfnc[sx][sy][sz])*(  wfnc[sx+1][sy][sz] + wfnc[sx-1][sy][sz] +
				                 wfnc[sx][sy+1][sz] + wfnc[sx][sy-1][sz] +
				                 wfnc[sx][sy][sz+1] + wfnc[sx][sy][sz-1] -
				                 ((dcomp) 6.)*wfnc[sx][sy][sz] )/(((dcomp) 2.)*A*A*MASS);
            }
	return res;
}

// a convenience method for getting energy of the current wave function
dcomp computeEnergy() 
{ return wfncEnergy(w);	}

// compute norm squared
dcomp wfncNorm2(dcomp*** wfnc) {
	dcomp norm=0;
    	for (int sx=1;sx<=NUM;sx++) 
      	  for (int sy=1;sy<=NUM;sy++) 
       	    for (int sz=1; sz<=NUMZ;sz++) { 
	      norm += conj(wfnc[sx][sy][sz])*wfnc[sx][sy][sz]; 
	    }
	return norm;
}

// compute expectation value of vinfinity
dcomp vInfinityExpectationValue(dcomp*** wfnc) {
	dcomp expectation=0;
    	for (int sx=1;sx<=NUM;sx++) 
      	  for (int sy=1;sy<=NUM;sy++) 
       	    for (int sz=1; sz<=NUMZ;sz++) { 
	      expectation += conj(wfnc[sx][sy][sz])*wfnc[sx][sy][sz]*potentialSub(sx,sy,sz); 
	    }
	return expectation;
}


// #ad.
// compute expectation value of r^2
dcomp r2ExpectationValue(dcomp*** wfnc) {
	dcomp expectation=0;
    	for (int sx=1;sx<=NUM;sx++) 
      	  for (int sy=1;sy<=NUM;sy++) 
       	    for (int sz=1; sz<=NUMZ;sz++) { 
	      expectation += conj(wfnc[sx][sy][sz])*wfnc[sx][sy][sz]*distsq(sx,sy,sz); 
	    }
	return expectation;
}



inline dcomp updateRule(int sx, int sy, int sz, double step) 
{
  return w[sx][sy][sz]*a[sx][sy][sz] + b[sx][sy][sz]*step*(
		  w[sx+1][sy][sz] + w[sx-1][sy][sz] +
		  w[sx][sy+1][sz] + w[sx][sy-1][sz] +
		  w[sx][sy][sz+1] + w[sx][sy][sz-1] -
		  ((dcomp) 6.)*w[sx][sy][sz])/(((dcomp) 2.)*A*A*MASS);
}

// load updated left and right boundaries into W
void updateBoundaries(double step) {

        for (int sx=0;sx<NUM+2;sx++)
      	  for (int sy=0;sy<NUM+2;sy++)
       	    for (int sz=1;sz<=NUMZ;sz+=NUMZ-1) {
				if (sx>=1 && sx<=NUM && sy>=1 && sy<=NUM) 
					W[sx][sy][sz] = updateRule(sx,sy,sz,step);
				else 
					W[sx][sy][sz] = w[sx][sy][sz]; 
			}
}

// update the grid; note you should always call updatedBondaries before calling this routine
void updateInterior(double step) {

    	for (int sx=0;sx<NUM+2;sx++) 
      	  for (int sy=0;sy<NUM+2;sy++)
       	    for (int sz=0;sz<NUMZ+2;sz++) {
		// no need to update the two slices which were already loaded by updatedBoundaries
		if (sz==1 || sz==NUMZ) continue;
		if (sx>=1 && sx<=NUM && sy>=1 && sy<=NUM && sz>=1 && sz<=NUMZ)
		  W[sx][sy][sz] = updateRule(sx,sy,sz,step);
		else 
		  W[sx][sy][sz] = w[sx][sy][sz]; 
	    }
}

void recordSnapshot(dcomp*** wfnc) { //, int step
	//int snap = (int) step/SNAPUPDATE;
  //int nsnaps = (int) step/SNAPUPDATE;
  for (int sx=0;sx<=NUM+1;sx++) 
      	  for (int sy=0;sy<=NUM+1;sy++)
       	    for (int sz=0; sz<=NUMZ+1;sz++) {
		wstore[1][sx][sy][sz] = wstore[0][sx][sy][sz];
		wstore[0][sx][sy][sz] = wfnc[sx][sy][sz];
            }
}

void normalizeWavefunction(dcomp*** wfnc) {
  dcomp norm = sqrt(normalizationCollect);
    	for (int sx=0;sx<=NUM+1;sx++) 
      	  for (int sy=0;sy<=NUM+1;sy++)
       	    for (int sz=0; sz<=NUMZ+1;sz++) 
		wfnc[sx][sy][sz] /= norm;
}
