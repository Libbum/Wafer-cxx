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
        
	w = new dcomp**[NUMX+6]; 
	for (int sx=0;sx<NUMX+6;sx++) w[sx] = new dcomp*[NUMY+6];
	for (int sx=0;sx<NUMX+6;sx++) for (int sy=0;sy<NUMY+6;sy++) w[sx][sy] = new dcomp[DISTNUMZ+6];

	W = new dcomp**[NUMX+6]; 
	for (int sx=0;sx<NUMX+6;sx++) W[sx] = new dcomp*[NUMY+6];
	for (int sx=0;sx<NUMX+6;sx++) for (int sy=0;sy<NUMY+6;sy++) W[sx][sy] = new dcomp[DISTNUMZ+6];

        W2 = new dcomp**[NUMX+6];
        for (int sx=0;sx<NUMX+6;sx++) W2[sx] = new dcomp*[NUMY+6];
        for (int sx=0;sx<NUMX+6;sx++) for (int sy=0;sy<NUMY+6;sy++) W2[sx][sy] = new dcomp[DISTNUMZ+6];

        v = new dcomp**[NUMX+6]; 
	for (int sx=0;sx<NUMX+6;sx++) v[sx] = new dcomp*[NUMY+6];
	for (int sx=0;sx<NUMX+6;sx++) for (int sy=0;sy<NUMY+6;sy++) v[sx][sy] = new dcomp[DISTNUMZ+6];

	a = new dcomp**[NUMX+6]; 
	for (int sx=0;sx<NUMX+6;sx++) a[sx] = new dcomp*[NUMY+6];
	for (int sx=0;sx<NUMX+6;sx++) for (int sy=0;sy<NUMY+6;sy++) a[sx][sy] = new dcomp[DISTNUMZ+6];

	b = new dcomp**[NUMX+6]; 
	for (int sx=0;sx<NUMX+6;sx++) b[sx] = new dcomp*[NUMY+6];
	for (int sx=0;sx<NUMX+6;sx++) for (int sy=0;sy<NUMY+6;sy++) b[sx][sy] = new dcomp[DISTNUMZ+6];

	int snaps = 2;
	wstore = new dcomp***[snaps]; 
	for (int n=0;n<snaps;n++) wstore[n] = new dcomp**[NUMX+6];
	for (int n=0;n<snaps;n++) for (int sx=0;sx<NUMX+6;sx++) wstore[n][sx] = new dcomp*[NUMY+6];
	for (int n=0;n<snaps;n++) for (int sx=0;sx<NUMX+6;sx++) for (int sy=0;sy<NUMY+6;sy++) wstore[n][sx][sy] = new dcomp[DISTNUMZ+6];
		
	return;
}

void allocateClusterMemory() {
   clust = (double *)calloc(CLUSTSIZE*3,sizeof(double));
   clustSpecies = (double *)calloc(CLUSTSIZE,sizeof(double));
}

void deallocateClusterMemory() {
   free(clust);
   free(clustSpecies);
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

        for (sx=0;sx<=NUMX+5;sx++)
        for (sy=0;sy<=NUMY+5;sy++)
        for (sz=0; sz<=DISTNUMZ+5;sz++) {
          v[sx][sy][sz] = potential(sx,sy,sz);
          b[sx][sy][sz] = 1./(1.+EPS*v[sx][sy][sz]/((dcomp) 2.));
          a[sx][sy][sz] = (1.-EPS*v[sx][sy][sz]/((dcomp) 2.))*b[sx][sy][sz];
        }  
}

// compute energy of a wave function
dcomp wfncEnergy(dcomp*** wfnc) {
	dcomp res=0;
    	for (int sx=3;sx<=2+NUMX;sx++) 
      	  for (int sy=3;sy<=2+NUMY;sy++) 
       	    for (int sz=3;sz<=2+DISTNUMZ;sz++) {
		    res += v[sx][sy][sz]*conj(wfnc[sx][sy][sz])*wfnc[sx][sy][sz] - 
			   conj(wfnc[sx][sy][sz])*(  ((dcomp) 2.)*wfnc[sx+3][sy][sz] - ((dcomp) 27.)*wfnc[sx+2][sy][sz] + ((dcomp) 270.)*wfnc[sx+1][sy][sz] + ((dcomp) 270.)*wfnc[sx-1][sy][sz] - ((dcomp) 27.)*wfnc[sx-2][sy][sz] + ((dcomp) 2.)*wfnc[sx-3][sy][sz]
			                           + ((dcomp) 2.)*wfnc[sx][sy+3][sz] - ((dcomp) 27.)*wfnc[sx][sy+2][sz] + ((dcomp) 270.)*wfnc[sx][sy+1][sz] + ((dcomp) 270.)*wfnc[sx][sy-1][sz] - ((dcomp) 27.)*wfnc[sx][sy-2][sz] + ((dcomp) 2.)*wfnc[sx][sy-3][sz]
			                           + ((dcomp) 2.)*wfnc[sx][sy][sz+3] - ((dcomp) 27.)*wfnc[sx][sy][sz+2] + ((dcomp) 270.)*wfnc[sx][sy][sz+1] + ((dcomp) 270.)*wfnc[sx][sy][sz-1] - ((dcomp) 27.)*wfnc[sx][sy][sz-2] + ((dcomp) 2.)*wfnc[sx][sy][sz-3]
				                       - ((dcomp) 1470.)*wfnc[sx][sy][sz] )/(((dcomp) 360.)*A*A*MASS);
            }
	return res;
}

// a convenience method for getting energy of the current wave function
dcomp computeEnergy() 
{ return wfncEnergy(w);	}

// compute norm squared
dcomp wfncNorm2(dcomp*** wfnc) {
	dcomp norm=0;
    	for (int sx=3;sx<=2+NUMX;sx++) 
      	  for (int sy=3;sy<=2+NUMY;sy++) 
       	    for (int sz=3; sz<=2+DISTNUMZ;sz++) { 
	      norm += conj(wfnc[sx][sy][sz])*wfnc[sx][sy][sz]; 
	    }
	return norm;
}

// compute expectation value of vinfinity
dcomp vInfinityExpectationValue(dcomp*** wfnc) {
	dcomp expectation=0;
    	for (int sx=3;sx<=2+NUMX;sx++) 
      	  for (int sy=3;sy<=2+NUMY;sy++) 
       	    for (int sz=3;sz<=2+DISTNUMZ;sz++) { 
	      expectation += conj(wfnc[sx][sy][sz])*wfnc[sx][sy][sz]*potentialSub(sx,sy,sz); 
	    }
	return expectation;
}

// compute expectation value of r^2
dcomp r2ExpectationValue(dcomp*** wfnc) {
	dcomp expectation=0;
    	for (int sx=3;sx<=2+NUMX;sx++) 
      	  for (int sy=3;sy<=2+NUMY;sy++) 
       	    for (int sz=3;sz<=2+DISTNUMZ;sz++) { 
	      expectation += conj(wfnc[sx][sy][sz])*wfnc[sx][sy][sz]*distsq(sx,sy,sz); 
	    }
	return expectation;
}

inline dcomp updateRule(int sx, int sy, int sz, double step) 
{
  return w[sx][sy][sz]*a[sx][sy][sz] + b[sx][sy][sz]*step*(
	        ((dcomp) 2.)*w[sx+3][sy][sz] - ((dcomp) 27.)*w[sx+2][sy][sz] + ((dcomp) 270.)*w[sx+1][sy][sz] + ((dcomp) 270.)*w[sx-1][sy][sz] - ((dcomp) 27.)*w[sx-2][sy][sz] + ((dcomp) 2.)*w[sx-3][sy][sz]
	      + ((dcomp) 2.)*w[sx][sy+3][sz] - ((dcomp) 27.)*w[sx][sy+2][sz] + ((dcomp) 270.)*w[sx][sy+1][sz] + ((dcomp) 270.)*w[sx][sy-1][sz] - ((dcomp) 27.)*w[sx][sy-2][sz] + ((dcomp) 2.)*w[sx][sy-3][sz]
	      + ((dcomp) 2.)*w[sx][sy][sz+3] - ((dcomp) 27.)*w[sx][sy][sz+2] + ((dcomp) 270.)*w[sx][sy][sz+1] + ((dcomp) 270.)*w[sx][sy][sz-1] - ((dcomp) 27.)*w[sx][sy][sz-2] + ((dcomp) 2.)*w[sx][sy][sz-3]
	      - ((dcomp) 1470.)*w[sx][sy][sz] )/(((dcomp) 360.)*A*A*MASS);
}

// load updated left and right boundaries into W
void updateBoundaries(double step) {
        for (int sx=0;sx<NUMX+6;sx++) {
      	  for (int sy=0;sy<NUMY+6;sy++) {
            //Left side
       	    for (int sz=3;sz<=5;sz++) {
				if (sx>=3 && sx<=2+NUMX && sy>=3 && sy<=2+NUMY) {
					W[sx][sy][sz] = updateRule(sx,sy,sz,step);
                } else { 
					W[sx][sy][sz] = w[sx][sy][sz]; 
                }
			}
            //Right side
       	    for (int sz=DISTNUMZ;sz<=DISTNUMZ+2;sz++) {
				if (sx>=3 && sx<=2+NUMX && sy>=3 && sy<=2+NUMY) {
					W[sx][sy][sz] = updateRule(sx,sy,sz,step);
                } else { 
					W[sx][sy][sz] = w[sx][sy][sz]; 
                }
			}
          }
        }
}

// update the grid; note you should always call updatedBondaries before calling this routine
void updateInterior(double step) {

    	for (int sx=0;sx<NUMX+6;sx++) 
      	  for (int sy=0;sy<NUMY+6;sy++)
       	    for (int sz=0;sz<DISTNUMZ+6;sz++) {
		// no need to update the slices which were already loaded by updatedBoundaries
		if (sz<3 || sz>2+DISTNUMZ) continue;
		if (sx>=3 && sx<=2+NUMX && sy>=3 && sy<=2+NUMY && sz>=3 && sz<=2+DISTNUMZ)
		  W[sx][sy][sz] = updateRule(sx,sy,sz,step);
		else 
		  W[sx][sy][sz] = w[sx][sy][sz]; 
	    }
}

void recordSnapshot(dcomp*** wfnc) { //, int step
	//int snap = (int) step/SNAPUPDATE;
  //int nsnaps = (int) step/SNAPUPDATE;
  for (int sx=0;sx<=NUMX+5;sx++) 
      	  for (int sy=0;sy<=NUMY+5;sy++)
       	    for (int sz=0; sz<=DISTNUMZ+5;sz++) {
		wstore[1][sx][sy][sz] = wstore[0][sx][sy][sz];
		wstore[0][sx][sy][sz] = wfnc[sx][sy][sz];
            }
}

void normalizeWavefunction(dcomp*** wfnc) {
  dcomp norm = sqrt(normalizationCollect);
    	for (int sx=0;sx<=NUMX+5;sx++) 
      	  for (int sy=0;sy<=NUMY+5;sy++)
       	    for (int sz=0; sz<=DISTNUMZ+5;sz++) 
		wfnc[sx][sy][sz] /= norm;
}
