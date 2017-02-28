/*

   grid.cpp

   Copyright (c) Michael Strickland
   Forked at v2.0; Additions by Tim DuBois

   GNU General Public License (GPLv3)
   See detailed text in license directory

*/

#include <cmath>
#include <cstring>
#include <stdlib.h>
#include <fstream>
	
#include "wafer.hpp"
#include "grid.hpp"
#include "outputroutines.hpp"
#include "initialconditions.hpp"
#include "potential.hpp"

#include <complex>

// FDTD 3d Schroedinger Eq Solver

using namespace std;

// this holds the current values of the wavefunction
dcomp ***w;

// this holds the updated values of the wavefunction
dcomp ***W;

// this holds the snapshots of the wavefunction
dcomp ****wstore;

// this holds the potential array
dcomp ***v;
dcomp ***vBase;

// these hold the alpha and beta arrays which are used during updates
dcomp ***a,***b;

// temp for point swap
dcomp ***tmp;

//nsnaps = STEPS/SNAPUPDATE;
int nsnaps = 0;

//Enegry penalty offset for beta
double epsilon = 0;

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
    
    v = new dcomp**[NUMX+6]; 
	for (int sx=0;sx<NUMX+6;sx++) v[sx] = new dcomp*[NUMY+6];
	for (int sx=0;sx<NUMX+6;sx++) for (int sy=0;sy<NUMY+6;sy++) v[sx][sy] = new dcomp[DISTNUMZ+6];

    vBase = new dcomp**[NUMX+6]; 
	for (int sx=0;sx<NUMX+6;sx++) vBase[sx] = new dcomp*[NUMY+6];
	for (int sx=0;sx<NUMX+6;sx++) for (int sy=0;sy<NUMY+6;sy++) vBase[sx][sy] = new dcomp[DISTNUMZ+6];

	a = new dcomp**[NUMX+6]; 
	for (int sx=0;sx<NUMX+6;sx++) a[sx] = new dcomp*[NUMY+6];
	for (int sx=0;sx<NUMX+6;sx++) for (int sy=0;sy<NUMY+6;sy++) a[sx][sy] = new dcomp[DISTNUMZ+6];

	b = new dcomp**[NUMX+6]; 
	for (int sx=0;sx<NUMX+6;sx++) b[sx] = new dcomp*[NUMY+6];
	for (int sx=0;sx<NUMX+6;sx++) for (int sy=0;sy<NUMY+6;sy++) b[sx][sy] = new dcomp[DISTNUMZ+6];

	wstore = new dcomp***[WAVEMAX+1]; 
	for (int n=0;n<=WAVEMAX;n++) wstore[n] = new dcomp**[NUMX+6];
	for (int n=0;n<=WAVEMAX;n++) for (int sx=0;sx<NUMX+6;sx++) wstore[n][sx] = new dcomp*[NUMY+6];
	for (int n=0;n<=WAVEMAX;n++) for (int sx=0;sx<NUMX+6;sx++) for (int sy=0;sy<NUMY+6;sy++) wstore[n][sx][sy] = new dcomp[DISTNUMZ+6];
		
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
    
    dcomp minima = potential(0,0,0);
    if (POTFILE) {
        potentialFromFile();
    }

    for (sx=0;sx<=NUMX+5;sx++)
    for (sy=0;sy<=NUMY+5;sy++)
    for (sz=0; sz<=DISTNUMZ+5;sz++) {
        if (!POTFILE) {
            v[sx][sy][sz] = potential(sx,sy,sz);
        }
        vBase[sx][sy][sz] = v[sx][sy][sz];
        b[sx][sy][sz] = 1./(1.+EPS*v[sx][sy][sz]/((dcomp) 2.));
        a[sx][sy][sz] = (1.-EPS*v[sx][sy][sz]/((dcomp) 2.))*b[sx][sy][sz];
        if (real(v[sx][sy][sz])<real(minima)) {
            minima = v[sx][sy][sz];
        }
    } 
    
    //Get 2*abs(min(potential)) for offset of beta
    epsilon = 2*abs(real(minima));
}

void potentialFromFile()
{
		// read from file
        string line;
	    char fname[32];
	    fstream input, debug_out;
	    int sx,sy,sz;
        double pot;

        sprintf(fname,"data/potential_%d.dat",nodeID);
	    input.open(fname, ios::in);
	    if (nodeID==1) cout << "==> Potential read from file" << endl;
	    //TODO: Make a failsafe for this
        //while( getline( input, line ) ) lines.push_back( line ) ;
       
        sprintf(fname,"debug/input_%d.txt",nodeID);
        debug_out.open(fname, ios::out);
       
		for (sx=0;sx<=NUMX+5;sx++)
			for (sy=0;sy<=NUMY+5;sy++)
				for (sz=0; sz<=DISTNUMZ+5;sz++) {
				    getline( input, line );	
                    pot = atof(line.c_str());
                    v[sx][sy][sz] = dcomp(pot,0.0);

                    debug_out << sx << " " << sy << " " << sz << " " << v[sx][sy][sz] << endl;
				}
        debug_out.close();
		input.close();

}

//Updates potintial with the current energy penalty
void updatePotential(dcomp beta)
{
    int sx,sy,sz;
        
    for (sx=0;sx<=NUMX+5;sx++)
    for (sy=0;sy<=NUMY+5;sy++)
    for (sz=0; sz<=DISTNUMZ+5;sz++) {
        v[sx][sy][sz] = vBase[sx][sy][sz] + epsilon*abs(beta*beta); //Add epsilon|beta|^2 offset
        b[sx][sy][sz] = 1./(1.+EPS*v[sx][sy][sz]/((dcomp) 2.));
        a[sx][sy][sz] = (1.-EPS*v[sx][sy][sz]/((dcomp) 2.))*v[sx][sy][sz];
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
        for (int sx=0;sx<=NUMX+5;sx++) {
      	  for (int sy=0;sy<=NUMY+5;sy++) {
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

    	for (int sx=0;sx<=NUMX+5;sx++) 
      	  for (int sy=0;sy<=NUMY+5;sy++)
       	    for (int sz=0;sz<=DISTNUMZ+5;sz++) {
		// no need to update the slices which were already loaded by updatedBoundaries
		if (sz<=5 || sz>=DISTNUMZ) continue;
		if (sx>=3 && sx<=2+NUMX && sy>=3 && sy<=2+NUMY && sz>=3 && sz<=2+DISTNUMZ)
		  W[sx][sy][sz] = updateRule(sx,sy,sz,step);
		else 
		  W[sx][sy][sz] = w[sx][sy][sz]; 
	    }
}

void storeConverged(dcomp*** wfnc,int num) { 
  for (int sx=0;sx<=NUMX+5;sx++) 
    for (int sy=0;sy<=NUMY+5;sy++)
        for (int sz=0; sz<=DISTNUMZ+5;sz++) {
		    wstore[num][sx][sy][sz] = wfnc[sx][sy][sz];
        }
}

void normalizeWavefunction(dcomp*** wfnc) {
  dcomp norm = sqrt(normalizationCollect);
    	for (int sx=0;sx<=NUMX+5;sx++) 
      	  for (int sy=0;sy<=NUMY+5;sy++)
       	    for (int sz=0; sz<=DISTNUMZ+5;sz++) 
		wfnc[sx][sy][sz] /= norm;
}
