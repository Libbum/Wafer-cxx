/*

   grid.hpp

   Copyright (c) Michael Strickland
   Forked at v2.0; Additions by Tim DuBois

   GNU General Public License (GPLv3)
   See detailed text in license directory

*/

#ifndef __grid_hpp__
#define __grid_hpp__

// this holds the current values of the wavefunction
extern dcomp ***w;

// this holds the updated values of the wavefunction
extern dcomp ***W;
// this holds the snapshots of the wavefunction
extern dcomp ****wstore;

// this is a temporary pointer used during array swap in copyDown
extern dcomp ***tmp;

// this holds the potential array
extern dcomp ***v;
extern dcomp ***vBase;

// these hold the alpha and beta arrays which are used during updates
extern dcomp ***a,***b;

// number of snapshots to capture in memory
extern int nsnaps;

//Energy penalty offset
extern double epsilon;

// return energy of the passed wavefnc
dcomp wfncEnergy(dcomp*** wfnc);

// return norm squared of the passed wavefnc
dcomp wfncNorm2(dcomp*** wfnc);

// return expectation value of vInfinity
dcomp vInfinityExpectationValue(dcomp*** wfnc);

// return expectation value of r^2     #ad.
dcomp r2ExpectationValue(dcomp*** wfnc);

// return energy of the current wavefnc
dcomp computeEnergy();

// other methods
void allocateMemory();
void allocateClusterMemory();
void deallocateClusterMemory();
void deallocateMemory();
void updateBoundaries(double eps);
void updateInterior(double eps);
void copyDown();
void loadPotentialArrays();
void potentialFromFile();
void updatePotential(dcomp beta);
void normalizeWavefunction(dcomp*** wfnc);
void storeConverged(dcomp*** wfnc, int num);

#endif /* __grid_hpp__ */
