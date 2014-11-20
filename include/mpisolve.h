/*

   mpisolve.h

   Copyright (c) Michael Strickland

   GNU General Public License (GPLv3)
   See detailed text in license directory

*/

#ifndef __mpisolve_h__
#define __mpisolve_h__

#include <iostream>
#include <iomanip>
#include <complex>
#include <fstream>
#include <limits>
#include <climits>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <cerrno>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/times.h>
#include <sys/types.h>
#include <unistd.h>

using namespace std;

typedef complex<double> dcomp;

extern "C"
{
    
  void dgemm_(char*,char*,int*,int*,int*,double*,double*,int*,
                         double*,int*,double*,double*,int*);
}

extern int DISTNUMZ,NUMX,NUMY,NUMZ,UPDATE,SNAPUPDATE,POTENTIAL,INITCONDTYPE,INITSYMMETRY,NF,SAVEWAVEFNCS,RUNTYPE,CLUSTSIZE,OUTPOT,EXCITEDSTATES;
extern double A,STEPS,EPS,MINTSTEP,SIG,MASS,T,TC,SIGMA,XI,TOLERANCE,ALX,ALY,ALZ,BOXSIZE;

extern int nodeID, numNodes, debug, waveNum;

extern dcomp energyCollect;
extern dcomp normalizationCollect;
extern dcomp vInfinityCollect;
extern dcomp rRMS2Collect;    // #ad.

extern double *clustSpecies;
extern double *clust;

/* debug values */
#define DEBUG_OFF		0
#define DEBUG_ON		1
#define DEBUG_FULL		2

/* message tags */
#define HELLO			1
#define DONE			2
#define SYNC_LEFT		3
#define SYNC_LEFT_MESSAGE	4
#define SYNC_RIGHT		5
#define SYNC_RIGHT_MESSAGE	6
#define NANERROR		7

// the main solve routines
void solveInitialize();
void solve();
void reInitSolver();
void solveRestart();
void solveFinalize();
void evolve(int);
void findExcitedStates();

// boundary sync stuff
void syncBoundaries(dcomp*** wfnc);
void sendLeftBoundary(dcomp*** wfnc);
void sendRightBoundary(dcomp*** wfnc);
void receiveLeftBoundary();
void receiveRightBoundary();
void loadRightBoundaryFromBuffer(dcomp*** wfnc);
void loadLeftBoundaryFromBuffer(dcomp*** wfnc);

#endif /* __mpisolve_h__ */
