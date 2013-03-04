/*

   ouputroutines.cpp

   Copyright (c) Michael Strickland

   GNU General Public License (GPLv3)
   See detailed text in license directory

*/

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <ctime>
#include <complex>

using namespace std;

#include "mpisolve.h"
#include "grid.h"
#include "outputroutines.h"
#include "potential.h"

typedef std::numeric_limits< double > dbl;

void outputMeasurements(const double time) {

	dcomp ener = energyCollect/normalizationCollect;
	dcomp binding = ener - vInfinityCollect/normalizationCollect;
	dcomp rRMS2 = rRMS2Collect/normalizationCollect;   // #ad.

	// output to screen

	cout.precision(12);
	cout.width(dwidth); cout << time;
	cout.width(dwidth); cout << setprecision (7) << ener;
	cout.width(dwidth); cout << setprecision (7) << binding;
	cout.width(dwidth); cout << setprecision (7) << sqrt(real(rRMS2));   // #ad.
	cout << endl;

	// output to files

	//energy_out << time;
	//energy_out << "\t" << ener;
	//energy_out << endl;
}

void outputSummaryData() {

      dcomp ener = energyCollect/normalizationCollect;
      dcomp binding = ener - vInfinityCollect/normalizationCollect;
      dcomp rRMS2 = rRMS2Collect/normalizationCollect;   // #ad.

      print_line();

      //cout << "==> Total energy : " << energyCollect << endl;
      //cout << "==> Normalization2 : " << normalizationCollect << endl;
      cout.precision(dbl::digits10);
      cout << "==> Ground State Energy : " << fixed << ener << endl;
      cout << "==> Ground State Binding Energy : " << binding << endl;
      cout << "==> Ground State r_RMS : " << sqrt(real(rRMS2)) << endl;  // #ad.
      cout << "==> Ground State L/r_RMS : " << float(NUM)/sqrt(real(rRMS2)) << endl;  // #ad.

}

void outputSnapshot(dcomp ***wfnc, char* label) {

  int x;
  static int h=NUM/2;
  static int hx=NUMX/2;

  fstream out;
  char fname[255];

  // dump wavefunction

  // output slices suitable for 2d viewing
  sprintf(fname,"data/snapshot/wavefunction_%s.dat",label);
  out.open(fname, ios::out);
  out.precision(10);
  for (int s=0;s<=NUMX+1;s++) {
    x=(nodeID-1)*NUMX + s;	  
    out << x << "\t";
    out << scientific << 0.5*(wfnc[s][h][h]+wfnc[s][h+1][h+1]) << "\t";
    out << endl;
  }
  out << "&&" << endl;
  for (int s=0;s<=NUM+1;s++) {
    out << s << "\t";
    out << scientific << 0.5*(wfnc[hx][s][h]+wfnc[hx+1][s][h+1]) << "\t";
    out << endl;
  }
  out << "&&" << endl;
  for (int s=0;s<=NUM+1;s++) {
    out << s << "\t";
    out << scientific << 0.5*(wfnc[hx][h][s]+wfnc[hx+1][h+1][s]) << "\t";
    out << endl;
  }
  out.close();
  
  return;
}

void outputWavefunction(dcomp ***wfnc, char* label) {

  int x;
  fstream out;
  char fname[255];

  // output full 3d wfnc
  sprintf(fname,"data/wavefunction_%s.dat",label);

  cout << "==> Dumping wave function to " << fname << endl;

  out.open(fname, ios::out);
  out.precision(12);
  for (int sx=1;sx<=NUMX;sx++) {
    x=(nodeID-1)*NUMX + sx;
    for (int sy=1;sy<=NUM;sy++) {
      for (int sz=1; sz<=NUM;sz++) {
                out << x  << "\t";
                out << sy << "\t";
                out << sz << "\t";
                out << real(wfnc[sx][sy][sz]) << "\t";
                out << imag(wfnc[sx][sy][sz]);
                out << endl;
  }}}
  out.close();

  return;

}

// output v 3d
void outputPotential(char* label) {

  int x;
  fstream out;
  char fname[255];

  // output full 3d wfnc
  sprintf(fname,"data/potential_%s.dat",label);

  cout << "==> Dumping potential to " << fname << endl;

  out.open(fname, ios::out);
  out.precision(12);
  for (int sx=1;sx<=NUMX;sx++) {
    x=(nodeID-1)*NUMX + sx;
    for (int sy=1;sy<=NUM;sy++) {
      for (int sz=1; sz<=NUM;sz++) {
                out << x  << "\t";
                out << sy << "\t";
                out << sz << "\t";
                out << real(v[sx][sy][sz]) << "\t";
                out << imag(v[sx][sy][sz]);
                out << endl;
  }}}
  out.close();

  return;

}

// output v along principal axes
void dumpPotential() {

  int h=NUM/2;
  fstream out;

  out.open("data/potential.dat", ios::out);
  for (int s=0;s<=NUM+1;s++) {
                out << s << "\t";
                out << v[s][h][h] << "\t";
                out << v[h][s][h] << "\t";
                out << v[h][h][s] << "\t";
                out << endl;
  }
  out.close();
  return;

}

void print_line() {
        for (int i=0;i<4*dwidth;i++) cout << "-"; cout << endl;
        return;
}

