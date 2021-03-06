/*
 *

 ouputroutines.cpp

 Copyright (c) Michael Strickland
 Forked at v2.0; Additions by Tim DuBois

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

#include "wafer.hpp"
#include "grid.hpp"
#include "outputroutines.hpp"
#include "potential.hpp"

typedef std::numeric_limits< double > dbl;

void outputMeasurements(const double time, dcomp laste) {

    dcomp ener = energyCollect/normalizationCollect;
    dcomp rRMS2 = rRMS2Collect/normalizationCollect;

    // output to screen

    cout.precision(12);
    cout.width(20); cout << time;
    if (POTENTIAL==22) {
        cout.width(35); cout << setprecision (dbl::digits10) << ener*1e6/239.2311;
    } else {
        cout.width(35); cout << setprecision (dbl::digits10) << ener;
    }
    cout.width(15); cout << setprecision (7) << sqrt(real(rRMS2));
    cout.width(25); cout << setprecision (7) << abs(real(laste)-real(ener))*1e6/239.2311; //ueV
    cout << endl;

}

void outputSummaryData(int WAVENUM) {

    dcomp ener = energyCollect/normalizationCollect;
    dcomp binding = ener - vInfinityCollect/normalizationCollect;
    dcomp rRMS2 = rRMS2Collect/normalizationCollect;   // #ad.

    print_line();

    //cout << "==> Total energy : " << energyCollect << endl;
    //cout << "==> Normalization2 : " << normalizationCollect << endl;
    if (POTENTIAL==22) {
        string state;
        switch (WAVENUM) {
            case 1:
                state = "First Excited ";
                break;
            case 2:
                state = "Second Excited ";
                break;
            case 3:
                state = "Third Excited ";
                break;
            case 4:
                state = "Fourth Excited ";
                break;
            case 5:
                state = "Fifth Excited ";
                break;
            default:
                state = "Ground ";
                break;
        }
        cout << "==> " << state << "State Energy : " << setprecision (dbl::digits10) << ener*1e6/239.2311 << " (ueV)" << endl;
        cout << "==> " << state << "State Binding Energy : " << binding*1e6/239.2311 << " (ueV)" << endl;
    } else {
        cout << "==> Ground State Energy.. : " << fixed << ener << endl;
        cout << "==> Ground State Binding Energy : " << binding << endl;
    }
    cout << "==> r_RMS : " << sqrt(real(rRMS2)) << endl;  // #ad.
    cout << "==> L/r_RMS : " << float(NUMX)/sqrt(real(rRMS2)) << endl;  // #ad.

}

void outputSnapshot(dcomp ***wfnc, char* label) {

    int z;
    static int hx=NUMX/2;
    static int hy=NUMY/2;
    static int hz=DISTNUMZ/2;

    fstream out;
    char fname[255];

    // dump wavefunction

    // output slices suitable for 2d viewing
    sprintf(fname,"data/snapshot/wavefunction_%s.dat",label);
    out.open(fname, ios::out);
    out.precision(10);
    for (int s=0;s<=NUMX+5;s++) {
        out << s << "\t";
        out << scientific << 0.5*(wfnc[s][hy][hz]+wfnc[s][hy+1][hz+1]) << "\t";
        out << endl;
    }
    out << "&&" << endl;
    for (int s=0;s<=NUMY+5;s++) {
        out << s << "\t";
        out << scientific << 0.5*(wfnc[hx][s][hz]+wfnc[hx+1][s][hz+1]) << "\t";
        out << endl;
    }
    out << "&&" << endl;
    for (int s=0;s<=DISTNUMZ+5;s++) {
        z=(nodeID-1)*DISTNUMZ + s;
        out << z << "\t";
        out << scientific << 0.5*(wfnc[hx][hy][s]+wfnc[hx+1][hy+1][s]) << "\t";
        out << endl;
    }
    out.close();

    return;
}

void outputWavefunction(dcomp ***wfnc, char* label) {

    int z;
    fstream out;
    char fname[255];

    // output full 3d wfnc
    sprintf(fname,"data/wavefunction_%s.dat",label);

    cout << "==> Dumping wave function to " << fname << endl;

    out.open(fname, ios::out);
    out.precision(12);
    for (int sx=3;sx<=2+NUMX;sx++) {
        for (int sy=3;sy<=2+NUMY;sy++) {
            for (int sz=3; sz<=2+DISTNUMZ;sz++) {
                z=(nodeID-1)*DISTNUMZ + sz;
                out << sx-2  << "\t";
                out << sy-2 << "\t";
                out << z-2 << "\t";
                out << real(wfnc[sx][sy][sz]) << "\t";
                out << imag(wfnc[sx][sy][sz]);
                out << endl;
            }}}
    out.close();

    return;

}

void outputWavefunctionBinary(dcomp ***wfnc, char* label) {

    int tx, ty, z;
    fstream out;
    char fname[255];
    double tmp;

    // output full 3d wfnc
    sprintf(fname,"data/wavefunction_%s.dat",label);

    cout << "==> Dumping wave function to " << fname << endl;

    out.open(fname, ios::out|ios::binary);
    out.precision(12);
    for (int sx=3;sx<=2+NUMX;sx++) {
        for (int sy=3;sy<=2+NUMY;sy++) {
            for (int sz=3; sz<=2+DISTNUMZ;sz++) {
                tx = sx-2;
                ty = sy-2;
                z=(nodeID-1)*DISTNUMZ + sz-2;

                out.write(reinterpret_cast<const char*>(&tx), sizeof(int));
                out.write(reinterpret_cast<const char*>(&ty), sizeof(int));
                out.write(reinterpret_cast<const char*>(&z), sizeof(int));
                tmp = real(wfnc[sx][sy][sz]);
                out.write(reinterpret_cast<const char*>(&tmp), sizeof(double));
                tmp = imag(wfnc[sx][sy][sz]);
                out.write(reinterpret_cast<const char*>(&tmp), sizeof(double));
            }}}
    out.close();

    return;

}

// output v 3d
void outputPotential(char* label) {

    int z;
    fstream out;
    char fname[255];

    double convert = 1e6/239.2311; //For POTENTIAL == 22

    // output full 3d wfnc
    sprintf(fname,"data/potential_%s.dat",label);

    if (POTENTIAL==22) {
        cout << "==> Dumping potential to " << fname << " (ueV)" << endl;
    } else {
        cout << "==> Dumping potential to " << fname << endl;
    }

    out.open(fname, ios::out);
    out.precision(12);
    for (int sx=3;sx<=2+NUMX;sx++) {
        for (int sy=3;sy<=2+NUMY;sy++) {
            for (int sz=3; sz<=2+DISTNUMZ;sz++) {
                z=(nodeID-1)*DISTNUMZ + sz;
                out << sx-2  << "\t";
                out << sy-2 << "\t";
                out << z-2 << "\t";
                if (POTENTIAL==22) {
                    out << real(v[sx][sy][sz])*convert << "\t";
                    out << imag(v[sx][sy][sz])*convert;
                } else {
                    out << real(v[sx][sy][sz]) << "\t";
                    out << imag(v[sx][sy][sz]);
                }
                out << endl;
            }}}
    out.close();

    return;

}

// output v 3d
void outputPotentialBinary(char* label) {

    int tx, ty, z;
    fstream out;
    char fname[255];
    double convert = 1e6/239.2311; //For POTENTIAL == 22
    double tmp;
    // output full 3d wfnc
    sprintf(fname,"data/potential_%s.dat",label);

    if (POTENTIAL==22) {
        cout << "==> Dumping potential to " << fname << " (ueV)" << endl;
    } else {
        cout << "==> Dumping potential to " << fname << endl;
    }

    out.open(fname, ios::out|ios::binary);
    out.precision(12);
    for (int sx=3;sx<=2+NUMX;sx++) {
        for (int sy=3;sy<=2+NUMY;sy++) {
            for (int sz=3; sz<=2+DISTNUMZ;sz++) {
                tx = sx-2;
                ty = sy-2;
                z=(nodeID-1)*DISTNUMZ + sz-2;
                out.write(reinterpret_cast<const char*>(&tx), sizeof(int));
                out.write(reinterpret_cast<const char*>(&ty), sizeof(int));
                out.write(reinterpret_cast<const char*>(&z), sizeof(int));
                if (POTENTIAL==22) {
                    tmp = real(v[sx][sy][sz])*convert;
                    out.write(reinterpret_cast<const char*>(&tmp), sizeof(double));
                    tmp = imag(v[sx][sy][sz])*convert;
                    out.write(reinterpret_cast<const char*>(&tmp), sizeof(double));
                } else {
                    tmp = real(v[sx][sy][sz]);
                    out.write(reinterpret_cast<const char*>(&tmp), sizeof(double));
                    tmp = imag(v[sx][sy][sz]);
                    out.write(reinterpret_cast<const char*>(&tmp), sizeof(double));
                }
            }}}
    out.close();

    return;

}


// output v along principal axes
void dumpPotential() {
    //WARNING: This is probably no good for arbitrary box. Currently not in use so it's not checked.
    int h=NUMZ/2;
    fstream out;

    out.open("data/potential.dat", ios::out);
    for (int s=0;s<=NUMZ+1;s++) {
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

