/*

   initialconditions.cpp

   Copyright (c) Michael Strickland
   Forked at v2.0; Additions by Tim DuBois

   GNU General Public License (GPLv3)
   See detailed text in license directory

*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <climits>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <ctime>
#include <cerrno>
#include <sys/stat.h>
#include <sys/time.h>
#include <string>
#include <vector>
#include <sstream>
#include <complex>

using namespace std;

#include "wafer.hpp"
#include "grid.hpp"
#include "initialconditions.hpp"
#include "random.hpp"

// initializes the variables 
void setInitialConditions(int seedMult)
{
  	double sig=SIG; // standard deviation
	int sx,sy,sz,tx,ty,tz,ii,oldnumy,minz,maxz,olddnumz,tmp;
	double dx,dy,dz,r,costheta,cosphi,temp,temp2;
	fstream input; //, debug_out;
	char fname[32];
	string line;
        string delim = "\t";
        string token = "";
	vector<string> lines;
	int fileSize,stridein=1,strideout=1,linenumber;
	
	
	//cout << "==> Initializing variables\n";

	//srand(((unsigned int)time(0))*seedMult);
	srand(19691123);

	switch (INITCONDTYPE) { 
	  case 0: 
		// read from file
	    sprintf(fname,"data/wavefunction_%d_%d.dat",WAVENUM,nodeID);
	    input.open(fname, ios::in);
	    if (nodeID==1) cout << "==> Initial wavefunction (state " << WAVENUM << ") : From file" << endl;
        //NOTE: numNodes must be the same for this run and the last! 
        //TODO: Check for this if at all possible.
	    if (!input) {
	        cout << "==> Error : Unable to open wavefunction file " << fname << ". Using random Gaussian instead." << endl;
		for (sx=0;sx<=NUMX+5;sx++)
	            for (sy=0;sy<=NUMY+5;sy++)
		        for (sz=0; sz<=DISTNUMZ+5;sz++)
                            w[sx][sy][sz] = dcomp(randGauss(sig),0.);
	    }
	    while( getline( input, line ) ) lines.push_back( line ) ;
       
        //sprintf(fname,"debug/input_%d.txt",nodeID);
        //debug_out.open(fname, ios::out);
       
        //Find parameters from input file so that readin works
		fileSize = lines.size();
        oldnumy = 0;
        minz = 100000;
        maxz = 0;
        for (ii=0; ii<fileSize; ii++) {
            //iterate through and find NUMY and DISTNUMZ of file
            line = lines.at(ii);
            
            line.erase(0, line.find(delim) + delim.length()); //remove numx
            
            token = line.substr(0, line.find(delim)); //numy
            tmp = atoi( token.c_str() ); //numy to int
            if (tmp > oldnumy) oldnumy = tmp;
            
            line.erase(0, line.find(delim) + delim.length()); //remove numy

            token = line.substr(0, line.find(delim)); //numz
            tmp = atoi( token.c_str() ); //numz to int

            if (tmp > maxz) maxz = tmp;
            if (tmp < minz) minz = tmp;
        }
       
        olddnumz = maxz-minz+1;

		if (nodeID==1) cout << "Previous run parameters: NUMY = " << oldnumy << ", DISTNUMZ = " << olddnumz << endl;

        if (DISTNUMZ > olddnumz) strideout = DISTNUMZ/olddnumz;
		if (DISTNUMZ < olddnumz) stridein = olddnumz/DISTNUMZ;
         
		for (sx=3;sx<=2+NUMX;sx++)
			for (sy=3;sy<=2+NUMY;sy++)
				for (sz=3; sz<=2+DISTNUMZ;sz++) {
					//if (debug && nodeID==1) cout << "Mark : " << sx << ", " << sy << ", " << sz << endl;
			        if (strideout==1 && strideout==1) {
						linenumber  = (sx-3)*oldnumy*olddnumz + (sy-3)*olddnumz + (sz-3);
					}					
			        if (strideout>1) { 
						// If input wavefunction has lower resolution, spread it out
						tx = ceil((sx-2)/((double)strideout));
						ty = ceil((sy-2)/((double)strideout));
						tz = ceil((sz-2)/((double)strideout));
						linenumber  = (tx-1)*oldnumy*olddnumz + (ty-1)*olddnumz + tz;
						//if (debug && nodeID==1) cout << "Respond : " << tx << ", " << ty << ", " << tz << endl;
					}
			        if (stridein>1) {
						// If input wavefunction has higher resolution, sample it						
						tx = (sx-2)*stridein;
						ty = (sy-2)*stridein;
						tz = (sz-2)*stridein;
						linenumber  = (tx-1)*oldnumy*olddnumz + (ty-1)*olddnumz + tz;
						//if (debug && nodeID==1) cout << "Respond : " << tx << ", " << ty << ", " << tz << endl;
					}					
					line = lines.at(linenumber);
                    //debug_out << "Line " << linenumber << ": " << line << endl;
					int space_index = line.find_last_of("\t");
					string linemod = line.substr(0,space_index-1);	
					int space_index2 = linemod.find_last_of("\t");
					std::istringstream stream;
					stream.str(line.substr(space_index,line.length()-space_index));
					stream >> temp; // imaginary part
					std::istringstream stream2;
					stream2.str(line.substr(space_index2,space_index-space_index2));
					stream2 >> temp2; // real part
                    //if (debug && nodeID==1) cout << "line #" << linenumber << " : " << line << endl;
					//if (debug && nodeID==1) cout << "real : " << temp2 << endl;
					//if (debug && nodeID==1) cout << "imag : " << temp << endl;
					w[sx][sy][sz] = dcomp(temp2,temp);
				}
        //debug_out.close();
		input.close();
		break;
	  case 1:
		// random
		if (nodeID==1) cout << "==> Initial wavefunction : Random" << endl;
        	for (sx=0;sx<=NUMX+5;sx++)
          	  for (sy=0;sy<=NUMY+5;sy++)
            	    for (sz=0; sz<=DISTNUMZ+5;sz++)
                      w[sx][sy][sz] = dcomp(randGauss(sig),0.);
		break;
	  case 2:
		// coulomb like
		if (nodeID==1) cout << "==> Initial wavefunction : Coulomb" << endl;
        	for (sx=0;sx<=NUMX+5;sx++)
          	  for (sy=0;sy<=NUMY+5;sy++)
            	    for (sz=0; sz<=DISTNUMZ+5;sz++) {
		      // coordinate system is centered in simulation volume
		      dx = sx - ((double)NUMX+5.)/2.;
		      dy = sy - ((double)NUMY+5.)/2.;
		      dz = sz - ((double)DISTNUMZ+5.)/2. + ( ((double)nodeID) - ((double)numNodes)/2. )*DISTNUMZ;
		      r = A*sqrt(dx*dx+dy*dy+dz*dz);
		      costheta = A*dz/r;
		      cosphi = A*dx/r;
                      w[sx][sy][sz]  = exp(-MASS*r);				// n=1
                      w[sx][sy][sz] += (2-MASS*r)*exp(-MASS*r/2);		// n=2,l=0
                      w[sx][sy][sz] += MASS*r*exp(-MASS*r/2)*costheta;		// n=2,l=1,m=0
                      w[sx][sy][sz] += MASS*r*exp(-MASS*r/2)*sqrt(1-costheta*costheta)*cosphi; //n=2,l=1,m=+-1
	    	    }
		break;
	  case 3:
		// constant
		if (nodeID==1) cout << "==> Initial wavefunction : Constant" << endl;
        	for (sx=0;sx<=NUMX+5;sx++)
          	  for (sy=0;sy<=NUMY+5;sy++)
            	    for (sz=0; sz<=DISTNUMZ+5;sz++)
                      w[sx][sy][sz] = 0.1;
		break;
	  case 4:
		// test grid
		if (nodeID==1) cout << "==> Initial wavefunction : Boolean Test" << endl;
        	for (sx=0;sx<=NUMX+5;sx++)
          	  for (sy=0;sy<=NUMY+5;sy++)
            	    for (sz=0; sz<=DISTNUMZ+5;sz++)
                      w[sx][sy][sz] = (sx%2)*(sy%2)*(sz%2);
		break;
	  default:
		cout << "Invalid initial condition type" << endl;
		exit(0);
		break;
	}

	// enforce BCs
        for (sx=0;sx<=NUMX+5;sx++)
          for (sy=0;sy<=NUMY+5;sy++) {
		w[sx][sy][0] = 0;
		w[sx][sy][1] = 0;
		w[sx][sy][2] = 0;
		w[sx][sy][2+DISTNUMZ+1] = 0;
		w[sx][sy][2+DISTNUMZ+2] = 0;
		w[sx][sy][2+DISTNUMZ+3] = 0;
	  }

        for (sz=0;sz<=DISTNUMZ+5;sz++)
          for (sy=0;sy<=NUMY+5;sy++) {
		w[0][sy][sz] = 0;
		w[1][sy][sz] = 0;
		w[2][sy][sz] = 0;
		w[2+NUMX+1][sy][sz] = 0;
		w[2+NUMX+2][sy][sz] = 0;
		w[2+NUMX+3][sy][sz] = 0;
	  }

        for (sx=0;sx<=NUMX+5;sx++)
          for (sz=0;sz<=DISTNUMZ+5;sz++) {
		w[sx][0][sz] = 0;
		w[sx][1][sz] = 0;
		w[sx][2][sz] = 0;
		w[sx][2+NUMY+1][sz] = 0;
		w[sx][2+NUMY+2][sz] = 0;
		w[sx][2+NUMY+3][sz] = 0;
	  }

	// zero out updated wavefnc for safety's sake
        for (sx=0;sx<=NUMX+5;sx++)
          for (sy=0;sy<=NUMY+5;sy++)
            for (sz=0;sz<=DISTNUMZ+5;sz++)
		W[sx][sy][sz] = 0;

	// symmetrize the intial condition
	symmetrizeWavefunction();
}

void symmetrizeWavefunction()
{
	int sx,sy,sz,x,y,z,sign=1;

	// set sign for inversions; convention is that odd numbers 
	// are symmetric and even numbers antisymmetric 
	sign = 2*(INITSYMMETRY%2)-1;

	switch (INITSYMMETRY) { 
	  case 0:
		// no symmetry
		break;
	  case 1:
		// symmetric in z
	  case 2:
		// antisymmetric in z
        	for (sx=0;sx<=NUMX+5;sx++)
		{
		  x=sx; 	
          	  for (sy=3;sy<=2+NUMY;sy++)
		  {
		    y=sy; 	
            	    for (sz=3;sz<=2+DISTNUMZ;sz++)
		    {
			z=sz;
		    	if (z>DISTNUMZ/2)
		      	  z = DISTNUMZ + 1 - z;
			if (sz>DISTNUMZ/2)
			  w[sx][sy][sz] = ((dcomp)sign)*w[x][y][z];
		    }
		  }
		}
		break;
	  case 3:
		// symmetric in y
	  case 4:
		// antisymmetric in y
        	for (sx=0;sx<=NUMX+5;sx++)
		{
		  x=sx; 	
          	  for (sy=3;sy<=2+NUMY;sy++)
		  {
		    y=sy; 	
		    if (y>NUMY/2)
		      y = NUMY + 1 - y;
            	    for (sz=3;sz<=2+DISTNUMZ;sz++)
		    {
			z=sz;
			if (sy>NUMY/2)
			  w[sx][sy][sz] = ((dcomp)sign)*w[x][y][z];
		    }
		  }
		}
		break;
	  default:
		cout << "Invalid symmetry type" << endl;
		exit(0);
		break;
	}
}


void readWavefunctionBinary(int waveNum) {
  //waveNum is the state that needs to be read
  //size of wavefunc is assumed to be in line with current params
  int tx, ty, z;
  int sx, sy, sz;
  fstream input; //, debug_out;
  char fname[255];
  double tmpre, tmpim;

  //sprintf(fname,"debug/input_%d.txt",nodeID);
  //debug_out.open(fname, ios::out);
  // input full 3d wfnc
  sprintf(fname,"data/wavefunction_%d_%d.dat",waveNum,nodeID);
  input.open(fname, ios::in|ios::binary);

  for (int sx=3;sx<=2+NUMX;sx++) {
    for (int sy=3;sy<=2+NUMY;sy++) {
      for (int sz=3; sz<=2+DISTNUMZ;sz++) {

                input.read((char*)&tx, sizeof(int));
                input.read((char*)&ty, sizeof(int));
                input.read((char*)&z, sizeof(int));
                input.read((char*)&tmpre, sizeof(double));
                input.read((char*)&tmpim, sizeof(double));
                //debug_out << endl << tx << "," << ty << "," << z << ";" << tmpre << "," << tmpim << endl;
                wstore[waveNum][sx][sy][sz] = dcomp(tmpre,tmpim);
  }}}
  input.close();
  //debug_out.close();

  return;

}

