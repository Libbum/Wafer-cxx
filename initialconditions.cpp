/*

   initialconditions.cpp

   Copyright (c) Michael Strickland

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

#include "mpisolve.h"
#include "grid.h"
#include "initialconditions.h"
#include "random.h"

// initializes the variables 
void setInitialConditions(int seedMult)
{
  	double sig=SIG; // standard deviation
	int sx,sy,sz,tx,ty,tz,ii,oldnumy,olddnumz,idx1,idx2,idx3,tmp;
	double dx,dy,dz,r,costheta,cosphi,temp,temp2;
	fstream input;
    fstream debug_out;
	char fname[32];
	string line;
	vector<string> lines;
	int inputLatticeSize,fileSize,stridein=1,strideout=1,linenumber;
	
	
	//cout << "==> Initializing variables\n";

	//srand(((unsigned int)time(0))*seedMult);
	srand(19691123);

	switch (INITCONDTYPE) { 
	  case 0: 
		// read from file
	    sprintf(fname,"data/wavefunction_0_%d.dat",nodeID);
	    input.open(fname, ios::in);
		if (nodeID==1) cout << "==> Initial wavefunction : From file" << endl;
        //NOTE: numNodes must be the same for this run and the last! 
        //TODO: Check for this if at all possible.
		if (!input) {
			cout << "==> Error : Unable to open wavefunction file " << fname << ". Using random Gaussian instead." << endl;
			for (sx=0;sx<NUMX+2;sx++)
				for (sy=0;sy<NUMY+2;sy++)
					for (sz=0; sz<DISTNUMZ+2;sz++)
                        w[sx][sy][sz] = dcomp(randGauss(sig),0.);
		}
		while( getline( input, line ) ) lines.push_back( line ) ;
		
        //OK, so this needs to be re-written. 
        //Lattice size per node
        inputLatticeSize = NUMX*NUMY*(NUMZ/(numNodes-1));  //round(pow((numNodes-1)*lines.size(),1/3.));
        //input data per node
		fileSize = lines.size();
        oldnumy = 0;
        olddnumz = 0;
        for (ii=0; ii<fileSize; ii++) {
            //iterate through and find NUMY and DISTNUMZ of file
            line = lines.at(ii);
            idx1 = line.find_first_of("\t"); //one before numy
            idx2 = line.substr(idx1+1,line.length()).find_first_of("\t"); //one before numz
            tmp = atoi( line.substr(idx1+1,idx1+idx2-1).c_str() ); //should be numy
            if (tmp > oldnumy) {
                oldnumy = tmp;
            }

            idx3 = line.substr(idx1+idx2+2,line.length()).find_first_of("\t"); //one after numz
            tmp = atoi( line.substr(idx1+idx2+2,idx1+idx2+idx3-1).c_str() ); //should be numz
        }
        sprintf(fname,"debug/debug_%d.txt",nodeID);
        debug_out.open(fname, ios::out);
        debug_out << "numy: " << oldnumy << ", numzrow = " << tmp << endl;
        if (inputLatticeSize > fileSize) strideout = inputLatticeSize/fileSize;
		if (inputLatticeSize < fileSize) stridein = fileSize/inputLatticeSize;
         
        //debug_out << "strideout: " << strideout << ", stridein: " << stridein << endl;
		for (sx=1;sx<=NUMX;sx++)
			for (sy=1;sy<=NUMY;sy++)
				for (sz=1; sz<=DISTNUMZ;sz++) {
					//if (debug && nodeID==1) cout << "Mark : " << sx << ", " << sy << ", " << sz << endl;
			        if (strideout==1 && strideout==1) {
						linenumber  = (sx-1)*inputLatticeSize*inputLatticeSize + (sy-1)*inputLatticeSize + (sz-1);
					}					
			        if (strideout>1) { 
						// If input wavefunction has lower resolution, spread it out
						tx = ceil(sx/((double)5.0));
						ty = ceil(sy/((double)5.0));
						tz = ceil(sz/((double)5.0));
						linenumber  = (tx-1)*143*2 + (ty-1)*2 + tz;
					    //if (linenumber > fileSize) {
                        //debug_out << linenumber << ", " << sx << ", " << sy << ", " << sz << "; " << tx << ", " << ty << ", " << tz << endl;
                        //}
						//tx = ceil(sx/((double)strideout));
						//ty = ceil(sy/((double)strideout));
						//tz = ceil(sz/((double)strideout));
						//linenumber  = (tx-1)*inputLatticeSize*inputLatticeSize + (ty-1)*inputLatticeSize + (tz-1);
						//if (debug && nodeID==1) cout << "Respond : " << tx << ", " << ty << ", " << tz << endl;
					}
			        if (stridein>1) {
						// If input wavefunction has higher resolution, sample it						
						tx = sx*stridein;
						ty = sy*stridein;
						tz = sz*stridein;
						linenumber  = (tx-1)*inputLatticeSize*inputLatticeSize + (ty-1)*inputLatticeSize + (tz-1);
						//if (debug && nodeID==1) cout << "Respond : " << tx << ", " << ty << ", " << tz << endl;
					}					
					line = lines.at(linenumber);
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
					//debug_out << "line #" << linenumber << " : " << line << endl;
					//debug_out << "real : " << temp2 << endl;
					//debug_out << "imag : " << temp << endl;
					w[sx][sy][sz] = dcomp(temp2,temp);
				}
        debug_out.close();  
		input.close();
		break;
	  case 1:
		// random
		if (nodeID==1) cout << "==> Initial wavefunction : Random" << endl;
        	for (sx=0;sx<NUMX+2;sx++)
          	  for (sy=0;sy<NUMY+2;sy++)
            	    for (sz=0; sz<DISTNUMZ+2;sz++)
                      w[sx][sy][sz] = dcomp(randGauss(sig),0.);
		break;
	  case 2:
		// coulomb like
		if (nodeID==1) cout << "==> Initial wavefunction : Coulomb" << endl;
        	for (sx=0;sx<=NUMX+1;sx++)
          	  for (sy=0;sy<=NUMY+1;sy++)
            	    for (sz=0; sz<=DISTNUMZ+1;sz++) {
		      // coordinate system is centered in simulation volume
		      dx = sx - ((double)NUMX+1.)/2.;
		      dy = sy - ((double)NUMY+1.)/2.;
		      dz = sz - ((double)DISTNUMZ+1.)/2. + ( ((double)nodeID) - ((double)numNodes)/2. )*DISTNUMZ;
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
        	for (sx=0;sx<=NUMX+1;sx++)
          	  for (sy=0;sy<=NUMY+1;sy++)
            	    for (sz=0; sz<=DISTNUMZ+1;sz++)
                      w[sx][sy][sz] = 0.1;
		break;
	  case 4:
		// test grid
		if (nodeID==1) cout << "==> Initial wavefunction : Boolean Test" << endl;
        	for (sx=0;sx<=NUMX+1;sx++)
          	  for (sy=0;sy<=NUMY+1;sy++)
            	    for (sz=0; sz<=DISTNUMZ+1;sz++)
                      w[sx][sy][sz] = (sx%2)*(sy%2)*(sz%2);
		break;
	  default:
		cout << "Invalid initial condition type" << endl;
		exit(0);
		break;
	}

	// enforce BCs
        for (sx=0;sx<=NUMX+1;sx++)
          for (sy=0;sy<=NUMY+1;sy++) {
		w[sx][sy][0] = 0;
		w[sx][sy][DISTNUMZ+1] = 0;
	  }

        for (sz=0;sz<=DISTNUMZ+1;sz++)
          for (sy=0;sy<=NUMY+1;sy++) {
		w[0][sy][sz] = 0;
		w[NUMX+1][sy][sz] = 0;
	  }

        for (sx=0;sx<=NUMX+1;sx++)
          for (sz=0;sz<=DISTNUMZ+1;sz++) {
		w[sx][0][sz] = 0;
		w[sx][NUMY+1][sz] = 0;
	  }

	// zero out updated wavefnc for safety's sake
        for (sx=0;sx<NUMX+2;sx++)
          for (sy=0;sy<NUMY+2;sy++)
            for (sz=0;sz<DISTNUMZ+2;sz++)
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
        	for (sx=0;sx<NUMX+2;sx++)
		{
		  x=sx; 	
          	  for (sy=1;sy<=NUMY;sy++)
		  {
		    y=sy; 	
            	    for (sz=1;sz<=DISTNUMZ;sz++)
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
        	for (sx=0;sx<NUMX+2;sx++)
		{
		  x=sx; 	
          	  for (sy=1;sy<=NUMY;sy++)
		  {
		    y=sy; 	
		    if (y>NUMY/2)
		      y = NUMY + 1 - y;
            	    for (sz=1;sz<=DISTNUMZ;sz++)
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
