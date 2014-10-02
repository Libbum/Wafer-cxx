/*
 
 Parallelized Finite Difference Time Domain Schrodinger Eq Solver
 
 mpisolve.c
 
 Copyright (c) Michael Strickland
 Forked at v2.0; Additions by Tim DuBois

 GNU General Public License (GPLv3)
 See detailed text in license directory
 
 */

#include "mpi.h" 
#include "mpisolve.h"
#include "grid.h"
#include "initialconditions.h"
#include "potential.h"
#include "outputroutines.h"
#include "paramreader.h"


typedef std::numeric_limits< double > dbl;

// these global vars are initialized from parameters file
// defaults set here are overridden by that file
int    DISTNUMZ=20,NUMX=20,NUMY=20,NUMZ=20,UPDATE=100,SNAPUPDATE=1000;
int    POTENTIAL=0,INITCONDTYPE=0,INITSYMMETRY=0,NF=2,SAVEWAVEFNCS=0,RUNTYPE=0,CLUSTSIZE=7,OUTPOT=0;
double  A=0.05,EPS=0.001,MINTSTEP=1.e-8,SIG=0.06,MASS=1.0,T=1.0,TC=0.2,SIGMA=1.0,XI=0.0,BOXSIZE=0.0,TOLERANCE=1.e-10,STEPS=40000;
double  ALX=4.7,ALY=4,ALZ=2.5788; //Aluminium Clusters & Grid Range

// species data if cluster is used
double *clustSpecies;

// cluster xyz coors if cluster is used
double *clust;

// mpi vars
int nodeID,numNodes;

// files
fstream debug_out;

// debug flag; options are DEBUG_{OFF,ON,FULL}
int debug = DEBUG_FULL;

// used for MPI non-blocking sends
double *leftSendBuffer,*rightSendBuffer;
MPI_Status leftSendStatus,rightSendStatus;
MPI_Status leftSendMessageStatus,rightSendMessageStatus;
MPI_Request leftSend,rightSend;
MPI_Request leftMessageSend,rightMessageSend;

// used for MPI non-blocking receives
double *leftReceiveBuffer,*rightReceiveBuffer;
MPI_Status leftReceiveStatus,rightReceiveStatus;
MPI_Status leftReceiveMessageStatus,rightReceiveMessageStatus;
MPI_Request leftReceive,rightReceive;
MPI_Request leftMessageReceive,rightMessageReceive;


// used for worker only intercommunication
MPI_Comm workers_comm;
MPI_Group workers_group;

// variables which will be loaded and reduced across nodes
dcomp energy=0;		// the local node energy
dcomp energyCollect=0;		// the total energy
dcomp normalization=0;		// the local node normalization squared
dcomp normalizationCollect=0;  // the total normalization squared
dcomp vInfinity=0;		// the local node expectation value of v_infty
dcomp vInfinityCollect=0;      // the total expectation value of v_infty
dcomp rRMS2=0;                 // the local node <r^2>  #ad.
dcomp rRMS2Collect=0;		// the total <r^2>  #ad.
int nanErrorCollect=0;

// #ad.: grnd-state energy and final time saved in global var after convergence
dcomp	EGrnd, EOne, timef;
int step = 0;

int main( int argc, char *argv[] ) 
{ 
    int done=0,checksum=0;
    char message[64]; 
    char fname[32];
    int exclude[1];
    struct tms starttime,endtime;
	
    MPI_Init(&argc,&argv); 
    MPI_Comm_size(MPI_COMM_WORLD,&numNodes); 
    MPI_Comm_rank(MPI_COMM_WORLD,&nodeID); 
    MPI_Status status;
    MPI_Group all_group;
   
    //char hostname[256];
    //gethostname(hostname,255);
    //cout << "Starting node " << nodeID << " of " << numNodes << " on " << hostname << endl;

    // setup group consisting of computational nodes only for internal communications
    MPI_Comm_group(MPI_COMM_WORLD, &all_group);
    exclude[0]=0;
    MPI_Group_excl(all_group, 1, exclude, &workers_group);
    MPI_Comm_create(MPI_COMM_WORLD, workers_group, &workers_comm);

    if (debug) {
		sprintf(fname,"debug/debug_%d.txt",nodeID);
		debug_out.open(fname, ios::out);
		debug_out << "==> Node " << nodeID << " is ready" << endl; 
    }
    // node 0 is the master
    if (nodeID == 0) {
        times(&starttime); // load start time into starttime structure
        print_line();
	cout << "Parameters from file" << endl;
	print_line();
	readParametersFromFile((char *)"params.txt",1);
	if (argc>1) {
		print_line();
		cout << "Parameters from commandline" << endl;
		print_line();
		readParametersFromCommandLine(argc,argv,1);
	}
        if (RUNTYPE == 1) {
           //Need to load cluster data and we really only want to do it once (per node).
           fstream input;
           string line;
           //just the size so we can allocate memory first 
           input.open("cluster.xyz", ios::in);
           if (!input) {
              cout << "==> Error : No cluster.xyz file present. Falling back to AL{X,Y,Z} values." << endl;
              RUNTYPE = 0;
           } else {
    
              getline(input,line);
              CLUSTSIZE = atoi(line.c_str())+1; //first line of XYZ file
              allocateClusterMemory();
              print_line();
              readClusterData((char *)"cluster.xyz", 1);
           }

           input.close();
           ALX = ceil((A*NUMX-A)/2);
           ALY = ceil((A*NUMY-A)/2);
           ALZ = ceil((A*NUMZ-A)/2);
           cout << "Using BoxSize of " << BOXSIZE << ". Grid values: NUMX = " << NUMX << ", NUMY = " << NUMY << ", (NUMZ = " << NUMZ << ")" << endl;
        } else if (RUNTYPE == 2) {
            //Forced NUMn values
             cout << "NUM values from params.txt: NUMX = " << NUMX << ", NUMY = " << NUMY << ", NUMZ = " << NUMZ << endl;
        } else {
            //Find NUMX and NUMY Values. These are overwritten in the next loop if it's invoked.
            if (ALX != 500) { //a value of 500 means unbound. Don't set the grid to separation length
                NUMX = ceil((2*ALX+A)/A);
            }
            if (ALY != 500) { //a value of 500 means unbound. Don't set the grid to separation length
                NUMY = ceil((2*ALY+A)/A);
            }

            if ((ALY == 500) || (ALX == 500)) {
                cout << "Unbound cluster. Grid values: NUMX = " << NUMX << ", NUMY = " << NUMY << ", (NUMZ = " << NUMZ << ")" << endl;
            } else {
                cout << "Calculated values for arbitrary grid: NUMX = " << NUMX << ", NUMY = " << NUMY << ", (NUMZ = " << NUMZ << ")" << endl;
            }
        }
    }
    else {
        readParametersFromFile((char *)"params.txt",0);
	    readParametersFromCommandLine(argc,argv,0);
        if (RUNTYPE == 1) {
           //Need to load cluster data and we really only want to do it once (per node).
           fstream input;
           string line;
           //just the size so we can allocate memory first 
           input.open("cluster.xyz", ios::in);
           if (!input) {
              RUNTYPE = 0;
           } else {
    
              getline(input,line);
              CLUSTSIZE = atoi(line.c_str())+1;
              allocateClusterMemory();
              readClusterData((char *)"cluster.xyz", 0);
           }

           input.close();
           ALX = ceil((A*NUMX-A)/2);
           ALY = ceil((A*NUMY-A)/2);
           ALZ = ceil((A*NUMZ-A)/2);
        } else if (RUNTYPE == 0) {
            //Find NUMX and NUMY Values. These are overwritten in the next loop if it's invoked.
            if (ALX != 500) { //a value of 500 means unbound. Don't set the grid to separation length
                NUMX = ceil((2*ALX+A)/A);
            }
            if (ALY != 500) { //a value of 500 means unbound. Don't set the grid to separation length
                NUMY = ceil((2*ALY+A)/A);
            }
        } //Do nothing for type 2
    }
	
    if (NUMZ%(numNodes-1)!=0) {
		if (nodeID==0)
			print_line();
        cout << "ERROR: Unable to partition lattice ... exiting" << endl;
        print_line();
		if (debug) {
			debug_out << "==> Goodbye from node " << nodeID << endl; 
			debug_out.close();
		}
        MPI_Finalize(); 
		exit(0);
    } else {
		DISTNUMZ = NUMZ/(numNodes-1);	
    }
	
    if (nodeID == 0) {
		
		//	
		// master node
		//	
		
		for (int node=1;node<numNodes;node++) {
			sprintf(message,"Hello to node %d",node); 
			MPI_Send(message, strlen(message), MPI_CHAR, node, HELLO, MPI_COMM_WORLD); 
		}
		
		// master loops and waits for children to report that they are ready to start
		checksum=0; 
		do {
			MPI_Recv(&done, 1, MPI_INT, MPI_ANY_SOURCE, HELLO, MPI_COMM_WORLD, &status); 
			checksum += done;
			if (debug) debug_out << "Received: hello from computational node" << endl;
		} while( checksum < numNodes-1 ); 
		
		// cluster is ready
		if (debug) debug_out << "==> Cluster ready" << endl;
	
                // Currently the master process does nothing.  
		// It simply starts and waits for the others
		
		// master loops and waits for children to report that they are done
		checksum=0;
		do {
			MPI_Recv(&done, 1, MPI_INT, MPI_ANY_SOURCE, DONE, MPI_COMM_WORLD, &status); 
			checksum += done;
                if (debug) debug_out << "Received: checkout from computational node" << endl;
			sleep(0.1); // sleep 0.1 seconds between checks in order to reduce CPU usage of master
		} while( checksum < numNodes-1 ); 
	        
		times(&endtime); // load end time into endtime structure
		cout << "==> User time: " << (endtime.tms_utime - starttime.tms_utime)/((double)sysconf(_SC_CLK_TCK)) << " seconds"<< endl;
		cout << "==> System time: " << (endtime.tms_stime - starttime.tms_stime)/((double)sysconf(_SC_CLK_TCK)) << " seconds" << endl;
		print_line();
		cout << "Done." << endl;
		print_line();
		
    } else {
		
		//	
		// computational node
		//	
		
		// set initial conditions and get ready for computation
		solveInitialize();

		// blocking receive to wait for master to fire up and say hello
		MPI_Recv(message, 20, MPI_CHAR, 0, HELLO, MPI_COMM_WORLD, &status); 
		if (debug == DEBUG_FULL) debug_out << "==> Received : "<< message << endl;

		// the master said hello so now let's get to work
		solve();
	        
        //If there's a nanError, extend EPS and resolve
        while (nanErrorCollect == 1) {
            solveRestart();

            nanErrorCollect = 0;
            MPI_Bcast(&nanErrorCollect, 1, MPI_INT, 0, workers_comm);
            if (EPS >= MINTSTEP) { //Don't loop forever
                solve();
            }
        }
                
		// done with main computation, now do any analysis required
		solveFinalize();
		
		// send message to master that we are done
		done = 1;
		MPI_Send( &done, 1, MPI_INT, 0, DONE, MPI_COMM_WORLD );
		
		cout.flush();
		
    }
	
    if (debug) { 
		debug_out << "==> Goodbye from node " << nodeID << endl; 
		debug_out.close();
    }
	
    MPI_Group_free(&workers_group);
    MPI_Finalize(); 
    return 0; 
} 

// solve initialize
void solveInitialize() {
	
		if (debug == DEBUG_FULL) debug_out << "Initialise" << endl;
    	// allocate memory
    	allocateMemory();
	
	// load the potential
	if (nodeID==1) {
		print_line();
		cout << "==> Loading Potential Arrays" << endl;
		flush(cout);
	}

		if (debug == DEBUG_FULL) debug_out << "Allocated memory" << endl;
		if (debug == DEBUG_FULL) debug_out << "NUMX = " << NUMX << ", NUMY = " << NUMY << endl;
    	loadPotentialArrays();
		if (debug == DEBUG_FULL) debug_out << "Loaded potential" << endl;
	if (OUTPOT) {
	   char label[64];
       sprintf(label,"0_%d",nodeID); 
       // output potential for debugging
       outputPotentialBinary(label);
    }
	
	if (nodeID==1) print_line();
	
	// set initial conditions
	setInitialConditions(nodeID+1);
	
	// output some summary information
	if (nodeID==1) { 
	print_line();
      	cout << "==> Number of computational nodes : " << numNodes-1 << endl; 
      	cout << "==> DISTNUMZ : " << DISTNUMZ << endl; 
      	print_line();
      	cout << "Spatial Step Size (A): " << A << endl;
      	cout << "Temporal Step Size (EPS): " << EPS << endl;
      	cout << "Standard Deviation of initial wavefunction noise (SIG): " << SIG << endl;
      	cout << "Box Dimensions: X = " << A*NUMX << ", Y = " << A*NUMY << ", Z = " << A*NUMZ << endl;
      	print_line();
      	cout.width(dwidth); cout << "Time";
      	cout.width(dwidth); cout << "Energy";
      	cout.width(dwidth); cout << "r_RMS";   
      	cout << endl;
      	print_line();
	}
}

void solveRestart() {
        double EPSold = EPS;
        EPS = EPSold / 2;

        if (nodeID==1) { 
            print_line();
            cout << "Temporal Step Size (EPS = " << EPSold << ") too large; rescaling to EPS = " << EPS << " and restarting..." << endl;
        }

    	loadPotentialArrays(); //to update a and b with new eps, wont need to reload cluster data if it's needed though.
        //Reinitialise w and W to ICs (usually rand gaussian)
        setInitialConditions(nodeID+1);

        if (nodeID==1) { 
            print_line();
      	    cout.width(dwidth); cout << "Time";
      	    cout.width(dwidth); cout << "Energy";
      	    cout.width(dwidth); cout << "r_RMS";   
      	    cout << endl;
      	    print_line();
        }
        
        //reset step
        step = 0;
}

// reduce observables across nodes to first worker node
void computeObservables(dcomp*** wfnc) {
	
	// sum energy across nodes
	double energy_re=0.,energy_im=0.;
	double energy_re_collect=0.,energy_im_collect=0.;
	energy = wfncEnergy(wfnc);
	energy_re = real(energy);
	energy_im = imag(energy);
	//cout << energy << endl;
	MPI_Reduce(&energy_re,&energy_re_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
	MPI_Reduce(&energy_im,&energy_im_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
	energyCollect = dcomp(energy_re_collect,energy_im_collect);
	
	// sum normalization squared across nodes
	double normalization_re=0.,normalization_im=0.;
	double normalization_re_collect=0.,normalization_im_collect=0.;
	normalization = wfncNorm2(wfnc);
	normalization_re = real(normalization);
	normalization_im = imag(normalization);
	MPI_Reduce(&normalization_re,&normalization_re_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
	MPI_Reduce(&normalization_im,&normalization_im_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
	normalizationCollect = dcomp(normalization_re_collect,normalization_im_collect);
	
	// sum expectation across nodes
	double vInfinity_re=0.,vInfinity_im=0.;
	double vInfinity_re_collect=0.,vInfinity_im_collect=0.;
	vInfinity = vInfinityExpectationValue(wfnc);
	vInfinity_re = real(vInfinity);
	vInfinity_im = imag(vInfinity);	
	MPI_Reduce(&vInfinity_re,&vInfinity_re_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
	MPI_Reduce(&vInfinity_im,&vInfinity_im_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
	vInfinityCollect = dcomp(vInfinity_re_collect,vInfinity_im_collect);
	
	// sum r-squared across nodes
	double rRMS2_re=0.,rRMS2_im=0.;
	double rRMS2_re_collect=0.,rRMS2_im_collect=0.;
	rRMS2 = r2ExpectationValue(wfnc);
	rRMS2_re = real(rRMS2);
	rRMS2_im = imag(rRMS2);
	MPI_Reduce(&rRMS2_re,&rRMS2_re_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
	MPI_Reduce(&rRMS2_im,&rRMS2_im_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
	rRMS2Collect = dcomp(rRMS2_re_collect,rRMS2_im_collect);
}

// main computational solve routine
void solve() {
	
	dcomp energytot,lastenergy = 1.0e50;
	int done=1; //step=0,
	char label[64]; 
	
	leftSendBuffer = (double *)malloc( 2 * (NUMX+6) * (NUMY+6) * sizeof(double) );
	rightSendBuffer = (double *)malloc( 2 * (NUMX+6) * (NUMY+6) * sizeof(double) );
	leftReceiveBuffer = (double *)malloc( 2 * (NUMX+6) * (NUMY+6) * sizeof(double) );
	rightReceiveBuffer = (double *)malloc( 2 * (NUMX+6) * (NUMY+6) * sizeof(double) );
	
	// send message to master that we are starting
	MPI_Send( &done, 1, MPI_INT, 0, HELLO, MPI_COMM_WORLD );
	
	// evolve lattice in steps of UPDATE and output data along the way
	do {
		
		// sync boundaries
		syncBoundaries(w);
		
		// reduce observables across nodes
		computeObservables(w);
                
		// output 2d snapshots of the wavefunction for inspection
		// and check convergence of ground state energy
		if (step%SNAPUPDATE==0) {
			// broadcast observables
			MPI_Bcast(&normalizationCollect, 1, MPI_DOUBLE_COMPLEX, 0, workers_comm);
			MPI_Bcast(&energyCollect, 1, MPI_DOUBLE_COMPLEX, 0, workers_comm);
			// force symmetry
			symmetrizeWavefunction();
            // normalize wavefunction
            normalizeWavefunction(w);
			energytot =  energyCollect/normalizationCollect;
            //break run if we have a nanError - use RMS for check as energy can have both nan and inf issues...
            if (!isfinite(real(energytot))) {
                nanErrorCollect = 1;
                if (debug) debug_out << "Nan Error Detected" << endl; 
                MPI_Bcast(&nanErrorCollect, 1, MPI_INT, 0, workers_comm);
                break;
            }
            if (abs(energytot-lastenergy)<TOLERANCE) {
	            if (nodeID==1) outputMeasurements(step*EPS);
                break;
	        } else {
		        lastenergy = energytot;
		        if (step!=STEPS) recordSnapshot(w);
	        }
        }
        if (nodeID==1) outputMeasurements(step*EPS);
        if (step<STEPS) evolve(UPDATE);
        step += UPDATE;
                
	} while ((step<=STEPS) && nanErrorCollect == 0);
	
	// save grd-state energy and tau_f in global variables
   if (nanErrorCollect == 0) {
        EGrnd = energytot;
        timef = step*EPS;
       
        if (debug) {
	        debug_out << "==> Unnormalized Energy : " << energy << endl;
	        debug_out << "==> Normalization2 : " << normalization << endl;
        }

        if (nodeID==1) {
            outputSummaryData();
        }	

    }


    free(rightSendBuffer);
    free(leftSendBuffer);
    free(rightReceiveBuffer);
    free(leftReceiveBuffer);
	
}

// solve finalize
void solveFinalize() {

    char label[64];
	// this routine currently computes the first excited state energy and wavefunction
	if (EPS >= MINTSTEP) {
       //Comment if higher order states are wanted
	    if (SAVEWAVEFNCS==1) {
		// save 3d wavefunction for states
	    	sprintf(label,"0_%d",nodeID); 
	    	outputWavefunctionBinary(w,label);
        // output potential for debugging
            //outputPotential(label);
	    }
       //Uncomment if higher order states are wanted
     //  findExcitedStates();
    } else {
       if (nodeID==1) cout << "ERROR: MINTSTEP value exceeded. Aborting; check memory conditions and alter input parameters" << endl;
    }
    deallocateClusterMemory();
}

// evolves solution nsteps
void evolve(int nsteps) {
	
	int leftTest,rightTest;
	
	for (int i=1;i<=nsteps;i++) {
		
		// receive boundary sync
		leftTest=1; 
		rightTest=1;
		if (i>1) {
			if (nodeID+1 < numNodes) { 
				receiveRightBoundary(); 
				rightTest = 0;
            }
			if (nodeID-1 >= 1) {
				receiveLeftBoundary();
				leftTest = 0;
            }
			while (!leftTest || !rightTest) {
				if (!rightTest) {
					MPI_Test(&rightReceive,&rightTest,&rightReceiveStatus);
					if (rightTest) loadRightBoundaryFromBuffer(w);
				}
				if (!leftTest) {
					MPI_Test(&leftReceive,&leftTest,&leftReceiveStatus);
					if (leftTest) loadLeftBoundaryFromBuffer(w);
				}
			}
		}
		
		// first update boundary so that the send can be happening while we update the interior
		updateBoundaries(EPS);

		// wait and make sure send buffers are ready
		if (i>1) {
			if (nodeID+1 < numNodes) MPI_Wait(&rightSend,&rightSendStatus);
			if (nodeID-1 >= 1) MPI_Wait(&leftSend,&leftSendStatus);
		}
		
		// send boundary sync 
		leftTest=1; 
		rightTest=1;
		if (i==1) {
			if (nodeID+1 < numNodes) sendRightBoundary(W);
			if (nodeID-1 >= 1) sendLeftBoundary(W);
		}
		else if (i!=nsteps && i>1) {
			if (nodeID+1 < numNodes) rightTest=0;
			if (nodeID-1 >= 1) leftTest=0;
			while (!leftTest || !rightTest) {
				if (!rightTest) {
					MPI_Test(&rightSend,&rightTest,&rightSendStatus);
					if (rightTest) sendRightBoundary(W);
				}
				if (!leftTest) {
					MPI_Test(&leftSend,&leftTest,&leftSendStatus);
					if (leftTest) sendLeftBoundary(W);
				}
			}
		}
		
		// evolve interior of the grid forward in time by EPS storing updated fields in capital vars
		updateInterior(EPS);
		
		// copy fields from capital vars (updated) down to lowercase vars (current)
		copyDown();
		
	}
	
}

/*---------------------------------------------------------------------------*/
/* Find excited states                                                       */
/*---------------------------------------------------------------------------*/

void findExcitedStates() {
	
	char label[64]; 
	dcomp overlap=0,overlapCollect=0;
	dcomp overlap2=0,overlapCollect2=0;

	// compute first excited state
	int snap = 0;  //most recently saved snapshot;
	
	// compute overlap
	for (int sx=1;sx<=NUMX;sx++) 
		for (int sy=1;sy<=NUMY;sy++)
       	    for (int sz=1; sz<=DISTNUMZ;sz++) 
				overlap += w[sx][sy][sz]*wstore[snap][sx][sy][sz];
	
	MPI_Reduce(&overlap,&overlapCollect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
	MPI_Bcast(&overlapCollect, 1, MPI_DOUBLE, 0, workers_comm);
	
	// subtract overlap
	for (int sx=0;sx<NUMX+6;sx++) 
		for (int sy=0;sy<NUMY+6;sy++)
       	    for (int sz=0; sz<DISTNUMZ+2;sz++) 
				W[sx][sy][sz] = wstore[snap][sx][sy][sz] - overlapCollect*w[sx][sy][sz];
	
	// compute observables
	computeObservables(W);
	
	// normalize
	MPI_Bcast(&normalizationCollect, 1, MPI_DOUBLE, 0, workers_comm);
	normalizeWavefunction(W);
	
	// output snapshot of excited states
	sprintf(label,"excited_1_%d",nodeID); 
	outputSnapshot(W,label);
	
	if (nodeID==1) {
		EOne = energyCollect/normalizationCollect;
		dcomp vinfty = vInfinityCollect/normalizationCollect;
		dcomp dEtau = (EOne-EGrnd)*timef;     // #ad.
		dcomp rRMS2 = rRMS2Collect/normalizationCollect;     // #ad.
		print_line();
        cout.precision(dbl::digits10);
		cout << "==> 1st excited state energy : " << fixed << EOne << endl; //setprecision (12) before << ener
		cout << "==> 1st excited state binding energy : " << EOne - vinfty << endl;
		// #ad.
		cout << "==> Ex. State r_RMS : " << sqrt(real(rRMS2)) << endl;
		cout << "==> Ex. State L/r_RMS : " << float(NUMX)/sqrt(real(rRMS2)) << endl;
		
		
		cout << "==> (E1-E0)*tau_f : " << dEtau << endl;
		if (real(dEtau)<4.) {
			cout << "==> WARNING: states nearly degenerate, tau_f too small!" << endl;
      cout << "E1-E0 " << EOne-EGrnd << " step " << step << " EPS " << EPS << endl;
		}
		print_line();
	}
	
	if (SAVEWAVEFNCS) {
		// save 3d wavefunction for states
		sprintf(label,"0_%d",nodeID); 
		outputWavefunctionBinary(w,label);
		sprintf(label,"1_%d",nodeID); 
		outputWavefunctionBinary(W,label);
	}
	
  //compute second excited state
  //copyDown();

  overlap = 0;
  overlapCollect = 0;
	
	// compute first excited state
	snap = 1; //second most recently saved snapshot
	
	// compute overlap
	for (int sx=1;sx<=NUMX;sx++) 
		for (int sy=1;sy<=NUMY;sy++)
       	    for (int sz=1; sz<=DISTNUMZ;sz++) {
				overlap += w[sx][sy][sz]*wstore[snap][sx][sy][sz];
                overlap2 += W[sx][sy][sz]*wstore[snap][sx][sy][sz];
            }
	
	MPI_Reduce(&overlap,&overlapCollect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
	MPI_Bcast(&overlapCollect, 1, MPI_DOUBLE, 0, workers_comm);
    MPI_Reduce(&overlap2,&overlapCollect2,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
    MPI_Bcast(&overlapCollect2, 1, MPI_DOUBLE, 0, workers_comm);

	// subtract overlap
	for (int sx=0;sx<NUMX+6;sx++) 
		for (int sy=0;sy<NUMY+6;sy++)
            for (int sz=0; sz<DISTNUMZ+2;sz++) {
                W2[sx][sy][sz] = wstore[snap][sx][sy][sz] - overlapCollect2*W[sx][sy][sz] - overlapCollect*w[sx][sy][sz];
      }
	
   // compute observables
  computeObservables(W2);
   
   // normalize
   MPI_Bcast(&normalizationCollect, 1, MPI_DOUBLE, 0, workers_comm);
   normalizeWavefunction(W2);
  


	// output snapshot of excited states
	sprintf(label,"excited_2_%d",nodeID); 
	outputSnapshot(W2,label);
	
	if (nodeID==1) {
		dcomp ener = energyCollect/normalizationCollect;
		dcomp snapdiff = ((dcomp) 1)/(ener-EOne);
        dcomp vinfty = vInfinityCollect/normalizationCollect;
		dcomp rRMS2 = rRMS2Collect/normalizationCollect;
		print_line();
        cout.precision(dbl::digits10);
		cout << "==> 2nd excited state energy : " << fixed << ener << endl; //setprecision (12) before << ener
		cout << "==> 2nd excited state binding energy : " << ener - vinfty << endl;
		cout << "==> Ex. State r_RMS : " << sqrt(real(rRMS2)) << endl;
		cout << "==> Ex. State L/r_RMS : " << float(NUMX)/sqrt(real(rRMS2)) << endl;
		
        cout << "==> 1/(E2-E1) : " << snapdiff << endl; 	
		if (real(snapdiff)*2>SNAPUPDATE) {
        cout << "==> WARNING: SNAPUPDATE is too small!" << endl;
    }
    if (real(snapdiff)<1) {
      cout << "==> WARNING: SNAPUPDATE is too large!" << endl;
    }

    print_line();
	}
	
	if (SAVEWAVEFNCS) {
		// save 3d wavefunction for states
		sprintf(label,"2_%d",nodeID); 
		outputWavefunctionBinary(W2,label);
  }
	return;
}

/*---------------------------------------------------------------------------*/
/* Boundary sync routines                                                    */
/*---------------------------------------------------------------------------*/

void syncBoundaries(dcomp ***wfnc) {
	
	// initiate sends and receives
	if (nodeID+1 < numNodes) { 
		sendRightBoundary(wfnc);
		receiveRightBoundary(); 
	}
	if (nodeID-1 >= 1) {
		sendLeftBoundary(wfnc);
		receiveLeftBoundary(); 
	}
	
	// now wait for communications to complete and sync wfnc when they do
	if (nodeID+1 < numNodes) { 
		MPI_Wait(&rightReceive,&rightReceiveStatus);
		loadRightBoundaryFromBuffer(wfnc);
		MPI_Wait(&rightSend,&rightSendStatus);
	}
	if (nodeID-1 >= 1) {
		MPI_Wait(&leftReceive,&leftReceiveStatus);
		loadLeftBoundaryFromBuffer(wfnc);
		MPI_Wait(&leftSend,&leftSendStatus);
    }
}

void sendRightBoundary(dcomp*** wfnc) {
	char message[64]; 
	if (debug > 1) {
		sprintf(message,"%d -> %d",nodeID,nodeID+1); 
		debug_out << "==> Sending : " << message << endl;
		MPI_Isend(message, strlen(message), MPI_CHAR, nodeID+1, SYNC_RIGHT_MESSAGE, MPI_COMM_WORLD, &rightMessageSend); 
	}
	for (int sx=0;sx<NUMX+6;sx++)
		for (int sy=0;sy<NUMY+6;sy++) {
			rightSendBuffer[sx*(NUMY+6)+sy] = real(wfnc[sx][sy][DISTNUMZ]);
			rightSendBuffer[sx*(NUMY+6)+sy + (NUMX+6)*(NUMY+6)] = imag(wfnc[sx][sy][DISTNUMZ]);
		}
	MPI_Isend(rightSendBuffer, 2*(NUMX+6)*(NUMY+6), MPI_DOUBLE, nodeID+1, SYNC_RIGHT, MPI_COMM_WORLD, &rightSend); 
}

void sendLeftBoundary(dcomp*** wfnc) {
	char message[64]; 
	if (debug > 1) {
		sprintf(message,"%d -> %d",nodeID,nodeID-1); 
		debug_out << "==> Sending : " << message << endl;
		MPI_Isend(message, strlen(message), MPI_CHAR, nodeID-1, SYNC_LEFT_MESSAGE, MPI_COMM_WORLD, &leftMessageSend); 
	}
	for (int sx=0;sx<NUMX+6;sx++)
		for (int sy=0;sy<NUMY+6;sy++) { 
			leftSendBuffer[sx*(NUMY+6)+sy] = real(wfnc[sx][sy][1]);
			leftSendBuffer[sx*(NUMY+6)+sy + (NUMX+6)*(NUMY+6)] = imag(wfnc[sx][sy][1]);
		}
	MPI_Isend(leftSendBuffer, 2*(NUMX+6)*(NUMY+6), MPI_DOUBLE, nodeID-1, SYNC_LEFT, MPI_COMM_WORLD, &leftSend);
}

void receiveRightBoundary() {
	char message[64]; 
	if (debug > 1) {
		MPI_Irecv(message, 255, MPI_CHAR, nodeID+1, SYNC_LEFT_MESSAGE, MPI_COMM_WORLD, &rightMessageReceive); 
		debug_out << "==> Received : " << message << endl;
	}
	MPI_Irecv(rightReceiveBuffer, 2*(NUMX+6)*(NUMY+6), MPI_DOUBLE, nodeID+1, SYNC_LEFT, MPI_COMM_WORLD, &rightReceive); 
}

inline void loadRightBoundaryFromBuffer(dcomp ***wfnc) {
	// update w array right boundary
	for (int sx=0;sx<NUMX+6;sx++)
		for (int sy=0;sy<NUMY+6;sy++)
			wfnc[sx][sy][DISTNUMZ+1] = dcomp(rightReceiveBuffer[sx*(NUMY+6)+sy],rightReceiveBuffer[sx*(NUMY+6)+sy+(NUMX+6)*(NUMY+6)]);
}

void receiveLeftBoundary() {
	char message[64];
	if (debug > 1) {
		MPI_Irecv(message, 255, MPI_CHAR, nodeID-1, SYNC_RIGHT_MESSAGE, MPI_COMM_WORLD, &leftMessageReceive); 
		debug_out << "==> Received : "<< message << endl;
	}
	MPI_Irecv(leftReceiveBuffer, 2*(NUMX+6)*(NUMY+6), MPI_DOUBLE, nodeID-1, SYNC_RIGHT, MPI_COMM_WORLD, &leftReceive); 
}

inline void loadLeftBoundaryFromBuffer(dcomp ***wfnc) {
	// update w array left boundary
	for (int sx=0;sx<NUMX+6;sx++)
		for (int sy=0;sy<NUMY+6;sy++)
			wfnc[sx][sy][0] = dcomp(leftReceiveBuffer[sx*(NUMY+6)+sy],leftReceiveBuffer[sx*(NUMY+6)+sy+(NUMX+6)*(NUMY+6)]);
}
