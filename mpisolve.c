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
int    DISTNUMZ=20,NUMX=20,NUMY=20,NUMZ=20,UPDATE=100,SNAPUPDATE=1000,WAVENUM=0,WAVEMAX=5;
int    POTENTIAL=0,INITCONDTYPE=0,INITSYMMETRY=0,NF=2,SAVEWAVEFNCS=0,RUNTYPE=0,CLUSTSIZE=7,OUTPOT=0,EXCITEDSTATES=0;
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
dcomp normalizationCollect=(1,1);  // the total normalization squared
dcomp beta=0;
dcomp betaCollect=0; 
dcomp vInfinity=0;		// the local node expectation value of v_infty
dcomp vInfinityCollect=0;      // the total expectation value of v_infty
dcomp rRMS2=0;                 // the local node <r^2>  #ad.
dcomp rRMS2Collect=0;		// the total <r^2>  #ad.
int nanErrorCollect=0;

//grnd-state energy and final time saved in global var after convergence
dcomp	EGrnd, EOne, timef;
int step = 0;

int main( int argc, char *argv[] ) 
{ 
    int done=0,checksum=0;
    char message[64]; 
    char label[64];
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
		if (debug == DEBUG_FULL) debug_out << "==> Rec Hello : "<< message << endl;

		// the master said hello so now let's get to work
        // Loops for each wavefunction to be calculated
		for (int ii=WAVENUM;ii<=WAVEMAX;ii++) { 
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
        
            if (debug == DEBUG_FULL) debug_out << "Solved" << endl;
            //solve() has found the ground state and is stored in w.
            storeConverged(w,WAVENUM);
            if (debug == DEBUG_FULL) debug_out << "Stored" << endl;
	        if (SAVEWAVEFNCS==1) {
		    // save 3d wavefunction for states
                sprintf(label,"%d_%d",WAVENUM, nodeID); 
                if (debug == DEBUG_FULL) debug_out << "Label" << endl;
	    	    outputWavefunctionBinary(w,label);
                if (debug == DEBUG_FULL) debug_out << "Output" << endl;
	        }
            
            if (ii<WAVEMAX) { //TODO: Remove hardcoded one
                //Set WAVENUM flag and go again.
                WAVENUM++;

                //Reset values and load new wavefunctions from matlab guess
                reInitSolver();
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

    if (WAVENUM>0) {
        //Starting from previously output states. Load lower ones to the store
        for (int ii=0; ii<WAVENUM; ii++) {
	        if (nodeID==1) { 
                cout << "==> Loading converged wavefunction (state " << ii << ") into memory" << endl;
            }
            readWavefunctionBinary(ii);
        }
    } 
	// set initial conditions
	setInitialConditions(nodeID+1);
	// output some summary information
	if (nodeID==1) { 
	    print_line();
      	cout << "==> Energy offset epsilon : " << epsilon << endl; 
      	cout << "==> Number of computational nodes : " << numNodes-1 << endl; 
      	cout << "==> DISTNUMZ : " << DISTNUMZ << endl; 
      	print_line();
      	cout << "Spatial Step Size (A): " << A << endl;
      	cout << "Temporal Step Size (EPS): " << EPS << endl;
      	cout << "Standard Deviation of initial wavefunction noise (SIG): " << SIG << endl;
      	cout << "Box Dimensions: X = " << A*NUMX << ", Y = " << A*NUMY << ", Z = " << A*NUMZ << endl;
      	print_line();
      	cout.width(20); cout << "Time";
      	cout.width(30); cout << "Energy";
      	cout.width(20); cout << "r_RMS";   
      	cout.width(20); cout << "Diff";   
      	cout << endl;
      	print_line();
	}
}

void reInitSolver() {
    //Reset Everything and load new wavefunction guess from disk
    if (nodeID==1) { 
        print_line();
        cout << "Energy converged to tolerance, wavefunction found. Loading higer state." << endl;
		flush(cout);
    }
   
    setInitialConditions(nodeID+1);
    
    //reset step
    step = 0;
	
    // output some summary information
	if (nodeID==1) { 
	print_line();
      	print_line();
      	cout.width(20); cout << "Time";
      	cout.width(30); cout << "Energy";
      	cout.width(20); cout << "r_RMS";   
      	cout.width(20); cout << "Diff";   
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
      	    cout.width(20); cout << "Time";
      	    cout.width(30); cout << "Energy";
      	    cout.width(20); cout << "r_RMS";   
      	    cout.width(20); cout << "Diff";   
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
	MPI_Reduce(&energy_re,&energy_re_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
	MPI_Reduce(&energy_im,&energy_im_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
	energyCollect = dcomp(energy_re_collect,energy_im_collect);
    
    getNormalization(wfnc);	

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

void getNormalization(dcomp*** wfnc) {
	// sum normalization squared across nodes
	double normalization_re=0.,normalization_im=0.;
	double normalization_re_collect=0.,normalization_im_collect=0.;
	normalization = wfncNorm2(wfnc);
	normalization_re = real(normalization);
	normalization_im = imag(normalization);
	MPI_Reduce(&normalization_re,&normalization_re_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
	MPI_Reduce(&normalization_im,&normalization_im_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
	normalizationCollect = dcomp(normalization_re_collect,normalization_im_collect);
}

void getOverlap(dcomp*** wfnc) {
    //Here for great justice, but may be removed if new GS works    

	MPI_Bcast(&normalizationCollect, 1, MPI_DOUBLE_COMPLEX, 0, workers_comm);
    dcomp norm = sqrt(normalizationCollect);

	// compute overlap
    for (int wnum=0; wnum<WAVENUM; wnum++)
	    for (int sx=3;sx<=2+NUMX;sx++) 
	        for (int sy=3;sy<=2+NUMY;sy++)
                for (int sz=3; sz<=2+DISTNUMZ;sz++) { 
                    beta += (wstore[wnum][sx][sy][sz]*(wfnc[sx][sy][sz] /= norm));
    }
	
	MPI_Reduce(&beta,&betaCollect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
	MPI_Bcast(&betaCollect, 1, MPI_DOUBLE, 0, workers_comm);
	
    //betaCollect is now beta, a number.
    updatePotential(betaCollect);
	
}

void gramSchmidt() {

    //w is normalised through snapupdate loop, no need to do it again
    
	dcomp overlapGS=0,overlapGSCollect=0;
	dcomp overlapES1=0,overlapES1Collect=0;
	dcomp overlapES2=0,overlapES2Collect=0;
	dcomp overlapES3=0,overlapES3Collect=0;
	dcomp overlapES4=0,overlapES4Collect=0;
	
	// compute overlaps
	for (int sx=3;sx<=2+NUMX;sx++) 
		for (int sy=3;sy<=2+NUMY;sy++)
       	    for (int sz=3; sz<=2+DISTNUMZ;sz++) {
                if (WAVENUM >= 1) {
                    overlapGS += wstore[0][sx][sy][sz]*w[sx][sy][sz];
                }
                if (WAVENUM >= 2) {
                    overlapES1 += wstore[1][sx][sy][sz]*w[sx][sy][sz];
                }
                if (WAVENUM >= 3) {
                    overlapES2 += wstore[2][sx][sy][sz]*w[sx][sy][sz];
                }
                if (WAVENUM >= 4) {
                    overlapES3 += wstore[3][sx][sy][sz]*w[sx][sy][sz];
                }
                if (WAVENUM >= 5) {
                    overlapES4 += wstore[4][sx][sy][sz]*w[sx][sy][sz];
                }
            }

    //collect overlaps
    if (WAVENUM >= 1) {
	    MPI_Reduce(&overlapGS,&overlapGSCollect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
	    MPI_Bcast(&overlapGSCollect, 1, MPI_DOUBLE, 0, workers_comm);
    }
    if (WAVENUM >= 2) {
	    MPI_Reduce(&overlapES1,&overlapES1Collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
	    MPI_Bcast(&overlapES1Collect, 1, MPI_DOUBLE, 0, workers_comm);
    }
    if (WAVENUM >= 3) {
	    MPI_Reduce(&overlapES2,&overlapES2Collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
	    MPI_Bcast(&overlapES2Collect, 1, MPI_DOUBLE, 0, workers_comm);
    }
    if (WAVENUM >= 4) {
	    MPI_Reduce(&overlapES3,&overlapES3Collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
	    MPI_Bcast(&overlapES3Collect, 1, MPI_DOUBLE, 0, workers_comm);
    }
    if (WAVENUM >= 5) {
	    MPI_Reduce(&overlapES4,&overlapES4Collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
	    MPI_Bcast(&overlapES4Collect, 1, MPI_DOUBLE, 0, workers_comm);
    }
	
	// subtract overlap
	for (int sx=0;sx<=NUMX+5;sx++) 
		for (int sy=0;sy<=NUMY+5;sy++)
       	    for (int sz=0; sz<=DISTNUMZ+5;sz++) {
                switch (WAVENUM) {
                    case 1:
                        W[sx][sy][sz] = w[sx][sy][sz] - wstore[0][sx][sy][sz]*overlapGSCollect;
                        break;
                    case 2:
                        W[sx][sy][sz] = w[sx][sy][sz] - wstore[0][sx][sy][sz]*overlapGSCollect - wstore[1][sx][sy][sz]*overlapES1Collect;
                        break;
                    case 3:
                        W[sx][sy][sz] = w[sx][sy][sz] - wstore[0][sx][sy][sz]*overlapGSCollect - wstore[1][sx][sy][sz]*overlapES1Collect - wstore[2][sx][sy][sz]*overlapES2Collect;
                        break;
                    case 4:
                        W[sx][sy][sz] = w[sx][sy][sz] - wstore[0][sx][sy][sz]*overlapGSCollect - wstore[1][sx][sy][sz]*overlapES1Collect - wstore[2][sx][sy][sz]*overlapES2Collect - wstore[3][sx][sy][sz]*overlapES3Collect;
                        break;
                    case 5:
                        W[sx][sy][sz] = w[sx][sy][sz] - wstore[0][sx][sy][sz]*overlapGSCollect - wstore[1][sx][sy][sz]*overlapES1Collect - wstore[2][sx][sy][sz]*overlapES2Collect - wstore[3][sx][sy][sz]*overlapES3Collect - wstore[4][sx][sy][sz]*overlapES4Collect;
                        break;
                }
            }
   
    copyDown(); 
	
}

// main computational solve routine
void solve() {
	
	dcomp energytot,lastenergy,dispenergy = 1.0e50;
	int done=1, gs=0;
	char label[64]; 
	
	leftSendBuffer = (double *)malloc( 6 * (NUMX+6) * (NUMY+6) * sizeof(double) );
	rightSendBuffer = (double *)malloc( 6 * (NUMX+6) * (NUMY+6) * sizeof(double) );
	leftReceiveBuffer = (double *)malloc( 6 * (NUMX+6) * (NUMY+6) * sizeof(double) );
	rightReceiveBuffer = (double *)malloc( 6 * (NUMX+6) * (NUMY+6) * sizeof(double) );
	
	// send message to master that we are starting
	MPI_Send( &done, 1, MPI_INT, 0, HELLO, MPI_COMM_WORLD );
	
	// evolve lattice in steps of UPDATE and output data along the way
	do {
	    gs = 0;	
		// sync boundaries
		syncBoundaries(w);
		
		// reduce observables across nodes
		computeObservables(w);
                
		// Check convergence of ground state energy
		if (step%SNAPUPDATE==0) { //TODO: I think we can do away with SNAPUPDATE now. Kill this if.
			// broadcast observables
			MPI_Bcast(&normalizationCollect, 1, MPI_DOUBLE_COMPLEX, 0, workers_comm);
			// force symmetry
			symmetrizeWavefunction();
            // normalize wavefunction
            normalizeWavefunction(w);
			
            // Orthoganalise wavefunction
            if (WAVENUM>0) {
                gramSchmidt();
                gs = 1;
		        syncBoundaries(w);
		        computeObservables(w);
			    MPI_Bcast(&normalizationCollect, 1, MPI_DOUBLE_COMPLEX, 0, workers_comm);
			    symmetrizeWavefunction();
                normalizeWavefunction(w);
            }
            
            MPI_Bcast(&energyCollect, 1, MPI_DOUBLE_COMPLEX, 0, workers_comm);
			energytot =  energyCollect/normalizationCollect;
            //break run if we have a nanError - use RMS for check as energy can have both nan and inf issues...
            if (!isfinite(real(energytot))) {
                nanErrorCollect = 1;
                if (debug) debug_out << "Nan Error Detected" << endl; 
                MPI_Bcast(&nanErrorCollect, 1, MPI_INT, 0, workers_comm);
                break;
            }
            if (abs(energytot-lastenergy)<TOLERANCE) {
	            if (nodeID==1) outputMeasurements(step*EPS, lastenergy);
                break;
	        } else {
                dispenergy = lastenergy;
		        lastenergy = energytot;
	        }
        }
        if (nodeID==1) {
            outputMeasurements(step*EPS, dispenergy);
        }
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
            outputSummaryData(WAVENUM);
        }	

    }


    free(rightSendBuffer);
    free(leftSendBuffer);
    free(rightReceiveBuffer);
    free(leftReceiveBuffer);
	
}

// solve finalize
void solveFinalize() {

	// this routine currently computes the first excited state energy and wavefunction
	if (EPS >= MINTSTEP) {
       //Comment if higher order states are wanted
       if (EXCITEDSTATES == 1) {
           //if higher order states are wanted
          findExcitedStates();
       }
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
		
        // Find overlap with lower level wavefunctions
        //if (WAVENUM>0) {
            //getNormalization(w);
            //getOverlap(w);
        //}
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
	int state = 0;  //Ground state;
	
	// compute overlap
	for (int sx=3;sx<=2+NUMX;sx++) 
		for (int sy=3;sy<=2+NUMY;sy++)
       	    for (int sz=3; sz<=2+DISTNUMZ;sz++) 
				overlap += w[sx][sy][sz]*wstore[state][sx][sy][sz];
	
	MPI_Reduce(&overlap,&overlapCollect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
	MPI_Bcast(&overlapCollect, 1, MPI_DOUBLE, 0, workers_comm);
	
	// subtract overlap
	for (int sx=0;sx<=NUMX+5;sx++) 
		for (int sy=0;sy<=NUMY+5;sy++)
       	    for (int sz=0; sz<=DISTNUMZ+5;sz++) 
				W[sx][sy][sz] = wstore[state][sx][sy][sz] - overlapCollect*w[sx][sy][sz];
	
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
        if (POTENTIAL == 22) {
		    cout << "==> 1st excited state energy : " << fixed << EOne*1e6/239.2311 << endl; //setprecision (12) before << ener
		    cout << "==> 1st excited state binding energy : " << EOne*1e6/239.2311 - vinfty << endl;
        } else {
		    cout << "==> 1st excited state energy : " << fixed << EOne << endl; //setprecision (12) before << ener
		    cout << "==> 1st excited state binding energy : " << EOne - vinfty << endl;
        }
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
		sprintf(label,"1_%d",nodeID); 
		outputWavefunctionBinary(W,label);
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
	for (int sx=0;sx<=NUMX+5;sx++)
		for (int sy=0;sy<=NUMY+5;sy++) {
			rightSendBuffer[sx*(NUMY+6)+sy] = real(wfnc[sx][sy][DISTNUMZ]);
			rightSendBuffer[sx*(NUMY+6)+sy + (NUMX+6)*(NUMY+6)] = imag(wfnc[sx][sy][DISTNUMZ]);
			rightSendBuffer[sx*(NUMY+6)+sy + (2*(NUMX+6)*(NUMY+6))] = real(wfnc[sx][sy][1+DISTNUMZ]);
			rightSendBuffer[sx*(NUMY+6)+sy + (3*(NUMX+6)*(NUMY+6))] = imag(wfnc[sx][sy][1+DISTNUMZ]);
			rightSendBuffer[sx*(NUMY+6)+sy + (4*(NUMX+6)*(NUMY+6))] = real(wfnc[sx][sy][2+DISTNUMZ]);
			rightSendBuffer[sx*(NUMY+6)+sy + (5*(NUMX+6)*(NUMY+6))] = imag(wfnc[sx][sy][2+DISTNUMZ]);
		}
	MPI_Isend(rightSendBuffer, 6*(NUMX+6)*(NUMY+6), MPI_DOUBLE, nodeID+1, SYNC_RIGHT, MPI_COMM_WORLD, &rightSend); 
}

void sendLeftBoundary(dcomp*** wfnc) {
	for (int sx=0;sx<=NUMX+5;sx++)
		for (int sy=0;sy<=NUMY+5;sy++) { 
			leftSendBuffer[sx*(NUMY+6)+sy] = real(wfnc[sx][sy][3]);
			leftSendBuffer[sx*(NUMY+6)+sy + (NUMX+6)*(NUMY+6)] = imag(wfnc[sx][sy][3]);
			leftSendBuffer[sx*(NUMY+6)+sy + (2*(NUMX+6)*(NUMY+6))] = real(wfnc[sx][sy][4]);
			leftSendBuffer[sx*(NUMY+6)+sy + (3*(NUMX+6)*(NUMY+6))] = imag(wfnc[sx][sy][4]);
			leftSendBuffer[sx*(NUMY+6)+sy + (4*(NUMX+6)*(NUMY+6))] = real(wfnc[sx][sy][5]);
			leftSendBuffer[sx*(NUMY+6)+sy + (5*(NUMX+6)*(NUMY+6))] = imag(wfnc[sx][sy][5]);
		}
	MPI_Isend(leftSendBuffer, 6*(NUMX+6)*(NUMY+6), MPI_DOUBLE, nodeID-1, SYNC_LEFT, MPI_COMM_WORLD, &leftSend);
}

void receiveRightBoundary() {
	MPI_Irecv(rightReceiveBuffer, 6*(NUMX+6)*(NUMY+6), MPI_DOUBLE, nodeID+1, SYNC_LEFT, MPI_COMM_WORLD, &rightReceive); 
}

inline void loadRightBoundaryFromBuffer(dcomp ***wfnc) {
	// update w array right boundary
	for (int sx=0;sx<=NUMX+5;sx++)
		for (int sy=0;sy<=NUMY+5;sy++) {
			wfnc[sx][sy][2+DISTNUMZ+1] = dcomp(rightReceiveBuffer[sx*(NUMY+6)+sy],rightReceiveBuffer[sx*(NUMY+6)+sy+(NUMX+6)*(NUMY+6)]);
			wfnc[sx][sy][2+DISTNUMZ+2] = dcomp(rightReceiveBuffer[sx*(NUMY+6)+sy + (2*(NUMX+6)*(NUMY+6))],rightReceiveBuffer[sx*(NUMY+6)+sy + (3*(NUMX+6)*(NUMY+6))]);
			wfnc[sx][sy][2+DISTNUMZ+3] = dcomp(rightReceiveBuffer[sx*(NUMY+6)+sy + (4*(NUMX+6)*(NUMY+6))],rightReceiveBuffer[sx*(NUMY+6)+sy + (5*(NUMX+6)*(NUMY+6))]);
        }
}

void receiveLeftBoundary() {
	MPI_Irecv(leftReceiveBuffer, 6*(NUMX+6)*(NUMY+6), MPI_DOUBLE, nodeID-1, SYNC_RIGHT, MPI_COMM_WORLD, &leftReceive); 
}

inline void loadLeftBoundaryFromBuffer(dcomp ***wfnc) {
	// update w array left boundary
	for (int sx=0;sx<=NUMX+5;sx++)
		for (int sy=0;sy<=NUMY+5;sy++) {
			wfnc[sx][sy][2] = dcomp(leftReceiveBuffer[sx*(NUMY+6)+sy + (4*(NUMX+6)*(NUMY+6))],leftReceiveBuffer[sx*(NUMY+6)+sy + (5*(NUMX+6)*(NUMY+6))]);
			wfnc[sx][sy][1] = dcomp(leftReceiveBuffer[sx*(NUMY+6)+sy + (2*(NUMX+6)*(NUMY+6))],leftReceiveBuffer[sx*(NUMY+6)+sy + (3*(NUMX+6)*(NUMY+6))]);
			wfnc[sx][sy][0] = dcomp(leftReceiveBuffer[sx*(NUMY+6)+sy],leftReceiveBuffer[sx*(NUMY+6)+sy+(NUMX+6)*(NUMY+6)]);
        }
}
