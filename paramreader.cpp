/*

   paramreader.cpp

   Copyright (c) Michael Strickland

   GNU General Public License (GPLv3)
   See detailed text in license directory

*/

#include <fstream>
#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <complex>

#include "mpisolve.h"

using namespace std;

// this workhorse examines a key to see if it corresponds to a var we are setting
// and then attempts to set the var corresponding to key by converting value to the
// appropriate type.  lots of hardcoding here
void setParameter(char *key, const char *value) {
	// integer params
	if (strcmp(key,"NUMZ")==0) NUMZ=atoi(value);
	if (strcmp(key,"UPDATE")==0) UPDATE=atoi(value);
	if (strcmp(key,"SNAPUPDATE")==0) SNAPUPDATE=atoi(value);
	if (strcmp(key,"POTENTIAL")==0) POTENTIAL=atoi(value);
	if (strcmp(key,"INITCONDTYPE")==0) INITCONDTYPE=atoi(value);
	if (strcmp(key,"INITSYMMETRY")==0) INITSYMMETRY=atoi(value);
	if (strcmp(key,"NF")==0) NF=atoi(value);
	if (strcmp(key,"SAVEWAVEFNCS")==0) SAVEWAVEFNCS=atoi(value);
        if (strcmp(key,"CLUSTER")==0) CLUSTER=atoi(value);
	// double/float params
	if (strcmp(key,"A")==0) A=atof(value);
	if (strcmp(key,"EPS")==0) EPS=atof(value);
	if (strcmp(key,"MINTSTEP")==0) MINTSTEP=atof(value);
	if (strcmp(key,"SIG")==0) SIG=atof(value);
	if (strcmp(key,"SIGMA")==0) SIGMA=atof(value);
	if (strcmp(key,"STEPS")==0) STEPS=atof(value);
	if (strcmp(key,"MASS")==0) MASS=atof(value);
	if (strcmp(key,"T")==0) T=atof(value);
	if (strcmp(key,"TC")==0) TC=atof(value);
	if (strcmp(key,"XI")==0) XI=atof(value);
	if (strcmp(key,"TOLERANCE")==0) TOLERANCE=atof(value);
    if (strcmp(key,"ALX")==0) ALX=atof(value);
    if (strcmp(key,"ALY")==0) ALY=atof(value);
    if (strcmp(key,"ALZ")==0) ALZ=atof(value);
  return;
}

//
// This routine assumes that paramters are in text file with
// each parameter on a new line in the format 
//
// PARAMKEY	PARAMVALUE
//
// The PARAMKEY must begin the line and only tabs and spaces
// can appear between the PARAMKEY and PARAMVALUE.
// 
// Lines which begin with 'commentmarker' defined below are ignored
//
void readParametersFromFile(char *filename, int echo) {
		
	string commentmarker = "//"; 
	char space = ' '; 
	char tab = '\t';

	int maxline = 128; // maximum line length used in the buffer for reading
	char buffer[maxline];
	ifstream paramFile(filename);
	
	while(!paramFile.eof()) {
		paramFile.getline(buffer,maxline,'\n');
		string line = buffer; int length = strlen(buffer);
		if (line.substr(0,commentmarker.length())!=commentmarker && line.length()>0) {
			char key[64]="",value[64]=""; int founddelim=0;
			for (int i=0;i<length;i++) {
				if (buffer[i]==space || buffer[i]==tab) founddelim=1;
				else {
					if (founddelim==0) key[strlen(key)] = buffer[i];
					else value[strlen(value)] = buffer[i];
				}
			}
			if (strlen(key)>0 && strlen(value)>0) {
				setParameter(key,value);
				if (echo) cout << key << " = " << value << endl;
			}
		}
	}
	
	return;	
}

//
// Read parameters from commandline
//
void readParametersFromCommandLine(int argc, char** argv, int echo) {
	int optind = 1;
	while (optind < argc) 
	{
	  if (argv[optind][0]=='-') {
	  	string key = argv[optind];
	  	key = key.substr(1,key.length()); // remove '-'
	  	string value = argv[optind+1]; // load value
  	  	if (echo) cout << key << " = " << value << endl;
	  	setParameter((char *)key.c_str(),value.c_str());
	  	optind++;
	  }
	  optind++;
	}
	return;
}

//
// Read cluster data from cluster.xyz
//
void readClusterData(char *filename, int clusterSize, int echo) {
    fstream input;
    string line, speciesstr, tmp;
    double xval, yval, zval, xmax=0, ymax=0;

    input.open(filename, ios::in);
    //first line is the number of atoms
    getline(input,line);
    getline(input,line); // ignore the comment line


    if (echo) cout << clusterSize << " atoms in cluster (these " << clusterSize-1 << " cage atoms and the delocalised oxygen):" << endl;

    for (int i = 0; i < clusterSize-1; i++) { //the minus one leaves room for the oxygen
        getline(input,line);
        speciesstr = line.substr(0, line.find(' ')); //species
        line.erase(0, line.find(' ') + 1); //remove numx
        line.erase(0, line.find_first_not_of(' ')); //remove whitespace padding
        tmp = line.substr(0, line.find(' '));
        xval = atof(tmp.c_str()); //x value
        line.erase(0, line.find(' ') + 1); //remove x
        line.erase(0, line.find_first_not_of(' ')); //remove whitespace padding
        tmp = line.substr(0, line.find(' '));
        yval = atof(tmp.c_str()); //y value
        line.erase(0, line.find(' ') + 1); //remove y
        line.erase(0, line.find_first_not_of(' ')); //remove whitespace padding
        tmp = line.substr(0, line.find(' '));
        zval = atof(tmp.c_str()); //z value
        
        if (speciesstr.compare("Al") == 0) {
            *(clustSpecies+i) = 1; //Aluminium
        } else {
            *(clustSpecies+i) = 2; //Oxygen
        }

        *(clust+(clusterSize*0)+i) = xval;
        *(clust+(clusterSize*1)+i) = yval;
        *(clust+(clusterSize*2)+i) = zval;
        
        if (xval > xmax) xmax = xval;
        if (yval > ymax) ymax = yval;
    
        if (echo) cout << speciesstr << "(" << *(clustSpecies+i) << "), " << *(clust+(clusterSize*0)+i) << ", " << *(clust+(clusterSize*1)+i) << ", " << *(clust+(clusterSize*2)+i) << endl;
    }
    *(clustSpecies+clusterSize-1) = 2; //final oxygen
    input.close();	
    
    //NUMZ (& DISTNUMZ) Should already be calculated and in params.txt
    //use max values of cluster to set new values for NUM{X,Y} 
    NUMX = ceil((2*xmax+A)/A);
    NUMY = ceil((2*ymax+A)/A);

    return;	
}
