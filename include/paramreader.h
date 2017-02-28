/*

   paramreader.h

   Copyright (c) Michael Strickland
   Forked at v2.0; Additions by Tim DuBois

   GNU General Public License (GPLv3)
   See detailed text in license directory

*/


#ifndef __paramreader_h__
#define __paramreader_h__

void readParametersFromFile(char *filename, int echo);
void readParametersFromCommandLine(int argc, char** argv, int echo);
void readClusterData(char *filename, int echo);

#endif /* __paramreader_h__ */
