/*

   grid.h

   Copyright (c) Michael Strickland
   Forked at v2.0; Additions by Tim DuBois

   GNU General Public License (GPLv3)
   See detailed text in license directory

*/

#ifndef __initialconditions_h__
#define __initialconditions_h__

void symmetrizeWavefunction();
void setInitialConditions(int seedMult);
void readWavefunctionBinary(int waveNum);

#endif /* __initialconditions_h__ */
