/* mexHatPotential.cpp
 * No frills potential builder designed for 3D calculations. dx independant.
 * Pulled from Matlab Mex routines and forced into OSX / gcc
 * Compile:
 * requrired flags for lapack:
 * -llapack -lblas
 * Written by Tim DuBois 24/09/12 */
//#include "mex.h"
#include <cmath>
#include <string>
#include <iostream>
#include <complex>

#if defined(VAYU)
#include <mkl.h>
#define LP_INT MKL_INT
#elif defined(TRIFID)
#define LP_INT int
#else
#error "One of VAYU or TRIFID must be defined"
#endif
//#include <Accelerate/Accelerate.h> //For Lapack & Blas

using namespace std;

#include "mexHatPotential.h"
#include "mpisolve.h"

/* Function is still cluster size dependant */
#define CLUSTER 7
#define PROD 49


/* Calculate matrix "\" division using LAPACK */
void matRightDivide(double *A, double *B, LP_INT M, LP_INT N)
{
    /* Create inputs for DGESV */
    LP_INT *iPivot;
    LP_INT info;
    iPivot = (LP_INT *)malloc((M+1)*sizeof(LP_INT));
    
    /* Call LAPACK */
    dgesv_(&M,&N,A,&M,iPivot,B,&M,&info);
    /* B now holds A\B */
    
    free(iPivot);
}

/* Calculate matrix multiplication using BLAS */
void matMultiply(double *A, double *B, double *C, int M, int N)
{
//    char *chn = "N";
    /* scalar values to use in dgemm */
    double one = 1.0, zero = 0.0;
    
    /* Pass arguments to Fortran by reference */
    dgemm_((char *)"N", (char *)"N", &M, &N, &M, &one, A, &M, B, &M, &zero, C, &M);
}

/* Calculate matrix inverse using LAPACK */
void matInverse(double *A, LP_INT N)
{
    double *WORK;  /* in/out arguments to dgetrf */
    LP_INT *iPivot;   /* inputs to dgetrf */
    LP_INT info;
    LP_INT prod;
    
    /* Create inputs for dgetrf */
    prod = N*N;
    WORK = (double *)malloc(prod*sizeof(double));
    iPivot = (LP_INT *)malloc((N+1)*sizeof(LP_INT));
    
    /* dgetrf(M,N,A,LDA,IPIV,INFO)
     * LU decomoposition of a general matrix, which means invert LDA columns of an M by N matrix
     * called A, sending the pivot indices to IPIV, and spitting error information to INFO. */
    dgetrf_(&N,&N,A,&N,iPivot,&info);
    
    /* generate inverse of a matrix given its LU decomposition */
    dgetri_(&N,A,&N,iPivot,WORK,&prod,&info);
    
    /* A is now the inverse */
    free(iPivot);
    free(WORK);    
}

/* Finds distance between two points */
double dist(double *pointsi, double *pointsj)
{
    return(sqrt(pow(pointsj[0]-pointsi[0],2)+pow(pointsj[1]-pointsi[1],2)+pow(pointsj[2]-pointsi[2],2)));
}

/* Calculates the integral over two 1s orbtials for Slater functions a
 * well as the integral between a point charge and each 1s orbital. */
void gammasm(double *za, double *zb, double *r, double *gambfa, double *gamfafb)
{
    if(r[0] > 1e-10) {
        /*Setup Constants*/
        double rr;
        double expzetaA,expzetaB;
        rr=1/r[0];
        
        /*fir*/
        expzetaA=exp(-2*za[0]*r[0]);
        expzetaB=exp(-2*zb[0]*r[0]);
        
        /*calculate*/
        gambfa[0] = (expzetaA*(za[0]+rr))+rr;
        
        if (fabs(za[0]-zb[0]) < 1e-8) {
            double zr;
            zr = za[0]*r[0];
            gamfafb[0] = (rr*(1+zr*(1.375+zr*(0.75+(zr/6))))*expzetaA);
        }
        else {
            double rzetasum,rzetadiff,zc1,zc2,zc3,zc4;
            rzetasum = 1/(za[0]+zb[0]);
            rzetadiff = 1/(za[0]-zb[0]);
            zc1 = pow(rzetasum,2)*pow(rzetadiff,2)*pow(zb[0],4)*za[0];
            zc2 = pow(rzetasum,2)*pow(rzetadiff,2)*pow(za[0],4)*zb[0];
            zc3 = pow(rzetasum,3)*pow(rzetadiff,3)*pow(zb[0],4)*(3*pow(za[0],2)-pow(zb[0],2));
            zc4 = pow(rzetasum,3)*pow(rzetadiff,3)*pow(za[0],4)*(3*pow(zb[0],2)-pow(za[0],2));
            
            gamfafb[0] = (-(zc1+zc3*rr)*expzetaA-(zc2+zc4*rr)*expzetaB);
        }
    }
    else {
        gamfafb[0] = 0;
        gambfa[0] = 0;
    }
}

/* Calculates the energy correction due to charge from the SM
 * electronegativity piece of the SM potential */
double electroneg(double *system)
{
    /* Constants, All in the form [Al O] */
    double J[] = {10.328655, 14.035715};
    double Xi[] = {0.968438, 2.143957};
    double chi[] = {0, 5.484763}; /*eV*/
    double Z[] = {0.746759, 0}; /*a.u.*/
    double kc = 14.4; /*eV \AA e^-2*/
    double chij[CLUSTER] = {0};
    LP_INT clust = CLUSTER, one = 1; /* For sending to LAPACK routines */
    
    /* Variables */
    double X[CLUSTER] = {0}, q[CLUSTER] = {0};
    double V[CLUSTER][CLUSTER] = {0}, invV[CLUSTER][CLUSTER] = {0}, V1[CLUSTER][CLUSTER] = {0};
    double Vchij[CLUSTER] = {0}, Vq[CLUSTER] = {0};
    double tmplhs = 0, tmprhs = 0, nmu = 0;
    double rij, gamjfi, gamfifj, jloop;
    int i, j, iatom, jatom, delta;
    
    
    /* Set chij */
    for (i = 0; i<CLUSTER; i++)
    {
        if (i != CLUSTER-1)
        {
            /* Aluminium */
            chij[i] = chi[1]; /* THIS IS ACTALLY BACKWARDS TO FIT WITH MATLAB VERSION. CHECK E.N. NOTES. */
        }
        else
        {
            /* Oxygen */
            chij[i] = chi[0];
        }
    }
    
    /* Solve B12, B13 and system neutral mu */
    for (i = 0; i<CLUSTER; i++)
    {
        if (i != CLUSTER-1)
        {
            /* Aluminium */
            iatom = 0;
        }
        else
        {
            /* Oxygen */
            iatom = 1;
        }
        jloop = 0;
        for (j = 0; j<CLUSTER; j++)
        {
            if (j != CLUSTER-1)
            {
                /* Aluminium */
                jatom = 0;
            }
            else
            {
                /* Oxygen */
                jatom = 1;
            }
            
            /* For V */
            if (i == j)
            {
                delta = 1;
            }
            else
            {
                delta = 0;
            }
            rij = dist(system+(3*i), system+(3*j)); /* Calculate distance */
            gammasm(&Xi[iatom],&Xi[jatom],&rij,&gamjfi,&gamfifj); /* call slater orbital solver */
            V[j][i] = J[iatom]*delta+kc*gamfifj; /*V was being built the wrong way, so indexes have been flipped (C vs MATLAB indexing)*/
            
            /* For X */
            if (j != i)
            {
                jloop += kc*Z[jatom]*(gamjfi-gamfifj);
            }
        }
        X[i] = chi[iatom]+jloop;
    }
    
    /* build copies (Some LAPACK functions use input variables as temp vars and will thus change them) */
    memcpy(&Vchij, &chij, CLUSTER*sizeof(double));
    memcpy(&invV, &V, PROD*sizeof(double));
    memcpy(&V1, &V, PROD*sizeof(double));
    /* V\chij */
    matRightDivide(&V1[0][0],&Vchij[0],clust,one); 
    /* inverse(V) */
    matInverse(&invV[0][0], clust);
    
    /* Calulate nuetmu */
    /*neutmu = sum(V\chij)/sum(sum(inv(V)));*/
   
    for (i = 0; i<CLUSTER; i++)
    {
        tmplhs += Vchij[i]; /*sum(V\chij)*/
        for (j = 0; j<CLUSTER; j++)
        {
            tmprhs += invV[j][i]; /*sum(sum(inv(V)))*/
        }
    }
    nmu = tmplhs/tmprhs;

    /* Calculate q */
    /* q = V\(neutmu-chij); */
    
    for (i = 0; i<CLUSTER; i++)
    {
        q[i] = nmu-chij[i]; /*nmu-chij*/
    }
    
    memcpy(&V1, &V, PROD*sizeof(double));
    matRightDivide(&V1[0][0],&q[0],clust,one);
     
    /* Calculate Ess */
    /* Ees = sum(q.*X)+0.5.*sum(V*q.*q); */
    
    /* Find V*q*/
    matMultiply(&V[0][0], &q[0], &Vq[0], CLUSTER, 1); //BLAS
    
    tmplhs = 0;
    tmprhs = 0;
    for (i = 0; i<CLUSTER; i++)
    {
        tmplhs += q[i]*X[i]; /* sum(q.*X) */
        tmprhs += Vq[i]*q[i]; /* sum(V*q.*q) */
    }
    
    return(tmplhs+0.5*tmprhs); /* Ees */
}

/* Builds SM Potential */
dcomp makepot(double *points)
{
    /* Variables for cluster */
    int i, j, l;
     
    /* Variables for SM */
    double r;
    double UBuck[CLUSTER][CLUSTER][3] = {0}, URyd[CLUSTER][CLUSTER][3] = {0}, rhoEAM[CLUSTER][CLUSTER][2] = {0}, UEAM[CLUSTER][CLUSTER][2] = {0};
    double UBuckTot = 0, URydTot = 0, UEAMTot = 0, UElecNegTot = 0;
    
    /* SM Constants
     * everything in ev or AA
     * From Al2O3 SM potential
     * Buckingham. Values correspond to:
     * Al	core	Al	core
     * Al	core	O	core
     * O	core	O	core */
    double ABuck[] = {4.474755,62.933909,3322.784218};
    double rhoBuck[] = {0.991317,0.443658,0.291065};
    double CBuck[] = {0,0,0};
    double RangeBuck[] = {0,20};
    
    /* Rydberg. Values correspond to:
     * Al	core	Al	core
     * Al	core	O	core
     * O	core	O	core */
    double ARyd[] = {0.159472,0.094594,1.865072};
    double BRyd[] = {5.949143672000,9.985407051900,16.822405075464};
    double r0Ryd[] = {3.365875,2.358570,2.005092};
    double RangeRyd[] = {0,20};
    
    /* Many Body EAM. Values correspond to:
     * Al	core	Al	core
     * Al	core	O	core
     * O	core	O	core */
    double RangeEAM[] = {0, 20};
    /* Values correspond to Al core; O core */
    double AEAM[] = {0.147699,1.000000};
    double BEAM[] = {2.017519,6.871329};
    double r0EAM[] = {3.365875,2.005092};
    int nEAM = 0;
    
    /* Values correspond to Al core; O core */
    double AfuncEAM[] = {1.987699,2.116850};
    
    
    for (i = 0; i < CLUSTER; i++) /* over mesh in x */
    {
        for (j = 0; j < CLUSTER; j++) /* over mesh in y */
        {
            if (i != j)
            {
                /* find r */
                r = dist(points+(3*i), points+(3*j));
                
                /* if statements check range. If r is not in the range the
                 * potential is not calculated */
                if ((r >= RangeBuck[0]) && (r <= RangeBuck[1]))
                {
                    UBuck[i][j][0] = ABuck[0]*exp(-r/rhoBuck[0])-(CBuck[0]/pow(r,6));
                    UBuck[i][j][1] = ABuck[1]*exp(-r/rhoBuck[1])-(CBuck[1]/pow(r,6));
                    UBuck[i][j][2] = ABuck[2]*exp(-r/rhoBuck[2])-(CBuck[2]/pow(r,6));
                }
                if ((r >= RangeRyd[0]) && (r <= RangeRyd[1]))
                {
                    URyd[i][j][0] = -ARyd[0]*(1+BRyd[0]*((r/r0Ryd[0])-1))*exp(-BRyd[0]*((r/r0Ryd[0])-1));
                    URyd[i][j][1] = -ARyd[1]*(1+BRyd[1]*((r/r0Ryd[1])-1))*exp(-BRyd[1]*((r/r0Ryd[1])-1));
                    URyd[i][j][2] = -ARyd[2]*(1+BRyd[2]*((r/r0Ryd[2])-1))*exp(-BRyd[2]*((r/r0Ryd[2])-1));
                }
                if ((r >= RangeEAM[0]) && (r <= RangeEAM[1]))
                {
                    rhoEAM[i][j][0] = AEAM[0]*pow(r,nEAM)*exp(-BEAM[0]*(r-r0EAM[0]));
                    rhoEAM[i][j][1] = AEAM[1]*pow(r,nEAM)*exp(-BEAM[1]*(r-r0EAM[1]));
                    /* check species of i and j to find what EAM params to use. */
                    
                    /* UEAM values: 0=Al Al, 1=Al O, 2=O O;
                     * rhoEAM/AfuncEAM values: 0=Al, 1=O; */
                    if ((i != CLUSTER-1) && (j != CLUSTER-1))
                    {
                        /*Al Al*/
                        UEAM[i][j][0] = -(AfuncEAM[0]*2)*sqrt(rhoEAM[i][j][0]);
                    }
                    else
                    {
                        /*Al O*/
                        UEAM[i][j][1] = -(AfuncEAM[0]*sqrt(rhoEAM[i][j][0])+AfuncEAM[1]*sqrt(rhoEAM[i][j][1]));
                    }
                    
                }
            }
        }
    }
    
    /* Sum totals */
    for (i = 0; i<CLUSTER; i++)
    {
        for (j = 0; j<CLUSTER; j++)
        {
            for (l = 0; l<3; l++)
            {
                if (!(UBuck[i][j][l] != UBuck[i][j][l])) /* NaN Check */
                {
                    UBuckTot += UBuck[i][j][l];
                }
                URydTot += URyd[i][j][l];
                if (l<2)
                {
                    UEAMTot += UEAM[i][j][l]; /* EAM only 2 dimensions */
                }
            }
        }
    }
    
    UElecNegTot = electroneg(points);

    //Curently potential is all in eV. Check to see if total is > 0 in \mueV, 
    //truncate if so. Return total in GeV for scaling reasons.
    if (((UBuckTot + URydTot + UEAMTot + UElecNegTot)*1e6) > 0.0)
    {
      return (dcomp)0.0;
    }
    else
    {
      return (dcomp)(UBuckTot + URydTot + UEAMTot + UElecNegTot); // eV *1e6; //return result in \mueV
    }

}

dcomp mexHatPotential(double dx, double dy, double dz)
  //dx,dy,dz are the coords of the system centered on the simulation volume
  //Essentially oxycube. Assume NUM is a global. Man need to #include something though
{
    //Top Al, Bottom Al, 3rd nn Al (X left), 4th nn Al (X right), 5th nn Al (Y left), 6th nn Al (Y right), Oxygen Position
    //double cluster[][3] = { { 0.0, 0.0, ALZ }, { 0.0, 0.0, -ALZ }, { ALX, 0.0, 0.0 }, { -ALX, 0.0, 0.0 }, { 0.0, ALY, 0.0 }, { 0.0, -ALY, 0.0 }, { ALX*1.5, ALY*1.5, ALZ*1.5 }, { -ALX*1.5, ALY*1.5, ALZ*1.5 }, { ALX*1.5, -ALY*1.5, ALZ*1.5 }, { -ALX*1.5, -ALY*1.5, ALZ*1.5 }, { ALX*1.5, ALY*1.5, -ALZ*1.5 }, { -ALX*1.5, ALY*1.5, -ALZ*1.5 }, { ALX*1.5, -ALY*1.5, -ALZ*1.5 }, { -ALX*1.5, -ALY*1.5, -ALZ*1.5 }, { dx, dy, dz } };
    //NO CAGE 7 49
    double cluster[][3] = { { 0.0, 0.0, ALZ }, { 0.0, 0.0, -ALZ }, { ALX, 0.0, 0.0 }, { -ALX, 0.0, 0.0 }, { 0.0, ALY, 0.0 }, { 0.0, -ALY, 0.0 }, { dx, dy, dz } };
    //GR2
    //double cluster[][3] = { { 0.0, 0.0, ALZ }, { 0.0, 0.0, -ALZ }, { ALX, 0.0, 0.0 }, { -ALX, 0.0, 0.0 }, { 0.0, ALY, 0.0 }, { 0.0, -ALY, 0.0 }, { 0, 0, -12.3 }, { -6.15, 0, -10.652 }, { -10.652, 0, -6.15 }, { -12.3, 0, 0 }, { -10.652, 0, 6.15 }, { -6.15, 0, 10.652 }, { 0, 0, 12.3 }, { 0, 0, -12.3 }, { -3.075, -5.3261, -10.652 }, { -5.3261, -9.225, -6.15 }, { -6.15, -10.652, 0 }, { -5.3261, -9.225, 6.15 }, { -3.075, -5.3261, 10.652 }, { 0, 0, 12.3 }, { 0, 0, -12.3 }, { 3.075, -5.3261, -10.652 }, { 5.3261, -9.225, -6.15 }, { 6.15, -10.652, 0 }, { 5.3261, -9.225, 6.15 }, { 3.075, -5.3261, 10.652 }, { 0, 0, 12.3 }, { 0, 0, -12.3 }, { 6.15, 0, -10.652 }, { 10.652, 0, -6.15 }, { 12.3, 0, 0 }, { 10.652, 0, 6.15 }, { 6.15, 0, 10.652 }, { 0, 0, 12.3 }, { 0, 0, -12.3 }, { 3.075, 5.3261, -10.652 }, { 5.3261, 9.225, -6.15 }, { 6.15, 10.652, 0 }, { 5.3261, 9.225, 6.15 }, { 3.075, 5.3261, 10.652 }, { 0, 0, 12.3 }, { 0, 0, -12.3 }, { -3.075, 5.3261, -10.652 }, { -5.3261, 9.225, -6.15 }, { -6.15, 10.652, 0 }, { -5.3261, 9.225, 6.15 }, { -3.075, 5.3261, 10.652 }, { 0, 0, 12.3 }, { 0, 0, -12.3 }, { -6.15, 0, -10.652 }, { -10.652, 0, -6.15 }, { -12.3, 0, 0 }, { -10.652, 0, 6.15 }, { -6.15, 0, 10.652 }, { 0, 0, 12.3 }, { dx, dy, dz } };
    //GR1.5 SPHERE 56, 3136
    //double cluster[][3] = { { 0.0, 0.0, ALZ }, { 0.0, 0.0, -ALZ }, { ALX, 0.0, 0.0 }, { -ALX, 0.0, 0.0 }, { 0.0, ALY, 0.0 }, { 0.0, -ALY, 0.0 }, {0, 0, -9.225 }, { -4.6125, 0, -7.9891 }, { -7.9891, 0, -4.6125 }, { -9.225, 0, 0 }, { -7.9891, 0, 4.6125 }, { -4.6125, 0, 7.9891 }, { 0, 0, 9.225 }, { 0, 0, -9.225 }, { -2.3062, -3.9945, -7.9891 }, { -3.9945, -6.9188, -4.6125 }, { -4.6125, -7.9891, 0 }, { -3.9945, -6.9188, 4.6125 }, { -2.3062, -3.9945, 7.9891 }, { 0, 0, 9.225 }, { 0, 0, -9.225 }, { 2.3063, -3.9945, -7.9891 }, { 3.9945, -6.9188, -4.6125 }, { 4.6125, -7.9891, 0 }, { 3.9945, -6.9188, 4.6125 }, { 2.3063, -3.9945, 7.9891 }, { 0, 0, 9.225 }, { 0, 0, -9.225 }, { 4.6125, 0, -7.9891 }, { 7.9891, 0, -4.6125 }, { 9.225, 0, 0 }, { 7.9891, 0, 4.6125 }, { 4.6125, 0, 7.9891 }, { 0, 0, 9.225 }, { 0, 0, -9.225 }, { 2.3063, 3.9945, -7.9891 }, { 3.9945, 6.9188, -4.6125 }, { 4.6125, 7.9891, 0 }, { 3.9945, 6.9188, 4.6125 }, { 2.3063, 3.9945, 7.9891 }, { 0, 0, 9.225 }, { 0, 0, -9.225 }, { -2.3062, 3.9945, -7.9891 }, { -3.9945, 6.9188, -4.6125 }, { -4.6125, 7.9891, 0 }, { -3.9945, 6.9188, 4.6125 }, { -2.3062, 3.9945, 7.9891 }, { 0, 0, 9.225 }, { 0, 0, -9.225 }, { -4.6125, 0, -7.9891 }, { -7.9891, 0, -4.6125 }, { -9.225, 0, 0 }, { -7.9891, 0, 4.6125 }, { -4.6125, 0, 7.9891 }, { 0, 0, 9.225 }, { dx, dy, dz } };
    /*Run potential builder*/
    return makepot(&cluster[0][0]);
}
