/* mexHatPotential.cpp
 * No frills potential builder designed for 3D calculations. dx independant, N independent.
 * Pulled from Matlab Mex routines 
 * Compile:
 * requrired flags for lapack:
 * -llapack -lblas
 * Written by Tim DuBois 24/09/12 
 * Updated to N independent version 21/10/13*/
#include <cmath>
#include <string>
#include <iostream>
#include <complex>

#if defined(RAIJIN)
#include <mkl.h>
#define LP_INT MKL_INT
#elif defined(TRIFID)
#define LP_INT int
#else
#error "One of RAIJIN or TRIFID must be defined"
#endif
//#include <Accelerate/Accelerate.h> //For Lapack & Blas

using namespace std;

#include "mexHatPotential.h"
#include "mpisolve.h"


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
double dist(double *pointsi, double *pointsj, int sizeC)
{
    return(sqrt(pow(*pointsj-*pointsi,2)+pow(*(pointsj+(sizeC*1))-*(pointsi+(sizeC*1)),2)+pow(*(pointsj+(sizeC*2))-*(pointsi+(sizeC*2)),2)));
}

/* Calculates the integral over two 1s orbtials for Slater functions as
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
double electroneg(double *system, double *species, int sizeC)
{
    /* Constants, All in the form [Al O] */
    double J[] = {10.328655, 14.035715};
    double Xi[] = {0.968438, 2.143957};
    double chi[] = {0, 5.484763}; /*eV*/
    double Z[] = {0.746759, 0}; /*a.u.*/
    double kc = 14.4; /*eV \AA e^-2*/
    LP_INT clust = sizeC, one = 1; /* For sending to LAPACK routines */

    /*Virtually allocated Constant*/
    double *chij;
    
    /* Variables */
    double tmplhs = 0, tmprhs = 0, nmu = 0;
    double rij, gamjfi, gamfifj, jloop;
    int i, j, iatom, jatom, delta;
    
    /*Virtually allocated Variables*/
    double *X, *q, *Vchij, *Vq, *V, *invV, *V1;
    
    /*Allocate memory*/    
    chij = (double *)calloc(sizeC,sizeof(double));
    X = (double *)calloc(sizeC,sizeof(double));
    q = (double *)calloc(sizeC,sizeof(double));
    Vchij = (double *)calloc(sizeC,sizeof(double));
    Vq = (double *)calloc(sizeC,sizeof(double));
    /*must use contiguous arrays for memcpy operations (cannot use V[x][y] types)*/
    V = (double *)calloc(sizeC*sizeC,sizeof(double));
    invV = (double *)calloc(sizeC*sizeC,sizeof(double));
    V1 = (double *)calloc(sizeC*sizeC,sizeof(double));
    
    /* Set chij */
    for (i = 0; i<sizeC; i++)
    {
        if (species[i] == 1) 
        {
            /* Aluminium */
            chij[i] = chi[0]; 
        }
        else
        {
            /* Oxygen */
            chij[i] = chi[1];
        }
    }
    
    /* Solve B12, B13 and system neutral mu */
    for (i = 0; i<sizeC; i++)
    {
        if (species[i] == 1)
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
        for (j = 0; j<sizeC; j++)
        {
            if (species[j] == 1) 
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
            rij = dist(system+i, system+j, sizeC); /* Calculate distance */
            gammasm(&Xi[iatom],&Xi[jatom],&rij,&gamjfi,&gamfifj); /* call slater orbital solver */
            *(V+(sizeC*j)+i) = J[iatom]*delta+kc*gamfifj; /*V[j][i] V was being built the wrong way, so indexes have been flipped (C vs MATLAB indexing)*/
            
            /* For X */
            if (j != i)
            {
                jloop += kc*Z[jatom]*(gamjfi-gamfifj);
            }
        }
        X[i] = chi[iatom]+jloop;
    }
    
    /* build copies (Some LAPACK functions use input variables as temp vars and will thus change them) */
    memcpy(Vchij, chij, sizeC*sizeof(double));
    memcpy(invV, V, (sizeC*sizeC)*sizeof(double));
    memcpy(V1, V, (sizeC*sizeC)*sizeof(double));
    /* V\chij */
    matRightDivide(V1,Vchij,clust,one); 
    /* inverse(V) */
    matInverse(invV, clust);
    
    /* Calulate nuetmu */
    /*neutmu = sum(V\chij)/sum(sum(inv(V)));*/
   
    for (i = 0; i<sizeC; i++)
    {
        tmplhs += Vchij[i]; /*sum(V\chij)*/
        for (j = 0; j<sizeC; j++)
        {
            tmprhs += *(invV+(sizeC*j)+i); /* invV[j][i] sum(sum(inv(V)))*/
        }
    }
    nmu = tmplhs/tmprhs;

    /* Calculate q */
    /* q = V\(neutmu-chij); */
    
    for (i = 0; i<sizeC; i++)
    {
        q[i] = nmu-chij[i]; /*nmu-chij*/
    }
    
    memcpy(V1, V, (sizeC*sizeC)*sizeof(double));
    matRightDivide(V1,q,clust,one);
     
    /* Calculate Ess */
    /* Ees = sum(q.*X)+0.5.*sum(V*q.*q); */
    
    /* Find V*q*/
    matMultiply(V, q, Vq, clust, one);
    
    tmplhs = 0;
    tmprhs = 0;
    for (i = 0; i<sizeC; i++)
    {
        tmplhs += q[i]*X[i]; /* sum(q.*X) */
        tmprhs += Vq[i]*q[i]; /* sum(V*q.*q) */
    }
    
    /*Free allocated variables*/
    free(chij); 
    free(X); 
    free(q); 
    free(Vchij); 
    free(Vq); 
    free(V);
    free(invV);
    free(V1);
        
    return(tmplhs+0.5*tmprhs); /* Ees */
}

/* Builds SM Potential */
dcomp makepot(double *points, double *species, int sizeC)
{
    /* Variables for cluster */
    int i, j, l;
     
    /* Variables for SM */
    double r;
    double UBuckTot = 0, URydTot = 0, UEAMTot = 0, UElecNegTot = 0;
    double rhoEAM[] = {0, 0};

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
    double RangeBuck[] = {0,12}; /*20*/
    
    /* Rydberg. Values correspond to:
     * Al	core	Al	core
     * Al	core	O	core
     * O	core	O	core */
    double ARyd[] = {0.159472,0.094594,1.865072};
    double BRyd[] = {5.949143672000,9.985407051900,16.822405075464};
    double r0Ryd[] = {3.365875,2.358570,2.005092};
    double RangeRyd[] = {0,12}; /*20*/
    
    /* Many Body EAM. Values correspond to:
     * Al	core	Al	core
     * Al	core	O	core
     * O	core	O	core */
    double RangeEAM[] = {0, 8}; /*20*/
    /* Values correspond to Al core; O core */
    double AEAM[] = {0.147699,1.000000};
    double BEAM[] = {2.017519,6.871329};
    double r0EAM[] = {3.365875,2.005092};
    int nEAM = 0;
    
    /* Values correspond to Al core; O core */
    double AfuncEAM[] = {1.987699,2.116850};
    
    for (i = 0; i < sizeC; i++) /* over mesh in x */
    {
        for (j = 0; j < sizeC; j++) /* over mesh in y*/
        {
            if (j > i)
            {
                /* find r */
                r = dist(points+i, points+j, sizeC);
                
//if ((j == 1) && (i == 0)) { 
    //double *pointsi, *pointsj;

    //pointsi = points+i;
    //pointsj = points+j;
    //cout << *pointsj << "," << *(pointsj+8) << "," << *(pointsj+9)  << "," << *(pointsj+10) << "," << *(pointsj+11) << "," << *(pointsj+12) << "," << *(pointsj+13) << "," << *(pointsj+14) << endl;
    //
    //return(sqrt(pow(*pointsj-*pointsi,2)+pow(*(pointsj+(sizeC*1))-*(pointsi+(sizeC*1)),2)+pow(*(pointsj+(sizeC*2))-*(pointsi+(sizeC*2)),2)));
    //cout << sizeC << ", " << species[j] << ", " << *(points) << ", " << *(points+1) << ", " << *(points+2) << endl;
//}
                /* if statements check range. If r is not in the range the
                 * potential is not calculated */
                if ((r >= RangeBuck[0]) && (r <= RangeBuck[1]))
                {
                    if ((species[i] == 1) && (species[j] == 1))
                    {
                        /*Al Al*/
                        UBuckTot += ABuck[0]*exp(-r/rhoBuck[0])-(CBuck[0]/pow(r,6));
                    }
                    else if ((species[i] == 2) && (species[j] == 2))
                    {
                        /*O O*/
                        UBuckTot += ABuck[2]*exp(-r/rhoBuck[2])-(CBuck[2]/pow(r,6));
                    }
                    else
                    {
                        /*Al O*/
                        UBuckTot += ABuck[1]*exp(-r/rhoBuck[1])-(CBuck[1]/pow(r,6));
                    }
                }
                if ((r >= RangeRyd[0]) && (r <= RangeRyd[1]))
                {
                    if ((species[i] == 1) && (species[j] == 1))
                    {
                        /*Al Al*/
                        URydTot += -ARyd[0]*(1+BRyd[0]*((r/r0Ryd[0])-1))*exp(-BRyd[0]*((r/r0Ryd[0])-1));
                    }
                    else if ((species[i] == 2) && (species[j] == 2))
                    {
                        /*O O*/
                        URydTot += -ARyd[2]*(1+BRyd[2]*((r/r0Ryd[2])-1))*exp(-BRyd[2]*((r/r0Ryd[2])-1));
                    }
                    else
                    {
                        /*Al O*/
                        URydTot += -ARyd[1]*(1+BRyd[1]*((r/r0Ryd[1])-1))*exp(-BRyd[1]*((r/r0Ryd[1])-1));
                    }   
                }
                if ((r >= RangeEAM[0]) && (r <= RangeEAM[1]))
                {
                    rhoEAM[0] = AEAM[0]*pow(r,nEAM)*exp(-BEAM[0]*(r-r0EAM[0])); 
                    rhoEAM[1] = AEAM[1]*pow(r,nEAM)*exp(-BEAM[1]*(r-r0EAM[1]));
                    /* check species of i and j to find what EAM params to use. */
                    
                    /* UEAM values: 0=Al Al, 1=Al O, 2=O O;
                     * rhoEAM/AfuncEAM values: 0=Al, 1=O; */
                    if ((species[i] == 1) && (species[j] == 1))
                    {
                        /*Al Al*/
                        UEAMTot += -(AfuncEAM[0]*2)*sqrt(rhoEAM[0]); 
                    }
                    else if ((species[i] == 2) && (species[j] == 2))
                    {
                        /*O O*/
                        UEAMTot += -(AfuncEAM[1]*2)*sqrt(rhoEAM[1]);
                    }
                    else
                    {
                        /*Al O*/
                        UEAMTot += -(AfuncEAM[0]*sqrt(rhoEAM[0])+AfuncEAM[1]*sqrt(rhoEAM[1]));
                    }
                }
            }
        }
    }
   
    UElecNegTot = electroneg(points,species,sizeC);

    //Curently potential is all in eV. Check to see if total is > 0, truncate if so.
    if ((UBuckTot + URydTot + UEAMTot + UElecNegTot) > 0.0) 
    {
      return (dcomp)0.0;
    }
    else
    {
      return (dcomp)(UBuckTot + URydTot + UEAMTot + UElecNegTot)*239.2311; //return result in eV*Scale Offset
    }

}

dcomp mexHatPotential(double dx, double dy, double dz)
  //dx,dy,dz are the coords of the system centered on the simulation volume
  //Essentially oxycube. Assume NUM is a global. 
{
    //Top Al, Bottom Al, 3rd nn Al (X left), 4th nn Al (X right), 5th nn Al (Y left), 6th nn Al (Y right), Oxygen Position
    //double cluster[][3] = { { 0.0, 0.0, ALZ }, { 0.0, 0.0, -ALZ }, { ALX, 0.0, 0.0 }, { -ALX, 0.0, 0.0 }, { 0.0, ALY, 0.0 }, { 0.0, -ALY, 0.0 }, { ALX*1.5, ALY*1.5, ALZ*1.5 }, { -ALX*1.5, ALY*1.5, ALZ*1.5 }, { ALX*1.5, -ALY*1.5, ALZ*1.5 }, { -ALX*1.5, -ALY*1.5, ALZ*1.5 }, { ALX*1.5, ALY*1.5, -ALZ*1.5 }, { -ALX*1.5, ALY*1.5, -ALZ*1.5 }, { ALX*1.5, -ALY*1.5, -ALZ*1.5 }, { -ALX*1.5, -ALY*1.5, -ALZ*1.5 }, { dx, dy, dz } };
    //NO CAGE 7 49
    
    
   // cluster = (double *)calloc(sizeC*3,sizeof(double));
    //GR2
    if (CLUSTER == 1) {
        //generate from cluster data
        dcomp V;
        
        double *cluster = (double *)calloc(CLUSTSIZE*3,sizeof(double)); 
        for (int j = 0; j<CLUSTSIZE; j++)
        {
            if (j < CLUSTSIZE-1)
            {
                // Add atoms 
                *(cluster+(CLUSTSIZE*0)+j) = *(clust+(CLUSTSIZE*0)+j); 
                *(cluster+(CLUSTSIZE*1)+j) = *(clust+(CLUSTSIZE*1)+j); 
                *(cluster+(CLUSTSIZE*2)+j) = *(clust+(CLUSTSIZE*2)+j); 
            }
            else
            {
                // Delocalised Oxygen 
                *(cluster+(CLUSTSIZE*0)+j) = dx; 
                *(cluster+(CLUSTSIZE*1)+j) = dy; 
                *(cluster+(CLUSTSIZE*2)+j) = dz; 
            }                
        }

        V = makepot(cluster,clustSpecies,CLUSTSIZE);
        free(cluster);
        return V;
    } else {
        //generate from AL{X,Y,Z} parameters, assuming species is all aluminium
        double species[] = {1, 1, 1, 1, 1, 1, 2};
        double cluster[] = {0, 0,  ALX,  -ALX, 0, 0,  dx, 0, 0, 0, 0,  ALY,  -ALY, dy,  ALZ,  -ALZ, 0, 0, 0, 0, dz};
        int sizeC = sizeof(species)/sizeof(double);
        return makepot(&cluster[0],&species[0],sizeC);
    }
    //double cluster[][3] = { { 0.0, 0.0, ALZ }, { 0.0, 0.0, -ALZ }, { ALX, 0.0, 0.0 }, { -ALX, 0.0, 0.0 }, { 0.0, ALY, 0.0 }, { 0.0, -ALY, 0.0 }, { 0, 0, -12.3 }, { -6.15, 0, -10.652 }, { -10.652, 0, -6.15 }, { -12.3, 0, 0 }, { -10.652, 0, 6.15 }, { -6.15, 0, 10.652 }, { 0, 0, 12.3 }, { 0, 0, -12.3 }, { -3.075, -5.3261, -10.652 }, { -5.3261, -9.225, -6.15 }, { -6.15, -10.652, 0 }, { -5.3261, -9.225, 6.15 }, { -3.075, -5.3261, 10.652 }, { 0, 0, 12.3 }, { 0, 0, -12.3 }, { 3.075, -5.3261, -10.652 }, { 5.3261, -9.225, -6.15 }, { 6.15, -10.652, 0 }, { 5.3261, -9.225, 6.15 }, { 3.075, -5.3261, 10.652 }, { 0, 0, 12.3 }, { 0, 0, -12.3 }, { 6.15, 0, -10.652 }, { 10.652, 0, -6.15 }, { 12.3, 0, 0 }, { 10.652, 0, 6.15 }, { 6.15, 0, 10.652 }, { 0, 0, 12.3 }, { 0, 0, -12.3 }, { 3.075, 5.3261, -10.652 }, { 5.3261, 9.225, -6.15 }, { 6.15, 10.652, 0 }, { 5.3261, 9.225, 6.15 }, { 3.075, 5.3261, 10.652 }, { 0, 0, 12.3 }, { 0, 0, -12.3 }, { -3.075, 5.3261, -10.652 }, { -5.3261, 9.225, -6.15 }, { -6.15, 10.652, 0 }, { -5.3261, 9.225, 6.15 }, { -3.075, 5.3261, 10.652 }, { 0, 0, 12.3 }, { 0, 0, -12.3 }, { -6.15, 0, -10.652 }, { -10.652, 0, -6.15 }, { -12.3, 0, 0 }, { -10.652, 0, 6.15 }, { -6.15, 0, 10.652 }, { 0, 0, 12.3 }, { dx, dy, dz } };
    //GR1.5 SPHERE 56, 3136
    //double cluster[][3] = { { 0.0, 0.0, ALZ }, { 0.0, 0.0, -ALZ }, { ALX, 0.0, 0.0 }, { -ALX, 0.0, 0.0 }, { 0.0, ALY, 0.0 }, { 0.0, -ALY, 0.0 }, {0, 0, -9.225 }, { -4.6125, 0, -7.9891 }, { -7.9891, 0, -4.6125 }, { -9.225, 0, 0 }, { -7.9891, 0, 4.6125 }, { -4.6125, 0, 7.9891 }, { 0, 0, 9.225 }, { 0, 0, -9.225 }, { -2.3062, -3.9945, -7.9891 }, { -3.9945, -6.9188, -4.6125 }, { -4.6125, -7.9891, 0 }, { -3.9945, -6.9188, 4.6125 }, { -2.3062, -3.9945, 7.9891 }, { 0, 0, 9.225 }, { 0, 0, -9.225 }, { 2.3063, -3.9945, -7.9891 }, { 3.9945, -6.9188, -4.6125 }, { 4.6125, -7.9891, 0 }, { 3.9945, -6.9188, 4.6125 }, { 2.3063, -3.9945, 7.9891 }, { 0, 0, 9.225 }, { 0, 0, -9.225 }, { 4.6125, 0, -7.9891 }, { 7.9891, 0, -4.6125 }, { 9.225, 0, 0 }, { 7.9891, 0, 4.6125 }, { 4.6125, 0, 7.9891 }, { 0, 0, 9.225 }, { 0, 0, -9.225 }, { 2.3063, 3.9945, -7.9891 }, { 3.9945, 6.9188, -4.6125 }, { 4.6125, 7.9891, 0 }, { 3.9945, 6.9188, 4.6125 }, { 2.3063, 3.9945, 7.9891 }, { 0, 0, 9.225 }, { 0, 0, -9.225 }, { -2.3062, 3.9945, -7.9891 }, { -3.9945, 6.9188, -4.6125 }, { -4.6125, 7.9891, 0 }, { -3.9945, 6.9188, 4.6125 }, { -2.3062, 3.9945, 7.9891 }, { 0, 0, 9.225 }, { 0, 0, -9.225 }, { -4.6125, 0, -7.9891 }, { -7.9891, 0, -4.6125 }, { -9.225, 0, 0 }, { -7.9891, 0, 4.6125 }, { -4.6125, 0, 7.9891 }, { 0, 0, 9.225 }, { dx, dy, dz } };
    /*Run potential builder*/
}
