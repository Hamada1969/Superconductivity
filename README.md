# Superconductivity
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// In this code, We use Lanczos Routine (Lanczos.c) for Finding the Eigenvalueof a Superconductor in 2D
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/resource.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Define external BLAS library functions which are used in this program ///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double dcopy_(unsigned long long int *N,double *SX, int *INCX, double *SY, int *INCY); // finds SY=SX 
double ddot_(unsigned long long int *N, double *SX, int *INCX, double *SY, int *INCY); // finds SX.SY
double dnrm2_(unsigned long long int *N, double *SX, int *INCX); // finds SX.SX 
double dscal(unsigned long long int *N, double *SA, double *SX, int *INCX); // finds SX=SA*SX 
double daxpy_(unsigned long long int *N, double *SA, double *SX, int *INCX, double *SY, int *INCY); // finds SY=SA*SX+SY
double dsymv_( char *uplo, unsigned long long int *n, double *alpha, double *m, int *llda, double *w, int *incx, double *beta, double *y, int *incy); // This finds y=m*x+alpha*y
double dstevr_(char *jobz, char *range, unsigned long long int *N, double *d, double *e, double *vl, double *vu, int *il, int *iu, double *abstol, int *M, double *w, double *z, int *lldz, int *isuppz, double * work, int *lwork, int * iwork, int *liwork, int *info); 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Start main() Program: which finds the eigenvalues of a symmetric tridiagonal matrix     /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main ()
{
int Ne, Nsi, Nehalf;                             // size and state parameters 
unsigned long long int Nst,Nst1;       //  number of states 
double t, J, J1, dist;                            // Hamiltonian model parameters (here: t-J Hubbard)
int L;                                                  // length of lattice (here hexagonal lattice)
int i,j,k,l;                                            // general counters for loops 
unsigned long long int k1,k2,k3,k4;  // input parameters 
char line[100];                                    // this stores the user input strings
printf("What is the length of the lattice?\n");fgets(line, sizeof(line),stdin);sscanf(line, "%d", &L);
printf("How many electrons?\n");fgets(line, sizeof(line),stdin);sscanf(line, "%d", &Ne);
printf("What is t?\n");fgets(line, sizeof(line),stdin);sscanf(line, "%lf", &t);
printf("What is J?\n");fgets(line, sizeof(line),stdin);sscanf(line, "%lf", &J);
/* find number of sites (to be changed for larger dimensions*/
Nsi = 2*L*L;                  // the number of sites for a honeycomb hexagonal lattice
Nehalf = Ne/2;
 /* find the number of states */
unsigned long long int Nsifac, Nsifac1,Nehalffac;
// find Nsi!
j=Nsi-1; Nsifac=Nsi;
while(j>Nsi-Nehalf) { Nsifac=Nsifac*j; j--; }
// find Nsi1!
j=Nsi-1; Nsifac1=Nsi;
while(j>Nsi-Ne) { Nsifac1=Nsifac1*j; j--; }   // find (Ne/2)!
j=Nehalf-1;
Nehalffac=Nehalf;
while(j>0) {Nehalffac=Nehalffac*j; j--;}
Nst1= Nsifac/Nehalffac;
Nst = Nsifac1/(Nehalffac*Nehalffac);
unsigned long long int num, num1, num2, num3;    // create state numbers and look up tables used in binary labeling 
num = pow(2,Nsi); // find the maximum number up to which binary numbers will be created 
// define maximum number of lanczos steps, to avoid using too much memory. 
int steps;
printf("How many Lanczos Steps?\n"); fgets(line, sizeof(line),stdin); sscanf(line, "%d", &steps);
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Create the binary states 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
unsigned long long int * statenumhalf, * statenum;
statenumhalf = (unsigned long long int *)
malloc(Nst1*sizeof(unsigned long long int));
statenum = (unsigned long long int *)
malloc(Nst*sizeof(unsigned long long int));
int * bin,* bin_change, elecct;
unsigned long long int tempnum, statenumrest
#include "bin_create_nodoubocc.c"     // the states have been created, 
////////////////////////////////////////////////////////////////////////////////////////////////// 
// Lanczos routine begins here 
///////////////////////////////////////////////////////////////////////////////////////////
int lda = Nst; int INCX = 1; int INCY = 1; /* variables used in the BLAS & LAPACK routines */
int il, iu, M, lwork, liwork, *isuppz, *iwork, info, ldz1, lwork1, liwork1;
unsigned long long int ldz;
double vl, vu, *work, abstol, *ZT, WT;
char jobz, range, uplo;
jobz = ’V’;range = ’I’;uplo = ’U’;
il = 1; iu = 1; M = iu - il + 1;
abstol = 0.000001;
////////////////////////////////////////////////////////////////////
// Define vectors used in Lanczos
////////////////////////////////////////////////////////////////////
double *w;   // the initial normalised trial vector
double *a;    // the diagonal matrix element vector of the tridiagonal matrix
double *b;   // the diagonal matrix element vector of the tridiagonal matrix
double *q;   // the initial "eigenvector"
double *y;   // an extra vector used in the
// matrix multiplication
double *t1; // an extra vector used in the code to copy vectors
double * cons; // a vector containing some constants used in Lanczos routine
double *Hamvec; // the row matrix used in the Hamiltonian generation
// The parameters used in the Hamiltonian code
int oldsite, newsite, i1, i2, j1, phase1,phase2,phaset;
unsigned long long int upper,lower,search;
int * bin_sign,  * dist_arr;
// Other parameters
unsigned long long int lan; // the lanczos step count
double * grndst; // the groundstate vector
double diff; // the difference between the eigenvalues
double check, norm, norm1;
int comp;
// Define  the allowable error
printf("What is difference in the eigenvalues?\n");fgets(line, sizeof(line),stdin); sscanf(line, "%lf", &diff);
// specify the step at which to start diagonalization
printf("When to start diagonalizing?\n");fgets(line, sizeof(line),stdin);sscanf(line, "%d", &comp);
// initialize the random number generator
gsl_rng *NUM= gsl_rng_alloc(gsl_rng_taus);
gsl_rng *PN= gsl_rng_alloc(gsl_rng_taus2);
// seed the rng
int SEED, SEED1; // the seeds for the generators
printf("Enter an integer\n");fgets(line, sizeof(line),stdin);sscanf(line, "%d", &SEED);
printf("Enter another integer\n");fgets(line, sizeof(line),stdin);sscanf(line, "%d", &SEED1);
// open a file to save debugging information
char data[100]; char data1[100]; char data2[100];
FILE *fp;
sprintf(data, "t_J_Honeycomb_LanczosRunInfo_L=%d_Ne=%d_J=%2.3lf_tol=%2.1e.txt",L,Ne,J,diff);
fp=fopen(data, "w");

#include "rand.c"                             // gsl will now generate a random normalised trial vector
#include "distance_honeycomb.c"  // This file defines the lattice
////////////////////////////////////////////////////////////////////////////////////////////
// Lanczos routine begins here
/////////////////////////////////////////////////////////////////////////////////////////////
// allocate space for the diagonal and subdiagonal matrix elements
a = (double *) malloc (steps*sizeof(double));
b = (double *) malloc (steps*sizeof(double));
// vectors which will be constantly rewritten over and over again
cons = (double *) calloc (2, sizeof(double));
y = (double *) calloc (Nst, sizeof(double));
q = (double *) calloc (Nst, sizeof(double));
#include "tJ_mat_unref_vec.c"       // Hamiltonian generation, 
cons[0] = 1.0;
daxpy_(&Nst, &cons[0], y, &INCX, q, &INCY);          //  q[i]=q[i]+y[i](completing q=q+M*w)
free(y);
a[0] = ddot_(&Nst, w, &INCX, q, &INCY);                 //  a[k]=w[i]*q[i]
cons[1] = -a[0];
daxpy_(&Nst, &cons[1], w, &INCX, q, &INCY);        //  q[i]=q[i]-a[k]*w[i]
b[0] = dnrm2_(&Nst, q, &INCX);                                 //  b[k]=q[i]*q[i]
free(cons);
for ( lan=1; lan<steps; lan++)
{
cons = (double *) calloc (4, sizeof(double));
t1 = (double *) calloc (Nst, sizeof(double));
dcopy_(&Nst, w, &INCX, t1, &INCY);                     //t=w[i]
free(w);
w = (double *) calloc (Nst, sizeof(double));
cons[0] = (1/b[lan-1]);
daxpy_(&Nst, &cons[0], q, &INCX, w, &INCY);    //w[i]=q[i]/b[k-1]
free(q);
q = (double *) calloc (Nst, sizeof(double));
cons[1] = -b[lan-1];
daxpy_(&Nst, &cons[1], t1, &INCX, q, &INCY);    //q[i]=-b[k-1]*t
free(t1);
y = (double *) calloc (Nst, sizeof(double));

#include "tJ_mat_unref_vec.c"       // the Hamiltonian generation 
cons[2] = 1.0;
daxpy_(&Nst, &cons[2], y, &INCX, q, &INCY);       //  q[i]=q[i]+y[i]
free(y);
a[lan] = ddot_(&Nst, w, &INCX, q, &INCY);            //  a[k]=w[i]*q[i]
cons[3] = -a[lan];
daxpy_(&Nst, &cons[3], w, &INCX, q, &INCY);     //  q[i]=q[i]-a[k]*w[i]
b[lan] = dnrm2_(&Nst, q, &INCX);                           //  b[k]=q[i]*q[i]
free(cons);
////////////////////////////////////////////////////////////////////////////////////////////
// after a given number of lanczos loops begin finding the groundstate T[i,j]
////////////////////////////////////////////////////////////////////////////////////////////
if (lan>comp)        
{
double *d, *e;                    // diag and sub diag components
d = (double *) calloc (lan, sizeof(double));
e = (double *) calloc (lan-1, sizeof(double));
dcopy_(&lan, a, &INCX, d, &INCY); // Copy diagonal vectors as they will be over written upon diagonalization
k1=lan-1;
dcopy_(&k1, b, &INCX, e, &INCY); // Copy subdiagonal vectors as they will be over written upon diagonalization
lwork1 = 33*lan;
liwork1 = 10*lan;
ldz1= lan;
ZT = (double *) malloc(ldz1*1*sizeof(double));   // This is the Lanczos groundstate 
isuppz = (int *) malloc(2*lan*sizeof(int));
work = (double *) malloc(33*lan*sizeof(double));
iwork = (int *) malloc(10*lan*sizeof(int));
dstevr_(&jobz, &range, &lan, d, e, &vl, &vu, &il, &iu, &abstol, &M, &WT, ZT, &ldz1, isuppz, work, &lwork1, iwork, &liwork1, &info);
free(d); free(e);
//////////////////////////////////////////////////////////////////////////////////////////////////////
// Find the tolerance condition and  save the results to the diagonstic file
///////////////////////////////////////////////////////////////////////////////////////////////////////
check = fabs(b[lan-1]*ZT[lan-1]);        
fprintf(fp,"%lld\t",lan); fprintf(fp,"%+2.10e\t",check); fprintf(fp,"%+2.10e\n",WT);  fflush(fp);
if (check <= diff) { steps = lan; break; }  // set the max steps for the vector code // test the tolerance condition
if(lan%10==0)  // save important info every tenth step in case something happens and code stopped prematurely
{
FILE *fp4;
sprintf(data1, "Z_t_J_Honeycomb_LanczosVector_L=%d_Ne=%d_J=%2.3lf.txt",L,Ne,J);
fp4=fopen(data1, "w"); for(i=0; i<lan; i++) { fprintf(fp4,"%+2.10e\n", ZT[i]); } fclose(fp4);
FILE *fp5;
sprintf(data2, "t_J_Honeycomb_RestartInfo_L=%d_Ne=%d_J=%2.3lf.txt",L,Ne,J);
fp5=fopen(data2, "w");
for(i=0; i<Nst; i++) { fprintf(fp5,"%+2.10e\n", w[i]); }
for(i=0; i<Nst; i++) { fprintf(fp5,"%+2.10e\n", q[i]); }
for(i=0; i<lan; i++) { fprintf(fp5,"%+2.10e\n", a[i]); }
for(i=0; i<lan; i++) { fprintf(fp5,"%+2.10e\n", b[i]); }
fclose(fp5); } // end file save
} } fclose(fp);  // end lan
////////////////////////////////////////////////////////////////////////////////////////////
// Print the groundstate energy to a file
////////////////////////////////////////////////////////////////////////////////////////////
FILE *fp2; 
sprintf(data, "t_J_Honeycomb_GSNRG_L=%d_Ne=%dJ=%2.3lf _tol=%2.1e.txt", L, Ne, J, diff);
fp2=fopen(data, "w"); fprintf(fp2,"%+2.10e\n\n", WT); fclose(fp2);
FILE *fp4;
sprintf(data1, "Z_t_J_Honeycomb_LanczosVector_L=%d_Ne=%d_J=%2.3lf.txt",L,Ne,J);
fp4=fopen(data1, "w"); for(i=0; i<lan; i++) { fprintf(fp4,"%+2.10e\n", ZT[i]); } fclose(fp4);
FILE *fpv; // create input file for vector code
char datav[100];
sprintf(datav, "in1_%d%d%2.1lf",L,Ne,J);
fpv=fopen(datav, "w"); fprintf(fpv, "%d\n%d\n%2.1lf\n%2.1lf\n%d\n%d\n%lld\n", L,Ne,t,J,SEED,SEED1,lan);
for(i=0; i<lan; i++) { fprintf(fpv, "%+2.10e\n",ZT[i]); } fclose(fpv);
FILE *fp5;
sprintf(data2, "t_J_Honeycomb_RestartInfo _L=%d_Ne=%d_J=%2.3lf.txt",L,Ne,J);
fp5=fopen(data2, "w");
for(i=0; i<Nst; i++) { fprintf(fp5,"%+2.10e\n", w[i]); }
for(i=0; i<Nst; i++) { fprintf(fp5,"%+2.10e\n", q[i]); }
for(i=0; i<lan; i++) { fprintf(fp5,"%+2.10e\n", a[i]); }
for(i=0; i<lan; i++) { fprintf(fp5,"%+2.10e\n", b[i]); }  fclose(fp5);
// free vectors used during Lanczos routine
free(w); free(a); free(b); free(q); free(ZT);
return 0;
}////////////////////////////End main()///////////////////////////////

Listing F-2. Generation of the Hamiltonian Matrix

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// tJ_mat_unref_vec.c /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

for (k1=0; k1<Nst; k1++)
{
Hamvec = (double *) calloc (Nst, sizeof(double));
bin = (int *)calloc(2*Nsi,sizeof(int));           // convert state to binary form
num=statenum[k1];
for (i1=0; i1<2*Nsi; i1++) {  if (num % 2 == 0) { num=num/2; } else {bin[i1]=1;num = (num-1)/2;} }
J1=0;   // determine the possibility of a diagonal entry  by analysing the state itself
for (i1=0; i1<Nsi; i1++) {l=(i1+1)%Nsi; J1=J1+(bin[i1]-bin[i1+Nsi])*(bin[l]-bin[l+Nsi]);}

Hamvec[k1] = J*J1/4;    // find the matrix element
// find all off-diagonal elements connected to the state,  consider first the up spins
// search for t-terms in the up spin state
// sort through the binaries to find the other possible arrangments
elecct = 0;
for (i1=0; i1<Nsi; i1++)
{
if (elecct==Nehalf){break;} // first find a site with an electron on it
if (bin[i1]==1) {
// find all possible nearest neighbours by sorting through the distance matrix for non-zero entries
elecct++; // count electrons to avoid sorting through empty sites
oldsite = i1;
for (i2=0; i2<Nsi; i2++) {if(dist_arr[i1*Nsi+i2]==1)
{
// search through the binary state to make sure  there is no electron on the site
if((bin[i2]==0)&&(bin[i2+Nsi]==0))
{
// Form the new state number
newsite =i2; tempnum = statenum[k1] -pow(2,oldsite)+pow(2,newsite);
// Find the state this is connected to
upper = Nst; lower = 0; search =lower+upper;
for (k3=0; k3<Nst/2; k3++)
{
if(search%2==0) {search = search/2; }
else{ search = (search-1)/2; }
if(tempnum == statenum[search]){ num = statenum[search]; k2=search;break;}
else if(tempnum < statenum[search]) { upper = search;search = lower + upper; }
else { lower = search; search = lower + upper; }
}
// Find the sign on the hopping parameter, t
bin_sign = (int *)calloc(2*Nsi,sizeof(int));
for (j1=0; j1<2*Nsi; j1++)
{
if (num % 2 == 0){num=num/2;}
else{bin_sign[j1]=1;num = (num-1)/2;}
}
// Calculate the sign on the hopping term
phase1=0; phase2=0;
for(j1=oldsite; j1<Nsi; j1++){phase1 = phase1 + bin[j1];}
for(j1=newsite; j1<Nsi; j1++){ phase2 = phase2 + bin_sign[j1]; }
free(bin_sign);
phaset = pow(-1,phase1+phase2); Hamvec[k2]=-phaset*t;}
// Search for J-terms
else if ((bin[i2]==0)&&(bin[i2+Nsi]==1))
{
// Look if there is a down spin on nearest neighbour site, make the site it is moving to empty
if(bin[i1+Nsi]==1){continue;}
newsite =i2;
tempnum = statenum[k1] - pow(2,oldsite) + pow(2,newsite) - pow(2,Nsi+newsite) + pow(2,Nsi+oldsite);
// Find the state this is connected to
upper = Nst; lower = 0;
search =lower+upper;
for (k3=0; k3<Nst/2; k3++)
{
if(search%2==0) { search = search/2;}
else { search = (search-1)/2; }
if(tempnum == statenum[search]) { num = statenum[search]; k2=search; break; }
else if(tempnum < statenum[search]) { upper = search; search = lower + upper;}
else { lower = search; search = lower + upper; }
}
Hamvec[k2]=-J/2; } // end down spin find if
else { continue;}
} // end distance if
} // end i2-loop
} // end electron search if bin[i1]
} // end i1-for
// Do this over again for the down spin,  you only need to find t-term
elecct=0;
for (i1=0; i1<Nsi; i1++) {
// find a site with an electron on it
if(elecct==Nehalf){break;}
if (bin[Nsi+i1]==1)
{
// Find all possible nearest neighbours by sorting the distance matrix for non-zero entries
elecct++; oldsite = i1;
for (i2=0; i2<Nsi;i2++) { if(dist_arr[i1*Nsi+i2]==1)
{
// Make sure there is no electron on the site
if((bin[i2]==1)||(bin[i2+Nsi]==1)){continue;}
// Form the new state number
newsite =i2;
tempnum = statenum[k1]-pow(2,Nsi+oldsite)+pow(2,Nsi+newsite);
// Find the new state this is connected to
upper = Nst; lower = 0;
search =lower+upper;
for (k3=0; k3<Nst/2; k3++)
{
if(search%2==0) { search = search/2;}
else { search = (search-1)/2; }
if(tempnum == statenum[search]) { num = statenum[search]; k2=search; break;}
else if(tempnum < statenum[search]) { upper = search; search = lower + upper;}
else { lower = search;search = lower + upper;}
}
// Find the sign on t
bin_sign = (int *)calloc(2*Nsi,sizeof(int));
for (j1=0; j1<2*Nsi; j1++) {
if (num % 2 == 0) {num=num/2; }
else {bin_sign[j1]=1; num = (num-1)/2;}
}
// Calculate the sign on the hopping term
phase1=0; phase2=0;
for(j1=(Nsi+oldsite); j1<2*Nsi; j1++) {phase1 = phase1 + bin[j1];}
for(j1=(Nsi+newsite); j1<2*Nsi; j1++){phase2 = phase2 + bin_sign[j1];}
free(bin_sign);
phaset = pow(-1,phase1+phase2);
Hamvec[k2]=-phaset*t;
} // End distance if
} // End i2-loop
} // End electron search if bin[i1]
} // End i1-for
// Multiply the Hamiltonian row by the current Lanczos vector
y[k1] = ddot_(&Nst, Hamvec, &INCX, w, &INCY);
free(Hamvec); // Free the Hamvec 
free(bin); // Free the bin 
}///////////////////////////////// end main (i - for)/////////////////////////////////////////////////////////////////

Listing F-3. Generation of the Binary States

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// bin_create_nodoubocc.c /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
for (num1=1; num1<num; num1++) { elecct = 0; num2=num1; 
for (i=0; i<Nsi; i++) { 
if (num2 % 2 == 0) { num2=num2/2; } 
else { elecct++; num2 = (num2-1)/2;}
}
// Decide if the statenum will be recorded based on if the number of electrons is proper
k1=0; 
if (elecct == Nehalf) { statenumhalf[k1]=num1; k1++; if(k1==Nst){break;} }
}
// Create the full state by multiplying every statenumhalf  by each other 
k1=0; for(k2=0; k2<Nst1; k2++) { num=statenumhalf[k2]; statenumrest=0; for (k=0; k<Nsi; k++) { if (num % 2 == 0) { num=num/2; } else { statenumrest = statenumrest+pow(2,Nsi+k); num = (num-1)/2; } }
for(k3=0; k3<Nst1; k3++) { num3=statenumhalf[k3]+statenumrest;
// Convert the temp_state to a binary to count the ones
bin = (int *)calloc(2*Nsi,sizeof(int)); 
for (k=0; k<2*Nsi; k++) { if (num3 % 2 == 0) { bin[k]=0; num3=num3/2; } else { bin[k]=1; num3 = (num3-1)/2; } }
// Count the ones
tempnum = 0; for (k4=0; k4<Nsi; k4++) { tempnum = tempnum + bin[k4]*bin[k4+Nsi]; } free(bin); if(tempnum==0) { statenum[k1]=statenumhalf[k3]+statenumrest; k1++; }}}
 free(statenumhalf);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

