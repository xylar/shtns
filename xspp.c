  ////////////////////////////////////
 // XSPP : XShells post-processing //
////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>
// FFTW : spatial derivative is d/dx = ik	(no minus sign !)
#include <fftw3.h>
// GSL for Legendre functions
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>

// parameters for SHT.c
#define NLAT 256
#define LMAX 255
#define NPHI 2
#define MMAX 0
#define MRES 1
// SHT on equal spaced grid + polar points.
#define SHT_EQUAL
#include "SHT.c"

// number of radial grid points.
long int NR,NU;		//  NR: total radial grid points. NU:for velocity field.
long int NG=0;		//  NG: grid points for inner core.

#include "grid.c"

double nu, eta;		// viscosity and magnetic diffusivity.
double Omega0;		// global rotation rate (of outer boundary) => Coriolis force .
double DeltaOmega;	// differential rotation (of inner core)

//#define DEB printf("%s:%u pass\n", __FILE__, __LINE__)
#define DEB (0)

#include "xshells_fields.c"
#include "xshells_io.c"

struct VectField B;
struct PolTor Blm;




void write_vect(char *fn, double *vec, long int N)
{
	FILE *fp; 
	long int i;
	
	fp = fopen(fn,"w");
	for (i=0;i<N;i++) {
		fprintf(fp,"%.6g ",vec[i]);
	}
	fclose(fp);
}

void write_mx(char *fn, double *mx, int N1, int N2)
{
	FILE *fp;
	int i,j;
	
	fp = fopen(fn,"w");
	for (i=0;i<N1;i++) {
		for(j=0;j<N2;j++) {
			fprintf(fp,"%.6g ",mx[i*N2+j]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}

void write_HS(char *fn, complex double **HS)
{
	FILE *fp;
	int ir,lm;
	
	fp = fopen(fn,"w");
	for (ir=0;ir<NR;ir++) {
		if (HS[ir] != NULL) {
			for(lm=0;lm<NLM;lm++) {
				fprintf(fp,"%.6g %.6g  ",creal(HS[ir][lm]),cimag(HS[ir][lm]));
			}
			fprintf(fp,"\n");
		}
	}
	fclose(fp);
}

void write_slice(char *fn, double **v, int im)
{
	FILE *fp;
	int i,j;

	fp = fopen(fn,"w");
		fprintf(fp,"0 ");			// first row = radius
		for(j=0;j<NLAT/2;j++) {
			fprintf(fp,"%.6g ",ct[j]);	// first line = cos(theta)
		}
		for(j=1;j<=NLAT/2;j++) {
			fprintf(fp,"-%.6g ",ct[NLAT/2-j]);	// first line = cos(theta)
		}
	for (i=0;i<NR;i++) {
		if (v[i] != NULL) {
			fprintf(fp,"\n%.6g ",r[i]);		// first row = radius
			for(j=0;j<NLAT;j++) {
				fprintf(fp,"%.6g ",v[i][im*NLAT + j]);		// data
			}
		}
	}
	fclose(fp);
}


void rot_PolTor(struct PolTor *PT, long int istart, long int iend)
{
	complex double tmp[NLM*2];
	complex double *lapl, *lapd, *lapt;
	complex double **p;
	long int ir,lm;

	lapl = tmp;	lapd = tmp + NLM;
	ir = istart;
		for (lm=0; lm<NLM; lm++) {		// Solenoidal deduced from radial derivative of Poloidal
			lapl[lm] = (Lr[ir].d -l2[lm]*r_2[ir])*PT->P[ir][lm] + Lr[ir].u*PT->P[ir+1][lm];
		}
	for (ir=istart+1; ir < iend; ir++) {
		for (lm=0; lm<NLM; lm++) {		// Solenoidal deduced from radial derivative of Poloidal
			lapd[lm] = Lr[ir].l*PT->P[ir-1][lm] + (Lr[ir].d -l2[lm]*r_2[ir])*PT->P[ir][lm] + Lr[ir].u*PT->P[ir+1][lm];
			PT->P[ir-1][lm] = -lapl[lm];
		}
		lapt = lapl;	lapl = lapd;	lapd = lapt;	// rotate buffers.
	}
	ir = iend;
		for (lm=0; lm<NLM; lm++) {		// Solenoidal deduced from radial derivative of Poloidal
			lapd[lm] = Lr[ir].l*PT->P[ir-1][lm] + (Lr[ir].d -l2[lm]*r_2[ir])*PT->P[ir][lm];
			PT->P[ir-1][lm] = -lapl[lm];
			PT->P[ir][lm] = -lapd[lm];
		}
	// switch Pol <> Tor
	p = PT->P;	PT->P = PT->T;	PT->T = p;
}

void usage()
{
	printf("\nUsage: xspp <poltor-file-saved-by-xshells> [op] command [args [...]]\n");
	printf("list of available optional ops :\n");
	printf(" curl : compute curl of field\n");
	printf("list of available commands :\n");
	printf(" axi  : write meridional slice of axisymetric component (m=0)\n");
	printf(" slice [angle]  : write meridional slice at phi=angle in degrees (default=0°) of vector field in spherical coordinates\n");
	printf(" HS  : write spherical harmonics decomposition\n");
}

int main (int argc, char *argv[])
{
	double t0, tmp;
	long int BC, irs,ire;
	long int ic;
	long int i,im,m,l,lm;

	printf(" [XSPP] Xshells Post-Processing   by N. Schaeffer / LGIT, build %s, %s\n",__DATE__,__TIME__);
	printf("  => compiled with: nlat=%d, nphi=%d,  lmax=%d, mmax=%d (mres=%d)\n",NLAT,NPHI,LMAX,MMAX,MRES);
	if (argc <= 2) { usage(); exit(1); }

// init
	fftw_plan_mode = FFTW_ESTIMATE;		// fast FFTW init.
	init_SH(0.);

//load
	Blm.P = NULL;		// require allocation by load_PolTor
	load_PolTor(argv[1], &Blm, &t0, &irs, &ire, &BC);
	alloc_VectField(&B, irs, ire);

// parse optional op...
	ic = 2;		// current argument count.
	if (strcmp(argv[ic],"curl") == 0)
	{
		rot_PolTor(&Blm, irs, ire);	// compute rotational of field.
		printf("> curl computed\n");
		ic++;	// go to next command line argument.
	}

	if (argc <= ic) runerr("missing command.");

// write radial grid
	write_vect("o_r",r,NR);
	printf("> radial grid points written to file : o_r\n");
	write_vect("o_ct",ct,NLAT);
	printf("> angular grid cos(theta) written to file : o_ct\n");

// parse commands ...
	if (strcmp(argv[ic],"axi") == 0)
	{
		for (i=irs; i<=ire; i++) {
			for (lm=LiM(MRES,1);lm<NLM;lm++) {	// zero out all non-axisymmetric modes.
				Blm.P[i][lm] = 0.0;	Blm.T[i][lm] = 0.0;
			}
		}
		PolTor_to_spat(&Blm, &B, irs, ire, BC);
		write_slice("o_Vp", B.p, 0);		// write phi component
		for (i=irs; i<=ire; i++) {
			SHtor_to_spat(Blm.P[i], B.t[i], B.p[i]);
			for (l=0;l<NLAT;l++) B.p[i][l] *= -r[i]*st[l];	// stream function
		}
		write_slice("o_Vpol", B.p, 0);
		printf("> axisymmetric component written to files : o_Vp (phi component) and o_Vpol (poloidal potential)\n");
		exit(0);
	}
	if (strcmp(argv[ic],"slice") == 0)
	{
		double phi = 0.0;
		ic++;
		if (argc > ic) sscanf(argv[ic],"%lf",&phi);
		PolTor_to_rot_spat(&Blm, 0.0, &B, irs+1, ire-1, BC);
		i = lround(phi/360. *NPHI) % NPHI;
		write_slice("o_Vr",B.r,i);	write_slice("o_Vt",B.t,i);	write_slice("o_Vp",B.p,i);
		printf("> meridional slice #%d (phi=%.1f°) written to files : o_Vr, o_Vt, o_Vp (spherical vector components)\n",i,i*360./NPHI);
		exit(0);
	}
	if (strcmp(argv[ic],"HS") == 0)
	{
		write_HS("o_Plm",Blm.P);	write_HS("o_Tlm",Blm.T);
		printf("> spherical harmonics decomposition written to files : o_Plm, o_Tlm (poloidal/toroidal components)\n");
		exit(0);
	}


	printf("!!! warning: command \"%s\" was not understood !!!\n",argv[ic]);
	exit(1);

/*
	fp = fopen("Pprof","w");
	for (i=0; i<NR; i++) {
		if (Blm.P[i] != NULL) fprintf(fp,"%.6g ",Blm.P[i][LM(2,0)]);
	}
	fclose(fp);
	fp = fopen("Tprof","w");
	for (i=0; i<NR; i++) {
		if (Blm.T[i] != NULL) fprintf(fp,"%.6g ",Blm.T[i][LM(1,0)]);
	}
	fclose(fp);
*/

}

