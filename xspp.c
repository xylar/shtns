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





int main (int argc, char *argv[])
{
	double t0;
	long int BC, irs,ire;
	long int i,im,m,l,lm;

	printf(" [XSPP] Xshells Post-Processing   by N. Schaeffer / LGIT, build %s, %s\n",__DATE__,__TIME__);
	if (argc < 3) {
		printf("\nUsage: xspp <poltor-file-saved-by-xshells> command [args [...]]\n");
		printf("list of available commands :\n");
		printf(" axisym   write meridional slice of axisymetric component (m=0)\n");
		printf(" slice    write meridional slice of vector field in spherical (default) or cylindrical coordinates\n");
		printf(" HS       write spherical harmonics decomposition\n");
		exit(1);
	}

// init
	fftw_plan_mode = FFTW_ESTIMATE;		// fast FFTW init.
	init_SH(0.);

//load
	Blm.P = NULL;		// require allocation by load_PolTor
	load_PolTor(argv[1], &Blm, &t0, &irs, &ire, &BC);
	alloc_VectField(&B, irs, ire);

// write radial grid
	write_vect("o_r",r,NR);
	printf("> radial grid points written to file : o_r\n");
	write_vect("o_ct",ct,NLAT);
	printf("> angular grid cos(theta) written to file : o_ct\n");

// parse command line...
	if (strcmp(argv[2],"axisym") == 0)
	{
		for (i=irs; i<=ire; i++) {
			for (lm=LiM(MRES,1);lm<NLM;lm++) {	// zero out all non-axisymmetric modes.
				Blm.P[i][lm] = 0.0;
				Blm.T[i][lm] = 0.0;
			}
		}
		PolTor_to_spat(&Blm, &B, irs, ire, BC);
		write_slice("o_Vp", B.p, 0);		// write phi component
		for (i=irs; i<=ire; i++) SH_to_spat(Blm.P[i], B.r[i]);		// Pol scalar to spatial domain
		write_slice("o_Vpol", B.r, 0);
		printf("> axisymmetric component written to files : o_Vp (phi component) and o_Vpol (spatial poloidal scalar)\n");
		exit(0);
	}
	if (strcmp(argv[2],"slice") == 0)
	{
		PolTor_to_spat(&Blm, &B, irs, ire, BC);
		write_slice("o_Vr",B.r,0);	write_slice("o_Vt",B.t,0);	write_slice("o_Vp",B.p,0);
		printf("> meridional slices written to files : o_Vr, o_Vt, o_Vp (spherical vector components)\n");
		exit(0);
	}
	if (strcmp(argv[2],"HS") == 0)
	{
		write_HS("o_Plm",Blm.P);	write_HS("o_Tlm",Blm.T);
		printf("> spherical harmonics decomposition written to files : o_Plm, o_Tlm (poloidal/toroidal components)\n");
		exit(0);
	}

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

