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
#define NLAT 201
#define LMAX 200
#define NPHI 1
#define MMAX 0
#define MRES 1
// SHT on equal spaced grid + polar points.
#define SHT_EQUAL
#include "SHT.c"

// number of radial grid points.
//long int NR,NU;		//  NR: total radial grid points. NU:for velocity field.
//long int NG=0;		//  NG: grid points for inner core.

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

long int BC;
long int irs, ire;	// loaded from file.
long int ips=0, ipe=NPHI-1;
long int its=0, ite=NLAT-1;


// find closest index to phi angle.
inline long int phi_to_idx(double phi) {
	return lround(phi/360. *NPHI*MRES) % NPHI;
}

// get phi angle from index
inline double phi_deg(long int i) {	return i*360./(NPHI*MRES);	}
inline double phi_rad(long int i) {	return i*2.*pi/(NPHI*MRES);	}


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

void write_HS(char *fn, complex double **HS)
{
	FILE *fp;
	long int ir,lm;
	
	fp = fopen(fn,"w");
	for (ir=irs;ir<=ire;ir++) {
		if (HS[ir] != NULL) {
			for(lm=0;lm<NLM;lm++) {
				fprintf(fp,"%.6g %.6g  ",creal(HS[ir][lm]),cimag(HS[ir][lm]));
			}
			fprintf(fp,"\n");
		}
	}
	fclose(fp);
}

void write_shell(char *fn, double *mx)
{
	FILE *fp;
	long int i,j;	//phi, theta
	
	fp = fopen(fn,"w");
		fprintf(fp,"0 ");			// first row = phi
		for(j=its;j<=ite;j++) {
			fprintf(fp,"%.6g ",acos(ct[j]));	// first line = theta (radians)
		} 
	for (i=ips; i<=ipe; i++) {
		fprintf(fp,"\n%.6g ",phi_rad(i));		// first row = phi (radians)
		for(j=its; j<=ite; j++) {
			fprintf(fp,"%.6g ",mx[i*NLAT+j]);
		}
	}
	fclose(fp);
}

void write_slice(char *fn, double **v, long int im)
{
	FILE *fp;
	long int i,j;

	fp = fopen(fn,"w");
		fprintf(fp,"0 ");			// first row = radius
		for(j=its; j<=ite; j++) {
			fprintf(fp,"%.6g ",ct[j]);	// first line = cos(theta)
		}
	for (i=irs;i<=ire;i++) {
		if (v[i] != NULL) {
			fprintf(fp,"\n%.6g ",r[i]);		// first row = radius
			for(j=its; j<=ite; j++) {
				fprintf(fp,"%.6g ",v[i][im*NLAT + j]);		// data
			}
		}
	}
	fclose(fp);
}

// equatorial cut : fac is used to change the sign of theta component to get z-component.
void write_equat(char *fn, double **v, double fac)
{
	FILE *fp;
	long int i,j;
	long int it = NLAT/2;	// equator is at NLAT/2.

	fp = fopen(fn,"w");
		fprintf(fp,"0 ");			// first row = radius
		for(j=ips; j<=ipe; j++) {
			fprintf(fp,"%.6g ",phi_rad(j));	// first line = phi angle (radians)
		}
	for (i=irs;i<=ire;i++) {
		if (v[i] != NULL) {
			fprintf(fp,"\n%.6g ",r[i]);		// first row = radius
			for(j=ips; j<=ipe; j++) {
				fprintf(fp,"%.6g ",v[i][j*NLAT + it] * fac);		// data
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
		for (lm=0; lm<NLM; lm++) {
			lapl[lm] = (Lr[ir].d -l2[lm]*r_2[ir])*PT->P[ir][lm] + Lr[ir].u*PT->P[ir+1][lm];
		}
	for (ir=istart+1; ir < iend; ir++) {
		for (lm=0; lm<NLM; lm++) {
			lapd[lm] = Lr[ir].l*PT->P[ir-1][lm] + (Lr[ir].d -l2[lm]*r_2[ir])*PT->P[ir][lm] + Lr[ir].u*PT->P[ir+1][lm];
			PT->P[ir-1][lm] = -lapl[lm];
		}
		lapt = lapl;	lapl = lapd;	lapd = lapt;	// rotate buffers.
	}
	ir = iend;
		for (lm=0; lm<NLM; lm++) {
			lapd[lm] = Lr[ir].l*PT->P[ir-1][lm] + (Lr[ir].d -l2[lm]*r_2[ir])*PT->P[ir][lm];
			PT->P[ir-1][lm] = -lapl[lm];
			PT->P[ir][lm] = -lapd[lm];
		}
	// switch Pol <> Tor
	p = PT->P;	PT->P = PT->T;	PT->T = p;
}

void usage()
{
	printf("\nUsage: xspp <poltor-file-saved-by-xshells> [op1] [op2] [...] command [args [...]]\n");
	printf("list of available optional ops :\n");
	printf(" curl : compute curl of field\n");
	printf(" rlim <rmin>:<rmax> : render only from rmin to rmax\n");
	printf(" philim <min>:<max> : render only from min to max azimutal degrees\n");
	printf("list of available commands :\n");
	printf(" axi  : write meridional slice of axisymetric component (m=0)\n");
	printf(" equat  : write equatorial cut of vector field in cylindrical coordinates (r, phi)\n");
	printf(" merid [angle]  : write meridional slice at phi=angle in degrees (default=0°) of vector field in spherical coordinates\n");
//	printf(" zcut [z]  : write slice at height=z (-1 to 1) of vector field in cylindrical coordinates (s, phi, z)\n");
	printf(" surf [r]  : write surface data at r (0 to 1) or ir (1 to NR-1) of vector field in spherical coordinates\n");
	printf(" HS  : write spherical harmonics decomposition\n");
	printf(" 3D  : write three-dimensional vector in cartesian coordinates on spherical grid.\n");
}

int main (int argc, char *argv[])
{
	double t0, tmp;
	long int ic;
	long int i,im,m,l,lm;

    int parse_op(int ic)	// returns number of parsed command line argument
    {
	double min, max;
	long int i;

	if (argc <= ic) return 0;
	if (strcmp(argv[ic],"curl") == 0)
	{
		rot_PolTor(&Blm, irs, ire);	// compute rotational of field.
		printf("> curl computed\n");
		return 1;	// 1 command line argument parsed.
	}
	if (strcmp(argv[ic],"rlim") == 0)
	{
		ic++;	// go to next command line argument.
		if (argc > ic) {
			sscanf(argv[ic],"%lf:%lf",&min,&max);
		} else runerr("rlim <rmin>:<rmax>  => missing arguments\n");
		i = r_to_idx(min);	if (i>irs) irs=i;
		i = r_to_idx(max);	if (i<ire) ire=i;
		printf("> restricting to shells #%d to #%d (that is r=[%f, %f])\n",irs,ire, r[irs], r[ire]);
		return 2;	// 2 command line argument parsed.
	}
	if (strcmp(argv[ic],"philim") == 0)
	{
		ic++;	// go to next command line argument.
		if (argc > ic) {
			sscanf(argv[ic],"%lf:%lf",&min,&max);
		} else runerr("philim <min>:<max>  => missing arguments\n");
		ips = phi_to_idx(min);	ipe = phi_to_idx(max);
		printf("> restricting to slice #%d to #%d (that is phi=[%f, %f])\n",ips,ipe, phi_deg(ips), phi_deg(ipe));
		return 2;	// 2 command line argument parsed.
	}
	return 0;
    }

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
	do {
		i = parse_op(ic);
		ic += i;
	} while(i!=0);

	if (argc <= ic) runerr("missing command.");

// write radial grid
	write_vect("o_r",&r[irs],ire-irs+1);
	printf("> radial grid points written to file : o_r\n");
	write_vect("o_cost",&ct[its],ite-its+1);
	printf("> angular grid cos(theta) written to file : o_cost\n");

// parse commands ...
	if (strcmp(argv[ic],"axi") == 0)
	{
		for (i=irs; i<=ire; i++) {
			for (lm=LiM(MRES,1);lm<NLM;lm++) {	// zero out all non-axisymmetric modes.
				Blm.P[i][lm] = 0.0;	Blm.T[i][lm] = 0.0;
			}
		}
		for (i=irs; i<=ire; i++) {
			SHtor_to_spat(Blm.T[i], (complex double *)B.t[i], (complex double *)B.p[i]);
			SHtor_to_spat(Blm.P[i], (complex double *)B.t[i], (complex double *)B.r[i]);
			for (l=0;l<NLAT;l++) B.r[i][l] *= -r[i]*st[l];	// stream function
		}
		write_slice("o_Vp", B.p, 0);		// write phi component
		write_slice("o_Vpol", B.r, 0);		// write stream function
		printf("> axisymmetric component written to files : o_Vp (phi component) and o_Vpol (poloidal stream function)\n");
		exit(0);
	}
	if (strcmp(argv[ic],"equat") == 0)
	{
		if (!(NLAT%2)) runerr("equatorial cut not supported for even NLAT");
		PolTor_to_spat(&Blm, &B, irs, ire, BC);
		write_equat("o_Vr",B.r, 1.);	write_equat("o_Vp",B.p, 1.);	write_equat("o_Vz",B.t, -1.);
		printf("> equatorial slice written to files : o_Vr, o_Vp, o_Vz (cylindrical vector components)\n");
		exit(0);
	}
	if (strcmp(argv[ic],"surf") == 0)
	{
		double rr = 1.0;
		ic++;
		if (argc > ic) sscanf(argv[ic],"%lf",&rr);
		i = r_to_idx(rr);
		if ((i<irs)||(i>ire)) runerr("requested r not available");
		PolTor_to_spat(&Blm, &B, i, i, BC);	// render just one shell.
		write_shell("o_Sr", B.r[i]);	write_shell("o_St", B.t[i]);	write_shell("o_Sp", B.p[i]);
		printf("> surface #%d (r=%.4f) written to files : o_Sr, o_St, o_Sp (spherical vector components)\n",i,r[i]);
		exit(0);
	}
	if (strcmp(argv[ic],"merid") == 0)
	{
		double phi = 0.0;
		ic++;
		if (argc > ic) sscanf(argv[ic],"%lf",&phi);
		PolTor_to_spat(&Blm, &B, irs, ire, BC);
		im = phi_to_idx(phi);
		write_slice("o_Vr",B.r,im);	write_slice("o_Vt",B.t,im);	write_slice("o_Vp",B.p,im);
		printf("> meridional slice #%d (phi=%.1f°) written to files : o_Vr, o_Vt, o_Vp (spherical vector components)\n",i,phi_deg(i));
		for (i=irs;i<=ire;i++) {
			for(l=0;l<NLAT;l++) {
				B.p[i][im*NLAT +l] = B.r[i][im*NLAT +l]*ct[l] - B.t[i][im*NLAT +l]*st[l];	// z
				B.r[i][im*NLAT +l] = B.r[i][im*NLAT +l]*st[l] + B.t[i][im*NLAT +l]*ct[l];	// s
			}
		}
		write_slice("o_Vs",B.r,im);	write_slice("o_Vz",B.p,im);
		printf("  and files o_Vs, o_Vp, o_Vz (cylindrical vector components)\n",i,phi_deg(i));
		exit(0);
	}
	if (strcmp(argv[ic],"HS") == 0)
	{
		write_HS("o_Plm",Blm.P);	write_HS("o_Tlm",Blm.T);
		printf("> spherical harmonics decomposition written to files : o_Plm, o_Tlm (poloidal/toroidal components)\n");
		exit(0);
	}
	if (strcmp(argv[ic],"3D") == 0)
	{
		PolTor_to_spat(&Blm, &B, irs, ire, BC);		//output o_X, o_Y, o_Z, o_Vx, o_Vy, o_Vz for direct use with matlab ?
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

