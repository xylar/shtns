  ////////////////////////////////////
 // XSPP : XShells post-processing //
////////////////////////////////////

#define DEB printf("%s:%u pass\n", __FILE__, __LINE__)
//#define DEB (0)

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
#define NLAT 301
#define LMAX 30
#define NPHI 64
#define MMAX 30
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

#include "xshells_fields.c"
#include "xshells_io.c"

struct VectField B;
struct PolTor Blm;

long int BC;
long int irs, ire;	// loaded from file.
long int ips=0, ipe=NPHI-1;
long int its=0, ite=NLAT-1;

long int mmin=0, mmax=MMAX;
long int lmin=0, lmax=LMAX;

double lspec[LMAX+1];	// spectra
double mspec[MMAX+1];

struct JobInfo jpar;	// parameters from loaded file

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
	fprintf(fp,"%% [XSHELLS] vector\n");
	for (i=0;i<N;i++) {
		fprintf(fp,"%.6g ",vec[i]);
	}
	fprintf(fp,"\n");	fclose(fp);
}

void write_HS(char *fn, complex double **HS)
{
	long int ir,lm,l,m;
	FILE *fp;
	
	fp = fopen(fn,"w");
	fprintf(fp,"%% [XSHELLS] Spherical Harmonics coefficients : MMAX=%d, LMAX=%d, MRES=%d. l-contiguous storing", mmax, lmax, jpar.mres);
	for (ir=irs;ir<=ire;ir++) {
		if (HS[ir] != NULL) {
			fprintf(fp,"\n%%  ir=%d, r=%f\n",ir,r[ir]);
			for (m=0; m<=mmax; m+=jpar.mres) {
				for (l=m; l<=lmax; l++) {
					lm = LM(l,m);
					fprintf(fp,"%.6g %.6g  ",creal(HS[ir][lm]),cimag(HS[ir][lm]));
				}
			}
		}
	}
	fprintf(fp,"\n");	fclose(fp);
}

void write_Spec_m(char *fn, complex double **HS)
{
	double sum, cr,ci;
	long int ir,lm,l,m;
	FILE *fp;

	fp = fopen(fn,"w");
	fprintf(fp,"%% [XSHELLS] Spherical Harmonics m-spectrum : MMAX=%d, MRES=%d. first row is r", mmax, jpar.mres);
	for (ir=irs;ir<=ire;ir++) {
		if (HS[ir] != NULL) {
			fprintf(fp,"\n%%  ir=%d, r=%f\n%.6g ",ir,r[ir],r[ir]);
			for (m=0; m<=mmax; m+=jpar.mres) {
				sum = 0.0;
				for (l=m; l<=lmax; l++) {
					lm = LM(l,m);
					cr = creal(HS[ir][lm]);	ci = cimag(HS[ir][lm]);
					sum += cr*cr + ci*ci;
				}
				fprintf(fp,"%.6g ",sum);
			}
		}
	}
	fprintf(fp,"\n");	fclose(fp);
}

void write_Spec_l(char *fn, complex double **HS)
{
	double sum, cr,ci;
	long int ir,lm,l,m;
	FILE *fp;

	fp = fopen(fn,"w");
	fprintf(fp,"%% [XSHELLS] Spherical Harmonics l-spectrum : LMAX=%d. first row is r", lmax);
	for (ir=irs;ir<=ire;ir++) {
		if (HS[ir] != NULL) {
			fprintf(fp,"\n%%  ir=%d, r=%f\n%.6g ",ir,r[ir],r[ir]);
			for (l=0; l<=lmax; l++) {
				sum = 0.0;
				for (m=0; (m<=mmax)&&(m<=l); m+=jpar.mres) {
					lm = LM(l,m);
					cr = creal(HS[ir][lm]);	ci = cimag(HS[ir][lm]);
					sum += cr*cr + ci*ci;
				}
				fprintf(fp,"%.6g ",sum);
			}
		}
	}
	fprintf(fp,"\n");	fclose(fp);
}

void write_shell(char *fn, double *mx)
{
	FILE *fp;
	long int i,j;	//phi, theta
	
	fp = fopen(fn,"w");
	fprintf(fp,"%% [XSHELLS] Surface data (sphere). first line is theta, first row is phi.\n0 ");
		for(j=its;j<=ite;j++) {
			fprintf(fp,"%.6g ",acos(ct[j]));	// first line = theta (radians)
		} 
	for (i=ips; i<=ipe; i++) {
		fprintf(fp,"\n%.6g ",phi_rad(i));		// first row = phi (radians)
		for(j=its; j<=ite; j++) {
			fprintf(fp,"%.6g ",mx[i*NLAT+j]);
		}
	}
	fprintf(fp,"\n");	fclose(fp);
}

void write_merid(char *fn, double **v, long int im)
{
	FILE *fp;
	long int i,j;

	fp = fopen(fn,"w");
	fprintf(fp,"%% [XSHELLS] Meridian slice. first line is cos(theta), first row is r.\n0 ");
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
	fprintf(fp,"\n");	fclose(fp);
}

// equatorial cut : fac is used to change the sign of theta component to get z-component.
void write_equat(char *fn, double **v, double fac)
{
	FILE *fp;
	long int i,j;
	long int it = NLAT/2;	// equator is at NLAT/2.

	fp = fopen(fn,"w");
	fprintf(fp,"%% [XSHELLS] Equatorial cut. first line is phi (in radian), first row is r.\n0 ");
		for(j=ips; j<=ipe; j++) {
			fprintf(fp,"%.6g ",phi_rad(j));	// first line = phi angle (radian)
		}
	for (i=irs;i<=ire;i++) {
		if (v[i] != NULL) {
			fprintf(fp,"\n%.6g ",r[i]);		// first row = radius
			for(j=ips; j<=ipe; j++) {
				fprintf(fp,"%.6g ",v[i][j*NLAT + it] * fac);		// data
			}
		}
	}
	fprintf(fp,"\n");	fclose(fp);
}

// vector vr, vt, vp is converted to cartesian coordinates.
void spher_to_cart(double ct,double p,double *vr,double *vt,double *vp)
{
	double st, cp,sp;
	double vx,vy,vz,vs;
	
	st = sqrt(1.0 -ct*ct);	cp = cos(p);	sp = sin(p);
	vz = *vr*ct - *vt*st;		vs = *vr*st + *vt*ct;
	vx = vs*cp - *vp*sp;		vy = vs*sp + *vp*cp;
	*vr = vx;	*vt = vy;	*vp = vz;
}

void write_line(char *fn,double x,double y,double z,double vx,double vy,double vz,int ni, struct PolTor *Blm)
{
	double rr,cost,phi;
	double br,bt,bp;
	int i;
	FILE *fp;
	
	fp = fopen(fn,"w");
	fprintf(fp,"%% [XSHELLS] line profile starting at %f,%f,%f with increment %f,%f,%f\n%% x y z\tr cos(theta) phi\tvr vt vp\n",x,y,z, vx,vy,vz);
	for (i=0; i<ni; i++) {
		rr = sqrt(x*x + y*y + z*z);		cost = z/rr;
		phi = atan(y/x);	if (x < 0.0) phi += pi;
		PolTor_to_point_interp(Blm, rr, cost, phi, &br, &bt, &bp);
		fprintf(fp,"%.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g\n",x,y,z, rr, cost, phi, br,bt,bp);
		x+=vx;	y+=vy;	z+=vz;
	}
	fclose(fp);
	printf("> line profile starting at at %.3f,%.3f,%.3f with increment %.4g,%.4g,%.4g written to %s\n",x,y,z, vx,vy,vz,fn);
}

void write_disc(char *fn,double x0,double y0,double z0, int nphi, struct PolTor *Blm)
{
	double r2, ct0, st0, phi0, s0;	// absolute coordinates of center.
	double s, ca, sa;			// disc coordinates
	double rr, cost,phi;		// spherical coordinates
	double x,y,z;				// cartesian coordinates
	double bx,by,bz;			// vector
	int ir, ip;
	FILE *fp;

	r2 = x0*x0 + y0*y0 + z0*z0;		rr = sqrt(r2);	// center position.
	ct0 = z0/rr;	st0 = sqrt(1. - ct0*ct0);
	phi0 = atan(y0/x0);	if (x0 < 0.) phi0 += pi;
	s0 = sqrt(x0*x0 + y0*y0);		// rotation of phi0.
	if (s0 == 0.) phi0 = 0.;

	fp = fopen(fn,"w");
	if (rr < r[NR-1]) {		// vector is normal and center of disc
		fprintf(fp,"%% [XSHELLS] disc slice centered at %f,%f,%f (first row is disc-radius)",x0,y0,z0);
	} else {		// vector is normal and center is 0,0,0
		fprintf(fp,"%% [XSHELLS] disc slice centered at 0 with normal %f,%f,%f (first row is disc-radius)",x0,y0,z0);
		x0 = 0.;	y0 = 0.;	z0 = 0.;	rr = 0.;	r2 = 0.;
	}
	ir = r_to_idx(rr);	if (r[ir] <= rr) ir++;	//r[ir] > rr;
	while((Blm->P[ir] == NULL)&&(ir<NR)) ir++;	// skip undef data.
	printf("> writing disc slice centered at x=%.3f, y=%.3f, z=%.3f to %s\n",x0,y0,z0,fn);
	printf("    normal : phi=%.1fdeg, cos(theta)=%.4f, sin(theta)=%.4f\n",phi0*180./pi,ct0,st0);
	printf("    start at ir=%d r=%.4f\n",ir, r[ir]);

	while ((ir<NR-1)&&(Blm->P[ir] != NULL)) {
		rr = r[ir];
		printf("ir=%d  r=%f\r",ir,rr);	fflush(stdout);
		s = sqrt(rr*rr - r2);		// disc radius.
		fprintf(fp,"\n%.6g ",s);	// first row = disc radius.
		for (ip=0; ip<nphi; ip++) {
			phi = ip*(2.*pi/nphi);		// disc angle
			ca = cos(phi);	sa = sin(phi);
			z = z0 + s*sa*st0;
			x = s0 - s*sa*ct0;
			y = s*ca;
			phi = atan(y/x);	if (x < 0.0) phi += pi;
			cost = z/rr;
			PolTor_to_point(Blm, ir, cost, phi+phi0, &bx, &by, &bz);
			spher_to_cart(cost, phi, &bx, &by, &bz);
			y = bz*st0 - bx*ct0;		// project into disc-relative coordinates.
			z = bz*ct0 + bx*st0;
			x = by;
//			fprintf(fp,"%.6g %.6g %.6g ",x*ca+y*sa, y*ca-x*sa, z);	// write data
			fprintf(fp,"%.6g %.6g %.6g ",x, y, z);	// write data
		}
		ir++;
	}
	fprintf(fp,"\n");	fclose(fp);
	printf("    end at ir=%d r=%.4f.\n",ir-1, rr);
}

// apply (l,m) filter
void filter_lm(struct PolTor *Blm, int lmin, int lmax, int mmin, int mmax)
{
	long int ir, m, l, lm;

	if (lmax > LMAX) lmax = LMAX;	if (mmax > MMAX) mmax = MMAX;

	for (ir=jpar.irs; ir<=jpar.ire; ir++) {		// filter all data, not only restricted data...
		for (m=0; m<mmin; m++) {
			for(l=m; l<=lmax; l++) {
				lm = LM(l,m);	Blm->P[ir][lm] = 0.0;	Blm->T[ir][lm] = 0.0;
			}
		}
		for (m=mmin; m<=mmax; m++) {
			for(l=m; l<lmin; l++) {
				lm = LM(l,m);	Blm->P[ir][lm] = 0.0;	Blm->T[ir][lm] = 0.0;
			}
			for(l=lmax+1; l<=LMAX; l++) {
				lm = LM(l,m);	Blm->P[ir][lm] = 0.0;	Blm->T[ir][lm] = 0.0;
			}
		}
		for (m=mmax+1; m<=MMAX; m++) {
			for(l=m; l<=lmax; l++) {
				lm = LM(l,m);	Blm->P[ir][lm] = 0.0;	Blm->T[ir][lm] = 0.0;
			}
		}
	}
}

void calc_spec(complex double **HS, double *spl, double *spm)
{
	double cr,ci;
	long int ir,lm,l,m;

	for (l=0; l<=LMAX; l++)	spl[l] = 0.0;
	for (m=0; m<=MMAX; m++)	spm[m] = 0.0;

	for (ir=irs;ir<=ire;ir++) {
		if (HS[ir] != NULL) {
			for (m=0; m<=mmax; m+=jpar.mres) {
				for (l=m; l<=lmax; l++) {
					lm = LM(l,m);
					cr = creal(HS[ir][lm]);	ci = cimag(HS[ir][lm]);
					spl[l] += cr*cr + ci*ci;
					spm[m] += cr*cr + ci*ci;
				}
			}
		}
	}
}

void usage()
{
	printf("\nUsage: xspp <poltor-file-saved-by-xshells> [op1] [op2] [...] command [args [...]]\n");
	printf("list of available optional ops :\n");
	printf(" curl : compute curl of field\n");
	printf(" rlim <rmin>:<rmax> : render only from rmin to rmax\n");
	printf(" philim <min>:<max> : render only from min to max azimutal degrees\n");
	printf(" llim <lmin>:<lmax> : use only spherical harmonic degree from lmin to lmax.\n");
	printf(" mlim <mmin>:<mmax> : use only spherical harmonic orders from mmin to mmax.\n");
	printf("list of available commands :\n");
	printf(" axi  : write meridional slice of axisymetric component (m=0)\n");
	printf(" equat  : write equatorial cut of vector field in cylindrical coordinates (r, phi)\n");
	printf(" merid [angle]  : write meridional slice at phi=angle in degrees (default=0°) of vector field in spherical coordinates\n");
//	printf(" zcut [z]  : write slice at height=z (-1 to 1) of vector field in cylindrical coordinates (s, phi, z)\n");
	printf(" surf [r]  : write surface data at r (0 to 1) or ir (1 to NR-1) of vector field in spherical coordinates\n");
	printf(" HS  : write full spherical harmonics decomposition\n");
	printf(" spec  : write spherical harmonic (l and m)-spectra for each shell\n");
	printf(" line ni x0,y0,z0 vx,vy,vz  : write ni points along line profile starting at (x0,y0,z0) along increment vector (vx,vy,vz)\n");
	printf(" disc nphi x0,y0,z0 : write disc slice centered at (x0,y0,z0) with nphi azimutal points.\n      if center is outside data domain, then it is taken as a normal vector, and the center is 0\n");
	printf(" 3D  : write three-dimensional vector in cartesian coordinates on spherical grid.\n");
}

int main (int argc, char *argv[])
{
	double tmp;
	long int ic;
	long int i,im,m,l,lm;
	int filter_req = 0;

    int parse_op(int ic)	// returns number of parsed command line argument
    {
	double min, max;
	long int i;

	if (argc <= ic) return 0;
	if (strcmp(argv[ic],"curl") == 0)
	{
		curl_PolTor(&Blm, irs, ire);	// compute rotational of field.
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
	if (strcmp(argv[ic],"mlim") == 0)
	{
		ic++;	// go to next command line argument.
		if (argc > ic) {
			sscanf(argv[ic],"%lf:%lf",&min,&max);
		} else runerr("mlim <min>:<max>  => missing arguments\n");
		mmin = min;	mmax = max;
		if (mmax > jpar.mmax) mmax = jpar.mmax;	if (mmin < 0) mmin = 0;
		printf("> using only m = #%d to #%d\n",mmin, mmax);
		filter_req = 1;
		return 2;	// 2 command line argument parsed.
	}
	if (strcmp(argv[ic],"llim") == 0)
	{
		ic++;	// go to next command line argument.
		if (argc > ic) {
			sscanf(argv[ic],"%lf:%lf",&min,&max);
		} else runerr("llim <min>:<max>  => missing arguments\n");
		lmin = min;	lmax = max;
		if (lmax > jpar.lmax) lmax = jpar.lmax;	if (lmin < 0) lmin = 0;
		printf("> using only l = #%d to #%d\n",lmin, lmax);
		filter_req = 1;
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
	load_PolTor(argv[1], &Blm, &jpar);
		irs = jpar.irs;	ire = jpar.ire;	BC = jpar.BC;
		if ((lmax < jpar.lmax)||(mmax < jpar.mmax))
			printf("  ! warning : data has higher resolution than compile time set up. => truncated.\n");
	alloc_VectField(&B, irs, ire);

// parse optional op...
	ic = 2;		// current argument count.
	do {
		i = parse_op(ic);
		ic += i;
	} while(i!=0);

	if (argc <= ic) runerr("missing command.");

// apply (l,m) restrictions
	if (lmax > LMAX) lmax = LMAX;	if (mmax > MMAX) mmax = MMAX;
	if (lmax > jpar.lmax) lmax=jpar.lmax;	if (mmax > jpar.mmax) mmax=jpar.mmax;
	if (filter_req) filter_lm(&Blm, lmin, lmax, mmin, mmax);

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
		write_merid("o_Vp", B.p, 0);		// write phi component
		write_merid("o_Vpol", B.r, 0);		// write stream function
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
		write_merid("o_Vr",B.r,im);	write_merid("o_Vt",B.t,im);	write_merid("o_Vp",B.p,im);
		printf("> meridional slice #%d (phi=%.1f°) written to files : o_Vr, o_Vt, o_Vp (spherical vector components)\n",im,phi_deg(im));
		for (i=irs;i<=ire;i++) {
			for(l=0;l<NLAT;l++) {
				B.p[i][im*NLAT +l] = B.r[i][im*NLAT +l]*ct[l] - B.t[i][im*NLAT +l]*st[l];	// z
				B.r[i][im*NLAT +l] = B.r[i][im*NLAT +l]*st[l] + B.t[i][im*NLAT +l]*ct[l];	// s
			}
		}
		write_merid("o_Vs",B.r,im);	write_merid("o_Vz",B.p,im);
		printf("  and files o_Vs, o_Vp, o_Vz (cylindrical vector components)\n",i,phi_deg(i));
		exit(0);
	}
	if (strcmp(argv[ic],"HS") == 0)
	{
		write_HS("o_Plm",Blm.P);	write_HS("o_Tlm",Blm.T);
		printf("> spherical harmonics decomposition written to files : o_Plm, o_Tlm (poloidal/toroidal components)\n");
		exit(0);
	}
	if (strcmp(argv[ic],"spec") == 0)
	{
		write_Spec_l("o_Plr",Blm.P);	write_Spec_l("o_Tlr",Blm.T);
		write_Spec_m("o_Pmr",Blm.P);	write_Spec_m("o_Tmr",Blm.T);
		printf("> spherical harmonics spectrum written to files : o_Plr, o_Tlr, o_Pmr, o_Tlr (poloidal/toroidal, l/m, at each r)\n");
		calc_spec(Blm.P, lspec, mspec);		write_vect("o_Pl",lspec, lmax+1);	write_vect("o_Pm",mspec, mmax+1);
		calc_spec(Blm.T, lspec, mspec);		write_vect("o_Tl",lspec, lmax+1);	write_vect("o_Tm",mspec, mmax+1);
		exit(0);
	}
	if (strcmp(argv[ic],"line") == 0)
	{
		double x,y,z, vx,vy,vz;
		long int ni;
		ic++;
		if (ic+2 >= argc) runerr("line definition is missing...\n");
		sscanf(argv[ic],"%lf",&x);	ni = x;
		sscanf(argv[ic+1],"%lf,%lf,%lf",&x, &y, &z);
		sscanf(argv[ic+2],"%lf,%lf,%lf",&vx, &vy, &vz);
		write_line("o_line",x,y,z,vx,vy,vz,ni,&Blm);
		exit(0);
	}
	if (strcmp(argv[ic],"disc") == 0)
	{
		double x,y,z;
		long int nphi;
		ic++;
		if (ic+1 >= argc) runerr("disc definition is missing...\n");
		sscanf(argv[ic],"%lf",&x);	nphi = x;
		sscanf(argv[ic+1],"%lf,%lf,%lf",&x, &y, &z);
		write_disc("o_disc",x,y,z,nphi,&Blm);
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

