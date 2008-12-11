/**
Compile-time Parmeter file for XSHELLS 
spherical harmonics : LMAX, NLAT, MMAX, MRES, NPHI
execution control : MASK, COIL, XS_LINEAR
initial fields : init_B0 and init_U0 functions.
*/
#ifndef LM	/* SHT size definition part (DO NOT REMOVE) */

///  SIZES  ///
// LMAX : maximum degree of Spherical Harmonic (even or odd)
#define LMAX 200
// NLAT : number of latitudinal (theta) gauss points, at least (LMAX+1)
// (and (LMAX+1)*2 for dealias) **can be EVEN or ODD**
#define NLAT 300
//#define NLAT 144

// MMAX : max number of fourrier decomposition (order = MMAX * MRES)
// ** MMAX*MRES <= LMAX **
#define MMAX 0
// MRES : azimutal symmetry
#define MRES 1
// NPHI : number of azimutal grid points, at least 2*MMAX+1, 3*MMAX for antialiasing and 4*MMAX for full dealias
// power of two is better (but not mandatory), 
#define NPHI 6

///  ADDITIONAL FINE TUNING  ///
// compute and print some debugging information...
//#define _SH_DEBUG_

#else	/* flow execution and imposed fields part (DO NOT REMOVE) */

///  EXECUTION CONTROL  ///
// use LINEAR computation : no u.grad(u), no J0xB0
//#define XS_LINEAR
// use non-linear computation, but removes J0xB0
#define NO_J0xB0
// use mask feature : allows computation in arbitrary fluid domains.
//#define MASK
// use coil : forcing by current.
//#define COIL

/// INITIAL OR IMPOSED FIELDS
// function init_B0 is called to set initial magnetic field value if b0 != 0. (see B0 in xshells.par)
// r is the radius, b0 and b1 are two parameters, Pol[NLM] and Tor[NLM] are poloidal/toroidal arrays at given radius.
// use the macros Set_Poloidal( val ) and Set_Toroidal( val ) to set Pol and Tor array for given r, l, m.
void init_B0(double r, double b0, double b1, complex double *Pol, complex double *Tor)
{
	long int l,m;
	double v;
	complex double z;
/*
#define B0_NAME "uniform vertical field (z-axis)"
	{
	l=1; m=0;
		Set_Poloidal( b0* r/2. * Y10_ct )	// Y10_ct  makes it unitary in physical space.
	}

#define B0_NAME "current-free dipole"
	{
	l=1; m=0;
		if (r < 0.1) r=0.1:		// avoid divergence at r=0
		z = 1.0/(2.*r*r);
		Set_Poloidal( b0* z ) 		// magnetic dipole like in E. Dormy's thesis (current free)
	}

#define B0_NAME "Jault 2008 (not current-free) normalized version"
	if (r != 0.0) {
	l=1; m=0;
		v = pi * r;
		z = ((sin(v)/v-cos(v))/v - 0.3*(sin(2*v)/(2*v)-cos(2*v))/(2*v)) *sqrt(4.*pi/(2*l+1));	// j1(pi*r) - 0.3*j1(2*pi*r)
		Set_Poloidal( b0* z )
	l=3; m=0;
		v = 5.7635 * r;
		z = -0.2*((sin(v)*(3./(v*v)-1.)-3.*cos(v)/v)/v) *sqrt(4.*pi/(2*l+1));	// -0.2*j3(k*r)
		Set_Poloidal( b0* z )
	}
*/
#define B0_NAME "3D based on Jault 2008 (dipole symmetry, not current-free, normalized version, no toroidal)"
	if (r != 0.0) {
	l=1; m=0;
		v = pi * r;
		z = ((sin(v)/v-cos(v))/v - 0.3*(sin(2*v)/(2*v)-cos(2*v))/(2*v));	// j1(pi*r) - 0.3*j1(2*pi*r)
		Set_Poloidal( b0* z *sqrt(4.*pi/(2*l+1)) )
//              Set_Toroidal( 10.0 * b0*(r-1.)*(r-0.35)* z *sqrt(4.*pi/(2*l+1)) )
	l=MRES*2+1; m=MRES;
		Set_Poloidal( 10.*r*r/(l*l) * b0* b1*z *sqrt(4.*pi/(2*l+1)) )
	l=3; m=0;
		v = 5.7635 * r;
		z = -0.2*((sin(v)*(3./(v*v)-1.)-3.*cos(v)/v)/v);	// -0.2*j3(k*r)
		Set_Poloidal( b0* z *sqrt(4.*pi/(2*l+1)) )
	l=MRES*2+4; m=MRES;
		Set_Poloidal( 10.*r*r/(l*l) * b0* b1*z *sqrt(4.*pi/(2*l+1)) )
	}
/*
#define B0_NAME "3D quadrupole symmetry (not current-free, normalized version)"
        if (r != 0.0) {
        l=2; m=0;
                v = pi * r;
                z = ((sin(v)/v-cos(v))/v - 0.3*(sin(2*v)/(2*v)-cos(2*v))/(2*v));        // j1(pi*r) - 0.3*j1(2*pi*r)
                Set_Poloidal( b0* z *sqrt(4.*pi/(2*l+1)) )
        l=MRES+2; m=MRES;
                Set_Poloidal( 10.*r*r/(l*l) * b0* b1*z *sqrt(4.*pi/(2*l+1)) )
        l=4; m=0;
                v = 5.7635 * r;
                z = -0.2*((sin(v)*(3./(v*v)-1.)-3.*cos(v)/v)/v);        // -0.2*j3(k*r)
                Set_Poloidal( b0* z *sqrt(4.*pi/(2*l+1)) )
        l=MRES+5; m=MRES;
                Set_Poloidal( 10.*r*r/(l*l) * b0* b1*z *sqrt(4.*pi/(2*l+1)) )
        }
*/
}

/*
void init_U0(double r, double u0, double u1, complex double *Pol, complex double *Tor)
{
	long int l,m;
	double v;
	complex double z;

#define U0_NAME "m=0 l=2 simple roll flow (eq 24 from Dudely & James 1989)"
	{
	l=2; m=0;
		z = r*sin(pi*r);
		Set_Toroidal( u0* z )
		Set_Poloidal( u0* z * 0.14 )
	}

#define U0_NAME "m=0 l=2 gubins flow (from Dudely & James 1989)"
	{
	l=2; m=0;
		z = -r*sin(2*pi*r) * tanh(2*pi*(1-r));
		Set_Toroidal( u0* z )
		Set_Poloidal( u0* z * 0.1 )
	}

#define U0_NAME "m=2 l=2 j2 pekeris flow (eq 20-21 Dudely & James 1989)"
	if (r != 0.0) {
	l=2; m=0;
		v =  5.7634591968447*r;
		z = 5.7634591968447 * ((3.0/(v*v*v) -1.0/v)*sin(v) - 3./(v*v) * cos(v));
		Set_Poloidal( u0 * z )
		Set_Toroidal( u0 * 5.7634591968447 * z )
	}
}
*/

#endif	/* DO NOT REMOVE */
