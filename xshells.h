/**
Parmeter file for XSHELLS 
spherical harmonics : LMAX, NLAT, MMAX, MRES, NPHI
execution control : MASK, COIL, XS_LINEAR
*/

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

///  EXECUTION CONTROL  ///
// use LINEAR computation : no u.grad(u), no J0xB0
//#define XS_LINEAR
// use non-linear computation, but removes J0xB0
#define NO_J0xB0
// use mask feature : allows computation in arbitrary fluid domains.
//#define MASK
// use coil : forcing by current.
//#define COIL

