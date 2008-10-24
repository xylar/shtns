/**
Parmeter file for Spherical Harmonics Transform.
LMAX, NLAT, MMAX, MRES, NPHI
*/

//  SIZES  //
// LMAX : maximum degree of Spherical Harmonic
#define LMAX 22
// NLAT : number of latitudinal (theta) gauss points, at least (LMAX+1)
// (and (LMAX+1)*2 for dealias) **can be EVEN or ODD**
#define NLAT 151
//#define NLAT 144

// MMAX : max number of fourrier decomposition (degree = MMAX * MRES)
// power of two is better (but not mandatory), MMAX*MRES <= LMAX
#define MMAX 3
// MRES : azimutal symmetry
#define MRES 3
// NPHI : number of azimutal grid points, at least 2*MMAX, 3*MMAX for antialiasing and 4*MMAX for full dealias
#define NPHI 8

//  ADDITIONAL FINE TUNING  //
// compute and print some debugging information...
//#define _SH_DEBUG_
