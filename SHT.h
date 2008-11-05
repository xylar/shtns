/**
Parmeter file for Spherical Harmonics Transform.
LMAX, NLAT, MMAX, MRES, NPHI
*/

//  SIZES  //
// LMAX : maximum degree of Spherical Harmonic
#define LMAX 6
// NLAT : number of latitudinal (theta) gauss points, at least (LMAX+1)
// (and (LMAX+1)*2 for dealias) **can be EVEN or ODD**
#define NLAT 16
//#define NLAT 144

// MMAX : max number of fourrier decomposition (degree = MMAX * MRES)
// ** MMAX*MRES <= LMAX **
#define MMAX 5
// MRES : azimutal symmetry
#define MRES 1
// NPHI : number of azimutal grid points, at least 2*MMAX+1, 3*MMAX for antialiasing and 4*MMAX for full dealias
// power of two is better (but not mandatory), 
#define NPHI 32

//  ADDITIONAL FINE TUNING  //
// compute and print some debugging information...
//#define _SH_DEBUG_
