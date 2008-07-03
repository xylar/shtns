/*
Parmeter file for Spherical Harmonics Transform.
LMAX, NLAT, MMAX, MRES, NPHI
*/

//  SIZES  //
// LMAX : maximum degree of Spherical Harmonic
#define LMAX 10
// NLAT : number of latitudinal (theta) gauss points, at least (LMAX+1), must be EVEN (and (LMAX+1)*2 for dealias)
#define NLAT 12
//#define NLAT 144

// MMAX : max number of fourrier decomposition (degree = MMAX * MRES)
// power of two is better (but not mandatory), MMAX*MRES <= LMAX
#define MMAX 8
// MRES : azimutal symmetry
#define MRES 1
// NPHI : number of azimutal grid points, at least 2*MMAX, 3*MMAX for antialiasing and 4*MMAX for full dealias
#define NPHI 16

//  ADDITIONAL FINE TUNING  //
// POLAR_OPT_THRESHOLD : value under wich the polar values of the Legendre Polynomials Plm are neglected,
// leading to increased performance (a few percent)
// 0 = no polar optimization;  1.e-14 = VERY safe;  1.e-10 = safe;  1.e-6 = aggresive.
//#define POLAR_OPT_THRESHOLD 1.e-6
#define POLAR_OPT_THRESHOLD 0
// compute and print some debugging information...
#define _SH_DEBUG_
