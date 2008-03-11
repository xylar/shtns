/*
Parmeter file for Spherical Harmonics Transform.
LMAX, NLAT, MMAX, MRES, NPHI
*/

//  SIZES  //
// LMAX : maximum degree of Spherical Harmonic
#define LMAX 79
// NLAT : number of latitudinal (theta) gauss points, at least (LMAX+1), must be EVEN
#define NLAT (LMAX+1)

// MMAX : max number of fourrier decomposition (degree = MMAX * MRES)
// power of two is better (but not mandatory), MMAX*MRES <= LMAX
#define MMAX 32
// MRES : azimutal symmetry
#define MRES 2
// NPHI : number of azimutal grid points, at least 2*MMAX or 3*MMAX for antialiasing.
#define NPHI (2*MMAX)

//  ADDITIONAL FINE TUNING  //
// POLAR_OPT_THRESHOLD : value under wich the polar values of the Legendre Polynomials Plm are neglected,
// leading to increased performance (a few percent)
// 0 = no polar optimization;  1.e-14 = very safe;  1.e-10 = safe;  1.e-6 = aggresive.
#define POLAR_OPT_THRESHOLD 1.e-10
// compute and print some debugging information...
#define _SH_DEBUG_
