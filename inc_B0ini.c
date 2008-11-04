/** Poloidal/Toroidal, Spherical Harmonics definition of magnetic field
 this is C code which will be included at compile time in xshells.c
 Set_Poloidal and Set_Toroidal are macros defined in xshells.c which
 automatically scale the result with the corresponding parameter in the .par file
*/

// double rr is the radius,
// set l and m before calling Set_Poloidal and/or Set_Toroidal macro with the desired value.

{
#undef INIT_FIELD_NAME
#define INIT_FIELD_NAME "uniform vertical field (z-axis)"
l=1; m=0;
	Set_Poloidal( rr/2. * Y10_ct)	// Y10_ct  makes it unitary in physical space.
}

/*{
#undef INIT_FIELD_NAME
#define INIT_FIELD_NAME "current-free dipole"
l=1; m=0;
	Set_Poloidal( 1.0/(2.*rr*rr) )		// magnetic dipole like in E. Dormy's thesis (current free)
	if (rr < 0.1) Set_Poloidal( 1.0/ (2.*0.1*0.1) )		// avoid divergence at r=0
}*/

/*{
#undef INIT_FIELD_NAME
#define INIT_FIELD_NAME "Jault 2008 (not current-free)"
l=1; m=0;
	double v = pi * rr;
	Set_Poloidal( ((sin(v)/v-cos(v))/v - 0.3*(sin(2*v)/(2*v)-cos(2*v))/(2*v)) )		// j1(pi*r) - 0.3*j1(2*pi*r)
	if (rr == 0.0) Set_Poloidal( 0.0 )
l=3; m=0;
	v = 5.7635 * rr;	
	Set_Poloidal( -0.2*((sin(v)*(3./(v*v)-1.)-3.*cos(v)/v)/v) )	// -0.2*j3(k*r)
	if (rr == 0.0) Set_Poloidal( 0.0 )
}*/

// do not remove the following lines
#undef Set_Poloidal
#undef Set_Toroidal
// end of file : do not add lines below this one.
