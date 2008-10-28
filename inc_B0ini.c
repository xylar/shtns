/** Poloidal/Toroidal, Spherical Harmonics definition of magnetic field
 this is C code which will be included at compile time in xshells.c
 Set_Poloidal and Set_Toroidal are macros defined in xshells.c which
 automatically scale the result with the corresponding parameter in the .par file
*/

// double rr is the radius,
// set l and m before calling Set_Poloidal and/or Set_Toroidal macro with the desired value.

{
#undef INIT_FIELD_NAME
#define INIT_FIELD_NAME "current free dipole"
l=1; m=0;
	Set_Poloidal( 1.0/(2.*rr*rr) )		// magnetic dipole like in E. Dormy's thesis (current free)
	if (rr < 0.1) Set_Poloidal( 1.0/ (2.*0.1*0.1) );	// avoid divergence at r=0
}

// do not remove the following lines
#undef Set_Poloidal
#undef Set_Toroidal
// end of file : do not add lines below this one.
