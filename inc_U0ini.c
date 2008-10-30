/** Poloidal/Toroidal, Spherical Harmonics definition of magnetic field
 this is C code which will be included at compile time in xshells.c
 Set_Poloidal and Set_Toroidal are macros defined in xshells.c which
 automatically scale the result with the corresponding parameter in the .par file
*/

// double rr is the radius,
// set l and m before calling Set_Poloidal and/or Set_Toroidal macro with the desired value.
// WARNING : compared to the article of Dudley & James, the Pol and Tor functions must be divided by r.

/*
{
#undef INIT_FIELD_NAME
#define INIT_FIELD_NAME "m=0 l=2 simple roll flow (eq 24 from Dudely & James 1989)"
l=2; m=0;
	Set_Toroidal( rr*sin(pi*rr) )
	Set_Poloidal( rr*sin(pi*rr) * 0.14 )
}
*/
{
#undef INIT_FIELD_NAME
#define INIT_FIELD_NAME "m=0 l=2 gubins flow (from Dudely & James 1989)"
l=2; m=0;
	Set_Toroidal( -rr*sin(2*pi*rr) * tanh(2*pi*(1-rr)) )
	Set_Poloidal( -rr*sin(2*pi*rr) * tanh(2*pi*(1-rr)) * 0.1 )
}
/*
{
#undef INIT_FIELD_NAME
#define INIT_FIELD_NAME "m=2 l=2 j2 pekeris flow (eq 20-21 Dudely & James 1989)"
l=2; m=0;
	if (rr != 0.0) {
		double z =  5.7634591968447*rr;
		double p = 5.7634591968447 * ((3.0/(z*z*z) -1.0/z)*sin(z) - 3./(z*z) * cos(z));
		Set_Poloidal( p )
		Set_Toroidal( 5.7634591968447 * p )
	}
}
*/

// do not remove the following lines
#undef Set_Poloidal
#undef Set_Toroidal
// end of file : do not add lines below this one.
