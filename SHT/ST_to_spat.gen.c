# This file is meta-code for SHT.c (spherical harmonic transform).
# it is intended for "make" to generates C code for 3 similar SHT functions,
# (namely SHsphtor_to_spat, SHsph_to_spat, SHtor_to_spat)
# from one generic function + tags.
# Basically, there are tags at the beginning of lines (S,T) that are information
# to keep or remove the line depending on the function to build.

/////////////////////////////////////////////////////
//   Spheroidal/Toroidal to (theta,phi) components inverse Spherical Harmonics Transform
// input  : Slm,Tlm = spherical harmonics coefficients of Spheroidal and Toroidal scalars : 
//          complex double array of size NLM [unmodified]
// output : BtF, BpF = theta, and phi vector components, spatial/fourrier data : 
//          complex double array of size NLAT*(NPHI/2+1) or double array of size NLAT*(NPHI/2+1)*2

#void SHsphtor_to_spat(complex double *Slm, complex double *Tlm, complex double *BtF, complex double *BpF)
#{
S	complex double se, so, dse, dso;	// spheroidal even and odd parts
T	complex double te, to, dte, dto;	// toroidal ...
S	complex double *Sl;
T	complex double *Tl;
	struct DtDp *dyl;
	long int i,im,m,l;

	im = 0;		// zonal part : d/dphi = 0;
		m = im*MRES;
S		Sl = &Slm[LiM(0,im)];	// virtual pointer for l=0 and im
T		Tl = &Tlm[LiM(0,im)];
		i=0;
		dyl = dylm[im] + i*(LMAX-m+1) -m;
		while (i<NLAT/2) {	// ops : NLAT/2 * [ (lmax-m+1)*2 + 4]	: almost twice as fast.
			l=m;
S			dse = 0.0;	dso = 0.0;
T			dte = 0.0;	dto = 0.0;
			while (l<LMAX) {	// compute even and odd parts
T				(double) dto += dyl[l].t * (double) Tl[l];	// m=0 : everything is real.
S				(double) dso += dyl[l].t * (double) Sl[l];
T				(double) dte += dyl[l+1].t * (double) Tl[l+1];
S				(double) dse += dyl[l+1].t * (double) Sl[l+1];
				l+=2;
			}
			if (l==LMAX) {
T				(double) dto += dyl[l].t * Tl[l];
S				(double) dso += dyl[l].t * Sl[l];
			}
			BtF[i] = 0.0
S				+ (dse+dso)			// Bt = dS/dt
				;
			BtF[NLAT-(i+1)] = 0.0
S		 		+ (dse-dso)
				;
			BpF[i] = 0.0
T		 		- (dte+dto)			// Bp = - dT/dt
				;
			BpF[NLAT-(i+1)] = 0.0
T				- (dte-dto)
				;
			i++;
			dyl += (LMAX-m+1);
		}
		BtF += NLAT;	BpF += NLAT;
	for (im=1; im<=MMAX; im++) {
		m = im*MRES;
S		Sl = &Slm[LiM(0,im)];	// virtual pointer for l=0 and im
T		Tl = &Tlm[LiM(0,im)];
		i=0;
		while (i<tm[im]) {	// polar optimization
			BtF[i] = 0.0;
			BtF[NLAT-tm[im] + i] = 0.0;	// south pole zeroes
			BpF[i] = 0.0;
			BpF[NLAT-tm[im] + i] = 0.0;	// south pole zeroes
			i++;
		}
		dyl = dylm[im] + i*(LMAX-m+1) -m;
		while (i<NLAT/2) {	// ops : NLAT/2 * [ (lmax-m+1)*2 + 4]	: almost twice as fast.
			l=m;
S			dse = 0.0;	dso = 0.0;	se = 0.0;	so = 0.0;
T			dte = 0.0;	dto = 0.0;	te = 0.0;	to = 0.0;
			while (l<LMAX) {	// compute even and odd parts
T				dto += dyl[l].t * Tl[l];
T				te += dyl[l].p * Tl[l];
S				dso += dyl[l].t * Sl[l];
S				se += dyl[l].p * Sl[l];
T				dte += dyl[l+1].t * Tl[l+1];
T				to += dyl[l+1].p * Tl[l+1];
S				dse += dyl[l+1].t * Sl[l+1];
S				so += dyl[l+1].p * Sl[l+1];
				l+=2;
			}
			if (l==LMAX) {
T				dto += dyl[l].t * Tl[l];
T				te += dyl[l].p * Tl[l];
S				dso += dyl[l].t * Sl[l];
S				se += dyl[l].p * Sl[l];
			}
			BtF[i] = 		// Bt = dS/dt       + I.m/sint *T
S					(dse+dso)
T					+ I*(te+to)
					;
			BtF[NLAT-(i+1)] =
S					(dse-dso)
T					+ I*(te-to)
					;
			BpF[i] = 		// Bp = I.m/sint * S - dT/dt
S					I*(se+so)
T					- (dte+dto)
					;
			BpF[NLAT-(i+1)] =
S					I*(se-so)
T					- (dte-dto)
					;
			i++;
			dyl += (LMAX-m+1);
		}
		BtF += NLAT;	BpF += NLAT;
	}
	for(im=MMAX+1; im<=NPHI/2; im++) {	// padding for high m's
		for (i=0;i<NLAT;i++) {
			BtF[i] = 0.0;	BpF[i] = 0.0;
		}
		BtF += NLAT;	BpF += NLAT;
	}

	BtF -= NLAT*(NPHI/2+1);		// restore original pointers
	BpF -= NLAT*(NPHI/2+1);
	fftw_execute_dft_c2r(ifft, BtF, (double *) BtF);
	fftw_execute_dft_c2r(ifft, BpF, (double *) BpF);
#}
