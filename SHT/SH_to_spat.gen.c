# This file is meta-code for SHT.c (spherical harmonic transform).
# it is intended for "make" to generate C code for similar SHT functions,
# from one generic function + tags.
# > See Makefile and SHT.c
# Basically, there are tags at the beginning of lines that are information
# to keep or remove the line depending on the function to build.
# tags :
# Q : line for scalar transform
# V : line for vector transform (spheroidal and toroidal)
# S : line for vector transfrom, spheroidal component
# T : line for vector transform, toroidal component.

/////////////////////////////////////////////////////
//   Inverse Spherical Harmonic Transform
// input  : Qlm,Slm,Tlm = spherical harmonics coefficients of Scalar, Spheroidal and Toroidal scalars : 
//          complex double array of size NLM [unmodified]
// output : BrF, BtF, BpF = theta, and phi vector components, spatial/fourrier data : 
//          complex double array of size NLAT*(NPHI/2+1) or double array of size NLAT*(NPHI/2+1)*2
#Q void SH_to_spat(complex double *Qlm, complex double *BrF)
#V void SHsphtor_to_spat(complex double *Slm, complex double *Tlm, complex double *BtF, complex double *BpF)
# {
Q	complex double fe, fo;		// even and odd parts
S	complex double se, so, dse, dso;	// spheroidal even and odd parts
T	complex double te, to, dte, dto;	// toroidal ...
Q	complex double *Ql;
S	complex double *Sl;
T	complex double *Tl;
Q	double *yl;
V	struct DtDp *dyl;
	long int i,im,m,l;

	im = 0;		// zonal part : d/dphi = 0;
		m = im*MRES;
Q		Ql = &Qlm[LiM(0,im)];	// virtual pointer for l=0 and im
S		Sl = &Slm[LiM(0,im)];	// virtual pointer for l=0 and im
T		Tl = &Tlm[LiM(0,im)];
		i=0;
Q		yl  = ylm[im] + i*(LMAX-m+1) -m;
V		dyl = dylm[im] + i*(LMAX-m+1) -m;
		while (i < NLAT_2) {	// ops : NLAT_2 * [ (lmax-m+1)*2 + 4]	: almost twice as fast.
			l=m;
Q			fe = 0.0;	fo = 0.0;
S			dse = 0.0;	dso = 0.0;
T			dte = 0.0;	dto = 0.0;
			while (l<LTR) {	// compute even and odd parts
QE				(double) fe += yl[l] * (double) Ql[l];		// fe += ylm[im][i*(LMAX-m+1) + (l-m)] * Qlm[LiM(l,im)];
QO				(double) fo += yl[l+1] * (double) Ql[l+1];	// fo += ylm[im][i*(LMAX-m+1) + (l+1-m)] * Qlm[LiM(l+1,im)];
TO				(double) dto += dyl[l].t * (double) Tl[l];	// m=0 : everything is real.
SE				(double) dso += dyl[l].t * (double) Sl[l];
TE				(double) dte += dyl[l+1].t * (double) Tl[l+1];
SO				(double) dse += dyl[l+1].t * (double) Sl[l+1];
				l+=2;
			}
			if (l==LTR) {
QE				(double) fe += yl[l] * (double) Ql[l];		// fe += ylm[im][i*(LMAX-m+1) + (l-m)] * Qlm[LiM(l,im)];
TO				(double) dto += dyl[l].t * Tl[l];
SE				(double) dso += dyl[l].t * Sl[l];
			}
Q			BrF[i] = fe + fo;
V			BtF[i] = 0.0
S				+ (dse+dso)			// Bt = dS/dt
V				;
V			BpF[i] = 0.0
T		 		- (dte+dto)			// Bp = - dT/dt
V				;
			i++;
QB			BrF[NLAT-i] = fe - fo;
VB			BtF[NLAT-i] = 0.0
SB		 		+ (dse-dso)
VB				;
VB			BpF[NLAT-i] = 0.0
TB				- (dte-dto)
VB				;
Q			yl  += (LMAX-m+1);
V			dyl += (LMAX-m+1);
		}
Q		BrF += NLAT;
V		BtF += NLAT;	BpF += NLAT;
	for (im=1; im<=MTR; im++) {
		m = im*MRES;
Q		Ql = &Qlm[LiM(0,im)];	// virtual pointer for l=0 and im
S		Sl = &Slm[LiM(0,im)];	// virtual pointer for l=0 and im
T		Tl = &Tlm[LiM(0,im)];
		i=0;
		while (i<tm[im]) {	// polar optimization
Q			BrF[i] = 0.0;
QB			BrF[NLAT-tm[im] + i] = 0.0;	// south pole zeroes <=> BrF[im*NLAT + NLAT-(i+1)] = 0.0;
V			BtF[i] = 0.0;
VB			BtF[NLAT-tm[im] + i] = 0.0;	// south pole zeroes
V			BpF[i] = 0.0;
VB			BpF[NLAT-tm[im] + i] = 0.0;	// south pole zeroes
			i++;
		}
Q		yl  = ylm[im] + i*(LMAX-m+1) -m;
V		dyl = dylm[im] + i*(LMAX-m+1) -m;
		while (i < NLAT_2) {	// ops : NLAT_2 * [ (lmax-m+1)*2 + 4]	: almost twice as fast.
			l=m;
Q			fe = 0.0;	fo = 0.0;
S			dse = 0.0;	dso = 0.0;	se = 0.0;	so = 0.0;
T			dte = 0.0;	dto = 0.0;	te = 0.0;	to = 0.0;
			while (l<LTR) {	// compute even and odd parts
QE				fe  += yl[l] * Ql[l];		// fe += ylm[im][i*(LMAX-m+1) + (l-m)] * Qlm[LiM(l,im)];
QO				fo  += yl[l+1] * Ql[l+1];	// fo += ylm[im][i*(LMAX-m+1) + (l+1-m)] * Qlm[LiM(l+1,im)];
TO				dto += dyl[l].t * Tl[l];
TO				te  += dyl[l].p * Tl[l];
SE				dso += dyl[l].t * Sl[l];
SE				se  += dyl[l].p * Sl[l];
TE				dte += dyl[l+1].t * Tl[l+1];
TE				to  += dyl[l+1].p * Tl[l+1];
SO				dse += dyl[l+1].t * Sl[l+1];
SO				so  += dyl[l+1].p * Sl[l+1];
				l+=2;
			}
			if (l==LTR) {
QE				fe  += yl[l] * Ql[l];		// fe += ylm[im][i*(LMAX-m+1) + (l-m)] * Qlm[LiM(l,im)];
TO				dto += dyl[l].t * Tl[l];
TO				te  += dyl[l].p * Tl[l];
SE				dso += dyl[l].t * Sl[l];
SE				se  += dyl[l].p * Sl[l];
			}
Q			BrF[i] = fe + fo;
V			BtF[i] = 		// Bt = dS/dt       + I.m/sint *T
S					(dse+dso)
T					+ I*(te+to)
V					;
V			BpF[i] = 		// Bp = I.m/sint * S - dT/dt
S					I*(se+so)
T					- (dte+dto)
V					;
			i++;
QB			BrF[NLAT-i] = fe - fo;
VB			BtF[NLAT-i] =
SB					(dse-dso)
TB					+ I*(te-to)
VB					;
VB			BpF[NLAT-i] =
SB					I*(se-so)
TB					- (dte-dto)
VB					;
Q			yl  += (LMAX-m+1);
V			dyl += (LMAX-m+1);
		}
Q		BrF += NLAT;
V		BtF += NLAT;	BpF += NLAT;
	}
	for (i=0; i < NLAT*(NPHI/2 -MTR); i++) {	// padding for high m's
	//for(im=MTR+1; im<=NPHI/2; im++) {	// padding for high m's
	//	for (i=0;i<NLAT;i++) {
Q			BrF[i] = 0.0;
V			BtF[i] = 0.0;	BpF[i] = 0.0;
	//	}
Q	//	BrF += NLAT;
V	//	BtF += NLAT;	BpF += NLAT;
	}

Q	//BrF -= NLAT*(NPHI/2+1);		// restore original pointer
V	//BtF -= NLAT*(NPHI/2+1);	BpF -= NLAT*(NPHI/2+1);		// restore original pointers
Q	BrF -= NLAT*(MTR+1);		// restore original pointer
V	BtF -= NLAT*(MTR+1);	BpF -= NLAT*(MTR+1);		// restore original pointers
 #if NPHI>1
Q	fftw_execute_dft_c2r(ifft, BrF, (double *) BrF);
V	fftw_execute_dft_c2r(ifft, BtF, (double *) BtF);
V	fftw_execute_dft_c2r(ifft, BpF, (double *) BpF);
 #else
Q	ifft_m0_c2r(BrF, (double *) BrF);
V	ifft_m0_c2r(BtF, (double *) BtF);	ifft_m0_c2r(BpF, (double *) BpF);
 #endif
# }
