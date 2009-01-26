# This file is meta-code for SHT.c (spherical harmonic transform).
# it is intended for "make" to generate C code for 3 similar SHT functions,
# (namely spat_to_SH [Q tag]), spat_to_SHsphtor [V tag], spat_to_SH3 [both Q&V tags])
# from one generic function + tags.
# Basically, there are tags at the beginning of lines (Q,V) that are information
# to keep or remove the line depending on the function to build. (Q for scalar, V for vector, # for comment)
#
//////////////////////////////////////////////////
//  Spherical Harmonics Transform
// input  : BrF, BtF, BpF = spatial/fourrier data : complex double array of size NLAT*(NPHI/2+1) or double array of size NLAT*(NPHI/2+1)*2
// output : Qlm, Slm, Tlm = spherical harmonics coefficients : complex double array of size NLM
#Q void spat_to_SH(complex double *BrF, complex double *Qlm)
#V void spat_to_SHsphtor(complex double *BtF, complex double *BpF, complex double *Slm, complex double *Tlm)
# {
  #if NPHI > 1
QB	complex double reo[2*NLAT_2];	// symmetric (even) and anti-symmetric (odd) parts, interleaved.
VB	complex double teo[2*NLAT_2], peo[2*NLAT_2];	// theta and phi even and odd parts
  #else
QB	double reo[2*NLAT_2];	// symmetric (even) and anti-symmetric (odd) parts, interleaved.
VB	double teo[2*NLAT_2], peo[2*NLAT_2];	// theta and phi even and odd parts
  #endif
Q	complex double *Ql;		// virtual pointers for given im
V	complex double *Sl, *Tl;	// virtual pointers for given im
Q	double *zl;
V	struct DtDp *dzl;
	long int i,im,m,l;

// defines how to access even and odd parts of data
Q1	#define re	BrF[i]
Q1	#define ro	BrF[i]
QB	#define re	reo[2*i]
QB	#define ro	reo[2*i+1]
V1	#define te	BtF[i]
V1	#define to	BtF[i]
VB	#define te	teo[2*i]
VB	#define to	teo[2*i+1]
V1	#define pe	BpF[i]
V1	#define po	BpF[i]
VB	#define pe	peo[2*i]
VB	#define po	peo[2*i+1]

 #if NPHI > 1
Q	fftw_execute_dft_r2c(fft,(double *) BrF, BrF);
V	fftw_execute_dft_r2c(fft,(double *) BtF, BtF);
V	fftw_execute_dft_r2c(fft,(double *) BpF, BpF);
 #else
Q	//fft_m0_r2c((double *) BrF, BrF);
V	//fft_m0_r2c((double *) BtF, BtF);	fft_m0_r2c((double *) BpF, BpF);
Q	fft_m0_r2eo((double *) BrF, reo);
V	fft_m0_r2eo((double *) BtF, teo);	fft_m0_r2eo((double *) BpF, peo);
 #endif
	im = 0;	m=0;	// dzl.p = 0.0 : and evrything is REAL
/*  #ifdef SHT_DCT
Q		fftw_execute_r2r(dctm0,(double *) BrF, (double *) BrF);		// DCT
		l=0;
Q		Ql = Qlm;		// virtual pointer for l=0 and im
V		Sl = Slm;	Tl = Tlm;		// virtual pointer for l=0 and im
Q		zl = zlm_dct[im];
		while (l<LTR) {		// l has parity of m
Q			Ql[l] = 0.0;	Ql[l+1] = 0.0;
			for (k=l; k<=NLAT; k+=2) {		// for m=0, zl coeff with k<l are zeros.
Q				(double) Ql[l]   += (double) BrF[k]   * zl[k];
Q				(double) Ql[l+1] += (double) BrF[k+1] * zl[k+1];
			}
			l+=2;
Q			zl += NLAT;
		}
		if ((LTR & 1) == 0) {	// if (l == LTR)  <=>  if ((LTR & 1) == 0) for m=0
Q			Ql[l] = 0.0;
			for (k=l; k<=NLAT; k+=2) {		// for m=0, DCT coeff with k<l are zeros.
Q				(double) Ql[l]   += (double) BrF[k]   * zl[k];
			}
		}
Q		BrF += NLAT;
  #else	*/
  #if NPHI > 1
 B		for (i=0;i<NLAT/2;i++) {	// compute symmetric and antisymmetric parts.
QB			(double) reo[2*i]   = (double) BrF[i] + (double) BrF[NLAT-(i+1)];
QB			(double) reo[2*i+1] = (double) BrF[i] - (double) BrF[NLAT-(i+1)];
VB			(double) teo[2*i]   = (double) BtF[i] + (double) BtF[NLAT-(i+1)];
VB			(double) teo[2*i+1] = (double) BtF[i] - (double) BtF[NLAT-(i+1)];
VB			(double) peo[2*i]   = (double) BpF[i] + (double) BpF[NLAT-(i+1)];
VB			(double) peo[2*i+1] = (double) BpF[i] - (double) BpF[NLAT-(i+1)];
 B		}
 B		if (i < NLAT_2) {		// NLAT is odd : special equator handling
QB			(double) reo[2*i] = (double) BrF[i];	(double) reo[2*i+1] = 0.0;
VB			(double) teo[2*i] = (double) BtF[i];	(double) teo[2*i+1] = 0.0;
VB			(double) peo[2*i] = (double) BpF[i];	(double) peo[2*i+1] = 0.0;
 B		}
  #endif
		l=m;
Q		Ql = Qlm;		// virtual pointer for l=0 and im
V		Sl = Slm;	Tl = Tlm;		// virtual pointer for l=0 and im
Q		zl = zlm[im];
V		dzl = dzlm[im];
QB		BrF += NLAT;
VB		BtF += NLAT;	BpF += NLAT;
		while (l<LTR) {		// ops : NLAT/2 * (2*(LMAX-m+1) + 4) : almost twice as fast.
QE			Ql[l] = 0.0;
QO			Ql[l+1] = 0.0;
VE			Sl[l] = 0.0;	Tl[l+1] = 0.0;
VO			Tl[l] = 0.0;	Sl[l+1] = 0.0;
			for (i=0; i < NLAT_2; i++) {
QE				(double) Ql[l]   += zl[2*i]   * (double) re;		// Qlm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
QO				(double) Ql[l+1] += zl[2*i+1] * (double) ro;	// Qlm[LiM(l+1,im)] += zlm[im][(l+1-m)*NLAT/2 + i] * fm[i];
VE				(double) Sl[l]   += dzl[2*i].t   * (double) to;
VO				(double) Tl[l]   -= dzl[2*i].t   * (double) po;
VO				(double) Sl[l+1] += dzl[2*i+1].t * (double) te;
VE				(double) Tl[l+1] -= dzl[2*i+1].t * (double) pe;
			}
			l+=2;
Q			zl += 2*NLAT_2;
V			dzl += 2*NLAT_2;
		}
		if (l==LMAX) {
QE			Ql[l] = 0.0;
VE			Sl[l] = 0.0;
VO			Tl[l] = 0.0;
			for (i=0;i<NLAT_2;i++) {
QE				(double) Ql[l] += zl[i] * (double) re;		// Qlm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
VE				(double) Sl[l] += dzl[i].t * (double) to;
VO				(double) Tl[l] -= dzl[i].t * (double) po;
			}
		} else if (l==LTR) {
QE			Ql[l] = 0.0;
VE			Sl[l] = 0.0;
VO			Tl[l] = 0.0;
			for (i=0; i < NLAT_2; i++) {
QE				(double) Ql[l] += (double) re * zl[2*i];		// Qlm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
VE				(double) Sl[l]   += dzl[2*i].t * (double) to;
VO				(double) Tl[l]   -= dzl[2*i].t * (double) po;
			}
		}
//  #endif
	for (im=1;im<=MTR;im++) {
		m=im*MRES;
 B		for (i=tm[im];i<NLAT/2;i++) {	// compute symmetric and antisymmetric parts.
QB			reo[2*i]   = BrF[i] + BrF[NLAT-(i+1)];
QB			reo[2*i+1] = BrF[i] - BrF[NLAT-(i+1)];
VB			teo[2*i]   = BtF[i] + BtF[NLAT-(i+1)];
VB			teo[2*i+1] = BtF[i] - BtF[NLAT-(i+1)];
VB			peo[2*i]   = BpF[i] + BpF[NLAT-(i+1)];
VB			peo[2*i+1] = BpF[i] - BpF[NLAT-(i+1)];
 B		}
 B		if (i<NLAT_2) {		// NLAT is odd : special equator handling
QB			reo[2*i] = BrF[i];		reo[2*i+1] = 0.0;
VB			teo[2*i] = BtF[i];		teo[2*i+1] = 0.0;
VB			peo[2*i] = BpF[i];		peo[2*i+1] = 0.0;
 B		}
		l=m;
Q		Ql = &Qlm[LiM(0,im)];		// virtual pointer for l=0 and im
V		Sl = &Slm[LiM(0,im)];	Tl = &Tlm[LiM(0,im)];		// virtual pointer for l=0 and im
Q		zl = zlm[im];
V		dzl = dzlm[im];
Q		BrF += NLAT;
V		BtF += NLAT;	BpF += NLAT;
		while (l<LTR) {		// ops : NLAT/2 * (2*(LMAX-m+1) + 4) : almost twice as fast.
QE			Ql[l] = 0.0;
QO			Ql[l+1] = 0.0;
VE			Sl[l] = 0.0;	Tl[l+1] = 0.0;		// Slm[LiM(l,im)] = 0.0;	Slm[LiM(l+1,im)] = 0.0;
VO			Tl[l] = 0.0;	Sl[l+1] = 0.0;
			for (i=tm[im]; i < NLAT_2; i++) {	// tm[im] : polar optimization
QE				Ql[l]   += re * zl[2*i];		// Qlm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
QO				Ql[l+1] += ro * zl[2*i+1];	// Qlm[LiM(l+1,im)] += zlm[im][(l+1-m)*NLAT/2 + i] * fm[i];
VE				Sl[l]   += dzl[2*i].t *to - dzl[2*i].p *pe*I;		// ref: these E. Dormy p 72.
VO				Tl[l]   -= dzl[2*i].t *po + dzl[2*i].p *te*I;
VO				Sl[l+1] += dzl[2*i+1].t *te - dzl[2*i+1].p *po*I;
VE				Tl[l+1] -= dzl[2*i+1].t *pe + dzl[2*i+1].p *to*I;
			}
			l+=2;
Q			zl += 2*NLAT_2;
V			dzl += 2*NLAT_2;
		}
		if (l==LMAX) {
QE			Ql[l] = 0.0;	// Qlm[LiM(l,im)] = 0.0;
VE			Sl[l] = 0.0;
VO			Tl[l] = 0.0;
			for (i=tm[im];i<NLAT_2;i++) {	// polar optimization
QE				Ql[l] += zl[i] * re;	// Qlm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
VE				Sl[l] += dzl[i].t *to - dzl[i].p *pe*I;
VO				Tl[l] -= dzl[i].t *po + dzl[i].p *te*I;
			}
		} else if (l==LTR) {
QE			Ql[l] = 0.0;
VE			Sl[l] = 0.0;
VO			Tl[l] = 0.0;
			for (i=tm[im]; i < NLAT_2; i++) {	// tm[im] : polar optimization
QE				Ql[l]   += re * zl[2*i];		// Qlm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
VE				Sl[l]   += dzl[2*i].t *to - dzl[2*i].p *pe*I;		// ref: these E. Dormy p 72.
VO				Tl[l]   -= dzl[2*i].t *po + dzl[2*i].p *te*I;
			}
		}
	}

Q	#undef re
Q	#undef ro
V	#undef te
V	#undef to
V	#undef pe
V	#undef po
# }