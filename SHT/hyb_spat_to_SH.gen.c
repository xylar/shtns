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
  #ifndef SHT_AXISYM
QB	complex double reo[2*NLAT_2];	// symmetric (even) and anti-symmetric (odd) parts, interleaved.
VB	complex double teo[2*NLAT_2], peo[2*NLAT_2];	// theta and phi even and odd parts
Q	complex double q0,q1;
V	complex double s0,s1,t0,t1;
  #else
QB	double reo[2*NLAT_2];	// symmetric (even) and anti-symmetric (odd) parts, interleaved.
VB	double teo[2*NLAT_2], peo[2*NLAT_2];	// theta and phi even and odd parts
Q	double q0,q1;
V	double s0,s1,t0,t1;
  #endif
Q	complex double *Ql;		// virtual pointers for given im
V	complex double *Sl, *Tl;	// virtual pointers for given im
Q	double *zl;
V	double *dzl0;
V	struct DtDp *dzl;
	long int ni;
	long int i,im,l;

// defines how to access even and odd parts of data
QB	#define re	reo[2*i]
QB	#define ro	reo[2*i+1]
VB	#define te	teo[2*i]
VB	#define to	teo[2*i+1]
VB	#define pe	peo[2*i]
VB	#define po	peo[2*i+1]

 #ifndef SHT_AXISYM
Q1	#define re	BrF[i]
Q1	#define ro	BrF[i]
V1	#define te	BtF[i]
V1	#define to	BtF[i]
V1	#define pe	BpF[i]
V1	#define po	BpF[i]
Q	fftw_execute_dft_r2c(fft,(double *) BrF, BrF);
V	fftw_execute_dft_r2c(fft,(double *) BtF, BtF);
V	fftw_execute_dft_r2c(fft,(double *) BpF, BpF);
Q	#define BR0 BrF
 #else
Q1	#define re	((double *)BrF)[i]
Q1	#define ro	((double *)BrF)[i]
V1	#define te	((double *)BtF)[i]
V1	#define to	((double *)BtF)[i]
V1	#define pe	((double *)BpF)[i]
V1	#define po	((double *)BpF)[i]
Q	#define BR0 ((double *)BrF)
 #endif

	ni = NLAT_2;	// copy NLAT_2 to a local variable for faster access (inner loop limit)

	im = 0;		// dzl.p = 0.0 : and evrything is REAL
  #ifndef SHT_NO_DCT
Q	if (MTR_DCT >= 0) {		// unfortunately, only scalar SHT can be faster with DCT.
Q		fftw_execute_r2r(dctm0,(double *) BrF, (double *) BrF);		// DCT
Q		l=0;
Q		Ql = Qlm;		// virtual pointer for l=0 and im
Q		zl = zlm_dct0;
Q		while (l<LTR) {		// l has parity of m
Q			q0 = 0.0;	q1 = 0.0;
Q			for (i=l; i<NLAT; i+=2) {		// for m=0, zl coeff with i<l are zeros.
Q				(double) q0 += BR0[i]   * zl[0];
Q				(double) q1 += BR0[i+1] * zl[1];
Q				zl+=2;
Q			}
Q			Ql[l] = q0;	Ql[l+1] = q1;
Q			l+=2;
Q		}
Q		if ((LTR & 1) == 0) {	// if (l == LTR)  <=>  if ((LTR & 1) == 0) for m=0
Q			q0 = 0.0;
Q			for (i=l; i<NLAT; i+=2) {		// for m=0, DCT coeff with it<l are zeros.
Q				(double) q0   += BR0[i]   * zl[0];
Q				zl+=2;
Q			}
Q			Ql[l] = q0;
Q			l++;
Q		}
Q  #ifdef SHT_VAR_LTR
Q		while( l<=LMAX ) {
Q			Ql[l] = 0.0;
Q			l++;
Q		}
Q  #endif
Q		BrF += NLAT;
Q	} else {
  #endif
  #ifndef SHT_AXISYM
 B		for (i=0;i<NLAT/2;i++) {	// compute symmetric and antisymmetric parts.
QB			(double) reo[2*i]   = (double) BrF[i] + (double) BrF[NLAT-1-i];
QB			(double) reo[2*i+1] = (double) BrF[i] - (double) BrF[NLAT-1-i];
VB			(double) teo[2*i]   = (double) BtF[i] + (double) BtF[NLAT-1-i];
VB			(double) teo[2*i+1] = (double) BtF[i] - (double) BtF[NLAT-1-i];
VB			(double) peo[2*i]   = (double) BpF[i] + (double) BpF[NLAT-1-i];
VB			(double) peo[2*i+1] = (double) BpF[i] - (double) BpF[NLAT-1-i];
 B		}
    #ifndef SHT_NLAT_EVEN
 B		if (NLAT & 1) {		// NLAT is odd : special equator handling
QB			(double) reo[2*i] = (double) BrF[i];	(double) reo[2*i+1] = 0.0;
VB			(double) teo[2*i] = (double) BtF[i];	(double) teo[2*i+1] = 0.0;
VB			(double) peo[2*i] = (double) BpF[i];	(double) peo[2*i+1] = 0.0;
 B		}
    #endif
  #else
QB	fft_m0_r2eo((double *) BrF, reo);
VB	fft_m0_r2eo((double *) BtF, teo);	fft_m0_r2eo((double *) BpF, peo);
  #endif
		l=0;
Q		Ql = Qlm;		// virtual pointer for l=0 and im
V		Sl = Slm;	Tl = Tlm;		// virtual pointer for l=0 and im
Q		zl = zlm[im];
V		dzl0 = (double *) dzlm[im];		// only theta derivative (d/dphi = 0 for m=0)
QB		BrF += NLAT;
VB		BtF += NLAT;	BpF += NLAT;
		while (l<LTR) {		// ops : NLAT/2 * (2*(LMAX-m+1) + 4) : almost twice as fast.
QE			q0 = 0.0;
QO			q1 = 0.0;
VE			s0 = 0.0;	t1 = 0.0;
VO			t0 = 0.0;	s1 = 0.0;
			for (i=0; i < ni; i++) {
QE				(double) q0 += zl[0] * (double) re;	// Qlm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
QO				(double) q1 += zl[1] * (double) ro;	// Qlm[LiM(l+1,im)] += zlm[im][(l+1-m)*NLAT/2 + i] * fm[i];
VE				(double) s0 += dzl0[0] * (double) to;
VO				(double) t0 -= dzl0[0] * (double) po;
VO				(double) s1 += dzl0[1] * (double) te;
VE				(double) t1 -= dzl0[1] * (double) pe;
Q				zl +=2;
V				dzl0 +=2;
			}
QE			Ql[l] = q0;
QO			Ql[l+1] = q1;
VE			Sl[l] = s0;	Tl[l+1] = t1;
VO			Tl[l] = t0;	Sl[l+1] = s1; 
			l+=2;
		}
		if (l==LMAX) {
QE			q0 = 0.0;
VE			s0 = 0.0;
VO			t0 = 0.0;
			for (i=0;i<ni;i++) {
QE				(double) q0 += zl[0] * (double) re;		// Qlm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
VE				(double) s0 += dzl0[0] * (double) to;
VO				(double) t0 -= dzl0[0] * (double) po;
Q				zl ++;
V				dzl0 ++;
			}
QE			Ql[l] = q0;
VE			Sl[l] = s0;
VO			Tl[l] = t0;
  #ifdef SHT_VAR_LTR
		} else {
		    if (l==LTR) {
QE			q0 = 0.0;
VE			s0 = 0.0;
VO			t0 = 0.0;
			for (i=0; i < ni; i++) {
QE				(double) q0 += (double) re * zl[0];		// Qlm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
VE				(double) s0 += dzl0[0] * (double) to;
VO				(double) t0 -= dzl0[0] * (double) po;
Q				zl +=2;
V				dzl0 +=2;
			}
QE			Ql[l] = q0;
VE			Sl[l] = s0;
VO			Tl[l] = t0;
			l++;
		    }
		    while( l<=LMAX ) {
Q			Ql[l] = 0.0;
V			Sl[l] = 0.0;	Tl[l] = 0.0;
			l++;
		    }
  #endif
		}
  #ifndef SHT_NO_DCT
Q	}
  #endif
  #ifndef SHT_AXISYM
	for (im=1;im<=MTR;im++) {
		l=im*MRES;
 B		for (i=tm[im];i<NLAT/2;i++) {	// compute symmetric and antisymmetric parts.
QB			reo[2*i]   = BrF[i] + BrF[NLAT-1-i];
QB			reo[2*i+1] = BrF[i] - BrF[NLAT-1-i];
VB			teo[2*i]   = BtF[i] + BtF[NLAT-1-i];
VB			teo[2*i+1] = BtF[i] - BtF[NLAT-1-i];
VB			peo[2*i]   = BpF[i] + BpF[NLAT-1-i];
VB			peo[2*i+1] = BpF[i] - BpF[NLAT-1-i];
 B		}
    #ifndef SHT_NLAT_EVEN
 B		if (NLAT & 1) {		// NLAT is odd : special equator handling
QB			reo[2*i] = BrF[i];		reo[2*i+1] = 0.0;
VB			teo[2*i] = BtF[i];		teo[2*i+1] = 0.0;
VB			peo[2*i] = BpF[i];		peo[2*i+1] = 0.0;
 B		}
    #endif
Q		Ql = &Qlm[LiM(0,im)];	// virtual pointer for l=0 and im
V		Sl = &Slm[LiM(0,im)];	Tl = &Tlm[LiM(0,im)];
Q		zl = zlm[im];
V		dzl = dzlm[im];
Q		BrF += NLAT;
V		BtF += NLAT;	BpF += NLAT;
		while (l<LTR) {		// ops : NLAT/2 * (2*(LMAX-m+1) + 4) : almost twice as fast.
Q			zl += 2*tm[im];
V			dzl += 2*tm[im];
QE			q0 = 0.0;
QO			q1 = 0.0;
VE			s0 = 0.0;	t1 = 0.0;		// Slm[LiM(l,im)] = 0.0;	Slm[LiM(l+1,im)] = 0.0;
VO			t0 = 0.0;	s1 = 0.0;
			for (i=tm[im]; i < ni; i++) {	// tm[im] : polar optimization
QE				q0 += re * zl[0];		// Qlm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
QO				q1 += ro * zl[1];	// Qlm[LiM(l+1,im)] += zlm[im][(l+1-m)*NLAT/2 + i] * fm[i];
VE				s0 += dzl[0].t *to - dzl[0].p *pe*I;		// ref: these E. Dormy p 72.
VO				t0 -= dzl[0].t *po + dzl[0].p *te*I;
VO				s1 += dzl[1].t *te - dzl[1].p *po*I;
VE				t1 -= dzl[1].t *pe + dzl[1].p *to*I;
Q				zl +=2;
V				dzl +=2;
			}
QE			Ql[l] = q0;
QO			Ql[l+1] = q1;
VE			Sl[l] = s0;	Tl[l+1] = t1;
VO			Tl[l] = t0;	Sl[l+1] = s1; 
			l+=2;
		}
		if (l==LMAX) {
Q			zl += tm[im];
V			dzl += tm[im];
QE			q0 = 0.0;	// Qlm[LiM(l,im)] = 0.0;
VE			s0 = 0.0;
VO			s0 = 0.0;
			for (i=tm[im];i<ni;i++) {	// polar optimization
QE				q0 += zl[0] * re;	// Qlm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
VE				s0 += dzl[0].t *to - dzl[0].p *pe*I;
VO				t0 -= dzl[0].t *po + dzl[0].p *te*I;
Q				zl++;
V				dzl++;
			}
QE			Ql[l] = q0;
VE			Sl[l] = s0;
VO			Tl[l] = t0;
    #ifdef SHT_VAR_LTR
		} else {
		    if (l==LTR) {
Q			zl += 2*tm[im];
V			dzl += 2*tm[im];
QE			q0 = 0.0;
VE			s0 = 0.0;
VO			t0 = 0.0;
			for (i=tm[im]; i < ni; i++) {	// tm[im] : polar optimization
QE				q0 += re * zl[0];		// Qlm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
VE				s0 += dzl[0].t *to - dzl[0].p *pe*I;		// ref: these E. Dormy p 72.
VO				t0 -= dzl[0].t *po + dzl[0].p *te*I;
Q				zl +=2;
V				dzl +=2;
			}
QE			Ql[l] = q0;
VE			Sl[l] = s0;
VO			Tl[l] = t0;
			l++;
		    }
		    while( l<=LMAX ) {
Q			Ql[l] = 0.0;
V			Sl[l] = 0.0;	Tl[l] = 0.0;
			l++;
		    }
    #endif
		}
	}
  #endif
Q	#undef re
Q	#undef ro
V	#undef te
V	#undef to
V	#undef pe
V	#undef po
# }
