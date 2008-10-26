# This file is meta-code for SHT.c (spherical harmonic transform).
# it is intended for "make" to generate C code for 3 similar SHT functions,
# (namely spat_to_SH [X tag]), spat_to_SHsphtor [V tag], spat_to_SH3 [both X&V tags])
# from one generic function + tags.
# Basically, there are tags at the beginning of lines (X,V) that are information
# to keep or remove the line depending on the function to build. (X for scalar, V for vector, # for comment)
#
//////////////////////////////////////////////////
//  Spherical Harmonics Transform
// input  : BrF, BtF, BpF = spatial/fourrier data : complex double array of size NLAT*(NPHI/2+1) or double array of size NLAT*(NPHI/2+1)*2
// output : Qlm, Slm, Tlm = spherical harmonics coefficients : complex double array of size NLM
#X void spat_to_SH(complex double *BrF, complex double *Qlm)
#V void spat_to_SHsphtor(complex double *BtF, complex double *BpF, complex double *Slm, complex double *Tlm)
# {
X	complex double reo[2*NLAT_2];	// symmetric (even) and anti-symmetric (odd) parts, interleaved.
V	complex double teo[2*NLAT_2], peo[2*NLAT_2];	// theta and phi even and odd parts
X	complex double *Ql;		// virtual pointers for given im
V	complex double *Sl, *Tl;		// virtual pointers for given im
X	double *zl;
V	struct DtDp *dzl;
	long int i,im,m,l;

 #if NPHI > 1
X	fftw_execute_dft_r2c(fft,(double *) BrF, BrF);
V	fftw_execute_dft_r2c(fft,(double *) BtF, BtF);
V	fftw_execute_dft_r2c(fft,(double *) BpF, BpF);
 #else
X	fft_m0_r2c((double *) BrF, BrF);
V	fft_m0_r2c((double *) BtF, BtF);	fft_m0_r2c((double *) BpF, BpF);
 #endif

	im = 0;		// dzl.p = 0.0 : and evrything is REAL
		m=im*MRES;
		for (i=0;i<NLAT/2;i++) {	// compute symmetric and antisymmetric parts.
X			(double) reo[2*i]   = (double) BrF[i] + (double) BrF[NLAT-(i+1)];
X			(double) reo[2*i+1] = (double) BrF[i] - (double) BrF[NLAT-(i+1)];
V			(double) teo[2*i]   = (double) BtF[i] + (double) BtF[NLAT-(i+1)];
V			(double) teo[2*i+1] = (double) BtF[i] - (double) BtF[NLAT-(i+1)];
V			(double) peo[2*i]   = (double) BpF[i] + (double) BpF[NLAT-(i+1)];
V			(double) peo[2*i+1] = (double) BpF[i] - (double) BpF[NLAT-(i+1)];
		}
		if (i < NLAT_2) {		// NLAT is odd : special equator handling
X			(double) reo[2*i] = (double) BrF[i];	(double) reo[2*i+1] = 0.0;
V			(double) teo[2*i] = (double) BtF[i];	(double) teo[2*i+1] = 0.0;
V			(double) peo[2*i] = (double) BpF[i];	(double) peo[2*i+1] = 0.0;
		}
		l=m;
X		Ql = &Qlm[LiM(0,im)];		// virtual pointer for l=0 and im
V		Sl = &Slm[LiM(0,im)];	Tl = &Tlm[LiM(0,im)];		// virtual pointer for l=0 and im
X		zl = zlm[im];
V		dzl = dzlm[im];
X		BrF += NLAT;
V		BtF += NLAT;	BpF += NLAT;
		while (l<LMAX) {		// ops : NLAT/2 * (2*(LMAX-m+1) + 4) : almost twice as fast.
X			Ql[l] = 0.0;	Ql[l+1] = 0.0;		// Qlm[LiM(l,im)] = 0.0;	Qlm[LiM(l+1,im)] = 0.0;
V			Sl[l] = 0.0;	Sl[l+1] = 0.0;
V			Tl[l] = 0.0;	Tl[l+1] = 0.0;
			for (i=0; i < 2*NLAT_2; i+=2) {
X				(double) Ql[l] += (double) reo[i] * zl[i];		// Qlm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
X				(double) Ql[l+1] += (double) reo[i+1] * zl[i+1];	// Qlm[LiM(l+1,im)] += zlm[im][(l+1-m)*NLAT/2 + i] * fm[i];				
V				(double) Sl[l]   += dzl[i].t * (double) teo[i+1];
V				(double) Tl[l]   -= dzl[i].t * (double) peo[i+1];
V				(double) Sl[l+1] += dzl[i+1].t * (double) teo[i];
V				(double) Tl[l+1] -= dzl[i+1].t * (double) peo[i];
			}
			l+=2;
X			zl += 2*NLAT_2;
V			dzl += 2*NLAT_2;
		}
		if (l==LMAX) {
X			Ql[l] = 0.0;
V			Sl[l] = 0.0;	Tl[l] = 0.0;
			for (i=0;i<NLAT_2;i++) {
X				(double) Ql[l] += zl[i] * (double) reo[2*i];	// Qlm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
V				(double) Sl[l] += dzl[i].t * (double) teo[2*i+1];
V				(double) Tl[l] -= dzl[i].t * (double) peo[2*i+1];
			}
		}
	for (im=1;im<=MMAX;im++) {
		m=im*MRES;
		for (i=tm[im];i<NLAT/2;i++) {	// compute symmetric and antisymmetric parts.
X			reo[2*i]   = BrF[i] + BrF[NLAT-(i+1)];
X			reo[2*i+1] = BrF[i] - BrF[NLAT-(i+1)];
V			teo[2*i]   = BtF[i] + BtF[NLAT-(i+1)];
V			teo[2*i+1] = BtF[i] - BtF[NLAT-(i+1)];
V			peo[2*i]   = BpF[i] + BpF[NLAT-(i+1)];
V			peo[2*i+1] = BpF[i] - BpF[NLAT-(i+1)];
		}
		if (i<NLAT_2) {		// NLAT is odd : special equator handling
X			reo[2*i] = BrF[i];		reo[2*i+1] = 0.0;
V			teo[2*i] = BtF[i];		teo[2*i+1] = 0.0;
V			peo[2*i] = BpF[i];		peo[2*i+1] = 0.0;
		}
		l=m;
X		Ql = &Qlm[LiM(0,im)];		// virtual pointer for l=0 and im
V		Sl = &Slm[LiM(0,im)];	Tl = &Tlm[LiM(0,im)];		// virtual pointer for l=0 and im		
X		zl = zlm[im];
V		dzl = dzlm[im];
X		BrF += NLAT;
V		BtF += NLAT;	BpF += NLAT;
		while (l<LMAX) {		// ops : NLAT/2 * (2*(LMAX-m+1) + 4) : almost twice as fast.
X			Ql[l] = 0.0;	Ql[l+1] = 0.0;		// Qlm[LiM(l,im)] = 0.0;	Qlm[LiM(l+1,im)] = 0.0;
V			Sl[l] = 0.0;	Sl[l+1] = 0.0;		// Slm[LiM(l,im)] = 0.0;	Slm[LiM(l+1,im)] = 0.0;
V			Tl[l] = 0.0;	Tl[l+1] = 0.0;
			for (i=tm[im]*2; i < 2*NLAT_2; i+=2) {	// tm[im] : polar optimization
X				Ql[l] += reo[i] * zl[i];		// Qlm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
X				Ql[l+1] += reo[i+1] * zl[i+1];	// Qlm[LiM(l+1,im)] += zlm[im][(l+1-m)*NLAT/2 + i] * fm[i];
V				Sl[l]   += dzl[i].t *teo[i+1] - dzl[i].p *peo[i]*I;		// ref: these E. Dormy p 72.
V				Tl[l]   -= dzl[i].t *peo[i+1] + dzl[i].p *teo[i]*I;				
V				Sl[l+1] += dzl[i+1].t *teo[i] - dzl[i+1].p *peo[i+1]*I;
V				Tl[l+1] -= dzl[i+1].t *peo[i] + dzl[i+1].p *teo[i+1]*I;
			}
			l+=2;
X			zl += 2*NLAT_2;
V			dzl += 2*NLAT_2;
		}
		if (l==LMAX) {
X			Ql[l] = 0.0;	// Qlm[LiM(l,im)] = 0.0;
V			Sl[l] = 0.0;	Tl[l] = 0.0;
			for (i=tm[im];i<NLAT_2;i++) {	// polar optimization
X				Ql[l] += zl[i] * reo[2*i];	// Qlm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
V				Sl[l] += dzl[i].t *teo[2*i+1] - dzl[i].p *peo[2*i]*I;
V				Tl[l] -= dzl[i].t *peo[2*i+1] + dzl[i].p *teo[2*i]*I;
			}
		}
	}
# }