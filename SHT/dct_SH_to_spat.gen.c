

void SHsphtor_to_spat_dct(complex double *Slm, complex double *Tlm, complex double *BtF, complex double *BpF)
{
//	complex double *Ql;
	complex double *Sl, *Tl;
//	double *yl;
	struct DtDp *dyl;

	long int k,im,m,l;

	im=0;	m=0;
//		Ql = &Qlm[LiM(0,im)];		// virtual pointer for l=0 and im
		Sl = &Slm[LiM(0,im)];
		Tl = &Tlm[LiM(0,im)];

//		yl = ylm_dct[im];
		dyl = dylm_dct[im];
		for (k=0; k<NLAT; k++) {
//			BrF[k] = 0.0;
			BtF[k] = 0.0;		// zero out array (includes DCT padding)
			BpF[k] = 0.0;
		}
		for (l=m; l<LMAX; l+=2) {
			for (k=0; k<=l; k+=2) {
//				BrF[k]   += yl[k]   * Ql[l];
//				BrF[k+1] += yl[k+1] * Ql[l+1];
				BtF[k]   += dyl[k].t   * Sl[l+1];
				BpF[k]   -= dyl[k].t   * Tl[l+1];
				BtF[k+1] += dyl[k+1].t * Sl[l];
				BpF[k+1] -= dyl[k+1].t * Tl[l];
			}
//			yl += (l+2 - (m&1));
			dyl += l+3;	//(l+1 + (m&1));
		}
		if (l==LMAX) {
			for (k=0; k<=l; k+=2) {
//				BrF[k]   += yl[k]   * Ql[l];
				BtF[k+1] += dyl[k+1].t * Sl[l];
				BpF[k+1] -= dyl[k+1].t * Tl[l];
			}
		}
//		BrF += NLAT;
		BtF += NLAT;	BpF += NLAT;

	for (im=1;im<=MMAX;im++) {
		m=im*MRES;
//		Ql = &Qlm[LiM(0,im)];		// virtual pointer for l=0 and im
		Sl = &Slm[LiM(0,im)];
		Tl = &Tlm[LiM(0,im)];

//		yl = ylm_dct[im];
		dyl = dylm_dct[im];
		for (k=0; k<NLAT; k++) {
//			BrF[k] = 0.0;
			BtF[k] = 0.0;		// zero out array (includes DCT padding)
			BpF[k] = 0.0;
		}
		for (l=m; l<LMAX; l+=2) {
//			for (k=0; k<=l; k+=2) {
			for (k=0; k<=l+1; k+=2) {
//				BrF[k]   += yl[k]   * Ql[l];
//				BrF[k+1] += yl[k+1] * Ql[l+1];
				BtF[k]   += dyl[k].p   * (I*Tl[l])   + dyl[k].t   * Sl[l+1];
				BpF[k]   += dyl[k].p   * (I*Sl[l])   - dyl[k].t   * Tl[l+1];
				BtF[k+1] += dyl[k+1].p * (I*Tl[l+1]) + dyl[k+1].t * Sl[l];
				BpF[k+1] += dyl[k+1].p * (I*Sl[l+1]) - dyl[k+1].t * Tl[l];
			}
//			yl += (l+2 - (m&1));
			dyl += l+3;	//(l+1 + (m&1));
		}
		if (l==LMAX) {
			for (k=0; k<=l; k+=2) {
//				BrF[k]   += yl[k]   * Ql[l];
				BtF[k]   += dyl[k].p   * (I*Tl[l]);
				BpF[k]   += dyl[k].p   * (I*Sl[l]);
				BtF[k+1] += dyl[k+1].t * Sl[l];
				BpF[k+1] -= dyl[k+1].t * Tl[l];
			}
		}
//		BrF += NLAT;
		BtF += NLAT;	BpF += NLAT;
	}
	for (k=0; k < NLAT*(NPHI/2 -MMAX); k++) {	// FFT padding for high m's
//		BrF[k] = 0.0;
		BtF[k] = 0.0;	BpF[k] = 0.0;
	}

//	BrF -= NLAT*(MMAX+1);		// restore original pointer
	BtF -= NLAT*(MMAX+1);	BpF -= NLAT*(MMAX+1);	// restore original pointer
//	fftw_execute_r2r(idct,(double *) BrF, (double *) BrF);		// iDCT
	fftw_execute_r2r(idct,(double *) BtF, (double *) BtF);		// iDCT
	fftw_execute_r2r(idct,(double *) BpF, (double *) BpF);		// iDCT
  #if NPHI>1
	if (MRES & 1) {		// odd m's must be multiplied by sin(theta) which was removed from ylm's
		for (im=1; im<=MMAX; im+=2) {	// odd m's
//			for (k=0; k<NLAT; k++) BrF[im*NLAT + k] *= st[k];
		}
		for (im=0; im<=MMAX; im+=2) {	//even m's
			for (k=0; k<NLAT; k++) {
				BtF[im*NLAT + k] *= st[k];
				BpF[im*NLAT + k] *= st[k];
			}
		}
	} else {	// only even m's
		for (im=0; im<=MMAX; im++) {
			for (k=0; k<NLAT; k++) {
				BtF[im*NLAT + k] *= st[k];
				BpF[im*NLAT + k] *= st[k];
			}
		}
	}
//	fftw_execute_dft_c2r(ifft, BrF, (double *) BrF);
	fftw_execute_dft_c2r(ifft, BtF, (double *) BtF);
	fftw_execute_dft_c2r(ifft, BpF, (double *) BpF);
  #else
	im=0;
		for (k=0; k<NLAT; k++) {
			BtF[im*NLAT + k] *= st[k];
			BpF[im*NLAT + k] *= st[k];
		}
//	ifft_m0_c2r(BrF, (double *) BrF);
	ifft_m0_c2r(BtF, (double *) BtF);	ifft_m0_c2r(BpF, (double *) BpF);
  #endif
}

