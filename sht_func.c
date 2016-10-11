/*
 * Copyright (c) 2010-2016 Centre National de la Recherche Scientifique.
 * written by Nathanael Schaeffer (CNRS, ISTerre, Grenoble, France).
 * 
 * nathanael.schaeffer@univ-grenoble-alpes.fr
 * 
 * This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software. You can use,
 * modify and/or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 * 
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 * 
 */

/** \internal \file sht_func.c
 * \brief Rotation of Spherical Harmonics.
 */


/** \addtogroup rotation Rotation of SH fields.
Rotation around axis other than Z should be considered of beta quality (they have been tested but may still contain bugs).
They also require \c mmax = \c lmax. They use an Algorithm inspired by the pseudospectral rotation described in
Gimbutas Z. and Greengard L. 2009 "A fast and stable method for rotating spherical harmonic expansions" <i>Journal of Computational Physics</i>.
doi:<a href="http://dx.doi.org/10.1016/j.jcp.2009.05.014">10.1016/j.jcp.2009.05.014</a>

These functions do only require a call to \ref shtns_create, but not to \ref shtns_set_grid.
*/
//@{

/// Rotate a SH representation Qlm around the z-axis by angle alpha (in radians),
/// which is the same as rotating the reference frame by angle -alpha.
/// Result is stored in Rlm (which can be the same array as Qlm).
void SH_Zrotate(shtns_cfg shtns, cplx *Qlm, double alpha, cplx *Rlm)
{
	int im, l, lmax, mmax, mres;

	lmax = shtns->lmax;		mmax = shtns->mmax;		mres = shtns->mres;

	if (Rlm != Qlm) {		// copy m=0 which does not change.
		l=0;	do { Rlm[l] = Qlm[l]; } while(++l <= lmax);
	}
	for (int im=1; im<=mmax; im++) {
		cplx eima = cos(im*mres*alpha) - I*sin(im*mres*alpha);		// rotate reference frame by angle -alpha
		for (l=im*mres; l<=lmax; ++l)	Rlm[LiM(shtns, l, im)] = Qlm[LiM(shtns, l, im)] * eima;
	}
}

//@}

/// \internal initialize pseudo-spectral rotations
static void SH_rotK90_init(shtns_cfg shtns)
{
	cplx *q;
	double *q0;
	int nfft, nrembed, ncembed;
	
//	if ((shtns->mres != 1) || (shtns->mmax != shtns->lmax)) runerr("Arbitrary rotations require lmax=mmax and mres=1");

#define NWAY 4

	const int lmax = shtns->lmax;
	const int fac = 2*VSIZE2*NWAY;		// we need a multiple of 2*VSIZE2*NWAY ...
	const int ntheta = fft_int( ((lmax+fac)/fac) , 7) * fac;		// ... and also an fft-friendly value

	// generate the equispaced grid for synthesis
	shtns->ct_rot = malloc( sizeof(double)*ntheta );
	shtns->st_rot = shtns->ct_rot + (ntheta/2);
	for (int k=0; k<ntheta/2; ++k) {
		double cost = cos(((0.5*M_PI)*(2*k+1))/ntheta);
		double sint = sqrt((1.0-cost)*(1.0+cost));
		shtns->ct_rot[k] = cost;
		shtns->st_rot[k] = sint;
	}

	// plan FFT
	size_t sze = sizeof(double)*(2*ntheta+2)*lmax;
	q0 = VMALLOC(sze);		// alloc.
	#ifdef OMP_FFTW
		int k = (lmax < 63) ? 1 : shtns->nthreads;
		fftw_plan_with_nthreads(k);
	#endif
	q = (cplx*) q0;		// in-place FFT
	nfft = 2*ntheta;	nrembed = nfft+2;		ncembed = nrembed/2;
	shtns->fft_rot = fftw_plan_many_dft_r2c(1, &nfft, lmax, q0, &nrembed, lmax, 1, q, &ncembed, lmax, 1, FFTW_MEASURE);

	VFREE(q0);
	shtns->npts_rot = ntheta;		// save ntheta, and mark as initialized.
}

/** \internal rotation kernel used by SH_Yrotate90(), SH_Xrotate90() and SH_rotate().
 Algorithm based on the pseudospectral rotation[1] :
 - rotate around Z by angle dphi0.
 - synthetize for each l the spatial description for phi=0 and phi=pi on an equispaced latitudinal grid.
 - Fourier ananlyze as data on the equator to recover the m in the 90 degrees rotated frame.
 - rotate around new Z by angle dphi1.
 [1] Gimbutas Z. and Greengard L. 2009 "A fast and stable method for rotating spherical harmonic expansions" Journal of Computational Physics. **/
static void SH_rotK90(shtns_cfg shtns, cplx *Qlm, cplx *Rlm, double dphi0, double dphi1)
{
//	if (shtns->npts_rot == 0)	SH_rotK90_init(shtns);

//	ticks tik0, tik1, tik2, tik3;

	const int lmax = shtns->lmax;
	const int ntheta = shtns->npts_rot;
	size_t sze = sizeof(double)*(2*ntheta+2)*lmax;
	double* const q0 = VMALLOC(sze);		// alloc.

	// rotate around Z by dphi0,  and also pre-multiply imaginary parts by m
	if (Rlm != Qlm) {		// copy m=0 which does not change.
		long l=0;	do { Rlm[l] = Qlm[l]; } while(++l <= lmax);
	}
	for (int m=1; m<=lmax; m++) {
		cplx eima = cos(m*dphi0) - I*sin(m*dphi0);		// rotate reference frame by angle -dphi0
		long lm = LiM(shtns,m,m);
		double em = m;
		for (long l=m; l<=lmax; ++l) {
			cplx qrot = Qlm[lm] * eima;
			((double*)Rlm)[2*lm]   = creal(qrot);
			((double*)Rlm)[2*lm+1] = cimag(qrot) * em;		// multiply every imaginary part by m  (part of im/sin(theta)*Ylm)
			lm++;
		}
	}
	Qlm = Rlm;

//	tik0 = getticks();

		rnd* const qve = (rnd*) VMALLOC( sizeof(rnd)*NWAY*4*lmax );	// vector buffer
		rnd* const qvo = qve + NWAY*2*lmax;		// for odd m
		double* const ct = shtns->ct_rot;
		double* const st = shtns->st_rot;
		double* const alm = shtns->alm;
		const long nk = ntheta/(2*VSIZE2);		// ntheta is a multiple of (2*VSIZE2)
		long k = 0;
		do {
			rnd cost[NWAY], y0[NWAY], y1[NWAY];
			long l=0;
			// m=0
			double*	al = alm;
			for (int j=0; j<NWAY; ++j) {
				y0[j] = vall(al[0]) / vread(st, j+k);		// l=0  (discarded) DIVIDED by sin(theta) [to be consistent with m>0]
				cost[j] = vread(ct, j+k);
			}
			for (int j=0; j<NWAY; ++j) {
				y1[j]  = vall(al[1]) * y0[j] * cost[j];
			}
			al += 2;	l+=2;
			while(l<=lmax) {
				for (int j=0; j<NWAY; ++j) {
					qve[ (l-2)*2*NWAY + 2*j]   = y1[j] * vall(creal(Qlm[l-1]));	// l-1
					qve[ (l-2)*2*NWAY + 2*j+1] = vall(0.0);
					qvo[ (l-2)*2*NWAY + 2*j]   = vall(0.0);
					qvo[ (l-2)*2*NWAY + 2*j+1] = vall(0.0);
					y0[j]  = vall(al[1])*(cost[j]*y1[j]) + vall(al[0])*y0[j];
				}
				for (int j=0; j<NWAY; ++j) {
					qve[ (l-1)*2*NWAY + 2*j]   = y0[j] * vall(creal(Qlm[l]));	// l
					qve[ (l-1)*2*NWAY + 2*j+1] = vall(0.0);
					qvo[ (l-1)*2*NWAY + 2*j]   = vall(0.0);
					qvo[ (l-1)*2*NWAY + 2*j+1] = vall(0.0);
					y1[j]  = vall(al[3])*(cost[j]*y0[j]) + vall(al[2])*y1[j];
				}
				al+=4;	l+=2;
			}
			if (l==lmax+1) {
				for (int j=0; j<NWAY; ++j) {
					qve[ (l-2)*2*NWAY + 2*j]   = y1[j] * vall(creal(Qlm[l-1]));	// l-1
					qve[ (l-2)*2*NWAY + 2*j+1] = vall(0.0);
					qvo[ (l-2)*2*NWAY + 2*j]   = vall(0.0);
					qvo[ (l-2)*2*NWAY + 2*j+1] = vall(0.0);
				}
			}
			// m > 0
			for (long m=1; m<=lmax; ++m) {
				rnd* qv = qve;
				if (m&1) qv = qvo;		// store even and odd m separately.
				double*	al = shtns->alm + m*(2*(lmax+1) -m+1);
				cplx* Ql = &Qlm[LiM(shtns, 0,m)];	// virtual pointer for l=0 and m
				rnd cost[NWAY], y0[NWAY], y1[NWAY];
				for (int j=0; j<NWAY; ++j) {
					cost[j] = vread(st, k+j);
					y0[j] = vall(2.0);		// *2 for m>0
				}
				long l=m-1;
				long int ny = 0;
				  if ((int)lmax <= SHT_L_RESCALE_FLY) {
					do {		// sin(theta)^(m-1)
						if (l&1) for (int j=0; j<NWAY; ++j) y0[j] *= cost[j];
						for (int j=0; j<NWAY; ++j) cost[j] *= cost[j];
					} while(l >>= 1);
				  } else {
					long int nsint = 0;
					do {		// sin(theta)^(m-1)		(use rescaling to avoid underflow)
						if (l&1) {
							for (int j=0; j<NWAY; ++j) y0[j] *= cost[j];
							ny += nsint;
							if (vlo(y0[0]) < (SHT_ACCURACY+1.0/SHT_SCALE_FACTOR)) {
								ny--;
								for (int j=0; j<NWAY; ++j) y0[j] *= vall(SHT_SCALE_FACTOR);
							}
						}
						for (int j=0; j<NWAY; ++j) cost[j] *= cost[j];
						nsint += nsint;
						if (vlo(cost[0]) < 1.0/SHT_SCALE_FACTOR) {
							nsint--;
							for (int j=0; j<NWAY; ++j) cost[j] *= vall(SHT_SCALE_FACTOR);
						}
					} while(l >>= 1);
				  }
				for (int j=0; j<NWAY; ++j) {
					y0[j] *= vall(al[0]);
					cost[j] = vread(ct, j+k);
				}
				for (int j=0; j<NWAY; ++j) {
					y1[j]  = (vall(al[1])*y0[j]) *cost[j];
				}
				l=m;		al+=2;
				while ((ny<0) && (l<lmax)) {		// ylm treated as zero and ignored if ny < 0
					for (int j=0; j<NWAY; ++j) {
						y0[j] = vall(al[1])*(cost[j]*y1[j]) + vall(al[0])*y0[j];
					}
					for (int j=0; j<NWAY; ++j) {
						y1[j] = vall(al[3])*(cost[j]*y0[j]) + vall(al[2])*y1[j];
					}
					l+=2;	al+=4;
					if (fabs(vlo(y0[NWAY-1])) > SHT_ACCURACY*SHT_SCALE_FACTOR + 1.0) {		// rescale when value is significant
						++ny;
						for (int j=0; j<NWAY; ++j) {
							y0[j] *= vall(1.0/SHT_SCALE_FACTOR);		y1[j] *= vall(1.0/SHT_SCALE_FACTOR);
						}
					}
				}
			  if (ny == 0) {
				while (l<lmax) {
					rnd qr = vall(creal(Ql[l]));		rnd qi = vall(cimag(Ql[l]));
					for (int j=0; j<NWAY; ++j) {
						qv[ (l-1)*2*NWAY + 2*j]   += y0[j] * qr;	// l
						qv[ (l-1)*2*NWAY + 2*j+1] += y0[j] * qi;
					}
					qr = vall(creal(Ql[l+1]));		qi = vall(cimag(Ql[l+1]));
					for (int j=0; j<NWAY; ++j) {
						qv[ (l)*2*NWAY + 2*j]   += y1[j] * qr;	// l+1
						qv[ (l)*2*NWAY + 2*j+1] += y1[j] * qi;
					}
					for (int j=0; j<NWAY; ++j) {
						y0[j] = vall(al[1])*(cost[j]*y1[j]) + vall(al[0])*y0[j];	// l+2
					}
					for (int j=0; j<NWAY; ++j) {
						y1[j] = vall(al[3])*(cost[j]*y0[j]) + vall(al[2])*y1[j];	// l+3
					}
					l+=2;	al+=4;
				}
				if (l==lmax) {
					rnd qr = vall(creal(Ql[l]));		rnd qi = vall(cimag(Ql[l]));
					for (int j=0; j<NWAY; ++j) {
						qv[ (l-1)*2*NWAY + 2*j]   += y0[j] * qr;	// l
						qv[ (l-1)*2*NWAY + 2*j+1] += y0[j] * qi;
					}
				}
			  }
			}
			// construct ring using symmetries + transpose...
			double signl = -1.0;
			double* qse = (double*) qve;
			double* qso = (double*) qvo;
			for (long l=1; l<=lmax; ++l) {
				for (int j=0; j<NWAY; j++) {
					for (int i=0; i<VSIZE2; i++) {
						double qre = qse[(l-1)*2*NWAY*VSIZE2 + 2*j*VSIZE2 + i];		// m even
						double qie = qse[(l-1)*2*NWAY*VSIZE2 + (2*j+1)*VSIZE2 + i];
						double qro = qso[(l-1)*2*NWAY*VSIZE2 + 2*j*VSIZE2 + i];		// m odd
						double qio = qso[(l-1)*2*NWAY*VSIZE2 + (2*j+1)*VSIZE2 + i];
						long ijk = (k+j)*VSIZE2 + i;
						qre *= st[ijk];			qro *= st[ijk];		// multiply real part by sin(theta)  [to get Ylm from Ylm/sin(theta)]
						// because qr and qi map on different parities with respect to the future Fourier tranform, we can add them !!
						// note that this may result in leak between even and odd m's if their amplitude is widely different.
						q0[ijk*lmax +(l-1)]              =  (qre + qro) - (qie + qio);
						q0[(ntheta+ijk)*lmax +(l-1)]     = ((qre + qro) + (qie + qio)) * signl;				// * (-1)^l
						q0[(2*ntheta-1-ijk)*lmax +(l-1)] =  (qre - qro) + (qie - qio);
						q0[(ntheta-1-ijk)*lmax +(l-1)]   = ((qre - qro) - (qie - qio)) * signl;				// (-1)^(l-m)
					}
				}
				signl *= -1.0;
			}
			k += NWAY;
		} while (k<nk);
	VFREE(qve);
#undef NWAY

//	tik1 = getticks();

	// perform FFT
	cplx* q = (cplx*) q0;		// in-place FFT
	fftw_execute_dft_r2c(shtns->fft_rot, q0, q);

//	tik2 = getticks();

	const int nphi = 2*ntheta;
	double ydyl[lmax+1];
	long m=0;		long lm=1;		// start at l=1,m=0
	long l;
		legendre_sphPlm_deriv_array_equ(shtns, lmax, m, ydyl+m);
		for (l=1; l<lmax; l+=2) {
			Rlm[lm] =  -creal(q[m*lmax +(l-1)])/(ydyl[l]*nphi);
			Rlm[lm+1] =  creal(q[m*lmax +l])/(ydyl[l+1]*nphi);
			lm+=2;
		}
		if (l==lmax) {
			Rlm[lm] =  -creal(q[m*lmax +(l-1)])/(ydyl[l]*nphi);
			lm+=1;
		}
	dphi1 += M_PI/nphi;	// shift rotation angle by angle of first synthesis latitude.
	for (m=1; m<=lmax; ++m) {
		legendre_sphPlm_deriv_array_equ(shtns, lmax, m, ydyl+m);
		cplx eimdp = (cos(m*dphi1) - I*sin(m*dphi1))/nphi;
		for (l=m; l<lmax; l+=2) {
			Rlm[lm] =  eimdp*q[m*lmax +(l-1)]*(1./ydyl[l]);
			Rlm[lm+1] =  eimdp*q[m*lmax +l]*(-1./ydyl[l+1]);
			lm+=2;
		}
		if (l==lmax) {
			Rlm[lm] =  eimdp*q[m*lmax +(l-1)]*(1./ydyl[l]);
			lm++;
		}
	}
	VFREE(q0);

//	tik3 = getticks();
//	printf("    tick ratio : %.3f  %.3f  %.3f\n", elapsed(tik1,tik0)/elapsed(tik3,tik0), elapsed(tik2,tik1)/elapsed(tik3,tik0), elapsed(tik3,tik2)/elapsed(tik3,tik0));

}


/// \addtogroup rotation
//@{

/// rotate Qlm by 90 degrees around X axis and store the result in Rlm.
/// shtns->mres MUST be 1, and lmax=mmax.
void SH_Xrotate90(shtns_cfg shtns, cplx *Qlm, cplx *Rlm)
{
	int lmax= shtns->lmax;
	if ((shtns->mres != 1) || (shtns->mmax < lmax)) shtns_runerr("truncature makes rotation not closed.");

	if (lmax == 1) {
		Rlm[0] = Qlm[0];	// l=0 is invariant.
		int l=1;													// rotation matrix for rotX(90), l=1 : m=[0, 1r, 1i]
			double q0 = creal(Qlm[LiM(shtns, l, 0)]);
			Rlm[LiM(shtns, l, 0)] = sqrt(2.0) * cimag(Qlm[LiM(shtns, l, 1)]);			//[m=0]     0        0    sqrt(2)
			Rlm[LiM(shtns, l ,1)] = creal(Qlm[LiM(shtns, l, 1)]) - I*(sqrt(0.5)*q0);	//[m=1r]    0        1      0
		return;																			//[m=1i] -sqrt(2)/2  0      0
	}

	SH_rotK90(shtns, Qlm, Rlm, 0.0,  -M_PI/2);
}

/// rotate Qlm by 90 degrees around Y axis and store the result in Rlm.
/// shtns->mres MUST be 1, and lmax=mmax.
void SH_Yrotate90(shtns_cfg shtns, cplx *Qlm, cplx *Rlm)
{
	int lmax= shtns->lmax;
	if ((shtns->mres != 1) || (shtns->mmax < lmax)) shtns_runerr("truncature makes rotation not closed.");

	if (lmax == 1) {
		Rlm[0] = Qlm[0];	// l=0 is invariant.
		int l=1;											// rotation matrix for rotY(90), l=1 : m=[0, 1r, 1i]
			double q0 = creal(Qlm[LiM(shtns, l, 0)]);									//[m=0]       0       0    sqrt(2)
			Rlm[LiM(shtns, l, 0)] = sqrt(2.0) * creal(Qlm[LiM(shtns, l, 1)]);			//[m=1r] -sqrt(2)/2   0      0
			Rlm[LiM(shtns, l ,1)] = I*cimag(Qlm[LiM(shtns, l, 1)]) - sqrt(0.5) * q0;	//[m=1i]      0       0      1
		return;
	}

	SH_rotK90(shtns, Qlm, Rlm, -M_PI/2, 0.0);
}

/// rotate Qlm around Y axis by arbitrary angle, using composition of rotations. Store the result in Rlm.
void SH_Yrotate(shtns_cfg shtns, cplx *Qlm, double alpha, cplx *Rlm)
{
	if ((shtns->mres != 1) || (shtns->mmax < shtns->lmax)) shtns_runerr("truncature makes rotation not closed.");

	SH_rotK90(shtns, Qlm, Rlm, 0.0, M_PI/2 + alpha);	// Zrotate(pi/2) + Yrotate90 + Zrotate(pi+alpha)
	SH_rotK90(shtns, Rlm, Rlm, 0.0, M_PI/2);			// Yrotate90 + Zrotate(pi/2)
}

//@}



/** \addtogroup operators Special operators
 * Apply special operators in spectral space: multiplication by cos(theta), sin(theta).d/dtheta.
*/
//@{


/// \internal generates the cos(theta) matrix up to lmax+1
/// \param mx : an array of 2*NLM double that will be filled with the matrix coefficients.
/// xq[lm] = mx[2*lm-1] * q[lm-1] + mx[2*lm] * q[lm+1];			[note the shift in indices compared to the public functions]
static void mul_ct_matrix_shifted(shtns_cfg shtns, double* mx)
{
	long int im,l,lm;
	double a_1;

	if (SHT_NORM == sht_schmidt) {
		lm=0;
		for (im=0; im<=MMAX; im++) {
			double* al = alm_im(shtns,im);
			long int m=im*MRES;
			a_1 = 1.0 / al[1];
			l=m;
			while(++l <= LMAX) {
				al+=2;				
				mx[2*lm+1] = a_1;
				a_1 = 1.0 / al[1];
				mx[2*lm] = -a_1*al[0];        // = -al[2*(lm+1)] / al[2*(lm+1)+1];
				lm++;
			}
			if (l == LMAX+1) {	// the last one needs to be computed (used in vector to scalar transform)
				mx[2*lm+1] = a_1;
				mx[2*lm] = sqrt((l+m)*(l-m))/(2*l+1);		// LMAX+1
				lm++;
			}
		}
	} else {
		lm=0;
		for (im=0; im<=MMAX; im++) {
			double* al = alm_im(shtns, im);
			l=im*MRES;
			while(++l <= LMAX+1) {	// compute coeff up to LMAX+1, it fits into the 2*NLM bloc, and is needed for vector <> scalar conversions.
				a_1 = 1.0 / al[1];
				mx[2*lm] = a_1;		// specific to orthonormal.
				mx[2*lm+1] = a_1;
				lm++;	al+=2;
			}
		}
	}
}

static void st_dt_matrix_shifted(shtns_cfg shtns, double* mx)
{
	mul_ct_matrix_shifted(shtns, mx);
	for (int lm=0; lm<NLM; lm++) {
		mx[2*lm]   *= -(shtns->li[lm] + 2);		// coeff (l+1)
		mx[2*lm+1] *=   shtns->li[lm];			// coeff (l-1)
	}
}


/// fill mx with the coefficients for multiplication by cos(theta)
/// \param mx : an array of 2*NLM double that will be filled with the matrix coefficients.
/// xq[lm] = mx[2*lm] * q[lm-1] + mx[2*lm+1] * q[lm+1];
void mul_ct_matrix(shtns_cfg shtns, double* mx)
{
	long int im,l,lm;
	double a_1;
	
	mx[0] = 0.0;
	mul_ct_matrix_shifted(shtns, mx+1);			// shift indices
	for (int im=1; im<=MMAX; im++) {				// remove the coeff for lmax+1 (for backward compatibility)
		int lm = LiM(shtns, im*MRES, im);
		mx[2*lm-1] = 0.0;		mx[2*lm] = 0.0;
	}
	mx[NLM-1] = 0.0;
}

/// fill mx with the coefficients of operator sin(theta).d/dtheta
/// \param mx : an array of 2*NLM double that will be filled with the matrix coefficients.
/// stdq[lm] = mx[2*lm] * q[lm-1] + mx[2*lm+1] * q[lm+1];
void st_dt_matrix(shtns_cfg shtns, double* mx)
{
	mul_ct_matrix(shtns, mx);
	for (int lm=0; lm<NLM; lm++) {
		mx[2*lm]   *=   shtns->li[lm] - 1;		// coeff (l-1)
		mx[2*lm+1] *= -(shtns->li[lm] + 2);		// coeff (l+1)
	}
}

/// Multiplication of Qlm by a matrix involving l+1 and l-1 only.
/// The result is stored in Rlm, which MUST be different from Qlm.
/// mx is an array of 2*NLM values as returned by \ref mul_ct_matrix or \ref st_dt_matrix
/// compute: Rlm[lm] = mx[2*lm] * Qlm[lm-1] + mx[2*lm+1] * Qlm[lm+1];
void SH_mul_mx(shtns_cfg shtns, double* mx, cplx *Qlm, cplx *Rlm)
{
	long int nlmlim, lm;
	v2d* vq = (v2d*) Qlm;
	v2d* vr = (v2d*) Rlm;
	nlmlim = NLM-1;
	lm = 0;
		s2d mxu = vdup(mx[1]);
		vr[0] = mxu*vq[1];
	for (lm=1; lm<nlmlim; lm++) {
		s2d mxl = vdup(mx[2*lm]);		s2d mxu = vdup(mx[2*lm+1]);
		vr[lm] = mxl*vq[lm-1] + mxu*vq[lm+1];
	}
	lm = nlmlim;
		s2d mxl = vdup(mx[2*lm]);
		vr[lm] = mxl*vq[lm-1];
}

//@}

// truncation at LMAX and MMAX
#define LTR LMAX
#define MTR MMAX

/** \addtogroup local Local and partial evaluation of SH fields.
 * These do only require a call to \ref shtns_create, but not to \ref shtns_set_grid.
 * These functions are not optimized and can be relatively slow, but they provide good
 * reference implemenation for the transforms.
*/
//@{

/// Evaluate scalar SH representation \b Qlm at physical point defined by \b cost = cos(theta) and \b phi
double SH_to_point(shtns_cfg shtns, cplx *Qlm, double cost, double phi)
{
	double yl[LMAX+1];
	double vr0, vr1;
	long int l,m,im;

	vr0 = 0.0;		vr1 = 0.0;
	m=0;	im=0;
		legendre_sphPlm_array(shtns, LTR, im, cost, &yl[m]);
		for (l=m; l<LTR; l+=2) {
			vr0 += yl[l] * creal( Qlm[l] );
			vr1 += yl[l+1] * creal( Qlm[l+1] );
		}
		if (l==LTR) {
			vr0 += yl[l] * creal( Qlm[l] );
		}
		vr0 += vr1;
	if (MTR>0) {
		im = 1;  do {
			m = im*MRES;
			legendre_sphPlm_array(shtns, LTR, im, cost, &yl[m]);
			v2d* Ql = (v2d*) &Qlm[LiM(shtns, 0,im)];	// virtual pointer for l=0 and im
			v2d vrm0 = vdup(0.0);		v2d vrm1 = vdup(0.0);
			for (l=m; l<LTR; l+=2) {
				vrm0 += vdup(yl[l]) * Ql[l];
				vrm1 += vdup(yl[l+1]) * Ql[l+1];
			}
			cplx eimp = 2.*(cos(m*phi) + I*sin(m*phi));		// we need something accurate here.
			vrm0 += vrm1;
			if (l==LTR) {
				vrm0 += vdup(yl[l]) * Ql[l];
			}
			vr0 += vcplx_real(vrm0)*creal(eimp) - vcplx_imag(vrm0)*cimag(eimp);
		} while(++im <= MTR);
	}
	return vr0;
}

void SH_to_grad_point(shtns_cfg shtns, cplx *DrSlm, cplx *Slm, double cost, double phi,
					   double *gr, double *gt, double *gp)
{
	double yl[LMAX+1];
	double dtyl[LMAX+1];
	double vtt, vpp, vr0, vrm;
	long int l,m,im;

	const double sint = sqrt((1.-cost)*(1.+cost));
	vtt = 0.;  vpp = 0.;  vr0 = 0.;  vrm = 0.;
	m=0;	im=0;
		legendre_sphPlm_deriv_array(shtns, LTR, im, cost, sint, &yl[m], &dtyl[m]);
		for (l=m; l<=LTR; ++l) {
			vr0 += yl[l] * creal( DrSlm[l] );
			vtt += dtyl[l] * creal( Slm[l] );
		}
	if (MTR>0) {
		im=1;  do {
			m = im*MRES;
			legendre_sphPlm_deriv_array(shtns, LTR, im, cost, sint, &yl[m], &dtyl[m]);
			cplx eimp = 2.*(cos(m*phi) + I*sin(m*phi));
			cplx imeimp = eimp*m*I;
			l = LiM(shtns, 0,im);
			v2d* Ql = (v2d*) &DrSlm[l];		v2d* Sl = (v2d*) &Slm[l];
			v2d qm = vdup(0.0);
			v2d dsdt = vdup(0.0);		v2d dsdp = vdup(0.0);
			for (l=m; l<=LTR; ++l) {
				qm += vdup(yl[l]) * Ql[l];
				dsdt += vdup(dtyl[l]) * Sl[l];
				dsdp += vdup(yl[l]) * Sl[l];
			}
			vrm += vcplx_real(qm)*creal(eimp) - vcplx_imag(qm)*cimag(eimp);			// dS/dr
			vtt += vcplx_real(dsdt)*creal(eimp) - vcplx_imag(dsdt)*cimag(eimp);		// dS/dt
			vpp += vcplx_real(dsdp)*creal(imeimp) - vcplx_imag(dsdp)*cimag(imeimp);	// + I.m/sint *S
		} while (++im <= MTR);
		vr0 += vrm*sint;
	}
	*gr = vr0;	// Gr = dS/dr
	*gt = vtt;	// Gt = dS/dt
	*gp = vpp;	// Gp = I.m/sint *S
}

/// Evaluate vector SH representation \b Qlm at physical point defined by \b cost = cos(theta) and \b phi
void SHqst_to_point(shtns_cfg shtns, cplx *Qlm, cplx *Slm, cplx *Tlm, double cost, double phi,
					   double *vr, double *vt, double *vp)
{
	double yl[LMAX+1];
	double dtyl[LMAX+1];
	double vtt, vpp, vr0, vrm;
	long int l,m,im;

	const double sint = sqrt((1.-cost)*(1.+cost));
	vtt = 0.;  vpp = 0.;  vr0 = 0.;  vrm = 0.;
	m=0;	im=0;
		legendre_sphPlm_deriv_array(shtns, LTR, im, cost, sint, &yl[m], &dtyl[m]);
		for (l=m; l<=LTR; ++l) {
			vr0 += yl[l] * creal( Qlm[l] );
			vtt += dtyl[l] * creal( Slm[l] );
			vpp -= dtyl[l] * creal( Tlm[l] );
		}
	if (MTR>0) {
		im=1;  do {
			m = im*MRES;
			legendre_sphPlm_deriv_array(shtns, LTR, im, cost, sint, &yl[m], &dtyl[m]);
			cplx eimp = 2.*(cos(m*phi) + I*sin(m*phi));
			cplx imeimp = eimp*m*I;
			l = LiM(shtns, 0,im);
			v2d* Ql = (v2d*) &Qlm[l];	v2d* Sl = (v2d*) &Slm[l];	v2d* Tl = (v2d*) &Tlm[l];
			v2d qm = vdup(0.0);
			v2d dsdt = vdup(0.0);		v2d dtdt = vdup(0.0);
			v2d dsdp = vdup(0.0);		v2d dtdp = vdup(0.0);
			for (l=m; l<=LTR; ++l) {
				qm += vdup(yl[l]) * Ql[l];
				dsdt += vdup(dtyl[l]) * Sl[l];
				dtdt += vdup(dtyl[l]) * Tl[l];
				dsdp += vdup(yl[l]) * Sl[l];
				dtdp += vdup(yl[l]) * Tl[l];
			}
			vrm += vcplx_real(qm)*creal(eimp) - vcplx_imag(qm)*cimag(eimp);
			vtt += (vcplx_real(dtdp)*creal(imeimp) - vcplx_imag(dtdp)*cimag(imeimp))	// + I.m/sint *T
					+ (vcplx_real(dsdt)*creal(eimp) - vcplx_imag(dsdt)*cimag(eimp));	// + dS/dt
			vpp += (vcplx_real(dsdp)*creal(imeimp) - vcplx_imag(dsdp)*cimag(imeimp))	// + I.m/sint *S
					- (vcplx_real(dtdt)*creal(eimp) - vcplx_imag(dtdt)*cimag(eimp));	// - dT/dt
		} while (++im <= MTR);
		vr0 += vrm*sint;
	}
	*vr = vr0;
	*vt = vtt;	// Bt = I.m/sint *T  + dS/dt
	*vp = vpp;	// Bp = I.m/sint *S  - dT/dt
}
//@}
	
#undef LTR
#undef MTR


/*
	SYNTHESIS AT A GIVEN LATITUDE
	(does not require a previous call to shtns_set_grid)
*/

/// synthesis at a given latitude, on nphi equispaced longitude points.
/// vr, vt, and vp arrays must have nphi doubles allocated.
/// It does not require a previous call to shtns_set_grid, but it is NOT thread-safe, 
/// unless called with a different shtns_cfg
/// \ingroup local
void SHqst_to_lat(shtns_cfg shtns, cplx *Qlm, cplx *Slm, cplx *Tlm, double cost,
					double *vr, double *vt, double *vp, int nphi, int ltr, int mtr)
{
	cplx vst, vtt, vsp, vtp, vrr;
	cplx *vrc, *vtc, *vpc;
	double* ylm_lat;
	double* dylm_lat;
	double st_lat;

	if (ltr > LMAX) ltr=LMAX;
	if (mtr > MMAX) mtr=MMAX;
	if (mtr*MRES > ltr) mtr=ltr/MRES;
	if (mtr*2*MRES >= nphi) mtr = (nphi-1)/(2*MRES);

	ylm_lat = shtns->ylm_lat;
	if (ylm_lat == NULL) {		// alloc memory for Legendre functions ?
		ylm_lat = (double *) malloc(sizeof(double)* 2*NLM);
		shtns->ylm_lat = ylm_lat;
	}
	dylm_lat = ylm_lat + NLM;

	st_lat = sqrt((1.-cost)*(1.+cost));	// sin(theta)
	if (cost != shtns->ct_lat) {		// compute Legendre functions ?
		shtns->ct_lat = cost;
		for (int m=0,j=0; m<=mtr; ++m) {
			legendre_sphPlm_deriv_array(shtns, ltr, m, cost, st_lat, &ylm_lat[j], &dylm_lat[j]);
			j += LMAX -m*MRES +1;
		}
	}

	vrc = (cplx*) fftw_malloc(sizeof(double) * 3*(nphi+2));
	vtc = vrc + (nphi/2+1);
	vpc = vtc + (nphi/2+1);

	if (nphi != shtns->nphi_lat) {		// compute FFTW plan ?
		if (shtns->ifft_lat) fftw_destroy_plan(shtns->ifft_lat);
		#ifdef OMP_FFTW
			fftw_plan_with_nthreads(1);
		#endif
		shtns->ifft_lat = fftw_plan_dft_c2r_1d(nphi, vrc, vr, FFTW_ESTIMATE);
		shtns->nphi_lat = nphi;
	}

	for (int m = 0; m<nphi/2+1; ++m) {	// init with zeros
		vrc[m] = 0.0;	vtc[m] = 0.0;	vpc[m] = 0.0;
	}
	long j=0;
	int m=0;
		vrr=0;	vtt=0;	vst=0;
		for(int l=m; l<=ltr; ++l, ++j) {
			vrr += ylm_lat[j] * creal(Qlm[j]);
			vst += dylm_lat[j] * creal(Slm[j]);
			vtt += dylm_lat[j] * creal(Tlm[j]);
		}
		j += (LMAX-ltr);
		vrc[m] = vrr;
		vtc[m] =  vst;	// Vt =   dS/dt
		vpc[m] = -vtt;	// Vp = - dT/dt
	for (int m=MRES; m<=mtr*MRES; m+=MRES) {
		vrr=0;	vtt=0;	vst=0;	vsp=0;	vtp=0;
		for(int l=m; l<=ltr; ++l, ++j) {
			vrr += ylm_lat[j] * Qlm[j];
			vst += dylm_lat[j] * Slm[j];
			vtt += dylm_lat[j] * Tlm[j];
			vsp += ylm_lat[j] * Slm[j];
			vtp += ylm_lat[j] * Tlm[j];
		}
		j+=(LMAX-ltr);
		vrc[m] = vrr*st_lat;
		vtc[m] = I*m*vtp + vst;	// Vt = I.m/sint *T  + dS/dt
		vpc[m] = I*m*vsp - vtt;	// Vp = I.m/sint *S  - dT/dt
	}
	fftw_execute_dft_c2r(shtns->ifft_lat,vrc,vr);
	fftw_execute_dft_c2r(shtns->ifft_lat,vtc,vt);
	fftw_execute_dft_c2r(shtns->ifft_lat,vpc,vp);
	fftw_free(vrc);
}

/// synthesis at a given latitude, on nphi equispaced longitude points.
/// vr arrays must have nphi doubles allocated.
/// It does not require a previous call to shtns_set_grid, but it is NOT thread-safe,
/// unless called with a different shtns_cfg
/// \ingroup local
void SH_to_lat(shtns_cfg shtns, cplx *Qlm, double cost,
					double *vr, int nphi, int ltr, int mtr)
{
	cplx vrr;
	cplx *vrc;
	double* ylm_lat;
	double* dylm_lat;
	double st_lat;

	if (ltr > LMAX) ltr=LMAX;
	if (mtr > MMAX) mtr=MMAX;
	if (mtr*MRES > ltr) mtr=ltr/MRES;
	if (mtr*2*MRES >= nphi) mtr = (nphi-1)/(2*MRES);

	ylm_lat = shtns->ylm_lat;
	if (ylm_lat == NULL) {
		ylm_lat = (double *) malloc(sizeof(double)* 2*NLM);
		shtns->ylm_lat = ylm_lat;
	}
	dylm_lat = ylm_lat + NLM;

	st_lat = sqrt((1.-cost)*(1.+cost));	// sin(theta)
	if (cost != shtns->ct_lat) {
		shtns->ct_lat = cost;
		for (int m=0,j=0; m<=mtr; ++m) {
			legendre_sphPlm_deriv_array(shtns, ltr, m, cost, st_lat, &ylm_lat[j], &dylm_lat[j]);
			j += LMAX -m*MRES +1;
		}
	}

	vrc = (cplx*) fftw_malloc(sizeof(double) * (nphi+2));

	if (nphi != shtns->nphi_lat) {
		if (shtns->ifft_lat) fftw_destroy_plan(shtns->ifft_lat);
		#ifdef OMP_FFTW
			fftw_plan_with_nthreads(1);
		#endif
		shtns->ifft_lat = fftw_plan_dft_c2r_1d(nphi, vrc, vr, FFTW_ESTIMATE);
		shtns->nphi_lat = nphi;
	}

	for (int m = 0; m<nphi/2+1; ++m) {	// init with zeros
		vrc[m] = 0.0;
	}
	long j=0;
	int m=0;
		vrr=0;
		for(int l=m; l<=ltr; ++l, ++j) {
			vrr += ylm_lat[j] * creal(Qlm[j]);
		}
		j += (LMAX-ltr);
		vrc[m] = vrr;
	for (int m=MRES; m<=mtr*MRES; m+=MRES) {
		vrr=0;
		for(int l=m; l<=ltr; ++l, ++j) {
			vrr += ylm_lat[j] * Qlm[j];
		}
		j+=(LMAX-ltr);
		vrc[m] = vrr*st_lat;
	}
	fftw_execute_dft_c2r(shtns->ifft_lat,vrc,vr);
	fftw_free(vrc);
}

// SPAT_CPLX transform indexing scheme:
// if (l<=MMAX) : l*(l+1) + m
// if (l>=MMAX) : l*(2*mmax+1) - mmax*mmax + m
///\internal
void SH_2real_to_cplx(shtns_cfg shtns, cplx* Rlm, cplx* Ilm, cplx* Zlm)
{
	// combine into complex coefficients
	unsigned ll = 0;
	unsigned lm = 0;
	for (unsigned l=0; l<=LMAX; l++) {
		ll += (l<=MMAX) ? 2*l : 2*MMAX+1;
		Zlm[ll] = creal(Rlm[lm]) + I*creal(Ilm[lm]);		// m=0
		lm++;
	}
	for (unsigned m=1; m<=MMAX; m++) {
		ll = (m-1)*m;
		for (unsigned l=m; l<=LMAX; l++) {
			ll += (l<=MMAX) ? 2*l : 2*MMAX+1;
			cplx rr = Rlm[lm];
			cplx ii = Ilm[lm];
			Zlm[ll+m] = rr + I*ii;			// m>0
			rr = conj(rr) + I*conj(ii);		// m<0, m even
			if (m&1) rr = -rr;				// m<0, m odd
			Zlm[ll-m] = rr;
			lm++;
		}
	}
}

///\internal
void SH_cplx_to_2real(shtns_cfg shtns, cplx* Zlm, cplx* Rlm, cplx* Ilm)
{
	// extract complex coefficients corresponding to real and imag
	unsigned ll = 0;
	unsigned lm = 0;
	for (unsigned l=0; l<=LMAX; l++) {
		ll += (l<=MMAX) ? 2*l : 2*MMAX+1;
		Rlm[lm] = creal(Zlm[ll]);		// m=0
		Ilm[lm] = cimag(Zlm[ll]);
		lm++;
	}
	double half_parity = 0.5;
	for (unsigned m=1; m<=MMAX; m++) {
		ll = (m-1)*m;
		half_parity = -half_parity;		// (-1)^m * 0.5
		for (unsigned l=m; l<=LMAX; l++) {
			ll += (l<=MMAX) ? 2*l : 2*MMAX+1;
			cplx b = Zlm[ll-m] * half_parity;		// (-1)^m for m negative.
			cplx a = Zlm[ll+m] * 0.5;
			Rlm[lm] = (conj(b) + a);		// real part
			Ilm[lm] = (conj(b) - a)*I;		// imag part
			lm++;
		}
	}
}

/// complex scalar transform.
/// in: complex spatial field.
/// out: alm[l*(l+1)+m] is the SH coefficients of order l and degree m (with -l <= m <= l)
/// for a total of (LMAX+1)^2 coefficients.
void spat_cplx_to_SH(shtns_cfg shtns, cplx *z, cplx *alm)
{
	long int nspat = shtns->nspat;
	double *re, *im;
	cplx *rlm, *ilm;

	if (MRES != 1) shtns_runerr("complex SH requires mres=1.");

	// alloc temporary fields
	re = (double*) VMALLOC( 2*(nspat + NLM*2)*sizeof(double) );
	im = re + nspat;
	rlm = (cplx*) (re + 2*nspat);
	ilm = rlm + NLM;

	// split z into real and imag parts.
	for (int k=0; k<nspat; k++) {
		re[k] = creal(z[k]);		im[k] = cimag(z[k]);
	}

	// perform two real transforms:
	spat_to_SH(shtns, re, rlm);
	spat_to_SH(shtns, im, ilm);

	// combine into complex coefficients
	SH_2real_to_cplx(shtns, rlm, ilm, alm);

	VFREE(re);
}

/// complex scalar transform.
/// in: alm[l*(l+1)+m] is the SH coefficients of order l and degree m (with -l <= m <= l)
/// for a total of (LMAX+1)^2 coefficients.
/// out: complex spatial field.
void SH_to_spat_cplx(shtns_cfg shtns, cplx *alm, cplx *z)
{
	long int nspat = shtns->nspat;
	double *re, *im;
	cplx *rlm, *ilm;

	if (MRES != 1) shtns_runerr("complex SH requires mres=1.");

	// alloc temporary fields
	re = (double*) VMALLOC( 2*(nspat + NLM*2)*sizeof(double) );
	im = re + nspat;
	rlm = (cplx*) (re + 2*nspat);
	ilm = rlm + NLM;

	// extract complex coefficients corresponding to real and imag
	SH_cplx_to_2real(shtns, alm, rlm, ilm);

	// perform two real transforms:
	SH_to_spat(shtns, rlm, re);
	SH_to_spat(shtns, ilm, im);

	// combine into z
	for (int k=0; k<nspat; k++)
		z[k] = re[k] + I*im[k];

	VFREE(re);
}

void SH_cplx_Xrotate90(shtns_cfg shtns, cplx *Qlm, cplx *Rlm)
{
	if (MRES != 1) shtns_runerr("complex SH requires mres=1.");

	// alloc temporary fields
	cplx* rlm = (cplx*) VMALLOC( NLM*2*sizeof(cplx) );
	cplx* ilm = rlm + NLM;

	// extract complex coefficients corresponding to real and imag
	SH_cplx_to_2real(shtns, Qlm, rlm, ilm);

	// perform two real rotations:
	SH_Xrotate90(shtns, rlm, rlm);
	SH_Xrotate90(shtns, ilm, ilm);

	// combine back into complex coefficients
	SH_2real_to_cplx(shtns, rlm, ilm, Rlm);

	VFREE(rlm);
}

void SH_cplx_Yrotate90(shtns_cfg shtns, cplx *Qlm, cplx *Rlm)
{
	if (MRES != 1) shtns_runerr("complex SH requires mres=1.");

	// alloc temporary fields
	cplx* rlm = (cplx*) VMALLOC( NLM*2*sizeof(cplx) );
	cplx* ilm = rlm + NLM;

	// extract complex coefficients corresponding to real and imag
	SH_cplx_to_2real(shtns, Qlm, rlm, ilm);

	// perform two real rotations:
	SH_Yrotate90(shtns, rlm, rlm);
	SH_Yrotate90(shtns, ilm, ilm);

	// combine back into complex coefficients
	SH_2real_to_cplx(shtns, rlm, ilm, Rlm);

	VFREE(rlm);
}

/// complex scalar rotation around Y
/// in: Qlm[l*(l+1)+m] is the SH coefficients of order l and degree m (with -l <= m <= l)
/// out: Qlm[l*(l+1)+m] is the rotated SH coefficients of order l and degree m (with -l <= m <= l)
void SH_cplx_Yrotate(shtns_cfg shtns, cplx *Qlm, double alpha, cplx *Rlm)
{
	if (MRES != 1) shtns_runerr("complex SH requires mres=1.");

	// alloc temporary fields
	cplx* rlm = (cplx*) VMALLOC( NLM*2*sizeof(cplx) );
	cplx* ilm = rlm + NLM;

	// extract complex coefficients corresponding to real and imag
	SH_cplx_to_2real(shtns, Qlm, rlm, ilm);

	// perform two real rotations:
	SH_Yrotate(shtns, rlm, alpha, rlm);
	SH_Yrotate(shtns, ilm, alpha, ilm);

	// combine back into complex coefficients
	SH_2real_to_cplx(shtns, rlm, ilm, Rlm);

	VFREE(rlm);
}

/// complex scalar rotation around Z
/// in: Qlm[l*(l+1)+m] is the SH coefficients of order l and degree m (with -l <= m <= l)
/// out: Qlm[l*(l+1)+m] is the rotated SH coefficients of order l and degree m (with -l <= m <= l)
void SH_cplx_Zrotate(shtns_cfg shtns, cplx *Qlm, double alpha, cplx *Rlm)
{
	if (MRES != 1) shtns_runerr("complex SH requires mres=1.");

	// alloc temporary fields
	cplx* rlm = (cplx*) VMALLOC( NLM*2*sizeof(cplx) );
	cplx* ilm = rlm + NLM;

	// extract complex coefficients corresponding to real and imag
	SH_cplx_to_2real(shtns, Qlm, rlm, ilm);

	// perform two real rotations:
	SH_Zrotate(shtns, rlm, alpha, rlm);
	SH_Zrotate(shtns, ilm, alpha, ilm);

	// combine back into complex coefficients
	SH_2real_to_cplx(shtns, rlm, ilm, Rlm);

	VFREE(rlm);
}

/*
/// Rotate a SH representation of complex field Qlm around the z-axis by angle alpha (in radians),
/// which is the same as rotating the reference frame by angle -alpha.
/// Result is stored in Rlm (which can be the same array as Qlm).
void SH_cplx_Zrotate(shtns_cfg shtns, cplx *Qlm, double alpha, cplx *Rlm)
{
	if (MRES != 1) shtns_runerr("complex SH requires mres=1.");

	cplx* eima = (cplx*) VMALLOC( (2*MMAX+1)*sizeof(cplx) );
	eima += MMAX;
	eima[0] = 1.0;
	for (int m=1; m<=MMAX; m++) {		// precompute the complex numbers
		double cma = cos(m*alpha);
		double sma = sin(m*alpha);
		eima[m]  = cma - I*sma;
		eima[-m] = cma + I*sma;
	}

	unsigned ll=0;
	for (unsigned l=0; l<=MMAX; l++) {
		for (int m=-l; m<=l; m++) {
			Rlm[ll] = Qlm[ll] * eima[m];
			ll++;
		}
	}
	for (unsigned l=MMAX+1; l<=LMAX; l++) {
		for (int m=-MMAX; m<=MMAX; m++) {
			Rlm[ll] = Qlm[ll] * eima[m];
			ll++;
		}
	}

	VFREE(eima);
}
*/

/*
void SH_to_spat_grad(shtns_cfg shtns, cplx *alm, double *gt, double *gp)
{
	double *mx;
	cplx *blm, *clm;
	
	blm = (cplx*) VMALLOC( 3*NLM*sizeof(cplx) );
	clm = blm + NLM;
	mx = (double*)(clm + NLM);

	st_dt_matrix(shtns, mx);
	SH_mul_mx(shtns, mx, alm, blm);
	int lm=0;
	for (int im=0; im<=MMAX; im++) {
		int m = im*MRES;
		for (int l=m; l<=LMAX; l++) {
			clm[lm] = alm[lm] * I*m;
			lm++;
		}
	}
	SH_to_spat(shtns, blm, gt);
	SH_to_spat(shtns, clm, gp);
	for (int ip=0; ip<NPHI; ip++) {
		for (int it=0; it<NLAT; it++) {
			gt[ip*NLAT+it] /= shtns->st[it];
			gp[ip*NLAT+it] /= shtns->st[it];
		}
	}
	VFREE(blm);
}
*/
