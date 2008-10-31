

#define FILE_VERSION 0
struct JobInfo {
	int version, lmax, mmax, mres, nlm, nr, irs, ire;	// 8 integers
	int BC, i2, i3, i4, i5, i6, i7, i8;			// 8 integers
	double Omega0, nu, eta, t, DeltaOmega, d6, d7, d8;		// 8 doubles
	char text[896];						// 896 chars => total 1024 bytes.
};


/// save spherical harmonic representation of Poloidal/Toroidal field PT.
void save_PolTor(char *fn, struct PolTor *PT, double time, int BC)
{
	long int ir;
	struct JobInfo fh;
	FILE *fp;

	ir = 0;
	while ((PT->P[ir] == NULL)&&(ir<NR)) ir++;	// find first populated shell.
	fh.irs = ir;
	while ((PT->P[ir] != NULL)&&(ir<NR)) ir++;	// find last populated shell.
	fh.ire = ir-1;

	fh.version = FILE_VERSION;	fh.BC = BC;
	fh.lmax = LMAX;		fh.mmax = MMAX;		fh.mres = MRES;		fh.nlm = NLM;
	fh.nr = NR;
	fh.Omega0 = Omega0;	fh.nu = nu;	fh.eta = eta;	fh.t = time;	fh.DeltaOmega = DeltaOmega;

	fp = fopen(fn,"w");
	fwrite(&fh, sizeof(fh), 1, fp);		// Header.
	fwrite(r, sizeof(double), NR, fp);	// grille radiale.

	for (ir= fh.irs; ir<= fh.ire; ir++) {
		fwrite(PT->P[ir], sizeof(complex double), NLM, fp);		// data
		fwrite(PT->T[ir], sizeof(complex double), NLM, fp);		// data
	}

	fclose(fp);
}

/// load spherical harmonic representation of Poloidal/Toroidal field PT.
void load_PolTor(char *fn, struct PolTor *PT, struct JobInfo *fh)
{
	long int ir,im,m,l,lm;
	long int imt,lmt;
	complex double **St, **Sp;
	FILE *fp;

	long int lim(long int l, long int im) {
		return im*(2*fh->lmax+3 -(im+2)*fh->mres)/2 + l;
	}
	
	fp = fopen(fn,"r");
	if (fp == NULL) runerr("[load] file not found !");

	fread(fh, sizeof(struct JobInfo), 1, fp);		// read Header
	printf("[load] from file '%s' : NR=%d, Lmax=%d, Mmax=%d, Mres=%d (Nlm=%d)\n",fn, fh->nr, fh->lmax, fh->mmax, fh->mres, fh->nlm);
	printf("       Omega0=%.3e, nu=%.3e, eta=%.3e, t=%.3e\n", fh->Omega0, fh->nu, fh->eta, fh->t );
	printf("       ir_start=%d, ir_end=%d, BC=%d\n", fh->irs, fh->ire, fh->BC);
//	*time = fh->t;	*BC = fh->BC;	*irs = fh->irs;	*ire = fh->ire;
	Omega0 = fh->Omega0;	DeltaOmega = fh->DeltaOmega;

	if (r == NULL) {
		r = (double *) malloc(fh->nr * sizeof(double));		// alloc radial grid (global scope)
		NR = fh->nr;
	}
	else if (fh->nr != NR) runerr("wrong NR : radial sizes must match.");

	fread(r, sizeof(double), fh->nr ,fp);			// load radial grid.
	init_Deriv_sph();					// init radial derivative matrices.

	if (PT->P == NULL)	{	// destination field not yet allocated ?
		alloc_DynamicField(PT, NULL, NULL, NULL, fh->irs, fh->ire);
	}

	Sp = (complex double **) malloc( fh->nr * sizeof(complex double *) );
	St = (complex double **) malloc( fh->nr * sizeof(complex double *) );

	for (ir=0; ir<fh->nr; ir++) {			// initialize at NULL pointer.
		Sp[ir] = NULL;	St[ir] = NULL;
	}
	for (ir= fh->irs; ir<= fh->ire; ir++) {
		Sp[ir] = (complex double *) malloc(fh->nlm * sizeof(complex double));	// alloc shell
		St[ir] = (complex double *) malloc(fh->nlm * sizeof(complex double));	// alloc shell
		fread(Sp[ir], sizeof(complex double), fh->nlm, fp);	// load data
		fread(St[ir], sizeof(complex double), fh->nlm, fp);	// load data
		for (im=0, lm=0, imt=0, lmt=0; im<=MMAX; im++) {
			m=im*MRES;
			while(m > imt*fh->mres) {		// skip data if not required.
				lmt+= fh->lmax +1 -imt*fh->mres;
				imt++;
			}
			if ((imt <= fh->mmax) && (m == imt*fh->mres)) {		// data present : copy
				for (l=m; l<= LMAX; l++) {
					if (l<=fh->lmax) {
						PT->P[ir][lm] = Sp[ir][lmt];
						PT->T[ir][lm] = St[ir][lmt];
						lmt++;
					} else {
						PT->P[ir][lm] = 0.0;
						PT->T[ir][lm] = 0.0;
					}
					lm++;
				}
				for (l=LMAX+1; l<fh->lmax; l++) {
					lmt++;
				}
			} else {	// data not present : set to zero
				for (l=m; l<= LMAX; l++) {
					PT->P[ir][lm] = 0.0;
					PT->T[ir][lm] = 0.0;
					lm++;
				}
			}
		}
	}

	for (ir=fh->nr-1; ir>=0; ir--) {	// free memory allocated for shells.
		if (St[ir] != NULL) { free(St[ir]);	free(Sp[ir]); }
	}
	free(St);	free(Sp);
	fclose(fp);
}
