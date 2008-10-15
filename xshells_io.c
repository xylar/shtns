
void write_vect(char *fn, double *vec, int N)
{
	FILE *fp; 
	int i;
	
	fp = fopen(fn,"w");
	for (i=0;i<N;i++) {
		fprintf(fp,"%.6g ",vec[i]);
	}
	fclose(fp);
}

void write_mx(char *fn, double *mx, int N1, int N2)
{
	FILE *fp;
	int i,j;
	
	fp = fopen(fn,"w");
	for (i=0;i<N1;i++) {
		for(j=0;j<N2;j++) {
			fprintf(fp,"%.6g ",mx[i*N2+j]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}

void write_HS(char *fn, complex double **HS)
{
	FILE *fp;
	int ir,lm;
	
	fp = fopen(fn,"w");
	for (ir=0;ir<NR;ir++) {
		if (HS[ir] != NULL) {
			for(lm=0;lm<NLM;lm++) {
				fprintf(fp,"%.6g %.6g  ",creal(HS[ir][lm]),cimag(HS[ir][lm]));
			}
			fprintf(fp,"\n");
		}
	}
	fclose(fp);
}

void write_slice(char *fn, double **v, int im)
{
	FILE *fp;
	int i,j;

	fp = fopen(fn,"w");
		fprintf(fp,"0 ");			// first row = radius
		for(j=0;j<NLAT/2;j++) {
			fprintf(fp,"%.6g ",ct[j]);	// first line = cos(theta)
		}
		for(j=1;j<=NLAT/2;j++) {
			fprintf(fp,"-%.6g ",ct[NLAT/2-j]);	// first line = cos(theta)
		}
	for (i=0;i<NR;i++) {
		if (v[i] != NULL) {
			fprintf(fp,"\n%.6g ",r[i]);		// first row = radius
			for(j=0;j<NLAT;j++) {
				fprintf(fp,"%.6g ",v[i][im*NLAT + j]);		// data
			}
		}
	}
	fclose(fp);
}





#define FILE_VERSION 0
struct Header {
	int version, lmax, mmax, mres, nlm, nr, irs, ire;	// 8 integers
	double Omega0, nu, eta, t, d5, d6, d7, d8;		// 8 doubles
	char text[928];						// 928 chars => total 1024 bytes.
};



/// save spherical harmonic representation of Poloidal/Toroidal field PT.
void save_PolTor(char *fn, struct PolTor *PT, double time)
{
	long int ir;
	struct Header fh;
	FILE *fp;

	ir = 0;
	while ((PT->P[ir] == NULL)&&(ir<NR)) ir++;	// find first populated shell.
	fh.irs = ir;
	while ((PT->P[ir] != NULL)&&(ir<NR)) ir++;	// find last populated shell.
	fh.ire = ir-1;

	fh.version = FILE_VERSION;
	fh.lmax = LMAX;		fh.mmax = MMAX;		fh.mres = MRES;		fh.nlm = NLM;
	fh.nr = NR;
	fh.Omega0 = Omega0;	fh.nu = nu;	fh.eta = eta;	fh.t = time;

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
void load_PolTor(char *fn, struct PolTor *PT, double time)
{
	long int ir,im,m,l,lm;
	long int imt,lmt;
	complex double **St, **Sp;
	struct Header fh;
	FILE *fp;

	long int lim(long int l, long int im) {
		return im*(2*fh.lmax+3 -(im+2)*fh.mres)/2 + l;
	}
	
	fp = fopen(fn,"r");
	if (fp == NULL) runerr("[load] file not found !");

	fread(&fh, sizeof(fh), 1, fp);		// read Header
	printf("[load] from file '%s' : NR=%d, Lmax=%d, Mmax=%d, Mres=%d (Nlm=%d)\n",fn, fh.nr, fh.lmax, fh.mmax, fh.mres, fh.nlm);
	printf("       Omega0=%.3e, nu=%.3e, eta=%.3e, t=%.3e\n", fh.Omega0, fh.nu, fh.eta, fh.t );

	if (r == NULL) {
		r = (double *) malloc(fh.nr * sizeof(double));		// alloc radial grid (global scope)
		NR = fh.nr;
	}
	else if (fh.nr != NR) runerr("wrong NR : radial sizes must match.");

	fread(r, sizeof(double), fh.nr ,fp);			// load radial grid.
	init_Deriv_sph();					// init radial derivative matrices.

	Sp = (complex double **) malloc( fh.nr * sizeof(complex double *) );
	St = (complex double **) malloc( fh.nr * sizeof(complex double *) );

	for (ir=0; ir<fh.nr; ir++) St[ir] = NULL;		// initialize at NULL pointer.
	for (ir= fh.irs; ir<= fh.ire; ir++) {
		Sp[ir] = (complex double *) malloc(fh.nlm * sizeof(complex double));	// alloc shell
		St[ir] = (complex double *) malloc(fh.nlm * sizeof(complex double));	// alloc shell
		fread(Sp[ir], sizeof(complex double), fh.nlm, fp);	// load data
		fread(St[ir], sizeof(complex double), fh.nlm, fp);	// load data
		for (im=0, lm=0, imt=0, lmt=0; im<=MMAX; im++) {
			m=im*MRES;
			while(m > imt*fh.mres) {		// skip data if not required.
				lmt+= fh.lmax +1 -imt*fh.mres;
				imt++;
			}
			if ((imt <= fh.mmax) && (m == imt*fh.mres)) {		// data present : copy
				for (l=m; l<= LMAX; l++) {
					if (l<=fh.lmax) {
						PT->P[ir][lm] = Sp[ir][lmt];
						PT->T[ir][lm] = St[ir][lmt];
						lmt++;
					} else {
						PT->P[ir][lm] = 0.0;
						PT->T[ir][lm] = 0.0;
					}
					lm++;
				}
				for (l=LMAX+1; l<fh.lmax; l++) {
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

	for (ir=fh.nr-1; ir>=0; ir--) {	// free memory allocated for shells.
		if (St[ir] != NULL) { free(St[ir]);	free(Sp[ir]); }
	}
	free(St);	free(Sp);
	fclose(fp);
}
