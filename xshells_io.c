

#define FILE_VERSION 0
struct Header {
	int version, lmax, mmax, mres, nlm, nr, ir0 , i8;	// 8 integers
	double Omega0, nu, eta, t, d5, d6, d7, d8;		// 8 doubles
	char text[928];						// 928 chars => total 1024 bytes.
};

/// save spherical harmonic representation of scalar field S.
void save_field(char *fn, complex double **S, double time)
{
	long int ir;
	struct Header fh;
	FILE *fp;

	ir = 0;
	while ((S[ir] == NULL)&&(ir<NR)) ir++;	// find first populated shell.

	fh.version = FILE_VERSION;
	fh.lmax = LMAX;		fh.mmax = MMAX;		fh.mres = MRES;		fh.nlm = NLM;
	fh.nr = NR;	fh.ir0 = ir;
	fh.Omega0 = Omega0;	fh.nu = nu;	fh.eta = eta;	fh.t = time;

	fp = fopen(fn,"w");
	fwrite(&fh, sizeof(fh), 1, fp);		// Header.
	fwrite(r, sizeof(double), NR, fp);	// grille radiale.

	while( (ir<NR) && (S[ir] != NULL) )
		fwrite(S[ir], sizeof(complex double), NLM, fp);		// data

	fclose(fp);
}

/// load spherical harmonic representation of scalar field S.
void load_field(char *fn, complex double **S, double time)
{
	long int ir;
	struct Header fh;
	FILE *fp;
	
	fp = fopen(fn,"r");
	if (fp == NULL) runerr("[load] file not found !");

	fread(&fh, sizeof(fh), 1, fp);		// read Header
	printf("[load] from file '%s' : NR=%d, Lmax=%d, Mmax=%d, Mres=%d (Nlm=%d)\n",fn, fh.nr, fh.lmax, fh.mmax, fh.mres, fh.nlm);
	printf("       Omega0=%.3e, nu=%.3e, eta=%.3e, t=%.3e\n", fh.Omega0, fh.nu, fh.eta, fh.time );

	if (


	DXY = fh.dxy;
	JET_W0 = fh.W0;

	DKZ = fh.dkz;

	Ut = (struct CplxVect *) malloc((2*KMAX+1)*(2*KMAX+1)*(fh.Kz_max+1) * sizeof(struct CplxVect));
	fread(Ut, sizeof(struct CplxVect), (2*KMAX+1)*(2*KMAX+1)*(fh.Kz_max+1), fp);

	for (ix=0; ix<=2*KMAX; ix++) {
		for (iy=0; iy<=2*KMAX; iy++) {
			iz=0;
			id = (KzMAX+1) * (iy + ix*(2*KMAX+1));
			idt = (fh.Kz_max+1) * (iy + ix*(2*KMAX+1));
			while ((iz <= fh.Kz_max) && (iz <= KzMAX)) {
				U[id].u = Ut[idt].u;
				U[id].v = Ut[idt].v;
				U[id].w = Ut[idt].w;
				iz++;	idt++;	id++;
			}
			while (iz<=KzMAX) {
				U[id].u = 0.0;  U[id].v = 0.0;  U[id].w = 0.0;
				iz++;   id++;
			}
		}
	}

	free(Ut);	// free temp memory
	fclose(fp);
	
}