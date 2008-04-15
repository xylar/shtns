///////////////////////////////////////////////
// grid : Radial finite difference setup
//////////////////////////////////////////////

/*	Inversion de matrices x-Bandes. (x = 3,5)
	Les conditions aux limites sont prises en compte (ou pas !)
	=> pour l'�volution temporelle semi-implicite des potentiels et champs.
        => Les Matrices ont une taille NR*(LMAX+1), la decomposition se fait sur NR.
*/

#define RL(i,l) ( i*(LMAX+1) + l )
#define RLM(i,lm) ( i*NLM + lm )

// LM_LOOP : loop over all (l,im) and perform "action"  : l,im,lm are defined. (but NOT m )
#define LM_LOOP( action ) for (im=0, lm=0; im<=MMAX; im++) { for (l=im*MRES; l<=LMAX; l++, lm++) { action } }

struct TriDiag {
	double l,d,u;		// lower, diagonal, upper.
};

struct TriDiagL {
	double l;		// lower
	double d[LMAX+1];	// diagonal (depends on l)
	double u;		// upper
};

struct CinqDiag {
	double l2,l1, d ,u1,u2;
};



// discretisation spatiale
double *r, *r_1, *r_2, *dr;		// r = rayon; r_1 = 1/r; r_2 = 1/(r*r); dr[i] = r[i+1]-r[i];
struct TriDiag *Gr, *Wr, *D2r, *Lr;	// Gr = Gradient; Wr= 1/r Gr( r .); D2r = d2/dr2; Lr = Laplacien radial scalaire.

/*
void runerr(const char * error_text)
{
	printf("*** Run-time error : %s\n",error_text);
	exit(1);
}
*/

// Multiplication d'un vecteur complexe par une matrice Tribande dependant de l
// y = M.x    (y and x MUST be different)
inline void cTriMul(struct TriDiagL *M, complex double **x, complex double **y, int istart, int iend)
{
	long int j,l,im,lm;

	j=istart;
		for (im=0, lm=0; im<=MMAX; im++) {
			for (l=im*MRES; l<=LMAX; l++, lm++)
				y[j][lm] = M[j].d[l] * x[j][lm] + M[j].u * x[j+1][lm];
		}
	for (j=istart+1; j<iend; j++) {
		for (im=0, lm=0; im<=MMAX; im++) {
			for (l=im*MRES; l<=LMAX; l++, lm++)
				y[j][lm] = M[j].l * x[j-1][lm] + M[j].d[l] * x[j][lm] + M[j].u * x[j+1][lm];
		}
	}
	j = iend;
		for (im=0, lm=0; im<=MMAX; im++) {
			for (l=im*MRES; l<=LMAX; l++, lm++)
				y[j][lm] = M[j].l * x[j-1][lm] + M[j].d[l] * x[j][lm];
		}
}

// Multiplication d'un vecteur complexe par une matrice Tribande dependant de l, ajoute au resultat.
// y += M.x    (y and x MUST be different)
inline void cTriMulAdd(struct TriDiagL *M, complex double **x, complex double **y, int istart, int iend)
{
	long int j,l,im,lm;

	j=istart;
		for (im=0, lm=0; im<=MMAX; im++) {
			for (l=im*MRES; l<=LMAX; l++, lm++)
				y[j][lm] += M[j].d[l] * x[j][lm] + M[j].u * x[j+1][lm];
		}
	for (j=istart+1; j<iend; j++) {
		for (im=0, lm=0; im<=MMAX; im++) {
			for (l=im*MRES; l<=LMAX; l++, lm++)
				y[j][lm] += M[j].l * x[j-1][lm] + M[j].d[l] * x[j][lm] + M[j].u * x[j+1][lm];
		}
	}
	j = iend;
		for (im=0, lm=0; im<=MMAX; im++) {
			for (l=im*MRES; l<=LMAX; l++, lm++)
				y[j][lm] += M[j].l * x[j-1][lm] + M[j].d[l] * x[j][lm];
		}
}

// decomposition PARTIELLE d'une matrice tribande. Les divisions ont lieu dans cette etape.
// seul l'element diagonal est ecras� !!!
void TriDec(struct TriDiagL *M, int istart, int iend)
{
	double tmp;
	int j,l;

	j = istart;
		tmp = 1.0;
		for (l=0;l<=LMAX;l++)
			tmp *= M[j].d[l];
		if (tmp == 0.0) runerr("[TriDec] first pivot is zero.");

		for (l=0;l<=LMAX;l++)
			M[j].d[l] = 1.0/M[j].d[l];
	while(j < iend) {	//Decomposition.
		j++;
		for (l=0;l<=LMAX;l++) {
			tmp = M[j].d[l] - M[j].l * M[j-1].u * M[j-1].d[l];
			if (tmp == 0.0) runerr("[TriDec] zero pivot encountered.");	// Algorithm fails
			M[j].d[l] = 1.0/tmp;
		}
	}
}

// Solve M x = b, where b and x can be the same array.
//   M is a TriDiagL Matrix, previously decomposed by TriDec.
//   b and x are full Ylm : (NR*LMMAX), [b and x can point to the same data].
inline void cTriSolve(struct TriDiagL *M, complex double **b, complex double **x, int istart, int iend)
{
	long int j,l,im,lm;

	j = istart;
	if (x != b) {		// don't copy if b and x are the same vector.
		for (lm=0; lm<NLM; lm++)  x[j][lm] = b[j][lm];
	}
	while (j < iend-1) {	// Forward substitution.
		j++;
		for (im=0, lm=0; im<=MMAX; im++) {
			for (l=im*MRES; l<=LMAX; l++, lm++)
				x[j][lm] = b[j][lm] - (M[j].l * M[j-1].d[l] * x[j-1][lm]);
		}
	}	// j = iend-1;
	j = iend;
		for (im=0, lm=0; im<=MMAX; im++) {
			for (l=im*MRES; l<=LMAX; l++, lm++)
				x[j][lm] = (b[j][lm] - (M[j].l * M[j-1].d[l] * x[j-1][lm])) * M[j].d[l];
		}
	while (j > istart) {	// Back-substitution.
		j--;
		for (im=0, lm=0; im<=MMAX; im++) {
			for (l=im*MRES; l<=LMAX; l++, lm++)
				x[j][lm] = ( x[j][lm] - M[j].u *x[j+1][lm] ) * M[j].d[l];
		}
	}	// j=istart
}

void CinqDec(struct CinqDiag *M, int n)
{
	double alp, bet;
	int i;
	
	if (M[0].d == 0.0) runerr("[CinqDec] first pivot is zero.");
	
	i = 1;
		bet = M[i].l1 / M[i-1].d;
		M[i].d -= bet * M[i-1].u1;
		M[i].u1 -= bet * M[i-1].u2;
	for(i=2;i<n;i++) {
		if (M[i-1].d == 0.0) runerr("[CinqDec] zero pivot encountered.");
		alp = M[i].l2 / M[i-2].d;
		M[i].l1 -= alp * M[i-2].u1;
		bet = M[i].l1 / M[i-1].d;
		M[i].d -= bet * M[i-1].u1 + alp * M[i-2].u2;
		M[i].u1 -= bet * M[i-1].u2;
	}

	for(i=0;i<n;i++) {
		M[i].d = 1.0 / M[i].d;
		M[i].l2 *= M[i].d;
		M[i].l1 *= M[i].d;
		M[i].u1 *= M[i].d;
		M[i].u2 *= M[i].d;
	}
}

void CinqSolve(struct CinqDiag *M, complex double *b, complex double *x, int n)
/* Solves M x = b, where b and x can be the same array.
   M is a CinqDiag Matrix, previously decomposed by CinqDec.
*/
{
	long int i;

// forward avec condition limite : x[-1] != 0;
/*	i=0;
		x[i] = M[i].d * b[i] - M[i].l1*x[i-1];
	for(i=1;i<n;i++)
		x[i] = M[i].d * b[i] - M[i].l1 * x[i-1] - M[i].l2 * x[i-2];
*/
// forward (suppose CL : x[-1] = 0)
	i=0;
		x[i] = M[i].d * b[i];
	i=1;
		x[i] = M[i].d * b[i] - M[i].l1 * x[i-1];
	for(i=2;i<n;i++)
		x[i] = M[i].d * b[i] - M[i].l1 * x[i-1] - M[i].l2 * x[i-2];

// backward (suppose CL : x[n] = 0)
/*	i = n-2;
		x[i] -= M[i].u1 * x[i+1];
	for(i=n-3;i>=0;i--)
		x[i] -= M[i].u1 * x[i+1] + M[i].u2 * x[i+2];
*/
// backward avec CL : x[n] != 0		(peut etre le cas pour le strain impos� m=2)
	i = n-1;
		x[i] -= M[i].u1 * x[i+1];
	for(i=n-2;i>=0;i--)	// tient compte des conditions limites !!!
		x[i] -= M[i].u1 * x[i+1] + M[i].u2 * x[i+2];
}




/*	G�n�ration des Grilles et des op�rateurs de d�rivation spatiale.
*/


void Reg_Grid(double rmin, double rmax)		// Grille r�guli�re.
{
	double dx;
	int i;

	dx = (rmax - rmin)/(NR-1);
	for(i=0;i<NR;i++)
		r[i] = ((rmax-rmin)*i)/(NR-1);

	printf("[GRID:Reg] NR=%d, r=[%f, %f], dr=%f\n",NR,r[0],r[NR-1],dx);
}

/*
void Exp_Grid(double dx0, double alpha)		// Grille exponentielle.
{
	double dx;
	int i;

	r[0] = 0.0;
	dx = dx0;
	for (i=1;i<NR;i++)
	{
		dx = dx*alpha;
		r[i] = r[i-1] + dx;
	}

	printf("[GRID:Exp] NR=%d alpha=%f rmax=%f, dxmin=%f dxmax=%f\n",NR,alpha,r[NR-1],dx0,dx);
}

void Mixed_Grid(int Nreg, double rreg_max, double dxmax_mul)	// G�n�re une grille r�guliere au centre, puis exponentielle.
{
	double dx, dxmax, alpha;
	int i, Nirr;

	if (Nreg > NR)	Nreg = NR;
// Grille : r�guili�re dans le coeur ...
	r[0] = 0.0;
	dx = rreg_max/(Nreg-1);
	for(i=1;i<Nreg;i++)
	{
		r[i] = r[0] + dx*i;
	}

// ... puis decroit exponentiellement, d'un facteur alpha.
	Nirr = (NR-1)-(Nreg+1);
	dxmax = dxmax_mul * log(Nirr);
	alpha = exp(log(dxmax/dx)/(double)Nirr);
	for (i=Nreg;i<NR;i++)
	{
		dx = dx*alpha;
		r[i] = r[i-1] + dx;
	}

	dxmax = dx;
	dx = r[1]-r[0];
	printf("[GRID:Mixed] NR=%d alpha=%f rmax=%f, dxmin=%f dxmax=%f\n",NR,alpha,r[NR-1],dx,dxmax);
}
*/

// Gradients, Laplaciens .... en spherique.
void init_Deriv_sph()
{
	double t, rm, rp, grm, grp;
	int i;

	// creation des variables globales :
	r_1 = (double *) malloc(NR * sizeof(double));
	r_2 = (double *) malloc(NR * sizeof(double));
	dr = (double *) malloc(NR * sizeof(double));
	
	Gr = (struct TriDiag *) malloc(NR * sizeof(struct TriDiag));
	Wr = (struct TriDiag *) malloc(NR * sizeof(struct TriDiag));
	Lr = (struct TriDiag *) malloc(NR * sizeof(struct TriDiag));
	D2r = (struct TriDiag *) malloc(NR * sizeof(struct TriDiag));
	
	for(i=1;i<NR;i++)
	{
		r_1[i] = 1./r[i];
		r_2[i] = 1./(r[i]*r[i]);
	}
	i = 0;
		r_1[i] = 0.0;	r_2[i] = 0.0;

	// increments
	for(i=0;i<NR-1;i++)
		dr[i] = r[i+1]-r[i];
	dr[NR-1] = 0.0;

	// gradient et Laplacien
	for (i=1;i<NR-1;i++)
	{
		t = 1.0/((dr[i-1]+dr[i])*dr[i]*dr[i-1]);
		Gr[i].l = -dr[i]*dr[i]*t;
		Gr[i].d = (dr[i]*dr[i] - dr[i-1]*dr[i-1])*t;	// =0 en grille reguliere.
		Gr[i].u = dr[i-1]*dr[i-1]*t;

		Wr[i].l = r_1[i]* Gr[i].l * r[i-1];	// Wr = 1/r * d/dr(r .), pour la vorticit� m=0 !
		Wr[i].d =         Gr[i].d;
		Wr[i].u = r_1[i]* Gr[i].u * r[i+1];

		D2r[i].l =  2.0*dr[i]*t;
		D2r[i].d =  -2.0*(dr[i-1]+dr[i])*t;
		D2r[i].u =  2.0*dr[i-1]*t;
/*
		Lr[i].l = r_1[i]* D2r[i].l * r[i-1];	// Laplacien radial : 1/r . d2/dr2(r .)
		Lr[i].d =         D2r[i].d;
		Lr[i].u = r_1[i]* D2r[i].u * r[i+1];
*/
		Lr[i].l = D2r[i].l + 2.0*r_1[i]*Gr[i].l;	// Laplacien radial : d2/dr2 + 2/r.d/r
		Lr[i].d = D2r[i].d + 2.0*r_1[i]*Gr[i].d;
		Lr[i].u = D2r[i].u + 2.0*r_1[i]*Gr[i].u;

	}
	// les extremites douvent etre d�termin�es par les CL au cas par cas.
	i = 0;
		Gr[i].l = 0.0; Gr[i].d = 0.0; Gr[i].u = 0.0;
		Wr[i].l = 0.0;	Wr[i].d = 0.0;	Wr[i].u = 0.0;
		D2r[i].l = 0.0; D2r[i].d = 0.0; D2r[i].u = 0.0;
		Lr[i].l = 0.0;	Lr[i].d = 0.0;	Lr[i].u = 0.0;
	i = NR-1;
		Gr[i].l = 0.0; Gr[i].d = 0.0; Gr[i].u = 0.0;
		Wr[i].l = 0.0;	Wr[i].d = 0.0;	Wr[i].u = 0.0;
		D2r[i].l = 0.0; D2r[i].d = 0.0; D2r[i].u = 0.0;
		Lr[i].l = 0.0;	Lr[i].d = 0.0;	Lr[i].u = 0.0;
}

/*
// Gradients, Laplaciens .... en cylindrique.
void init_Deriv_cyl()
{
	double t, rm, rp, grm, grp;
	int i;

	// creation des variables globales :
	r_1 = (double *) malloc(NR * sizeof(double));
	r_2 = (double *) malloc(NR * sizeof(double));
	dr = (double *) malloc(NR * sizeof(double));
	
	Gr = (struct TriDiag *) malloc(NR * sizeof(struct TriDiag));
	Wr = (struct TriDiag *) malloc(NR * sizeof(struct TriDiag));
	Lr = (struct TriDiag *) malloc(NR * sizeof(struct TriDiag));
	D2r = (struct TriDiag *) malloc(NR * sizeof(struct TriDiag));
	
	for(i=1;i<NR;i++)
	{
		r_1[i] = 1./r[i];
		r_2[i] = 1./(r[i]*r[i]);
	}
	i = 0;
		r_1[i] = 0.0;	r_2[i] = 0.0;

	// increments
	for(i=0;i<NR-1;i++)
		dr[i] = r[i+1]-r[i];
	dr[NR-1] = 0.0;

	// gradient et Laplacien
	for (i=1;i<NR-1;i++)
	{
		t = 1.0/((dr[i-1]+dr[i])*dr[i]*dr[i-1]);
		Gr[i].l = -dr[i]*dr[i]*t;
		Gr[i].d = (dr[i]*dr[i] - dr[i-1]*dr[i-1])*t;	// =0 en grille reguliere.
		Gr[i].u = dr[i-1]*dr[i-1]*t;

		Wr[i].l = r_1[i]* Gr[i].l * r[i-1];	// Wr = 1/r * d/dr(r .), pour la vorticit� m=0 !
		Wr[i].d = Gr[i].d;
		Wr[i].u = r_1[i]* Gr[i].u * r[i+1];

		D2r[i].l =  2.0*dr[i]*t;
		D2r[i].d =  -2.0*(dr[i-1]+dr[i])*t;
		D2r[i].u =  2.0*dr[i-1]*t;

		Lr[i].l =  2.0*dr[i]*t            + r_1[i]*Gr[i].l;
		Lr[i].d =  -2.0*(dr[i-1]+dr[i])*t + r_1[i]*Gr[i].d;
		Lr[i].u =  2.0*dr[i-1]*t          + r_1[i]*Gr[i].u;
	}
	// les extremites douvent etre d�termin�es par les CL au cas par cas.
	i = 0;
		Gr[i].l = 0.0; Gr[i].d = 0.0; Gr[i].u = 0.0;
		Wr[i].l = 0.0;	Wr[i].d = 0.0;	Wr[i].u = 0.0;
		D2r[i].l = 0.0; D2r[i].d = 0.0; D2r[i].u = 0.0;
		Lr[i].l = 0.0;	Lr[i].d = 0.0;	Lr[i].u = 0.0;
	i = NR-1;
		Gr[i].l = 0.0; Gr[i].d = 0.0; Gr[i].u = 0.0;
		Wr[i].l = 0.0;	Wr[i].d = 0.0;	Wr[i].u = 0.0;
		D2r[i].l = 0.0; D2r[i].d = 0.0; D2r[i].u = 0.0;
		Lr[i].l = 0.0;	Lr[i].d = 0.0;	Lr[i].u = 0.0;

}
*/

// genere la grille radiale, et les matrices de derivees spatiales
void init_rad_sph(double rmin, double rmax)
{
	r = (double *) malloc(NR * sizeof(double));	// allocation de la grille radiale
	Reg_Grid(rmin, rmax);
	init_Deriv_sph();
}

