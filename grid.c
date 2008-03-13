
/*	Inversion de matrices x-Bandes. (x = 3,5,7)
	Les conditions aux limites sont prises en compte (ou pas !)
	=> pour l'évolution temporelle semi-implicite des potentiels et champs.
*/

struct TriDiag {
	double l, d ,u;
};

struct CinqDiag {
	double l2,l1, d ,u1,u2;
};

struct SeptDiag {
	double l3,l2,l1, d ,u1,u2,u3;
};


void runerr(const char * error_text)
{
	printf("*** Run-time error : %s\n",error_text);
	exit(1);
}

void TriDec(struct TriDiag *M, int n)
// decomposition d'une matrice tribande. Les divisions ont lieu dans cette etape.
// les elements b et c sont écrasés !
{
	double tmp;
	int j;

	if (M[0].d == 0.0) runerr("[TriDec] first pivot is zero.");
	
	M[0].d = 1.0/M[0].d;
	for (j=1; j<n; j++)	//Decomposition.
	{
		M[j-1].u = M[j-1].u * M[j-1].d;
		tmp = M[j].d - M[j].l * M[j-1].u;
		if (tmp == 0.0) runerr("[TriDec] zero pivot encountered.");	// Algorithm fails
		M[j].d = 1.0/tmp;
	}
}

void TriSolve(struct TriDiag *M, double *b, double *x, int n)
/* Solves M x = b, where b and x can be the same array.
   M is a TriDiag Matrix, previously decomposed by TriDec. */
{
	int j;

	x[0] = b[0] * M[0].d;
	for (j=1;j<n;j++)		// Forward substitution.
		x[j] = (b[j] - M[j].l * x[j-1]) * M[j].d;
	for (j=(n-2);j>=0;j--)		//Backsubstitution.
		x[j] -= M[j].u *x[j+1];
}

void cTriSolve(struct TriDiag *M, complex double *b, complex double *x, int n)
/* Solves M x = b, where b and x can be the same array.
   M is a TriDiag Matrix, previously decomposed by TriDec. */
{
	int j;

	x[0] = b[0] * M[0].d;
	for (j=1;j<n;j++)		// Forward substitution.
		x[j] = (b[j] - M[j].l * x[j-1]) * M[j].d;
	for (j=(n-2);j>=0;j--)		//Backsubstitution.
		x[j] -= M[j].u *x[j+1];
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
	int i;

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
// backward avec CL : x[n] != 0		(peut etre le cas pour le strain imposé m=2)
	i = n-1;
		x[i] -= M[i].u1 * x[i+1];
	for(i=n-2;i>=0;i--)	// tient compte des conditions limites !!!
		x[i] -= M[i].u1 * x[i+1] + M[i].u2 * x[i+2];
}


void SeptDec(struct SeptDiag *M, int n)
{
	double alp, bet, gam;
	int i;
	
	if (M[0].d == 0.0) runerr("[SeptDec] first pivot is zero.");
	
	i=1;
		gam = M[i].l1 / M[i-1].d;
		M[i].d -= gam * M[i-1].u1;
		M[i].u1 -= gam * M[i-1].u2;
		M[i].u2 -= gam * M[i-1].u3;
	i=2;
		if (M[i-1].d == 0.0) runerr("[SeptDec] second pivot turned out to be zero.");
		bet = M[i].l2 / M[i-2].d;
		M[i].l1 -= bet * M[i-2].u1;
		gam = M[i].l1 / M[i-1].d;
		M[i].d -= gam * M[i-1].u1 + bet * M[i-2].u2;
		M[i].u1 -= gam * M[i-1].u2 + bet * M[i-2].u3;
		M[i].u2 -= gam * M[i-1].u3;
	for(i=3;i<n;i++) {
		if (M[i-1].d == 0.0) runerr("[SeptDec] zero pivot encountered.");
		alp = M[i].l3 / M[i-3].d;
		M[i].l2 -= alp * M[i-3].u1;
		bet = M[i].l2 / M[i-2].d;
		M[i].l1 -= bet * M[i-2].u1 + alp * M[i-3].u2;
		gam = M[i].l1 / M[i-1].d;
		M[i].d -= gam * M[i-1].u1 + bet * M[i-2].u2 + alp * M[i-3].u3;
		M[i].u1 -= gam * M[i-1].u2 + bet * M[i-2].u3;
		M[i].u2 -= gam * M[i-1].u3;
	}

	for(i=0;i<n;i++) {
		M[i].d = 1.0 / M[i].d;
		M[i].l3 *= M[i].d;
		M[i].l2 *= M[i].d;
		M[i].l1 *= M[i].d;
		M[i].u1 *= M[i].d;
		M[i].u2 *= M[i].d;
		M[i].u3 *= M[i].d;
	}
}

void SeptSolve(struct SeptDiag *M, complex double *b, complex double *x, int n)
/* Solves M x = b, where b and x can be the same array.
   M is a SeptDiag Matrix, previously decomposed by SeptDec.
*/
{
	int i;

// forward (suppose CL : x[-1] = 0)
	i=0;
		x[i] = M[i].d * b[i];
	i=1;
		x[i] = M[i].d * b[i] - M[i].l1 * x[i-1];
	i=2;
		x[i] = M[i].d * b[i] - M[i].l1 * x[i-1] - M[i].l2 * x[i-2];
	for(i=3;i<n;i++)
		x[i] = M[i].d * b[i] - M[i].l1 * x[i-1] - M[i].l2 * x[i-2] - M[i].l3 * x[i-3];

// backward (suppose CL : x[n] = 0) ce qui est le cas pour les champs poloidaux consideres.
	i = n-2;
		x[i] -= M[i].u1 * x[i+1];
	i = n-3;
		x[i] -= M[i].u1 * x[i+1] + M[i].u2 * x[i+2];
	for(i=n-4;i>=0;i--)
		x[i] -= M[i].u1 * x[i+1] + M[i].u2 * x[i+2] + M[i].u3 * x[i+3];

}



/*	Génération des Grilles et des opérateurs de dérivation spatiale.
*/

// discretisation spatiale
double *r, *r_1, *r_2, *dr;
struct TriDiag *Gr, *Lr, *Wr, *D2r;		// Laplacien, Wr = 1/r Gr( r .), D2r = d2/dr2

void Reg_Grid(double dx)		// Grille régulière.
{
	int i;

	r[0] = 0.0;
	for(i=1;i<Nr;i++)
		r[i] = r[i-1] + dx;

	printf("[GRID:Reg] Nr=%d rmax=%f, dx=%f\n",Nr,r[Nr-1],dx);
}

/*
void Exp_Grid(double dx0, double alpha)		// Grille exponentielle.
{
	double dx;
	int i;

	r[0] = 0.0;
	dx = dx0;
	for (i=1;i<Nr;i++)
	{
		dx = dx*alpha;
		r[i] = r[i-1] + dx;
	}

	printf("[GRID:Exp] Nr=%d alpha=%f rmax=%f, dxmin=%f dxmax=%f\n",Nr,alpha,r[Nr-1],dx0,dx);
}

void Mixed_Grid(int Nreg, double rreg_max, double dxmax_mul)	// Génère une grille réguliere au centre, puis exponentielle.
{
	double dx, dxmax, alpha;
	int i, Nirr;

	if (Nreg > Nr)	Nreg = Nr;
// Grille : réguilière dans le coeur ...
	r[0] = 0.0;
	dx = rreg_max/(Nreg-1);
	for(i=1;i<Nreg;i++)
	{
		r[i] = r[0] + dx*i;
	}

// ... puis decroit exponentiellement, d'un facteur alpha.
	Nirr = (Nr-1)-(Nreg+1);
	dxmax = dxmax_mul * log(Nirr);
	alpha = exp(log(dxmax/dx)/(double)Nirr);
	for (i=Nreg;i<Nr;i++)
	{
		dx = dx*alpha;
		r[i] = r[i-1] + dx;
	}

	dxmax = dx;
	dx = r[1]-r[0];
	printf("[GRID:Mixed] Nr=%d alpha=%f rmax=%f, dxmin=%f dxmax=%f\n",Nr,alpha,r[Nr-1],dx,dxmax);
}
*/

void init_Deriv()
{
	double t, rm, rp, grm, grp;
	int i;

	// creation des variables globales :
	r_1 = (double *) malloc(Nr * sizeof(double));
	r_2 = (double *) malloc(Nr * sizeof(double));
	dr = (double *) malloc(Nr * sizeof(double));
	
	Gr = (struct TriDiag *) malloc(Nr * sizeof(struct TriDiag));
	Wr = (struct TriDiag *) malloc(Nr * sizeof(struct TriDiag));
	Lr = (struct TriDiag *) malloc(Nr * sizeof(struct TriDiag));
	D2r = (struct TriDiag *) malloc(Nr * sizeof(struct TriDiag));
	
	for(i=1;i<Nr;i++)
	{
		r_1[i] = 1./r[i];
		r_2[i] = 1./(r[i]*r[i]);
	}
	i = 0;
		r_1[i] = 0.0;	r_2[i] = 0.0;

	// increments
	for(i=0;i<Nr-1;i++)
		dr[i] = r[i+1]-r[i];
	dr[Nr-1] = 0.0;

	// gradient et Laplacien
	for (i=1;i<Nr-1;i++)
	{
		t = 1.0/((dr[i-1]+dr[i])*dr[i]*dr[i-1]);
		Gr[i].l = -dr[i]*dr[i]*t;
		Gr[i].d = (dr[i]*dr[i] - dr[i-1]*dr[i-1])*t;	// =0 en grille reguliere.
		Gr[i].u = dr[i-1]*dr[i-1]*t;

		Wr[i].l = r_1[i]* Gr[i].l * r[i-1];	// Wr = 1/r * d/dr(r .), pour la vorticité m=0 !
		Wr[i].d = Gr[i].d;
		Wr[i].u = r_1[i]* Gr[i].u * r[i+1];

		D2r[i].l =  2.0*dr[i]*t;
		D2r[i].d =  -2.0*(dr[i-1]+dr[i])*t;
		D2r[i].u =  2.0*dr[i-1]*t;

		Lr[i].l =  2.0*dr[i]*t            + r_1[i]*Gr[i].l;
		Lr[i].d =  -2.0*(dr[i-1]+dr[i])*t + r_1[i]*Gr[i].d;
		Lr[i].u =  2.0*dr[i-1]*t          + r_1[i]*Gr[i].u;
	}
	// les extremites douvent etre déterminées par les CL au cas par cas.
	i = 0;
		Gr[i].l = 0.0; Gr[i].d = 0.0; Gr[i].u = 0.0;
		Wr[i].l = 0.0;	Wr[i].d = 0.0;	Wr[i].u = 0.0;
		D2r[i].l = 0.0; D2r[i].d = 0.0; D2r[i].u = 0.0;
		Lr[i].l = 0.0;	Lr[i].d = 0.0;	Lr[i].u = 0.0;
	i = Nr-1;
		Gr[i].l = 0.0; Gr[i].d = 0.0; Gr[i].u = 0.0;
		Wr[i].l = 0.0;	Wr[i].d = 0.0;	Wr[i].u = 0.0;
		D2r[i].l = 0.0; D2r[i].d = 0.0; D2r[i].u = 0.0;
		Lr[i].l = 0.0;	Lr[i].d = 0.0;	Lr[i].u = 0.0;

}

// genere la grille radiale, et les matrices de derivees spatiales
void init_Rad(double dr)
{
	r = (double *) malloc(Nr * sizeof(double));	// allocation de la grille radiale

	Reg_Grid(dr);

	init_Deriv();
}

