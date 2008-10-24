///////////////////////////////////////////////
// grid : Radial finite difference setup
//////////////////////////////////////////////

/*	x-banded matrix inversion. (x = 3,5)
	Bondary condition are taken into acount (or not !)
	=> for temporal semi-implicit evolution of fields.
        => matrices have size NR*(LMAX+1), decomposition is done on NR.
*/

#ifndef NR
long int NR;
#endif

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

struct PentaDiag {
	double l2,l1, d ,u1,u2;
};

// discretisation spatiale
double* r = NULL;		// NULL at init time.
double *r_1, *r_2, *dr;		// r = rayon; r_1 = 1/r; r_2 = 1/(r*r); dr[i] = r[i+1]-r[i];
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
		LM_LOOP(  y[j][lm] = M[j].d[l] * x[j][lm] + M[j].u * x[j+1][lm];  )
	for (j=istart+1; j<iend; j++) {
		LM_LOOP(  y[j][lm] = M[j].l * x[j-1][lm] + M[j].d[l] * x[j][lm] + M[j].u * x[j+1][lm];  )
	}
	j = iend;
		LM_LOOP(  y[j][lm] = M[j].l * x[j-1][lm] + M[j].d[l] * x[j][lm];  )
}

// x has elements istart-1 and iend+1 set (if zero => same as cTriSolve)
inline void cTriMulBC(struct TriDiagL *M, complex double **x, complex double **y, int istart, int iend)
{
	long int j,l,im,lm;

	for (j=istart; j<=iend; j++) {
		for (im=0, lm=0; im<=MMAX; im++) {
			for (l=im*MRES; l<=LMAX; l++, lm++)
				y[j][lm] = M[j].l * x[j-1][lm] + M[j].d[l] * x[j][lm] + M[j].u * x[j+1][lm];
		}
	}
}

/*
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
*/

// decomposition PARTIELLE d'une matrice tribande. Les divisions ont lieu dans cette etape.
// seul l'element diagonal est ecrase !!!
void TriDec(struct TriDiagL *M, int istart, int iend)
{
	double tmp;
	long int j,l;

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

// x has elements istart-1 and iend+1 set (if zero => same as cTriSolve)
inline void cTriSolveBC(struct TriDiagL *M, complex double **b, complex double **x, int istart, int iend)
{
	long int j,l,im,lm;

	j = istart;	// handle lower boundary : b(j) -> b(j) - Ml(j).x(j-1)
		for (lm=0; lm<NLM; lm++)  x[j][lm] = b[j][lm] - M[j].l * x[j-1][lm];

	while (j < iend-1) {	// Forward substitution.
		j++;
		for (im=0, lm=0; im<=MMAX; im++) {
			for (l=im*MRES; l<=LMAX; l++, lm++)
				x[j][lm] = b[j][lm] - (M[j].l * M[j-1].d[l] * x[j-1][lm]);
		}
	}	// j = iend-1;
	j = iend;	// handle upper boundary : b(j) -> b(j) - Mu(j).x(j+1)
		for (im=0, lm=0; im<=MMAX; im++) {
			for (l=im*MRES; l<=LMAX; l++, lm++)
				x[j][lm] = (b[j][lm] - M[j].u*x[j+1][lm] - (M[j].l * M[j-1].d[l] * x[j-1][lm])) * M[j].d[l];
		}
	while (j > istart) {	// Back-substitution.
		j--;
		for (im=0, lm=0; im<=MMAX; im++) {
			for (l=im*MRES; l<=LMAX; l++, lm++)
				x[j][lm] = ( x[j][lm] - M[j].u *x[j+1][lm] ) * M[j].d[l];
		}
	}	// j=istart
}


// Multiply complex vector by a Penta-diagonal matrix (l-dependant)
// y = M.x    (y and x MUST be different)
inline void cPentaMul(struct PentaDiag **M, complex double **x, complex double **y, int istart, int iend)
{
	long int j,l,im,lm;

	j=istart;
		for (im=0, lm=0; im<=MMAX; im++) {
			for (l=im*MRES; l<=LMAX; l++, lm++)
				y[j][lm] = M[j][l].d * x[j][lm] + M[j][l].u1 * x[j+1][lm] + M[j][l].u2 * x[j+2][lm];
		}
	j=istart+1;
		for (im=0, lm=0; im<=MMAX; im++) {
			for (l=im*MRES; l<=LMAX; l++, lm++)
				y[j][lm] = M[j][l].l1 * x[j-1][lm] + M[j][l].d * x[j][lm] + M[j][l].u1 * x[j+1][lm] + M[j][l].u2 * x[j+2][lm];
		}
	for (j=istart+2; j<iend-1; j++) {
		for (im=0, lm=0; im<=MMAX; im++) {
			for (l=im*MRES; l<=LMAX; l++, lm++)
				y[j][lm] = M[j][l].l2 * x[j-2][lm] + M[j][l].l1 * x[j-1][lm] + M[j][l].d * x[j][lm] + M[j][l].u1 * x[j+1][lm] + M[j][l].u2 * x[j+2][lm];
		}
	}
	j = iend-1;
		for (im=0, lm=0; im<=MMAX; im++) {
			for (l=im*MRES; l<=LMAX; l++, lm++)
				y[j][lm] = M[j][l].l2 * x[j-2][lm] + M[j][l].l1 * x[j-1][lm] + M[j][l].d * x[j][lm] + M[j][l].u1 * x[j+1][lm];
		}
	j = iend;
		for (im=0, lm=0; im<=MMAX; im++) {
			for (l=im*MRES; l<=LMAX; l++, lm++)
				y[j][lm] = M[j][l].l2 * x[j-2][lm] + M[j][l].l1 * x[j-1][lm] + M[j][l].d * x[j][lm];
		}
}

void PentaDec(struct PentaDiag **M, long int istart, long int iend)
{
	double alp, bet;
	long int i,l;
/*
        u(:,2,4)=u(:,2,4)-u(:,2,2)*u(:,1,5)/u(:,1,3)
        u(:,2,3)=u(:,2,3)-u(:,2,2)*u(:,1,4)/u(:,1,3)

        do i=3,NS-1
           u(:,i,2)=u(:,i,2)-u(:,i,1)*u(:,i-2,4)/u(:,i-2,3)
           u(:,i,4)=u(:,i,4)-u(:,i,2)*u(:,i-1,5)/u(:,i-1,3)
           u(:,i,3)=u(:,i,3)-u(:,i,2)*u(:,i-1,4)/u(:,i-1,3)-u(:,i,1)*
     c          u(:,i-2,5)/u(:,i-2,3)
        enddo
c! invert d , scale u1,u2,l1,l2
        do i=1,NS-1
           u(:,i,3)=1.0/u(:,i,3)
           u(:,i,4)=u(:,i,4)*u(:,i,3)
           u(:,i,5)=u(:,i,5)*u(:,i,3)
           u(:,i,2)=u(:,i,2)*u(:,i,3)
           u(:,i,1)=u(:,i,1)*u(:,i,3)
        enddo
*/
	i = istart;
		bet = 1.0;
		for (l=0;l<=LMAX;l++)
			bet *= M[i][l].d;
		if (bet == 0.0) runerr("[PentaDec] first pivot is zero.");

	i = istart+1;
		for (l=0;l<=LMAX;l++) {
			bet = M[i][l].l1 / M[i-1][l].d;
			M[i][l].d -= bet * M[i-1][l].u1;
			M[i][l].u1 -= bet * M[i-1][l].u2;
		}
	for(i=istart+2; i<=iend; i++) {
		for (l=0;l<=LMAX;l++) {
			if (M[i-1][l].d == 0.0) runerr("[PentaDec] zero pivot encountered.");
			alp = M[i][l].l2 / M[i-2][l].d;
			M[i][l].l1 -= alp * M[i-2][l].u1;
			bet = M[i][l].l1 / M[i-1][l].d;
			M[i][l].d -= bet * M[i-1][l].u1 + alp * M[i-2][l].u2;
			M[i][l].u1 -= bet * M[i-1][l].u2;
		}
	}

	for(i=istart; i<=iend; i++) {
		for (l=0;l<=LMAX;l++) {
			M[i][l].d = 1.0 / M[i][l].d;
			M[i][l].l2 *= M[i][l].d;
			M[i][l].l1 *= M[i][l].d;
			M[i][l].u1 *= M[i][l].d;
			M[i][l].u2 *= M[i][l].d;
		}
	}
}


inline void cPentaSolve(struct PentaDiag **M, complex double **b, complex double **x, long int istart, long int iend)
/* Solves M x = b, where b and x can be the same array. (NR*NLM)
   M is a CinqDiag Matrix, previously decomposed by CinqDec. (NR*LMAX)
*/
{
	long int i,l,im,lm;

// forward
	i = istart;
		for (im=0, lm=0; im<=MMAX; im++) {
			for (l=im*MRES; l<=LMAX; l++, lm++)
				x[i][lm] = M[i][l].d * b[i][lm];
		}
	i = istart+1;
		for (im=0, lm=0; im<=MMAX; im++) {
			for (l=im*MRES; l<=LMAX; l++, lm++)
				x[i][lm] = M[i][l].d * b[i][lm] - M[i][l].l1 * x[i-1][lm];
		}
	for(i=istart+2; i<=iend; i++) {
		for (im=0, lm=0; im<=MMAX; im++) {
			for (l=im*MRES; l<=LMAX; l++, lm++)
				x[i][lm] = M[i][l].d * b[i][lm] - M[i][l].l1 * x[i-1][lm] - M[i][l].l2 * x[i-2][lm];
		}
	}
// backward
	i = iend-1;
		for (im=0, lm=0; im<=MMAX; im++) {
			for (l=im*MRES; l<=LMAX; l++, lm++)
				x[i][lm] -= M[i][l].u1 * x[i+1][lm];
		}
	for(i= iend-2; i>=istart; i--) {
		for (im=0, lm=0; im<=MMAX; im++) {
			for (l=im*MRES; l<=LMAX; l++, lm++)
				x[i][lm] -= M[i][l].u1 * x[i+1][lm] + M[i][l].u2 * x[i+2][lm];
		}
	}
}


/*	Grid Generation with radial derivation operators.
*/


long int Reg_Grid(double rmin, double rg, double rmax)		// Regular Grid
{
	double dx1,dx2;
	long int i;

	if (rmin > rg) runerr("rmin > rg");
	if (rg >= rmax) runerr("rg >= rmax");

	NG = lround( (rg-rmin)/(rmax-rmin) * NR );	// number of radial points in inner core.

	// inner core
	dx1 = (rg - rmin)/NG;
	for (i=0; i<NG; i++)
		r[i] = ((rg-rmin)*i)/NG +rmin;
	r[NG] = rg;

	// outer core
	dx2 = (rmax - rg)/(NR-NG-1);
	for (i=NG+1; i<NR; i++)
		r[i] = ((rmax-rg)*(i-NG))/(NR-NG-1) +rg;

	printf("[GRID:Reg] NR=%d, r=[%f, %f], NG=%d, Ric=%f, dx1=%f, dx2=%f\n", NR,r[0], r[NR-1], NG, r[NG], dx1, dx2);
}

// BL_Grid : grid with densification in Boundary Layers, based on code from D. Jault.
// in:  rmin, rg, rmax are respectively the minimal, inner core and outer core radius.
// out: r[NR] is set with grid points (has to be allocated before function call)
//      NG is set with inner core index.
void BL_Grid(double rmin, double rg, double rmax)
{
	// La somme des N premieres puissances de l'inconnue x est y
	// on suppose x<1
	double decre(int N,double y)
	{	double e,f,x;
		
		e=y/(y+1.);	f=N+1.;
		x=e;
		while(fabs(y+1-(1-pow(x,f))/(1-x))>0.0001)
			x=e+pow(x,f)/(y+1);
		return x;
	}

	int nh = 0;		// can be made a parameter if required : add nh points to the inner hartman layer.
	int nr1,nr2, i,j;
	double e,q,hu,h;

	if (rmin != 0.0) runerr("rmin != 0 not supported.");
	if (rmax != 1.0) runerr("rmax != 1 not supported.");

	NG = (NR-nh-1)/5 +nh/4;
	nr1 = NG + (NR-nh-1)/5 +3*nh/4;
//	nr2=nr1+(NR-nh-1)/3;
	nr2 = (NR-1) - (NR-nh-1)/5;
	
	// each BL has thickness h = 0.15 of the gap
	// nr1 index of external shell of internal BL
	// nr2 index of internal shell of external BL
	h = (rmax-rg)*0.15;
	hu=(1.0 -h*2. -rg)/(nr2-nr1);
	r[NG]=rg;
	r[nr1]=rg+h;	r[nr2]=rmax-h;
	r[0]=rmin;	r[NR-1]=rmax;

	// Uniform grid in the bulk
	for(i=nr1+1; i<=nr2-1; i++)
		r[i] = r[nr1] + (i-nr1)*hu;

	// Outer boundary layer
	q=decre(NR-1-nr2,h/hu);
	e=hu;
	for(i=nr2+1; i<=NR-2; i++) {
		e*=q;
		r[i]=r[i-1]+e;
	}

	// Inner boundary layer
	q=decre(nr1-NG,h/hu);
	e=hu;
	for(i=nr1-1; i>=NG+1; i--) {
		e*=q;
		r[i]=r[i+1]-e;
	}

	// Inner core
	r[NG-1] = 2*rg - r[NG+1];
	e = rg - r[NG-1];
	q = 1.14;
	for(i=NG-2; e<=r[i+1]/(i+1); i--) {
		e*=q;
		r[i]=r[i+1]-e;
	}
	hu=r[i+2]/(i+2);
	for(j=0; j<=i+1; j++)
		r[j]=j*hu;

	printf("[GRID:BL] NR=%d, r=[%f, %f], NG=%d, Ric=%f\n", NR,r[0], r[NR-1], NG, r[NG]);
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

		Wr[i].l = r_1[i]* Gr[i].l * r[i-1];	// Wr = 1/r * d/dr(r .)
		Wr[i].d =         Gr[i].d;
		Wr[i].u = r_1[i]* Gr[i].u * r[i+1];

		D2r[i].l =  2.0*dr[i]*t;
		D2r[i].d =  -2.0*(dr[i-1]+dr[i])*t;
		D2r[i].u =  2.0*dr[i-1]*t;

		Lr[i].l = r_1[i]* D2r[i].l * r[i-1];	// Laplacien radial : 1/r . d2/dr2(r .)
		Lr[i].d =         D2r[i].d;
		Lr[i].u = r_1[i]* D2r[i].u * r[i+1];
/*
		Lr[i].l = D2r[i].l + 2.0*r_1[i]*Gr[i].l;	// Laplacien radial : d2/dr2 + 2/r.d/r
		Lr[i].d = D2r[i].d + 2.0*r_1[i]*Gr[i].d;
		Lr[i].u = D2r[i].u + 2.0*r_1[i]*Gr[i].u;
*/
	}
	// start and end must be determined by BC on a case by case basis.
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
void init_rad_sph(double rmin, double Ric, double rmax)
{
	r = (double *) malloc(NR * sizeof(double));	// allocation de la grille radiale
//	Reg_Grid(rmin, Ric, rmax);
	BL_Grid(rmin, Ric, rmax);
	init_Deriv_sph();
	NU = NR-NG;
}

/// finds the index to the shell with closest radius to radius rr
long int r_to_idx(double rr)
{
	long int i;

	i = NR-2;
	while((r[i] > rr)&&(i>0)) i--;
	if ((rr-r[i]) > (r[i+1]-rr)) i++;
	return i;	// i is always between 0 and NR-1
}
