
// common sht functions and api calls

/* regular transforms */

void SH_to_spat(shtns_cfg shtns, cplx *Qlm, double *Vr) {
	shtns->ftable.sy1(shtns, Qlm, Vr, shtns->lmax);
}

void spat_to_SH(shtns_cfg shtns, double *Vr, cplx *Qlm) {
	shtns->ftable.an1(shtns, Vr, Qlm, shtns->lmax);
}

void SHsphtor_to_spat(shtns_cfg shtns, cplx *Slm, cplx *Tlm, double *Vt, double *Vp) {
	shtns->ftable.sy2(shtns, Slm, Tlm, Vt, Vp, shtns->lmax);
}

void spat_to_SHsphtor(shtns_cfg shtns, double *Vt, double *Vp, cplx *Slm, cplx *Tlm) {
	shtns->ftable.an2(shtns, Vt, Vp, Slm, Tlm, shtns->lmax);
}

void SHqst_to_spat(shtns_cfg shtns, cplx *Qlm, cplx *Slm, cplx *Tlm, double *Vr, double *Vt, double *Vp) {
	shtns->ftable.sy3(shtns, Qlm, Slm, Tlm, Vr, Vt, Vp, shtns->lmax);
}

void spat_to_SHqst(shtns_cfg shtns, double *Vr, double *Vt, double *Vp, cplx *Qlm, cplx *Slm, cplx *Tlm) {
	shtns->ftable.an3(shtns, Vr, Vt, Vp, Qlm, Slm, Tlm, shtns->lmax);
}

void SHtor_to_spat(shtns_cfg shtns, cplx *Tlm, double *Vt, double *Vp) {
	shtns->ftable.gto(shtns, Tlm, Vt, Vp, shtns->lmax);
}

void SHsph_to_spat(shtns_cfg shtns, cplx *Slm, double *Vt, double *Vp) {
	shtns->ftable.gsp(shtns, Slm, Vt, Vp, shtns->lmax);
}


/* variable l-truncated transforms */

void SH_to_spat_l(shtns_cfg shtns, cplx *Qlm, double *Vr, int ltr) {
	shtns->ftable_l.sy1(shtns, Qlm, Vr, ltr);
}

void spat_to_SH_l(shtns_cfg shtns, double *Vr, cplx *Qlm, int ltr) {
	shtns->ftable_l.an1(shtns, Vr, Qlm, ltr);
}

void SHsphtor_to_spat_l(shtns_cfg shtns, cplx *Slm, cplx *Tlm, double *Vt, double *Vp, int ltr) {
	shtns->ftable_l.sy2(shtns, Slm, Tlm, Vt, Vp, ltr);
}

void spat_to_SHsphtor_l(shtns_cfg shtns, double *Vt, double *Vp, cplx *Slm, cplx *Tlm, int ltr) {
	shtns->ftable_l.an2(shtns, Vt, Vp, Slm, Tlm, ltr);
}

void SHqst_to_spat_l(shtns_cfg shtns, cplx *Qlm, cplx *Slm, cplx *Tlm, double *Vr, double *Vt, double *Vp, int ltr) {
	shtns->ftable_l.sy3(shtns, Qlm, Slm, Tlm, Vr, Vt, Vp, ltr);
}

void spat_to_SHqst_l(shtns_cfg shtns, double *Vr, double *Vt, double *Vp, cplx *Qlm, cplx *Slm, cplx *Tlm, int ltr) {
	shtns->ftable_l.an3(shtns, Vr, Vt, Vp, Qlm, Slm, Tlm, ltr);
}

void SHtor_to_spat_l(shtns_cfg shtns, cplx *Tlm, double *Vt, double *Vp, int ltr) {
	shtns->ftable_l.gto(shtns, Tlm, Vt, Vp, ltr);
}

void SHsph_to_spat_l(shtns_cfg shtns, cplx *Slm, double *Vt, double *Vp, int ltr) {
	shtns->ftable_l.gsp(shtns, Slm, Vt, Vp, ltr);
}


/* successive scalar + vector for 3D transform (can be faster than simultaneous transform) */
static
void SHqst_to_spat_2(shtns_cfg shtns, cplx *Qlm, cplx *Slm, cplx *Tlm, double *Vr, double *Vt, double *Vp) {
	const int ltr = shtns->lmax;
	shtns->ftable.sy1(shtns, Qlm, Vr, ltr);
	shtns->ftable.sy2(shtns, Slm, Tlm, Vt, Vp, ltr);
}

static
void SHqst_to_spat_2l(shtns_cfg shtns, cplx *Qlm, cplx *Slm, cplx *Tlm, double *Vr, double *Vt, double *Vp, int ltr) {
	shtns->ftable_l.sy1(shtns, Qlm, Vr, ltr);
	shtns->ftable_l.sy2(shtns, Slm, Tlm, Vt, Vp, ltr);
}

static
void spat_to_SHqst_2(shtns_cfg shtns, double *Vr, double *Vt, double *Vp, cplx *Qlm, cplx *Slm, cplx *Tlm) {
	const int ltr = shtns->lmax;
	shtns->ftable.an1(shtns, Vr, Qlm, ltr);
	shtns->ftable.an2(shtns, Vt, Vp, Slm, Tlm, ltr);
}

static
void spat_to_SHqst_2l(shtns_cfg shtns, double *Vr, double *Vt, double *Vp, cplx *Qlm, cplx *Slm, cplx *Tlm, int ltr) {
	shtns->ftable_l.an1(shtns, Vr, Qlm, ltr);
	shtns->ftable_l.an2(shtns, Vt, Vp, Slm, Tlm, ltr);
}


#if defined(SHT_F77_API)

/*  Fortran 77 api  */

extern shtns_cfg sht_data;

// Fortran API : Call from fortran without the trailing '_'
//@{

	// regular
/// \ingroup fortapi
void shtns_sh_to_spat_(cplx *Qlm, double *Vr) {
	sht_data->ftable.sy1(sht_data, Qlm, Vr, sht_data->lmax);
}

/// \ingroup fortapi
void shtns_spat_to_sh_(double *Vr, cplx *Qlm) {
	sht_data->ftable.an1(sht_data, Vr, Qlm, sht_data->lmax);
}

/// \ingroup fortapi
void shtns_sphtor_to_spat_(cplx *Slm, cplx *Tlm, double *Vt, double *Vp) {
	sht_data->ftable.sy2(sht_data, Slm, Tlm, Vt, Vp, sht_data->lmax);
}

/// \ingroup fortapi
void shtns_spat_to_sphtor_(double *Vt, double *Vp, cplx *Slm, cplx *Tlm) {
	sht_data->ftable.an2(sht_data, Vt, Vp, Slm, Tlm, sht_data->lmax);
}

/// \ingroup fortapi
void shtns_qst_to_spat_(cplx *Qlm, cplx *Slm, cplx *Tlm, double *Vr, double *Vt, double *Vp) {
	sht_data->ftable.sy3(sht_data, Qlm, Slm, Tlm, Vr, Vt, Vp, sht_data->lmax);
}

/// \ingroup fortapi
void shtns_spat_to_qst_(double *Vr, double *Vt, double *Vp, cplx *Qlm, cplx *Slm, cplx *Tlm) {
	sht_data->ftable.an3(sht_data, Vr, Vt, Vp, Qlm, Slm, Tlm, sht_data->lmax);
}

/// \ingroup fortapi
void shtns_sph_to_spat_(cplx *Slm, double *Vt, double *Vp) {
	sht_data->ftable.gsp(sht_data, Slm, Vt, Vp, sht_data->lmax);
}

/// \ingroup fortapi
void shtns_tor_to_spat_(cplx *Tlm, double *Vt, double *Vp) {
	sht_data->ftable.gto(sht_data, Tlm, Vt, Vp, sht_data->lmax);
}


	// variable l-truncation
/// \ingroup fortapi
void shtns_sh_to_spat_l_(cplx *Qlm, double *Vr, int* ltr) {
	sht_data->ftable_l.sy1(sht_data, Qlm, Vr, *ltr);
}

/// \ingroup fortapi
void shtns_spat_to_sh_l_(double *Vr, cplx *Qlm, int* ltr) {
	sht_data->ftable_l.an1(sht_data, Vr, Qlm, *ltr);
}

/// \ingroup fortapi
void shtns_sphtor_to_spat_l_(cplx *Slm, cplx *Tlm, double *Vt, double *Vp, int* ltr) {
	sht_data->ftable_l.sy2(sht_data, Slm, Tlm, Vt, Vp, *ltr);
}

/// \ingroup fortapi
void shtns_spat_to_sphtor_l_(double *Vt, double *Vp, cplx *Slm, cplx *Tlm, int* ltr) {
	sht_data->ftable_l.an2(sht_data, Vt, Vp, Slm, Tlm, *ltr);
}

/// \ingroup fortapi
void shtns_qst_to_spat_l_(cplx *Qlm, cplx *Slm, cplx *Tlm, double *Vr, double *Vt, double *Vp, int* ltr) {
	sht_data->ftable_l.sy3(sht_data, Qlm, Slm, Tlm, Vr, Vt, Vp, *ltr);
}

/// \ingroup fortapi
void shtns_spat_to_qst_l_(double *Vr, double *Vt, double *Vp, cplx *Qlm, cplx *Slm, cplx *Tlm, int* ltr) {
	sht_data->ftable_l.an3(sht_data, Vr, Vt, Vp, Qlm, Slm, Tlm, *ltr);
}

/// \ingroup fortapi
void shtns_sph_to_spat_l_(cplx *Slm, double *Vt, double *Vp, int* ltr) {
	sht_data->ftable_l.gsp(sht_data, Slm, Vt, Vp, *ltr);
}

/// \ingroup fortapi
void shtns_tor_to_spat_l_(cplx *Tlm, double *Vt, double *Vp, int* ltr) {
	sht_data->ftable_l.gto(sht_data, Tlm, Vt, Vp, *ltr);
}



//@}
#endif
