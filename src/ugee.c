/* ugee.c -- unified c source code for gee for Splus */
/* ugee.c -- ugee.c 4.13 98/01/27 */
/* /proj/stdevs/stdev0f/SLIBS/gee.dev/GEE97/SCCS/s.ugee.c */
/* combines sources from the following v 4.9 files */
/* cgee.c          chanmat.c       clinluxxy.c     diag_as_vec.c   normlib.c */

#include <R.h>

/* undefine duplicated macros */
#undef min
#undef max
#undef PI
#undef VOID

/* gee support @(#) previously cgee.c 4.12 98/01/26 */
#include "ugee.h"

static MATRIX *VC_GEE_diag_as_vec();
static MATRIX *VC_GEE_matsqrt();
static MATRIX *VC_GEE_mat1over();

void Cgee(x, y, id, n, offset, nobs, p,
	  parmvec, M_parm, S_beta, S_naivvar, S_robvar,
	  S_phi, S_R, tol, maxsz, S_iter, silent, errorstate, scale_fix,
	  compatflag)
double *x, *y, *id, *offset, *n;
Sint *nobs, *p, *M_parm, *compatflag;
Sint *silent, *parmvec;
Sint *S_iter;
double *S_beta, *S_naivvar, *S_robvar, *S_phi, *S_R, *tol;
Sint *maxsz, *errorstate, *scale_fix;
{
/* MAIN DECLS */

    MATRIX *xin, *yin, *idin, *offin, *nin;
    MATRIX **X, **Y, **OFFSET;
    MATRIX **N;
    /* MATRIX  *xpx, *xpy, *tmp; */
    MATRIX *beta, *betasave;
    MATRIX *alpha;
    MATRIX *R, *tmpR, *Ri;
    MATRIX *One, *mui; /*, *Opmui; */
    MATRIX *Ai, *ei, *ete;
    MATRIX *S1, *S2, *Di, *this_R, *S5, *S2i;
    MATRIX *tempmat1, *tempmat2, *Aop, *Dop, *zi, *DRop;
    MATRIX *robvar, *naivvar, *lag_wts, *tmpeep, *scratch, *wt;
    double phi, dni, phiLZ;
    int iter, ini, i2, j2, k;
    int alpha_VC_GEE_bandwidth;
    int *onep, one, nclust, i; /*, j; */
    double alpha_scalar;
    double alpha_scalar_LZ, exdiv_LZ;
    int maxni, ni, *maxnip;
    int link, var_mean_rel, corstruct;
    int *maxiter;
    double nnsclust = 0.;  /* for counting non-singletons */

    maxiter = (int*) malloc(sizeof(int));
    *errorstate = UNKNOWN_FAILURE ;

    if (!(*silent)) Rprintf("@(#) ugee.c 98/01/26 Cgee: GEE C source version chanlib 4.12 \n");

/* Initialize data */

    alpha = NULL;
    alpha_scalar = 0.;
    alpha_scalar_LZ = 0.;
    exdiv_LZ = 0;
    iter = 0;
    one = 1.;

    *maxiter = *S_iter;
    link = (*parmvec) -1;
    var_mean_rel = (*(parmvec+1)) -1;
    corstruct = (*(parmvec+2)) -1;

    alpha_VC_GEE_bandwidth = *M_parm + 1;

    onep = &one;

    from_S(y, nobs, onep, yin)
    from_S(x, nobs, p, xin)
    from_S(id, nobs, onep, idin)
    from_S(offset, nobs, onep, offin)
    from_S(n, nobs, onep, nin)
    from_S(S_beta, p, onep, beta);

    nclust = VC_GEE_nchanges(idin);

#define set_matrix_array(arrname, nel)				\
    if (!(arrname = (MATRIX **)malloc((unsigned)(nel * sizeof(MATRIX *))))) \
	error("set_matrix_array (mac): out of memory, requesting %d elements",nel);

    set_matrix_array(X, nclust)
    set_matrix_array(Y, nclust)
    set_matrix_array(OFFSET, nclust)
    set_matrix_array(N, nclust)

    VC_GEE_split(xin, idin, X);
    VC_GEE_split(yin, idin, Y);
    VC_GEE_split(offin, idin, OFFSET);
    VC_GEE_split(nin, idin, N);

    VC_GEE_destroy_matrix(xin);
    VC_GEE_destroy_matrix(yin);
    VC_GEE_destroy_matrix(idin);
    VC_GEE_destroy_matrix(offin);
    VC_GEE_destroy_matrix(nin);

    maxni = Y[0]->nrows;
    for (i = 1 ; i < nclust ; i++) {
	ni = Y[i]->nrows;
	if (ni > maxni) maxni = ni;
    }
    *maxsz = maxni;
    maxnip = &maxni;

/* xpx = VC_GEE_create_matrix(*p, *p, EPHEMERAL);
   xpy = VC_GEE_create_matrix(*p, one, EPHEMERAL); */

    if (corstruct == (int) fixed)
    {
	from_S(S_R, maxnip, maxnip, R);
    }

    if (!(*silent))
    {
	Rprintf("Cgee will use: ");
	switch(link) {
	case VC_GEE_identity:
	    Rprintf("VC_GEE_identity link, ");
	    break;
	case logarithm:
	    Rprintf("log link, ");
	    break;
	case logit:
	    Rprintf("logit link, ");
	    break;
	case reciprocal:
	    Rprintf("recip link, ");
	    break;
	case probit:
	    Rprintf("probit link, ");
	    break;
	case cloglog:
	    Rprintf("cloglog link, ");
	    break;
	default:
	    error("unknown link");
	    break;
	}
	switch(var_mean_rel)
	{
	case Gaussian:
	    Rprintf("Gaussian var, ");
	    break;
	case Poisson:
	    Rprintf("Poisson var, ");
	    break;
	case Binomial:
	    Rprintf("Binomial var, ");
	    break;
	case Gamma:
	    Rprintf("Gamma var, ");
	    break;
	default:
	    error("Cgee: unknown var_mean_rel. Dies.");
	    break;
	}
	switch (corstruct)
	{
	case independence:
	    Rprintf("indep corstr . ");
	    break;
	case exchangeable:
	    Rprintf("exch corstr . ");
	    break;
	case stat_M_dep:
	    Rprintf("stat corstr . ");
	    break;
	case AR_M:
	    Rprintf("AR-M corstr . ");
	    break;
	case non_stat_M_dep:
	    Rprintf("non-statM corstr . ");
	    break;
	case unstructured:
	    Rprintf("unstr corstr . ");
	    break;
	case fixed:
	    Rprintf("fixed corstr . ");
	    break;
	default:
	    error("unknown corstr");
	    break;
	}
	Rprintf("\n");
    }


    scratch = VC_GEE_create_matrix(maxni, maxni, EPHEMERAL);

    do {
	betasave = VC_GEE_matcopy(beta);
	phi = 0.;
	phiLZ = 0.;
	nnsclust = 0.;

	tmpR = VC_GEE_create_matrix(maxni, maxni, EPHEMERAL);


	switch (corstruct)
	{
	case independence:
	    break;
	case exchangeable:
	    alpha_scalar = 0.;
	    alpha_scalar_LZ = 0.;
	    exdiv_LZ = 0.;
	    break;
	case stat_M_dep:
	case AR_M:
	    alpha = VC_GEE_create_matrix(1, alpha_VC_GEE_bandwidth, EPHEMERAL);
	    break;
	case non_stat_M_dep:
	case unstructured:
	    alpha = VC_GEE_create_matrix(maxni, maxni, EPHEMERAL);
	    break;
	case fixed:
	    break;
	default:
	    break;
	}

	make_permanent(beta);
	for (i = 0 ; i < nclust ; i++)
	{
	    ni = Y[i]->nrows;
	    dni = (double)ni;
	    One = VC_GEE_col_1s(ni);
	    make_permanent(One);
	    tempmat1 = VC_GEE_matadd(VC_GEE_matmult(X[i], beta), OFFSET[i]);

	    switch (link)
	    {
		double maxfitted;
	    case VC_GEE_identity:
		mui = VC_GEE_matcopy(tempmat1);
		break;
	    case logarithm:
		mui = VC_GEE_matexp(tempmat1);
		break;
	    case logit:
		tempmat1 = VC_GEE_matexp(tempmat1);
		make_permanent(tempmat1);
		tempmat2 = VC_GEE_matadd(One, tempmat1);
		mui = VC_GEE_pxq_divby_px1(VC_GEE_px1_times_pxq(N[i],
								tempmat1), tempmat2);
		make_permanent(mui);
		maxfitted = VC_GEE_matmax(VC_GEE_pxq_divby_px1(mui,N[i]));
		if ((maxfitted >= .9999 || maxfitted <= .0001 
		     || ISNAN(maxfitted))
		    && var_mean_rel == (int) Binomial)
		{
		    error("Cgee: error: logistic model for probability has fitted value very close to 1.\nestimates diverging; iteration terminated.");
		}
		VC_GEE_destroy_matrix(tempmat1);
		break;
	    case reciprocal:
		mui = VC_GEE_pxq_divby_px1(One, tempmat1);
		break;
/* probit case added by pj catalano*/  /* OK */
	    case probit:
		tempmat1 = VC_GEE_matncdf(tempmat1);
		mui = VC_GEE_px1_times_pxq(N[i], tempmat1);
		make_permanent(mui);
		maxfitted = VC_GEE_matmax(VC_GEE_pxq_divby_px1(mui, N[i]));
//		if ((maxfitted >= .9999) && var_mean_rel == Binomial)
		if ((maxfitted >= .9999 || maxfitted <= .0001 
		     || ISNAN(maxfitted))
		    && var_mean_rel == (int) Binomial)
		{
		    error("Cgee: estimates diverging; iteration terminated");
		}
		break;

		/* cloglog case added by jh maindonald */
	    case cloglog:
		tempmat1=VC_GEE_matanticlog(tempmat1);
		mui=VC_GEE_px1_times_pxq(N[i], tempmat1);
		make_permanent(mui);
		maxfitted = VC_GEE_matmax(VC_GEE_pxq_divby_px1(mui, N[i]));
//		if ((maxfitted >= .999999) && var_mean_rel == (int) Binomial)
		if ((maxfitted >= .9999 || maxfitted <= .0001 
		     || ISNAN(maxfitted))
		    && var_mean_rel == (int) Binomial)
		{
		    error("Cgee: error: cloglog model for probability has fit ted value very close to 1.\nestimates diverging; iteration terminated.");
		}
		break;

	    default:
		error("Cgee: unknown link. Dies.");
		break;
	    }
	    make_permanent(mui);
	    ei = VC_GEE_matsub(Y[i], mui);

	    switch(var_mean_rel)
	    {
	    case Gaussian:
		/*	Ai = VC_GEE_ident(ni); */
		break;
	    case Poisson:
		Ai = VC_GEE_form_diag(mui);
		break;
	    case Binomial:
		Ai = VC_GEE_form_diag(VC_GEE_px1_times_pxq(
					   mui, VC_GEE_matsub(One, VC_GEE_pxq_divby_px1(mui,N[i]))));
		break;
	    case Gamma:
		Ai = VC_GEE_form_diag(VC_GEE_px1_times_pxq(mui, mui));
		break;
	    default:
		error("Cgee: unknown var_mean_rel. Dies.");
		break;
	    }
	    VC_GEE_destroy_matrix(mui);

	    if (var_mean_rel != Gaussian)
		ei = VC_GEE_px1_times_pxq(VC_GEE_mat1over(VC_GEE_matsqrt(VC_GEE_diag_as_vec(Ai))), ei);

	    make_permanent(ei);

	    ete = VC_GEE_matmult(VC_GEE_transp(ei), ei);

	    phi += MEL(ete,0,0)/dni;
	    phiLZ += MEL(ete,0,0);
	    if (dni>1) nnsclust += 1.0;

	    if (corstruct == (int) stat_M_dep || corstruct == (int) AR_M)
	    {
		if (dni < (double)alpha_VC_GEE_bandwidth)
		    error("cgee: M-dependence, M=%d, but clustsize=%d\nfatal error for this model",(int)*M_parm,(int)dni);

		lag_wts = VC_GEE_create_matrix(1, alpha_VC_GEE_bandwidth, PERMANENT);
		MEL(lag_wts, 0, 0) = (double)1.;

		for (ini = 1 ; ini < alpha_VC_GEE_bandwidth ; ini++)
		{
		    MEL(lag_wts, 0, ini) = dni/(dni - (double)ini);
		}
	    }

	    switch (corstruct)
	    {
	    case independence:
		break;
	    case exchangeable:
		if (ni > 1) alpha_scalar +=
		    1./(dni*(dni-1.)) * (VC_GEE_elsum(VC_GEE_matmult(ei, VC_GEE_transp(ei))) - MEL(ete,0,0));
		if (ni > 1) alpha_scalar_LZ +=
		    (VC_GEE_elsum(VC_GEE_matmult(ei, VC_GEE_transp(ei))) - MEL(ete,0,0));
		if (ni > 1) exdiv_LZ += dni*(dni-1); /* suppress .5 because num
							  is redundant */
		break;
	    case stat_M_dep:
	    case AR_M:

		alpha = VC_GEE_matadd(alpha, VC_GEE_transp(VC_GEE_px1_times_pxq (VC_GEE_transp(lag_wts), VC_GEE_transp(VC_GEE_covlag(ei, alpha_VC_GEE_bandwidth, 0)))));
		break;
	    case non_stat_M_dep:
	    case unstructured:
		tmpeep = VC_GEE_matmult(ei, VC_GEE_transp(ei));
		VC_GEE_plug(tmpeep, scratch, 0, 0);
		make_permanent(scratch);
		alpha = VC_GEE_matadd(alpha, scratch);
		scratch = VC_GEE_matsub(scratch, scratch);
		make_ephemeral(scratch);
		break;
	    case fixed:
		break;
	    default:
		error("corstruct not implemented.");
		break;
	    }
	    VC_GEE_destroy_matrix(ei);
	    VC_GEE_destroy_matrix(One);
	}

        if  (alpha != NULL) alpha = VC_GEE_scalar_times_matrix(
	    (double)nclust/(phi*nnsclust), alpha);

	alpha_scalar /= phi;
        alpha_scalar *= (double)nclust/nnsclust;

	phi /= (double)nclust;
	phiLZ /= ((double)(*nobs - *p));

	alpha_scalar_LZ /= (phiLZ*(exdiv_LZ-2.*(double)*p));
	if (*compatflag == 0) {  /* abandon compatibility with early macro */
	    alpha_scalar = alpha_scalar_LZ;
	    phi = phiLZ;
	}

	switch(corstruct)
	{
	case independence:  /* this is ridiculous */
	    R = VC_GEE_ident(maxni);
	    break;
	case exchangeable:
	    R = VC_GEE_matadd(VC_GEE_scalar_times_matrix(alpha_scalar, VC_GEE_matmult(VC_GEE_col_1s(maxni), VC_GEE_transp(VC_GEE_col_1s(maxni)))), VC_GEE_scalar_times_matrix((double)1.-alpha_scalar,VC_GEE_ident(maxni)));
	    break;
	case stat_M_dep:
	    tmpR = VC_GEE_create_matrix(1, maxni, EPHEMERAL);
	    VC_GEE_plug(alpha, tmpR, 0, 0);
	    MEL(tmpR, 0, 0) = (double)1.;
	    R = VC_GEE_toeplitz(tmpR);
	    break;
	case AR_M:
	    tmpR = VC_GEE_create_matrix(1, maxni, EPHEMERAL);
	    MEL(alpha, 0, 0) = (double)1.;
	    make_permanent(alpha);
	    VC_GEE_plug(alpha, tmpR, 0,0);
	    tempmat1 = VC_GEE_extract_cols(alpha, 0, alpha_VC_GEE_bandwidth-2);
	    wt = VC_GEE_toeplitz(tempmat1);
	    tempmat1 = VC_GEE_extract_cols(alpha, 1, alpha_VC_GEE_bandwidth-1);
	    alpha = VC_GEE_matmult(tempmat1, VC_GEE_luinv(wt));
	    for (i2 = alpha_VC_GEE_bandwidth ; i2 < maxni ; i2++)
	    {
		for (j2 = 0 ; j2 < alpha_VC_GEE_bandwidth-1 ; j2++)
		{
		    MEL(tmpR, 0, i2) += MEL(alpha,0,j2) *
			MEL(tmpR,0,(i2-j2)-1);
		}
	    }
	    R = VC_GEE_toeplitz(tmpR);
	    VC_GEE_destroy_matrix(alpha);
	    break;
	case non_stat_M_dep:
	    R = VC_GEE_band(alpha, alpha_VC_GEE_bandwidth);
	    for (k = 0 ; k < R->ncols ; k++)
	    {
		MEL(R, k, k) = (double)1.;
	    }
	    break;
	case unstructured:
	    R = VC_GEE_matcopy(alpha);
	    for (k = 0 ; k < R->ncols ; k++)
	    {
		MEL(R, k, k) = (double)1.;
	    }
	    break;
	case fixed:
	    break;
	default:
	    error("corstruct not implemented.");
	    break;
	}
	make_permanent(R);

	S1 = VC_GEE_create_matrix((int)*p, 1, EPHEMERAL);
	S2 = VC_GEE_create_matrix((int)*p, (int)*p, EPHEMERAL);

	for (i = 0 ; i < nclust ; i++)
	{
	    ni = Y[i]->nrows;
	    dni = (double)ni;
	    One = VC_GEE_col_1s(ni);
	    make_permanent(One);

	    tempmat1 = VC_GEE_matadd(VC_GEE_matmult(X[i], beta), OFFSET[i]);
	    switch (link)
	    {
	    case VC_GEE_identity:
		Di = VC_GEE_matcopy(X[i]);
		mui = VC_GEE_matcopy(tempmat1);
		break;
	    case logarithm:
		tempmat1 = VC_GEE_matexp(tempmat1);
		mui = VC_GEE_matcopy(tempmat1);
		Di = VC_GEE_px1_times_pxq(tempmat1, X[i]);
		break;
	    case logit:
		tempmat1 = VC_GEE_matexp(tempmat1);
		make_permanent(tempmat1);
		tempmat2 = VC_GEE_matadd(One, tempmat1);
		tempmat2 = VC_GEE_pxq_divby_px1(
		    VC_GEE_px1_times_pxq(N[i], tempmat1), tempmat2);
		make_permanent(tempmat2);
		VC_GEE_destroy_matrix(tempmat1);
		tempmat1 = VC_GEE_matsub(One, VC_GEE_pxq_divby_px1(tempmat2,N[i]));
		tempmat1 = VC_GEE_px1_times_pxq(tempmat1,
						 tempmat2);
		Di = VC_GEE_px1_times_pxq(tempmat1, X[i]);
		mui = VC_GEE_matcopy(tempmat2);
		VC_GEE_destroy_matrix(tempmat2);
		break;
	    case reciprocal:
		tempmat1 = VC_GEE_pxq_divby_px1(
		    One, tempmat1);
		mui = VC_GEE_matcopy(tempmat1);
		tempmat2 = VC_GEE_matcopy(tempmat1);
		tempmat1 = VC_GEE_px1_times_pxq(tempmat1, tempmat2);
		tempmat1 = VC_GEE_px1_times_pxq(
		    tempmat1, X[i]);
		Di = VC_GEE_scalar_times_matrix(-1., tempmat1);
		break;
		/* probit case added by pj catalano */
	    case probit:
		tempmat2 = VC_GEE_matcopy(tempmat1);
		tempmat2 = VC_GEE_matncdf(tempmat2);
		mui = VC_GEE_px1_times_pxq(N[i],tempmat2);
		tempmat1 = VC_GEE_matnpdf(tempmat1);
		tempmat2 = VC_GEE_px1_times_pxq(N[i],tempmat1);
		Di = VC_GEE_px1_times_pxq(tempmat2, X[i]);
		break;
/* cloglog case added by jh maindonald */
	    case cloglog:
		tempmat2=VC_GEE_matcopy(tempmat1);
		tempmat1=VC_GEE_matanticlog(tempmat1);
		make_permanent(tempmat1);
		mui=VC_GEE_px1_times_pxq(N[i], tempmat1);
		tempmat2=VC_GEE_matexp(tempmat2);
		tempmat2=VC_GEE_px1_times_pxq(N[i],tempmat2);
		tempmat1=VC_GEE_px1_times_pxq(tempmat2,VC_GEE_matsub(One,tempmat1));
		Di=VC_GEE_px1_times_pxq(tempmat1,X[i]);
		VC_GEE_destroy_matrix(tempmat1);
		break;

	    default:
		error("Cgee: unknown link. Dies.");
		break;
	    }
	    make_permanent(mui);
	    make_permanent(Di);

	    zi = VC_GEE_matadd(VC_GEE_matmult(Di, beta), VC_GEE_matsub(Y[i],
									  mui));

	    switch(var_mean_rel)
	    {
	    case Gaussian:
		/* Ai = VC_GEE_ident(ni); */
		break;
	    case Poisson:
		Ai = VC_GEE_form_diag(mui);
		break;
	    case Binomial:
		Ai = VC_GEE_form_diag(VC_GEE_px1_times_pxq(
					   mui, VC_GEE_matsub(One, VC_GEE_pxq_divby_px1(mui,N[i]))));
		break;
	    case Gamma:
		Ai = VC_GEE_form_diag(VC_GEE_px1_times_pxq(
					   mui, mui));
		break;
	    default:
		error("Cgee: unknown var_mean_rel. Dies.");
		break;
	    }
	    VC_GEE_destroy_matrix(mui);

	    if (var_mean_rel != Gaussian)  /* else Ai is VC_GEE_identity */
		Aop = VC_GEE_mat1over(VC_GEE_matsqrt(VC_GEE_diag_as_vec(Ai)));
	    if (var_mean_rel != Gaussian)  /* else Ai is VC_GEE_identity */
		make_permanent(Aop);


	    make_ephemeral(Di);
	    if (var_mean_rel != Gaussian)
		Di = VC_GEE_px1_times_pxq(Aop, Di);
	    make_permanent(Di);
	    if (var_mean_rel != Gaussian)
		zi = VC_GEE_px1_times_pxq(Aop, zi);

	    if (corstruct != independence)
	    {
		this_R = VC_GEE_corner(R, ni, ni);
        	/* Rprintf("cluster %d g\n",i); */

		/* Rprintf("VC_GEE_start inver\n");
		   Rprintf("D: %d x %d, R: %d x %d\n",Di->nrows,
		   Di->ncols, this_R->nrows, this_R->ncols); */
		Ri = VC_GEE_luinv(this_R);
		/* Rprintf("D: %d x %d, Ri: %d x %d\n",Di->nrows,
		   Di->ncols, Ri->nrows, Ri->ncols); */
		Dop = VC_GEE_matmult(VC_GEE_transp(Di), Ri);
		/* Rprintf("end inver\n"); */
	    }
	    else Dop = VC_GEE_transp(Di);

	    make_permanent(Dop);

	    S1 = VC_GEE_matadd(S1, VC_GEE_matmult(Dop, zi));
	    S2 = VC_GEE_matadd(S2, VC_GEE_matmult(Dop, Di));

	    VC_GEE_destroy_matrix(Dop);
	    if (var_mean_rel != Gaussian)
		VC_GEE_destroy_matrix(Aop);
	    VC_GEE_destroy_matrix(Di);
	    VC_GEE_destroy_matrix(One);

	}

	VC_GEE_destroy_matrix(beta);
	beta = VC_GEE_matmult(VC_GEE_luinv(S2), S1);
	make_permanent(beta);

	if (!(*silent)) Rprintf("current parameter estimates:\n");

	if (!(*silent)) VC_GEE_matdump(beta);

	if (!(*silent)) Rprintf("completed iteration %d\n",iter);
	iter++;



    } while((VC_GEE_matmax(VC_GEE_matabs(VC_GEE_matsub(VC_GEE_pxq_divby_px1(betasave, beta), VC_GEE_col_1s((int)*p)))) > *tol) && (iter < *maxiter));

    if (iter >= *maxiter)
    {
	warning("Maximum number of iterations consumed");
	warning("Convergence not achieved; results suspect");
	*errorstate = MAXITER_EXCEEDED ;
    }

    S2 = VC_GEE_create_matrix((int)*p, (int)*p, EPHEMERAL);
    S5 = VC_GEE_create_matrix((int)*p, (int)*p, EPHEMERAL);
    phi = 0.;
    phiLZ = 0.;

    for (i = 0 ; i < nclust ; i++)
    {
	ni = Y[i]->nrows;
	dni = (double)ni;
	One = VC_GEE_col_1s(ni);
	make_permanent(One);

	tempmat1 = VC_GEE_matadd(VC_GEE_matmult(X[i], beta), OFFSET[i]);
	switch (link)
	{
	case VC_GEE_identity:
	    Di = VC_GEE_matcopy(X[i]);
	    mui = VC_GEE_matcopy(tempmat1);
	    break;
	case logarithm:
	    tempmat1 = VC_GEE_matexp(tempmat1);
	    mui = VC_GEE_matcopy(tempmat1);
	    Di = VC_GEE_px1_times_pxq(tempmat1, X[i]);
	    break;
	case logit:
	    tempmat1 = VC_GEE_matexp(tempmat1);
	    make_permanent(tempmat1);
	    tempmat2 = VC_GEE_matadd(One, tempmat1);
	    tempmat2 = VC_GEE_pxq_divby_px1(
		VC_GEE_px1_times_pxq(N[i], tempmat1), tempmat2);
	    make_permanent(tempmat2);
	    VC_GEE_destroy_matrix(tempmat1);
	    tempmat1 = VC_GEE_matsub(One, VC_GEE_pxq_divby_px1(tempmat2,N[i]));
	    tempmat1 = VC_GEE_px1_times_pxq(tempmat1,
					     tempmat2);
	    Di = VC_GEE_px1_times_pxq(tempmat1, X[i]);
	    mui = VC_GEE_matcopy(tempmat2);
	    VC_GEE_destroy_matrix(tempmat2);
	    break;
	case reciprocal:
	    tempmat1 = VC_GEE_pxq_divby_px1(
		One, tempmat1);
	    mui = VC_GEE_matcopy(tempmat1);
	    tempmat2 = VC_GEE_matcopy(tempmat1);
	    tempmat1 = VC_GEE_px1_times_pxq(tempmat1, tempmat2);
	    tempmat1 = VC_GEE_px1_times_pxq(
		tempmat1, X[i]);
	    Di = VC_GEE_scalar_times_matrix(-1., tempmat1);
	    break;
	    /* probit case added by pj catalano*/
	case probit:
	    tempmat2 = VC_GEE_matcopy(tempmat1);
	    tempmat2 = VC_GEE_matncdf(tempmat2);
	    mui = VC_GEE_px1_times_pxq(N[i],tempmat2);
	    tempmat1 = VC_GEE_matnpdf(tempmat1);
	    tempmat2 = VC_GEE_px1_times_pxq(N[i],tempmat1);
	    Di = VC_GEE_px1_times_pxq(tempmat2, X[i]);
	    break;
/* maind */
	case cloglog:
	    tempmat2=VC_GEE_matcopy(tempmat1);
	    tempmat1=VC_GEE_matanticlog(tempmat1);
	    make_permanent(tempmat1);
	    mui=VC_GEE_px1_times_pxq(N[i], tempmat1);
	    tempmat2=VC_GEE_matexp(tempmat2);
	    tempmat2=VC_GEE_px1_times_pxq(N[i],tempmat2);
	    tempmat1=VC_GEE_px1_times_pxq(tempmat2,VC_GEE_matsub(One,tempmat1));
	    Di=VC_GEE_px1_times_pxq(tempmat1,X[i]);
	    VC_GEE_destroy_matrix(tempmat1);
	    break;

	default:
	    error("Cgee: unknown link. Dies.");
	    break;
	}
	make_permanent(mui);

	ei = VC_GEE_matsub(Y[i],mui);
	switch(var_mean_rel)
	{
	case Gaussian:
	    /* Ai = VC_GEE_ident(ni); */
	    break;
	case Poisson:
	    Ai = VC_GEE_form_diag(mui);
	    break;
	case Binomial:
	    Ai = VC_GEE_form_diag(VC_GEE_px1_times_pxq(
				       mui, VC_GEE_matsub(One, VC_GEE_pxq_divby_px1(mui,N[i]))));
	    break;
	case Gamma:
	    Ai = VC_GEE_form_diag(VC_GEE_px1_times_pxq(
				       mui, mui));
	    break;
	default:
	    error("Cgee: unknown var_mean_rel. Dies.\n");
	    break;
	}
	VC_GEE_destroy_matrix(mui);
	VC_GEE_destroy_matrix(One);

	if (var_mean_rel != Gaussian)
	{
	    Aop = VC_GEE_mat1over(VC_GEE_matsqrt(VC_GEE_diag_as_vec(Ai)));
	    make_permanent(Aop);
	    Di = VC_GEE_px1_times_pxq(Aop, Di);
	}
	make_permanent(Di);

	if (var_mean_rel != Gaussian)
	    ei = VC_GEE_px1_times_pxq(Aop, ei);
	make_permanent(ei);

	ete = VC_GEE_matmult(VC_GEE_transp(ei), ei);
	phi += MEL(ete,0,0) / dni;
	phiLZ += MEL(ete,0,0) ;

	if (corstruct != independence)
	    DRop = VC_GEE_matmult(VC_GEE_transp(Di), VC_GEE_luinv(VC_GEE_corner(R, ni, ni)));
	else
	    DRop = VC_GEE_transp(Di);

	make_permanent(DRop);

	S2 = VC_GEE_matadd(S2, VC_GEE_matmult(DRop, Di));
/* Rprintf("VC_GEE_start S5\n"); */
	S5 = VC_GEE_matadd(S5, VC_GEE_matmult(VC_GEE_matmult(DRop, ei),
						VC_GEE_matmult(VC_GEE_transp(ei), VC_GEE_transp(DRop))));
	/* Rprintf("end S5\n"); */

	if (var_mean_rel != Gaussian)
	    VC_GEE_destroy_matrix(Aop);
	VC_GEE_destroy_matrix(Di);
	VC_GEE_destroy_matrix(ei);
	VC_GEE_destroy_matrix(DRop);
    }

    phi /= (double)nclust;
    phiLZ /= (double)(*nobs - *p);
    if (*compatflag == 0) {  /* abandon compatibility with early macro */
	phi = phiLZ;
    }

    S2i = VC_GEE_luinv(S2);
    make_permanent(S2i);
    if (*scale_fix) phi = *S_phi;
    if (*scale_fix && var_mean_rel == Gaussian)
	warning("Scale parameter fixed at %f with Gaussian variance function",
		*S_phi);
    naivvar = VC_GEE_scalar_times_matrix(phi, S2i);
    robvar = VC_GEE_matmult(VC_GEE_matmult (S2i, S5), S2i);

    to_S(beta, S_beta)
    to_S(naivvar, S_naivvar)
    to_S(robvar, S_robvar)
    to_S(R, S_R)

    *S_phi = phi;
    *S_iter = iter;

    VC_GEE_destroy_matrix(beta);
    VC_GEE_destroy_matrix(naivvar);
    VC_GEE_destroy_matrix(robvar);
    VC_GEE_destroy_matrix(R);

    for (i = 0 ; i < nclust ; i++)
    {
	VC_GEE_destroy_matrix(X[i]);
	VC_GEE_destroy_matrix(Y[i]);
	VC_GEE_destroy_matrix(N[i]);
	VC_GEE_destroy_matrix(OFFSET[i]);
    }

    if (iter < *maxiter) *errorstate = NO_ERROR;
}



/* gee support @(#) chanmat.c 4.12 98/01/26 */

/* @(#) chanmat.nw 1.3 94/03/09 */

static MATRIX *VC_GEE_create_matrix(nrows, ncols, permanence)
/* $Y = |VC_GEE_create_matrix|(r,c,\cdot) :\Rightarrow Y \in M_{r\times c} \wedge
Y = 0 $ */
int nrows, ncols, permanence;
{
    MATRIX *tmp;
    double *head;
    int i;

    tmp = (MATRIX *) calloc (1, sizeof (struct matrix));

    if (tmp == NULL)
    {
	error("VC_GEE_create_matrix: malloc failed %d",
	      sizeof(struct matrix));
    }

    tmp->nrows = nrows;
    tmp->ncols = ncols;
    tmp->permanence = permanence;

    tmp->data = (double *) calloc (1,  nrows * ncols * sizeof (double)) ;

    if (tmp->data == NULL)
    {
	error("VC_GEE_create_matrix: malloc failed, nrows=%d ncols=%d",
	      nrows, ncols);
    }

    head = tmp->data;
    for (i = 0 ; i < nrows*ncols ; i++)
    {
	*head = (double)0.;
	head++;
    }

    return tmp;
}

static void VC_GEE_destroy_matrix(mat)
MATRIX *mat;
{
    if (mat == (MATRIX *) NULL) return;
    mat->nrows = 0;
    mat->ncols = 0;
/*if (mat->data != (double *)NULL) free((char *) mat->data);*/
    mat->data = (double *)NULL;
/*if (mat != (MATRIX *)NULL) free((char *) mat);*/
    mat = (MATRIX *)NULL;
}



static MATRIX *VC_GEE_transp(mat)
/*
$Y = VC_GEE_transp(X_{r\times c}) :\Rightarrow Y \in M_{c\times r} \wedge
y_{ji} = x_{ij}. $
*/
MATRIX *mat;
{
    double *telem, *inelem, *tbase;
    int nelem;
    MATRIX *tmp;

    tmp = VC_GEE_create_matrix(mat->ncols, mat->nrows, EPHEMERAL);
    inelem = mat->data;
    tbase = tmp->data;
    telem = tbase;
    for (nelem = 0 ; nelem < (mat->ncols * mat->nrows) ; nelem++)
    {
	*telem = *(inelem++);
	telem += mat->nrows;
	if (nelem % mat->ncols == (mat->ncols)-1)
	    telem = ++tbase;
    }
    if (is_ephemeral(mat)) VC_GEE_destroy_matrix(mat);
    return tmp;

}

static MATRIX *VC_GEE_corner(mat, nr, nc)
/*
$ r \leq r^{\prime} \wedge c \leq c^{\prime}
\wedge X \in M_{r^{\prime} \times c^{\prime}}
\wedge Y = $|VC_GEE_corner(X,r,c)|$ :\Rightarrow $ */
/* $ Y \in M_{r \times c} \wedge y_{ij} = x_{ij},
\thinspace 1 \leq i \leq r, \thinspace 1 \leq j \leq c $ */

MATRIX *mat;
int nr, nc;
{
    MATRIX *tmp;
    double *load;
    int i,j,sr, sc;
    sr = mat->nrows;
    sc = mat->ncols;
    if ((nr > sr) || (nc > sc))
    {
	error("VC_GEE_corner: request not a submatrix.\nfatal error");
    }
    tmp = VC_GEE_create_matrix(nr, nc, EPHEMERAL);
    load = tmp->data;
    for (i = 0 ; i < nr ; i++)
    {
	for (j = 0 ; j < nc ; j++)
	{
	    *(load++) = MEL(mat, i, j);
	}
    }
    free_if_ephemeral(mat);
    return tmp;
}


static MATRIX *VC_GEE_extract_rows(Source,VC_GEE_start,end)
	/* purely zero-based */

MATRIX *Source;
int VC_GEE_start, end;
{
    MATRIX *temp;
    int rows_to_get, i, j;

    rows_to_get = end - VC_GEE_start + 1;

    temp = VC_GEE_create_matrix(rows_to_get,Source->ncols,EPHEMERAL);

    for (i = 0 ; i < rows_to_get ; i++)
    {
	for (j = 0 ; j < Source->ncols ; j++)
	{
	    *(ELREF(temp,i,j)) = *(ELREF(Source,VC_GEE_start,j));
	}
	VC_GEE_start++;
    }
    /* DOES NOT CLEAN */
    return temp;
}

static MATRIX *VC_GEE_extract_cols(x, VC_GEE_start, end)
MATRIX *x;
int VC_GEE_start, end;
{
    MATRIX *tmp;
    tmp = VC_GEE_transp(x);
    tmp = VC_GEE_extract_rows(tmp, VC_GEE_start, end);
    tmp = VC_GEE_transp(tmp);
    free_if_ephemeral(x);
    return tmp;
}

static MATRIX *VC_GEE_matcopy(inmat)
MATRIX *inmat;
{
    int i, j;
    MATRIX *outmat;

    outmat = VC_GEE_create_matrix(inmat->nrows,inmat->ncols,EPHEMERAL);
    for (i = 0 ; i < inmat->nrows ; i++)
    {
	for (j = 0 ; j < inmat->ncols ; j++)
	{
	    *(ELREF(outmat,i,j)) = *(ELREF(inmat,i,j));
	}
    }
/* DOES NOT CLEAN */
    return outmat;
}

static int VC_GEE_split(matptr, discptr, matarrptr)
MATRIX *matptr, *discptr, *matarrptr[];
{   /* discriminator VC_GEE_vector assumed to be integer-valued dbls */
    int i, iVC_GEE_start, k, VC_GEE_start, end;
    if (discptr->ncols != 1)
    {
	error("VC_GEE_split: discriminator must be column vec.\nVC_GEE_split: ncols = %d", discptr->ncols);
    }

    k = 0;

    iVC_GEE_start = (int)MEL(discptr, 0, 0);
    VC_GEE_start = 0;
    end = 0;
    for (i = 1 ; i <= discptr->nrows ; i++)
    {
	if (i == discptr->nrows || MEL(discptr, i, 0) != iVC_GEE_start)
	{
	    matarrptr[k] = VC_GEE_matcopy(VC_GEE_extract_rows(matptr, VC_GEE_start, end));
	    make_permanent(matarrptr[k]);
	    k++;
	    VC_GEE_start = end+1;
	    if (i < discptr->nrows) /* don't need iVC_GEE_start at end of loop */
		iVC_GEE_start = MEL(discptr, i, 0);
	}
	if (VC_GEE_start < discptr->nrows) end++ ;
    }
    /* DOES NOT CLEAN */
    return k;
}


static void VC_GEE_plug(VC_GEE_plugm, socket, row, col)
int row, col;
MATRIX *VC_GEE_plugm, *socket;  /* not a unix socket */
{
    int pcol, prow;
    double *sockload, *VC_GEE_plughead, *sockrow_VC_GEE_start;
    int i,j;

    pcol = VC_GEE_plugm->ncols;
    prow = VC_GEE_plugm->nrows;

    if (pcol+col > socket->ncols || prow+row > socket->nrows)
    {
	error("M+-: VC_GEE_plug: socket too small");
    }

    sockload = socket->data + col + row*(socket->ncols);
    VC_GEE_plughead = VC_GEE_plugm->data;
    sockrow_VC_GEE_start = sockload;

    for (i = 0 ; i < prow ; i++)
    {
	sockload = sockrow_VC_GEE_start;
	for (j = 0 ; j < pcol ; j++)
	{
	    *(sockload++) = *(VC_GEE_plughead++);
	}
	sockrow_VC_GEE_start += socket->ncols;
    }
    free_if_ephemeral(VC_GEE_plugm);
}

static MATRIX *VC_GEE_form_diag(vec)
MATRIX *vec;
{
    MATRIX *tmp;
    int i, ord;

    ord = vec->nrows;
    tmp = VC_GEE_create_matrix(ord,  ord, EPHEMERAL);
    for (i = 0 ; i < ord ; i++)
	*(ELREF(tmp,i,i)) = MEL(vec,i,0);
    free_if_ephemeral(vec);
    return tmp;
}

static MATRIX *VC_GEE_band(in, wid)
MATRIX *in;
int wid;
{
    MATRIX *tmp;
    int i, j;
    tmp = VC_GEE_matcopy(in);
    for (i = 0 ; i < in->nrows ; i++)
    {
	for (j = i+wid ; j < in->ncols ; j++)
	{
	    MEL(tmp, i, j) = (double)0.;
	    if ((i < in->ncols) && (j < in->nrows))
	    {
		MEL(tmp, j, i) = (double)0.;
	    }
	}
    }
    free_if_ephemeral(in);
    return tmp;
}

static MATRIX *VC_GEE_toeplitz(in)
MATRIX *in;
{
    MATRIX *toep, *tin, *tmp;
    int n, p, inrows, incols, i, j;

    inrows = in->nrows;
    incols = in->ncols;

    if ((inrows > incols) ? inrows % incols : incols % inrows)
    {
	error("M+-:VC_GEE_toeplitz: argument invalid");
    }

    if (inrows > incols)
    {
	p = incols;
	n = inrows/p;
	tin = VC_GEE_matcopy(in);
	free_if_ephemeral(in);
    }
    else
    {
	p = inrows;
	n = incols/p;
	tin = VC_GEE_transp(in);
    }

    toep = VC_GEE_create_matrix(n*p, n*p, EPHEMERAL);

    for (i = 0 ; i < n ; i ++)
    {
	tmp = VC_GEE_extract_rows(tin, i*p, (i*p)+p-1);
	make_permanent(tmp);
	if (i == 0)
	{
	    for (j = 0 ; j < n ; j++)
	    {
		if (inrows > incols)
		    VC_GEE_plug(tmp, toep, j*p, j*p);
		else
		    VC_GEE_plug(VC_GEE_transp(tmp), toep, j*p, j*p);
	    }
	}
	else
	{
	    for (j = 0 ; j < n-i ; j++)
	    {
		VC_GEE_plug(VC_GEE_transp(tmp), toep, j*p, (j+i)*p);
		VC_GEE_plug(tmp, toep, (j+i)*p, j*p);
	    }
	}
	VC_GEE_destroy_matrix(tmp);
    }
    VC_GEE_destroy_matrix(tin);
    return toep;
}





#define get_nelem(x) (((x)->nrows) * ((x)->ncols))

static double VC_GEE_elsum(x)
MATRIX *x;
{
    double t=0.;
    double *loc;
    int i, nelem;

    nelem = get_nelem(x);
    loc = x->data;
    for (i = 0 ; i < nelem ; i++)
	t += *(loc++);
    if (is_ephemeral(x)) VC_GEE_destroy_matrix(x);
    return t;
}

static MATRIX *VC_GEE_matabs(x)
MATRIX *x;
{
    double *load, *look;
    MATRIX *tmp;
    int nelem, i;

    nelem = get_nelem(x);
    tmp = VC_GEE_create_matrix(x->nrows, x->ncols, EPHEMERAL);
    load = tmp->data;
    look = x->data;
    for (i = 0 ; i < nelem ; i++)
	*(load++) = fabs(*look++);
    free_if_ephemeral(x);
    return tmp ;
}

static double VC_GEE_matmax(x)
MATRIX *x;
{
    double t;
    double *loc;
    int i, nelem;

    nelem = get_nelem(x);
    loc = x->data;
    t = MEL(x,0,0);
    for (i = 0 ; i < nelem ; i++)
    {
	if (*(loc) > t) t = *(loc);
	loc++;
    }
    free_if_ephemeral(x);
    return t;
}


static MATRIX *VC_GEE_matexp(x)
MATRIX *x;
{
    double *load, *look;
    MATRIX *tmp;
    int nelem, i;

    nelem = get_nelem(x);
    tmp = VC_GEE_create_matrix(x->nrows, x->ncols, EPHEMERAL);
    load = tmp->data;
    look = x->data;
    for (i = 0 ; i < nelem ; i++)
	*(load++) = exp(*look++);
    free_if_ephemeral(x);
    return tmp ;
}



static MATRIX *VC_GEE_matadd(mat1, mat2)
MATRIX *mat1, *mat2;
{
    MATRIX *result;
    double *mat1base, *mat2base, *resbase;
    int i, j;
    if ((mat1->ncols != mat2->ncols) || (mat1->nrows != mat2->nrows))
    {
	error("VC_GEE_matadd: args (%dx%d) + (%dx%d) don't conform.\nfatal error",
	      mat1->nrows, mat1->ncols, mat2->nrows, mat2->ncols);
    }
    result = VC_GEE_create_matrix(mat1->nrows, mat1->ncols, EPHEMERAL);
    resbase = result->data;
    mat1base = mat1->data;
    mat2base = mat2->data;
    for (j = 0 ; j < result->nrows ; j++)
    {
	for (i = 0 ; i < result->ncols ; i++)
	{
	    *resbase = *mat1base + *mat2base ;
	    resbase++ ; mat1base++ ; mat2base++ ;
	    /* *(resbase++) = *(mat1base++) + *(mat2base++); */
	}
    }
    if (is_ephemeral(mat1)) VC_GEE_destroy_matrix(mat1);
    if (is_ephemeral(mat2)) VC_GEE_destroy_matrix(mat2);
    return result;
}

static MATRIX *VC_GEE_matsub(mat1, mat2)
MATRIX *mat1, *mat2;
{
    MATRIX *result;
    double *mat1base, *mat2base, *resbase;
    int i, j;
    if ((mat1->ncols != mat2->ncols) || (mat1->nrows != mat2->nrows))
    {
	error("VC_GEE_matsub: args (%dx%d) + (%dx%d) don't conform.\n",
	      mat1->nrows, mat1->ncols, mat2->nrows, mat2->ncols);
    }
    result = VC_GEE_create_matrix(mat1->nrows, mat1->ncols, EPHEMERAL);
    resbase = result->data;
    mat1base = mat1->data;
    mat2base = mat2->data;
    for (j = 0 ; j < result->nrows ; j++)
    {
	for (i = 0 ; i < result->ncols ; i++)
	{
	    *resbase = *mat1base - *mat2base ;
	    resbase++ ; mat1base++ ; mat2base++ ;
	    /* *(resbase++) = *(mat1base++) - *(mat2base++); */
	}
    }
    if (is_ephemeral(mat1)) VC_GEE_destroy_matrix(mat1);
    if (is_ephemeral(mat2)) VC_GEE_destroy_matrix(mat2);
    return result;
}

static MATRIX *VC_GEE_matmult(mat1, mat2)
MATRIX *mat1, *mat2;
{
    double *mat1base, *mat1loc, *mat2base, *mat2loc, *resbase;
    MATRIX *result;
    int i, rows, j;

    if (mat1->ncols != mat2->nrows)
    {
	error("VC_GEE_matmult: args (%dx%d) * (%dx%d) don't conform.\n",
	      mat1->nrows, mat1->ncols, mat2->nrows, mat2->ncols);
    }

    result = VC_GEE_create_matrix(mat1->nrows, mat2->ncols, EPHEMERAL);

    resbase = result->data;
    mat1base = mat1->data;
    mat2base = mat2->data;

    for (j = 0 ; j < result->nrows ; j++)
    {
	for (i = 0 ; i < result->ncols ; i++)
	{
	    mat1loc = mat1base;
	    mat2loc = mat2base;
	    for (rows = 0 ; rows < mat2->nrows ; rows++)
	    {
		*resbase += *(mat1loc++) * *mat2loc;
		mat2loc += mat2->ncols;
	    }
	    ++resbase;
	    ++mat2base;
	}
	mat1base += mat1->ncols;
	mat2base = mat2->data;
    }
    if (is_ephemeral(mat1)) VC_GEE_destroy_matrix(mat1);
    if (is_ephemeral(mat2)) VC_GEE_destroy_matrix(mat2);
    return result;
}



static MATRIX *VC_GEE_px1_times_pxq(px1, pxq) /* mult elements of a colvec */
				/* across corresp row of mat */
MATRIX *px1, *pxq;
{
    MATRIX *tmp;
    double *load, colel;
    int i, j;

    if (px1->ncols != 1)
    {
	error("M+-: VC_GEE_px1_times_pxq: arg1 not a col-vec");
    }
    if (px1->nrows != pxq->nrows)
    {
	error("M+-: VC_GEE_px1_times_pxq: args not conforming");
    }
    tmp = VC_GEE_matcopy(pxq);
    load = tmp->data;
    for (i = 0 ; i < tmp->nrows ; i++)
    {
	colel = MEL(px1, i, 0);
	for (j = 0 ; j < tmp->ncols ; j++)
	{
	    *load *= colel ;
	    load++ ;
	}
    }
    free_if_ephemeral(px1);
    free_if_ephemeral(pxq);
    return tmp;
}

static MATRIX *VC_GEE_pxq_divby_px1(pxq, px1) /* divide elements of a colvec */
				/* into corresp row of mat */
MATRIX *px1, *pxq;
{
    MATRIX *tmp;
    double *load, colel;
    int i, j;
    if (px1->ncols != 1)
    {
	error("M+-: VC_GEE_pxq_divby_px1: arg2 not a col-vec");
    }
    if (px1->nrows != pxq->nrows)
    {
	error("M+-: VC_GEE_pxq_divby_px1: args not conforming");
    }

    tmp = VC_GEE_matcopy(pxq);
    load = tmp->data;
    for ( i = 0 ; i < tmp->nrows ; i++)
    {
	colel = MEL(px1, i, 0);
	for (j = 0 ; j < tmp->ncols ; j++)
	{
	    *load = (*load) / colel ;
	    load++ ;
	}
    }
    free_if_ephemeral(px1);
    free_if_ephemeral(pxq);
    return tmp;
}

static MATRIX *VC_GEE_scalar_times_matrix(a, X)
double a;
MATRIX *X;
{
    MATRIX *tmp;
    double *tbase;
    int i, nelem;
    tmp = VC_GEE_matcopy(X);
    nelem = get_nelem(tmp);
    tbase = tmp->data;
    for (i = 0 ; i < nelem ; i++) {
	*tbase *= a ;
	tbase++ ;
    }
    free_if_ephemeral(X);
    return tmp;
}





static void VC_GEE_matdump(mat)
MATRIX *mat;
{
    double *curel;
    int outtok = 0;
    int nel;

    nel = mat->nrows * mat->ncols;

    for (curel = mat->data ;  curel < mat->data + nel ; curel++)
    {
	Rprintf( ((fabs(*curel)<.00001) && (fabs(*curel)>0.)) ? "%.4e%c" : "%.4f%c", *curel,
		 (outtok++%mat->ncols == mat->ncols-1) ? '\n' : ' ');
    }
/* DOES NOT CLEAN */
}

static MATRIX *VC_GEE_luinv(X)
MATRIX *X;
{
    MATRIX *Y;
    double det[2], *work;
    int job[1];

    double *y;
    int info[1], *ipvt;
    int nrows, ncols;

    det[0] = 0.;
    det[1] = 0.;

    info[0] = 0;
    job[0] = 0;

    Y = VC_GEE_matcopy(X);  /* inversion in situ */

    nrows = X->nrows;
    ncols = X->ncols;

    /* really R_alloc */
    ipvt = (int *)malloc((unsigned)nrows*sizeof(int));
    work = (double *)malloc((unsigned)nrows*sizeof(double));

    y = (double *)Y->data;
    F77_NAME(dgefa)(y,&nrows,&ncols,ipvt,info);


    job[0] = 11;
    F77_NAME(dgedi)(y,&nrows,&ncols,ipvt,det,work,job);

    /* free(ipvt);
       free(work); */
    free_if_ephemeral(X);

    return Y;
}

static MATRIX *VC_GEE_covlag(inmat, lag, demean)
MATRIX *inmat;
int lag, demean;
{
    MATRIX *xrows[MAX_COVLAG], *res, *temp;
    int n, i, j, nv, q;
    double nrec;

    n = inmat->nrows;
    nrec = (double)1./(double)n;
    if (n > MAX_COVLAG)
    {
	error("VC_GEE_covlag: arg has > MAX_COVLAG rows");
    }

    nv = inmat->ncols;

    res = VC_GEE_create_matrix(nv, lag*nv, EPHEMERAL);

    for (q = 0 ; q < n ; q++)
    {
	xrows[q] = VC_GEE_extract_rows(inmat, q, q);
	make_permanent(xrows[q]);
    }


    for (i = 0 ; i < lag ; i++)
    {
	temp = VC_GEE_create_matrix(nv, nv, EPHEMERAL);
	for (j = i ; j < n ; j++)
	{
	    if ((j-i) < n) temp = VC_GEE_matadd(temp,
						   VC_GEE_matmult(VC_GEE_transp(xrows[j]),xrows[j-i]));
	}
	VC_GEE_plug(VC_GEE_scalar_times_matrix(nrec, temp), res, 0, i*nv);
    }

    for (q = 0 ; q < n ; q++)
    {
	VC_GEE_destroy_matrix(xrows[q]);
    }
    return res;
}




static MATRIX *VC_GEE_ident(ord)
int ord;
{
    MATRIX *I;
    int i;

    I = VC_GEE_create_matrix(ord, ord, EPHEMERAL);
    for (i = 0 ; i < ord ; i++)
	*(ELREF(I,i,i)) = (double)1.0;
    return I;
}

static MATRIX *VC_GEE_col_1s(k)
int k;
{
    MATRIX *tmp;
    int i;
    tmp = VC_GEE_create_matrix(k, 1, EPHEMERAL);
    for (i = 0 ; i < k ; i++)
    {
	MEL(tmp,i,0) = 1.;
    }
    return tmp;
}


static int VC_GEE_nchanges(X)
MATRIX *X;
{
/* returns integer telling how often the value of X */
/* changes from row to row.  X must be column VC_GEE_vector */

    int tmp = 1, iVC_GEE_start, i;

    if (X->ncols != 1)
    {
	error("VC_GEE_nchanges:  must be column VC_GEE_vector; ncols = %d",
	      X->ncols);
    }

    iVC_GEE_start = MEL(X, 0, 0);

    for (i = 1 ; i < X->nrows ; i++)
    {
	if (MEL (X, i, 0) != iVC_GEE_start)
	{
	    tmp++;
	    iVC_GEE_start = MEL (X, i, 0);
	}
    }
    return tmp;
}

static  MATRIX *VC_GEE_matanticlog(x)
MATRIX *x;
{
    double *load, *look;
    double exp();
    MATRIX *tmp;
    int nelem, i;

    nelem = get_nelem(x);
    tmp = VC_GEE_create_matrix(x->nrows, x->ncols, EPHEMERAL);
    load = tmp->data;
    look = x->data;
    for (i = 0 ; i < nelem ; i++)
	*(load++) = 1-exp(-exp(*look++));
    free_if_ephemeral(x);
    return tmp ;
}



/* gee support @(#) VC_GEE_diag_as_vec.c 3.3 94/03/09 */


static MATRIX *VC_GEE_diag_as_vec(inmat)
MATRIX *inmat;
{
    int i;
    MATRIX *outmat;

    if(inmat->ncols!=inmat->nrows)
    {
	error("M+-: VC_GEE_diag_as_vec: arg is not a square matrix");
    }

    outmat= VC_GEE_create_matrix(inmat->nrows,1,EPHEMERAL);
    for(i= 0;i<inmat->nrows;i++)
    {
	*(ELREF(outmat,i,0))=  *(ELREF(inmat,i,i));
    }
    return outmat;
}



static MATRIX *VC_GEE_matsqrt(x)
MATRIX *x;
{
    int i,j;
    MATRIX *tmp;
    tmp= VC_GEE_matcopy(x);
    for(i= 0;i<x->ncols;i++)
    {
	for(j= 0;j<x->nrows;j++)
	{
	    MEL(tmp,i,j)= sqrt(MEL(x,i,j));
	}
    }
    if(is_ephemeral(x))VC_GEE_destroy_matrix(x);
    return tmp;
}

static MATRIX *VC_GEE_mat1over(x)
MATRIX *x;
{
    int i,j;
    MATRIX *tmp;
    tmp = VC_GEE_matcopy(x);
    for(i=0;i<x->ncols;i++)
    {
	for(j=0;j<x->nrows;j++)
	{
	    MEL(tmp,i,j)= 1./(MEL(x,i,j));
	}
    }
    if(is_ephemeral(x))VC_GEE_destroy_matrix(x);
    return tmp;
}


/* gee support @(#) normlib.c 4.12 98/01/26 */

/* normlib.c */
/*-------------------------------------------------*/

/* library  of routines to support VC_GEE_npdf() and VC_GEE_ncdf() functions to
   compute the normal pdf and cdf used in adding the probit link to gee
   code added by pj catalano                                               */

/* dependence on Numerical Recipes code removed by Thomas Lumley */
/* normal pdf/cdf now taken from internal R routines (more accurate, too) */

#include <Rmath.h>
#define dnorm1(x) dnorm(x, 0.0, 1.0, 0)
#define pnorm1(x) pnorm(x, 0.0, 1.0, 1, 0)


/* following added by pj catalano to support probit link option in cgee.c */

static MATRIX *VC_GEE_matnpdf(x)
MATRIX *x;
{
    double *load, *look;
    double VC_GEE_npdf();
    MATRIX *tmp;
    int nelem, i;

    nelem = get_nelem(x);
    tmp = VC_GEE_create_matrix(x->nrows, x->ncols, EPHEMERAL);
    load = tmp->data;
    look = x->data;
    for (i = 0 ; i < nelem ; i++)
	*(load++) = dnorm1(*look++);
    free_if_ephemeral(x);
    return tmp ;
}

static MATRIX *VC_GEE_matncdf(x)
MATRIX *x;
{
    double *load, *look;
    double VC_GEE_ncdf();
    MATRIX *tmp;
    int nelem, i;

    nelem = get_nelem(x);
    tmp = VC_GEE_create_matrix(x->nrows, x->ncols, EPHEMERAL);
    load = tmp->data;
    look = x->data;
    for (i = 0 ; i < nelem ; i++)
	*(load++) = pnorm1(*look++);
    free_if_ephemeral(x);
    return tmp ;
}

/* end of pj catalano additions  */

#include <R_ext/Rdynload.h>

static const R_CMethodDef CEntries[] = {
    {"Cgee", (DL_FUNC) &Cgee, 21},
   {NULL, NULL, 0}
};

#include <Rversion.h>
void
R_init_gee(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
#if defined(R_VERSION) && R_VERSION >= R_Version(2, 16, 0)
    R_forceSymbols(dll, TRUE);
#endif
}
