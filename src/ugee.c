/* ugee.c -- unified c source code for gee for Splus */
/* ugee.c -- ugee.c 4.13 98/01/27 */
/* /proj/stdevs/stdev0f/SLIBS/gee.dev/GEE97/SCCS/s.ugee.c */
/* combines sources from the following v 4.9 files */
/* cgee.c          chanmat.c       clinluxxy.c     diag_as_vec.c   normlib.c */

/* needed for printf etc in S-PLUS 4.0 and later */
#include <S.h>
#if defined(SPLUS_VERSION) && SPLUS_VERSION >= 4000
#  include <newredef.h>
#endif
#ifdef USING_R
#include <Rversion.h>
typedef int Sint;
#else
typedef long Sint;
#endif

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

void Cgee( x, y, id, n, offset, nobs, p,
	   parmvec, M_parm, S_beta, S_naivvar, S_robvar,
	   S_phi, S_R, tol, maxsz, S_iter, silent, errorstate , scale_fix ,
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
MATRIX *One, *mui; /* , *Opmui; */
MATRIX *Ai, *ei , *ete;
MATRIX *S1, *S2, *Di, *this_R, *S5, *S2i;
MATRIX *tempmat1, *tempmat2, *Aop, *Dop, *zi, *DRop;
MATRIX *robvar, *naivvar, *lag_wts, *tmpeep, *scratch, *wt;
double phi, dni, phiLZ;
int iter, ini, i2, j2, k;
int alpha_VC_GEE_bandwidth;
int *onep, one, nclust, i; /* , j; */
double alpha_scalar;
double alpha_scalar_LZ, exdiv_LZ;
int maxni, ni, *maxnip;
int link, var_mean_rel, corstruct;
int *maxiter;
double nnsclust = 0.;  /* for counting non-singletons */

maxiter = (int* ) malloc(sizeof(int));
*errorstate = UNKNOWN_FAILURE ;
errorbranch( Bail ); 

if (!(*silent)) printf("@(#) ugee.c 98/01/26 Cgee: GEE C source version chanlib 4.12 \n");

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

from_S( y , nobs , onep , yin )
from_S( x , nobs , p , xin )
from_S( id , nobs , onep , idin )
from_S( offset, nobs, onep, offin )
from_S( n, nobs, onep, nin )
from_S( S_beta, p, onep, beta );

nclust = VC_GEE_nchanges( idin );

#define set_matrix_array( arrname, nel ) \
   if (!( arrname = (MATRIX **)malloc((unsigned)(nel * sizeof(MATRIX *))))) \
     { fprintf(stderr, \
            "set_matrix_array (mac): out of memory, requesting %d elements\n",nel); \
     exit(1); \
     }

set_matrix_array( X, nclust )
set_matrix_array( Y, nclust )
set_matrix_array( OFFSET, nclust )
set_matrix_array( N, nclust )

VC_GEE_split( xin , idin , X );
VC_GEE_split( yin , idin, Y );
VC_GEE_split( offin , idin , OFFSET );
VC_GEE_split( nin , idin , N );

VC_GEE_destroy_matrix( xin );
VC_GEE_destroy_matrix( yin );
VC_GEE_destroy_matrix( idin );
VC_GEE_destroy_matrix( offin );
VC_GEE_destroy_matrix( nin );

maxni = Y[0]->nrows;
for ( i = 1 ; i < nclust ; i++ )
	{
	ni = Y[i]->nrows;
	if ( ni > maxni ) maxni = ni;
	}
*maxsz = maxni;
maxnip = &maxni;

/* xpx = VC_GEE_create_matrix( *p, *p , EPHEMERAL );
xpy = VC_GEE_create_matrix( *p, one , EPHEMERAL ); */

if ( corstruct == (int) fixed ) 
	{
	from_S( S_R, maxnip, maxnip , R );
	}

if (!(*silent))
	{
	printf("Cgee will use: ");
	switch( link )
		{
		case VC_GEE_identity:
			printf("VC_GEE_identity link, ");
			break;
		case logarithm:
			printf("log link, ");
			break;
		case logit:
			printf("logit link, ");
			break;
		case reciprocal:
			printf("recip link, ");
			break;
		case probit:
			printf("probit link, ");
			break;
                case cloglog:
                        printf("cloglog link, ");
                        break;
		default:
			printf("unknown link.\n");
			exit( UNKNOWN_LINK );
			break;
		}
	switch( var_mean_rel )
		{
		case Gaussian:
			printf("Gaussian var, ");
			break;
		case Poisson:
			printf("Poisson var, ");
			break;
		case Binomial:
			printf("Binomial var, ");
			break;
		case Gamma:
			printf("Gamma var, ");
			break;
		default:
			fprintf( stderr,
			"Cgee: unknown var_mean_rel. Dies.\n");
			exit( UNKNOWN_VAR_MEAN_REL );
			break;
		}
	switch ( corstruct )
		{
		case independence:
			printf("indep corstr . ");
			break;
		case exchangeable:
			printf("exch corstr . ");
			break;
		case stat_M_dep:
			printf("stat corstr . ");
			break;
		case AR_M:
			printf("AR-M corstr . ");
			break;
		case non_stat_M_dep:
			printf("non-statM corstr . ");
			break;
		case unstructured:
			printf("unstr corstr . ");
			break;
		case fixed:
			printf("fixed corstr . ");
			break;
		default:
			printf("unknown corstr. ");
			exit(1);
			break;
		}
	printf("\n");
	}


scratch = VC_GEE_create_matrix( maxni, maxni, EPHEMERAL );

do
	{
	betasave = VC_GEE_matcopy(beta);
	phi = 0.;
	phiLZ = 0.;
	nnsclust = 0.;

	tmpR = VC_GEE_create_matrix( maxni, maxni, EPHEMERAL );


	switch ( corstruct )
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
			alpha = VC_GEE_create_matrix( 1 , alpha_VC_GEE_bandwidth , EPHEMERAL);
			break;
		case non_stat_M_dep:
		case unstructured:
			alpha = VC_GEE_create_matrix( maxni, maxni, EPHEMERAL );
			break;
		case fixed:
			break;
		default:
			break;
		}

	make_permanent( beta );
	for ( i = 0 ; i < nclust ; i++ ) 
		{
		ni = Y[i]->nrows;
		dni = (double)ni;
		One = VC_GEE_col_1s( ni );
		make_permanent( One );
		tempmat1 = VC_GEE_matadd( VC_GEE_matmult( X[i] , beta ), OFFSET[i] );

		switch ( link )
			{
			double maxfitted;
			case VC_GEE_identity:
				mui = VC_GEE_matcopy( tempmat1 );
				break;
			case logarithm:
				mui = VC_GEE_matexp( tempmat1 );
				break;
			case logit:
				tempmat1 = VC_GEE_matexp( tempmat1 );
				make_permanent( tempmat1 );
				tempmat2 = VC_GEE_matadd( One, tempmat1);
				mui = VC_GEE_pxq_divby_px1( VC_GEE_px1_times_pxq( N[i],
					tempmat1), tempmat2);
				make_permanent(mui);
				maxfitted = VC_GEE_matmax(VC_GEE_pxq_divby_px1(mui,N[i]));
				if ( (maxfitted >= .9999 ) && var_mean_rel == (int) Binomial ) 
					{
					fprintf(stderr,"Cgee: error: logistic model for probability has fitted value very close to 1.\n");
					fprintf(stderr,"Cgee: estimates diverging; iteration terminated.\n");
					fprintf(stderr,"Cgee: returning to S\n");
					return;
					}  
				VC_GEE_destroy_matrix(tempmat1);
				break;
			case reciprocal:
				mui = VC_GEE_pxq_divby_px1(
					One, tempmat1 );
				break;
/* probit case added by pj catalano*/  /* OK */
                        case probit:
                                tempmat1 = VC_GEE_matncdf(tempmat1);
                                mui = VC_GEE_px1_times_pxq( N[i], tempmat1);
                                make_permanent(mui);
                                maxfitted = VC_GEE_matmax(VC_GEE_pxq_divby_px1(mui, N[i]));
                                if ( (maxfitted >= .9999 ) && var_mean_rel == Binomial)
                                        {
		fprintf(stderr,"Cgee: estimates diverging; iteration terminated.\n");
                                        fprintf(stderr,"Cgee: returning to S\n");
                                        *errorstate = LOGISTIC_DIVERGENCE;
                                        return;
                                        }
                                 break;

 /* cloglog case added by jh maindonald */
                        case cloglog:
                                tempmat1=VC_GEE_matanticlog(tempmat1);
                                mui=VC_GEE_px1_times_pxq(N[i], tempmat1);
                                make_permanent(mui);
                                maxfitted = VC_GEE_matmax(VC_GEE_pxq_divby_px1(mui, N[i]));
                                if ( (maxfitted >= .999999 ) && var_mean_rel == (int) Binomial )
                                        {
                                        fprintf(stderr,"Cgee: error: cloglog model for probability has fit ted value very close to 1.\n");
                                        fprintf(stderr,"Cgee: estimates diverging; iteration terminated.\n ");
                                        fprintf(stderr,"Cgee: returning to S\n");
                                        return;
                                        }
                                break;

			default:
				fprintf( stderr,
				"Cgee: unknown link. Dies.\n");
				exit( UNKNOWN_LINK );
				break;
			}
		make_permanent( mui );
		ei = VC_GEE_matsub( Y[i], mui );
		
		switch( var_mean_rel )
			{
			case Gaussian:
			/*	Ai = VC_GEE_ident(ni); */
				break;
			case Poisson:
				Ai = VC_GEE_form_diag( mui );
				break;
			case Binomial:
				Ai = VC_GEE_form_diag( VC_GEE_px1_times_pxq(
				mui, VC_GEE_matsub( One, VC_GEE_pxq_divby_px1(mui,N[i]))));
				break;
			case Gamma:
				Ai = VC_GEE_form_diag( VC_GEE_px1_times_pxq(
				mui, mui ));
				break;
			default:
				fprintf( stderr,
				"Cgee: unknown var_mean_rel. Dies.\n");
				exit( UNKNOWN_VAR_MEAN_REL );
				break;
			}
		VC_GEE_destroy_matrix(mui);

		if (var_mean_rel != Gaussian )
		ei = VC_GEE_px1_times_pxq( VC_GEE_mat1over(VC_GEE_matsqrt(VC_GEE_diag_as_vec(Ai))), ei );

		make_permanent( ei );

		ete = VC_GEE_matmult( VC_GEE_transp(ei), ei );

		phi += MEL(ete,0,0)/dni;
		phiLZ += MEL(ete,0,0);
		if (dni>1) nnsclust += 1.0;

		if ( corstruct == (int) stat_M_dep || corstruct == (int) AR_M )
			{
			if ( dni < (double)alpha_VC_GEE_bandwidth )
				{
fprintf(stderr,"cgee: M-dependence, M=%d, but clustsize=%d\n",(int)*M_parm,(int)dni);
fprintf(stderr,"cgee: M-dependence, fatal error for this model.\n");
exit(M_DEP_SIZE_LT_M);
				}
			
			lag_wts = VC_GEE_create_matrix( 1, alpha_VC_GEE_bandwidth, PERMANENT );
			MEL( lag_wts, 0, 0 ) = (double)1.;

			for ( ini = 1 ; ini < alpha_VC_GEE_bandwidth ; ini++ )
				{
				MEL(lag_wts, 0, ini) = dni/(dni - (double)ini);
				}
			}
		
		switch ( corstruct )
			{
			case independence:
				break;
			case exchangeable:
				if ( ni > 1 ) alpha_scalar +=
1./(dni*(dni-1.)) * (VC_GEE_elsum( VC_GEE_matmult( ei, VC_GEE_transp(ei) ) ) - MEL(ete,0,0));
				if ( ni > 1 ) alpha_scalar_LZ +=
(VC_GEE_elsum( VC_GEE_matmult( ei, VC_GEE_transp(ei) ) ) - MEL(ete,0,0));
				if ( ni > 1 ) exdiv_LZ += dni*(dni-1); /* suppress .5 because num
is redundant */
				break;
			case stat_M_dep:
			case AR_M:

alpha = VC_GEE_matadd( alpha, VC_GEE_transp( VC_GEE_px1_times_pxq ( VC_GEE_transp(lag_wts), VC_GEE_transp(VC_GEE_covlag( ei , alpha_VC_GEE_bandwidth , 0 ) ) ) ) );
				break;
			case non_stat_M_dep:
			case unstructured:
				tmpeep = VC_GEE_matmult( ei, VC_GEE_transp(ei) );
				VC_GEE_plug( tmpeep , scratch , 0, 0);
				make_permanent(scratch);
				alpha = VC_GEE_matadd( alpha, scratch );
				scratch = VC_GEE_matsub( scratch, scratch );
				make_ephemeral(scratch);
				break;
			case fixed:
				break;
			default:
			fprintf(stderr,"corstruct not implemented.\n");
				exit(1);
				break;
			}
		VC_GEE_destroy_matrix(ei);
		VC_GEE_destroy_matrix(One);
	      }

        if  (alpha != NULL ) alpha = VC_GEE_scalar_times_matrix( 
             (double)nclust/(phi*nnsclust) , alpha );

	alpha_scalar /= phi;
        alpha_scalar *= (double)nclust/nnsclust;

	phi /= (double)nclust;
	phiLZ /= ((double)(*nobs - *p));

	alpha_scalar_LZ /= (phiLZ*(exdiv_LZ-2.*(double)*p));
	if (*compatflag == 0) {  /* abandon compatibility with early macro */
		alpha_scalar = alpha_scalar_LZ;
		phi = phiLZ;
		}

	switch( corstruct )
		{
		case independence:  /* this is ridiculous */
			R = VC_GEE_ident( maxni );
			break;
		case exchangeable:
R = VC_GEE_matadd( VC_GEE_scalar_times_matrix( alpha_scalar , VC_GEE_matmult( VC_GEE_col_1s(maxni), VC_GEE_transp( VC_GEE_col_1s(maxni) ))), VC_GEE_scalar_times_matrix( (double)1.-alpha_scalar,VC_GEE_ident(maxni)));
			break;
		case stat_M_dep:
			tmpR = VC_GEE_create_matrix( 1, maxni, EPHEMERAL );
			VC_GEE_plug( alpha, tmpR, 0, 0 );
			MEL( tmpR, 0, 0 ) = (double)1.;
			R = VC_GEE_toeplitz( tmpR );
			break;
		case AR_M:
			tmpR = VC_GEE_create_matrix( 1, maxni, EPHEMERAL );
			MEL( alpha, 0, 0) = (double)1.;
			make_permanent(alpha);
			VC_GEE_plug( alpha, tmpR, 0,0);
			tempmat1 = VC_GEE_extract_cols( alpha, 0, alpha_VC_GEE_bandwidth-2 );
			wt = VC_GEE_toeplitz( tempmat1 );
			tempmat1 = VC_GEE_extract_cols( alpha, 1, alpha_VC_GEE_bandwidth-1 );
			alpha = VC_GEE_matmult( tempmat1, VC_GEE_luinv( wt ) );
			for ( i2 = alpha_VC_GEE_bandwidth ; i2 < maxni ; i2++ )
				{
				for ( j2 = 0 ; j2 < alpha_VC_GEE_bandwidth-1 ; j2++ )
					{
					MEL( tmpR, 0, i2 ) += MEL(alpha,0,j2) *
					MEL(tmpR,0,(i2-j2)-1);
					}
				}
			R = VC_GEE_toeplitz( tmpR );
			VC_GEE_destroy_matrix(alpha);
			break;
		case non_stat_M_dep:
			R = VC_GEE_band( alpha, alpha_VC_GEE_bandwidth );
			for ( k = 0 ; k < R->ncols ; k++ )
				{
				MEL( R, k, k) = (double)1.;
				}
			break;
		case unstructured:
			R = VC_GEE_matcopy( alpha );
			for ( k = 0 ; k < R->ncols ; k++ )
				{
				MEL( R, k, k) = (double)1.;
				}
			break;
		case fixed:
			break;
		default:
			fprintf(stderr,"corstruct not implemented.\n");
			exit(1);
			break;
		}
	make_permanent(R);

	S1 = VC_GEE_create_matrix( (int)*p, 1, EPHEMERAL);
	S2 = VC_GEE_create_matrix( (int)*p, (int)*p, EPHEMERAL);

	for ( i = 0 ; i < nclust ; i++ )
		{
		ni = Y[i]->nrows;
		dni = (double)ni;
		One = VC_GEE_col_1s(ni);
		make_permanent( One );

		tempmat1 = VC_GEE_matadd( VC_GEE_matmult( X[i], beta ), OFFSET[i] );
		switch ( link )
			{
			case VC_GEE_identity:
				Di = VC_GEE_matcopy(X[i]);
				mui = VC_GEE_matcopy( tempmat1 );
				break;
			case logarithm:
				tempmat1 = VC_GEE_matexp( tempmat1 );
				mui = VC_GEE_matcopy( tempmat1 );
				Di = VC_GEE_px1_times_pxq( tempmat1, X[i]);
				break;
			case logit:
				tempmat1 = VC_GEE_matexp( tempmat1 );
				make_permanent( tempmat1 );
				tempmat2 = VC_GEE_matadd( One, tempmat1);
				tempmat2 = VC_GEE_pxq_divby_px1(
					VC_GEE_px1_times_pxq( N[i], tempmat1), tempmat2 );
				make_permanent( tempmat2 );
				VC_GEE_destroy_matrix( tempmat1 );
				tempmat1 = VC_GEE_matsub( One, VC_GEE_pxq_divby_px1(tempmat2,N[i]));
				tempmat1 = VC_GEE_px1_times_pxq( tempmat1, 
					tempmat2 );
				Di = VC_GEE_px1_times_pxq( tempmat1, X[i]);
				mui = VC_GEE_matcopy( tempmat2 );
				VC_GEE_destroy_matrix( tempmat2 );
				break;
			case reciprocal:
				tempmat1 = VC_GEE_pxq_divby_px1(
					One, tempmat1 );
				mui = VC_GEE_matcopy(tempmat1);
				tempmat2 = VC_GEE_matcopy(tempmat1);
				tempmat1 = VC_GEE_px1_times_pxq( tempmat1, tempmat2 );
				tempmat1 = VC_GEE_px1_times_pxq(
				tempmat1, X[i] );
				Di = VC_GEE_scalar_times_matrix( -1., tempmat1 );
				break;
 /* probit case added by pj catalano */
                        case probit:
                                tempmat2 = VC_GEE_matcopy(tempmat1);
                                tempmat2 = VC_GEE_matncdf(tempmat2);
                                mui = VC_GEE_px1_times_pxq( N[i],tempmat2 );
                                tempmat1 = VC_GEE_matnpdf(tempmat1);
                                tempmat2 = VC_GEE_px1_times_pxq( N[i],tempmat1 );
                                Di = VC_GEE_px1_times_pxq( tempmat2, X[i] );
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
				fprintf( stderr,
				"Cgee: unknown link. Dies.\n");
				exit( UNKNOWN_LINK );
				break;
			}
		make_permanent( mui );
		make_permanent( Di );

		zi = VC_GEE_matadd( VC_GEE_matmult( Di, beta), VC_GEE_matsub( Y[i],
		mui ) );

		switch( var_mean_rel )
			{
			case Gaussian:
				/* Ai = VC_GEE_ident(ni); */
				break;
			case Poisson:
				Ai = VC_GEE_form_diag( mui );
				break;
			case Binomial:
				Ai = VC_GEE_form_diag( VC_GEE_px1_times_pxq(
				mui, VC_GEE_matsub( One, VC_GEE_pxq_divby_px1(mui,N[i]))));
				break;
			case Gamma:
				Ai = VC_GEE_form_diag( VC_GEE_px1_times_pxq(
				mui, mui ));
				break;
			default:
				fprintf( stderr,
				"Cgee: unknown var_mean_rel. Dies.\n");
				exit( UNKNOWN_VAR_MEAN_REL );
				break;
			}
		VC_GEE_destroy_matrix(mui);

		if ( var_mean_rel != Gaussian )  /* else Ai is VC_GEE_identity */
		Aop = VC_GEE_mat1over( VC_GEE_matsqrt( VC_GEE_diag_as_vec( Ai )) );
		if ( var_mean_rel != Gaussian )  /* else Ai is VC_GEE_identity */
		make_permanent( Aop );


		make_ephemeral( Di );
		if ( var_mean_rel != Gaussian )
		  Di = VC_GEE_px1_times_pxq( Aop, Di );
		make_permanent( Di );
		if ( var_mean_rel != Gaussian )
		  zi = VC_GEE_px1_times_pxq( Aop, zi );

		if ( corstruct != independence )
		{
		this_R = VC_GEE_corner( R, ni, ni );
        	/* printf("cluster %d g\n",i); */

		/* printf("VC_GEE_start inver\n");
		printf("D: %d x %d, R: %d x %d\n",Di->nrows,
			   Di->ncols, this_R->nrows, this_R->ncols); */
		Ri = VC_GEE_luinv(this_R);
		/* printf("D: %d x %d, Ri: %d x %d\n",Di->nrows,
			   Di->ncols, Ri->nrows, Ri->ncols); */
		Dop = VC_GEE_matmult( VC_GEE_transp(Di), Ri );
		/* printf("end inver\n"); */
		}
		else Dop = VC_GEE_transp(Di);

		make_permanent( Dop );

		S1 = VC_GEE_matadd( S1, VC_GEE_matmult( Dop, zi ) );
		S2 = VC_GEE_matadd( S2, VC_GEE_matmult( Dop, Di ) );

		VC_GEE_destroy_matrix(Dop);
		if ( var_mean_rel != Gaussian )
		VC_GEE_destroy_matrix(Aop);
		VC_GEE_destroy_matrix(Di);
		VC_GEE_destroy_matrix(One);

		}

	VC_GEE_destroy_matrix(beta);
	beta = VC_GEE_matmult( VC_GEE_luinv( S2 ), S1 );
	make_permanent(beta);

	if (!(*silent)) printf("current parameter estimates:\n");

	if (!(*silent)) VC_GEE_matdump(beta);

	if (!(*silent)) printf( "completed iteration %d\n",iter);
	iter++;


		
	} while(( VC_GEE_matmax( VC_GEE_matabs( VC_GEE_matsub( VC_GEE_pxq_divby_px1( betasave, beta ), VC_GEE_col_1s((int)*p)))) > *tol) && (iter < *maxiter ));

if ( iter >= *maxiter ) 
	{
	fprintf(stderr,"Maximum number of iterations consumed.\n");
	fprintf(stderr,"Convergence not achieved; results suspect.\n");
	*errorstate = MAXITER_EXCEEDED ;
	}

S2 = VC_GEE_create_matrix( (int)*p, (int)*p, EPHEMERAL);
S5 = VC_GEE_create_matrix( (int)*p, (int)*p, EPHEMERAL);
phi = 0.;
phiLZ = 0.;

for ( i = 0 ; i < nclust ; i++ )
	{
	ni = Y[i]->nrows;
	dni = (double)ni;
	One = VC_GEE_col_1s(ni);
	make_permanent( One );

	tempmat1 = VC_GEE_matadd( VC_GEE_matmult( X[i], beta ), OFFSET[i] );
	switch ( link )
		{
		case VC_GEE_identity:
			Di = VC_GEE_matcopy(X[i]);
			mui = VC_GEE_matcopy( tempmat1 );
			break;
		case logarithm:
			tempmat1 = VC_GEE_matexp( tempmat1 );
			mui = VC_GEE_matcopy( tempmat1 );
			Di = VC_GEE_px1_times_pxq( tempmat1, X[i]);
			break;
		case logit:
				tempmat1 = VC_GEE_matexp( tempmat1 );
				make_permanent( tempmat1 );
				tempmat2 = VC_GEE_matadd( One, tempmat1);
				tempmat2 = VC_GEE_pxq_divby_px1(
					VC_GEE_px1_times_pxq( N[i], tempmat1), tempmat2 );
				make_permanent( tempmat2 );
				VC_GEE_destroy_matrix( tempmat1 );
				tempmat1 = VC_GEE_matsub( One, VC_GEE_pxq_divby_px1(tempmat2,N[i]));
				tempmat1 = VC_GEE_px1_times_pxq( tempmat1, 
					tempmat2 );
				Di = VC_GEE_px1_times_pxq( tempmat1, X[i]);
				mui = VC_GEE_matcopy( tempmat2 );
				VC_GEE_destroy_matrix( tempmat2 );
				break;
		case reciprocal:
			tempmat1 = VC_GEE_pxq_divby_px1(
				One, tempmat1 );
			mui = VC_GEE_matcopy(tempmat1);
			tempmat2 = VC_GEE_matcopy(tempmat1);
			tempmat1 = VC_GEE_px1_times_pxq( tempmat1, tempmat2 );
			tempmat1 = VC_GEE_px1_times_pxq(
			tempmat1, X[i] );
			Di = VC_GEE_scalar_times_matrix( -1., tempmat1 );
			break;
 /* probit case added by pj catalano*/
                case probit:
                        tempmat2 = VC_GEE_matcopy(tempmat1);
                        tempmat2 = VC_GEE_matncdf(tempmat2);
                        mui = VC_GEE_px1_times_pxq( N[i],tempmat2 );
                        tempmat1 = VC_GEE_matnpdf(tempmat1);
                        tempmat2 = VC_GEE_px1_times_pxq( N[i],tempmat1 );
                        Di = VC_GEE_px1_times_pxq( tempmat2, X[i] );
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
			fprintf( stderr,
			"Cgee: unknown link. Dies.\n");
			exit( UNKNOWN_LINK );
			break;
		}
	make_permanent( mui );

	ei = VC_GEE_matsub( Y[i] ,mui );
	switch( var_mean_rel )
		{
		case Gaussian:
			/* Ai = VC_GEE_ident(ni); */
			break;
		case Poisson:
			Ai = VC_GEE_form_diag( mui );
			break;
		case Binomial:
			Ai = VC_GEE_form_diag( VC_GEE_px1_times_pxq(
			mui, VC_GEE_matsub( One, VC_GEE_pxq_divby_px1(mui,N[i]))));
			break;
		case Gamma:
			Ai = VC_GEE_form_diag( VC_GEE_px1_times_pxq(
			mui, mui ));
			break;
		default:
			fprintf( stderr,
			"Cgee: unknown var_mean_rel. Dies.\n");
			exit( UNKNOWN_VAR_MEAN_REL );
			break;
		}
	VC_GEE_destroy_matrix(mui);
	VC_GEE_destroy_matrix(One);

	if ( var_mean_rel != Gaussian )
	  {
	    Aop = VC_GEE_mat1over( VC_GEE_matsqrt( VC_GEE_diag_as_vec( Ai ) ) );
	    make_permanent( Aop );
	    Di = VC_GEE_px1_times_pxq( Aop, Di );
	  }
	make_permanent( Di );

	if ( var_mean_rel != Gaussian )
	  ei = VC_GEE_px1_times_pxq( Aop, ei );
	make_permanent( ei );

	ete = VC_GEE_matmult( VC_GEE_transp(ei), ei);
	phi += MEL(ete,0,0) / dni;
	phiLZ += MEL(ete,0,0) ;

	if (corstruct != independence)
	DRop = VC_GEE_matmult( VC_GEE_transp(Di), VC_GEE_luinv( VC_GEE_corner( R, ni, ni ) ) );
	else
	DRop = VC_GEE_transp(Di);

	make_permanent( DRop );

	S2 = VC_GEE_matadd( S2, VC_GEE_matmult( DRop , Di ) );
/* printf("VC_GEE_start S5\n"); */
	S5 = VC_GEE_matadd( S5, VC_GEE_matmult( VC_GEE_matmult( DRop , ei ),
			VC_GEE_matmult(VC_GEE_transp(ei), VC_GEE_transp(DRop))));
	/* printf("end S5\n"); */

	if ( var_mean_rel != Gaussian )
	VC_GEE_destroy_matrix( Aop );
	VC_GEE_destroy_matrix( Di );
	VC_GEE_destroy_matrix( ei );
	VC_GEE_destroy_matrix( DRop );
	}

phi /= (double)nclust;
phiLZ /= (double)(*nobs - *p);
if (*compatflag == 0) {  /* abandon compatibility with early macro */
phi = phiLZ;
}

S2i = VC_GEE_luinv(S2);
make_permanent( S2i );
if ( *scale_fix ) phi = *S_phi;
if ( *scale_fix && var_mean_rel == Gaussian )
	fprintf(stderr,"Scale parameter fixed at %f with Gaussian variance function.\n",
			      *S_phi);
naivvar = VC_GEE_scalar_times_matrix( phi, S2i );
robvar = VC_GEE_matmult( VC_GEE_matmult ( S2i, S5 ), S2i );

to_S( beta, S_beta )
to_S( naivvar, S_naivvar )
to_S( robvar, S_robvar )
to_S( R, S_R )

*S_phi = phi;
*S_iter = iter;

VC_GEE_destroy_matrix( beta );
VC_GEE_destroy_matrix( naivvar );
VC_GEE_destroy_matrix( robvar );
VC_GEE_destroy_matrix( R );

for ( i = 0 ; i < nclust ; i++ )
	{
	VC_GEE_destroy_matrix( X[i] );
	VC_GEE_destroy_matrix( Y[i] );
	VC_GEE_destroy_matrix( N[i] );
	VC_GEE_destroy_matrix( OFFSET[i] );
	}

if ( iter < *maxiter ) *errorstate = NO_ERROR; 
Bail:
;
}



/* gee support @(#) chanmat.c 4.12 98/01/26 */

/* @(#) chanmat.nw 1.3 94/03/09 */

static MATRIX *VC_GEE_create_matrix( nrows, ncols, permanence )
/* $Y = |VC_GEE_create_matrix|(r,c,\cdot) :\Rightarrow Y \in M_{r\times c} \wedge
Y = 0 $ */
int nrows, ncols, permanence;
{
MATRIX *tmp;
double *head;
int i;

tmp = (MATRIX *) calloc ( 1, sizeof ( struct matrix ) );

if ( tmp == NULL )
	{
	fprintf( stderr , "VC_GEE_create_matrix: malloc attempt %d d.\n",
				sizeof( struct matrix ));
	Seterr_and_terminate( NO_MEM_MATSTRUCT );
	}

tmp->nrows = nrows;
tmp->ncols = ncols;
tmp->permanence = permanence;

tmp->data = ( double * ) calloc ( 1,  nrows * ncols * sizeof ( double ) ) ;

if ( tmp->data == NULL )
	{
	fprintf( stderr , "VC_GEE_create_matrix: malloc attempt %d d.\n",
				(unsigned)nrows*ncols);
	fprintf( stderr , "VC_GEE_create_matrix: nrows=%d ncols=%d.\n",
				nrows, ncols );
	Seterr_and_terminate( NO_MEM_MATDATA );
	}

head = tmp->data;
for ( i = 0 ; i < nrows*ncols ; i++ )
	{
	*head = (double)0.;
	head++;
	}

return tmp;
}

static void VC_GEE_destroy_matrix( mat )
MATRIX *mat;
{
if (mat == (MATRIX *) NULL ) return;
mat->nrows = 0;
mat->ncols = 0;
if (mat->data != (double *)NULL ) free((char *) mat->data );
mat->data = (double *)NULL;
if (mat != (MATRIX *)NULL) free((char *) mat );
mat = (MATRIX *)NULL;
}

 

static MATRIX *VC_GEE_transp( mat )
/*
$Y = VC_GEE_transp(X_{r\times c}) :\Rightarrow Y \in M_{c\times r} \wedge
y_{ji} = x_{ij}. $
*/
MATRIX *mat;
{
double *telem, *inelem, *tbase;
int nelem;
MATRIX *tmp;

tmp = VC_GEE_create_matrix( mat->ncols, mat->nrows , EPHEMERAL );
inelem = mat->data;
tbase = tmp->data;
telem = tbase;
for ( nelem = 0 ; nelem < ( mat->ncols * mat->nrows ) ; nelem++ )
	{
	*telem = *(inelem++);
	telem += mat->nrows;
	if ( nelem % mat->ncols == (mat->ncols)-1 )
		telem = ++tbase;
	}
if ( is_ephemeral( mat ) ) VC_GEE_destroy_matrix( mat );
return tmp;

}
 
static MATRIX *VC_GEE_corner( mat, nr, nc )
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
	fprintf( stderr, "VC_GEE_corner: request not a submatrix.\n");
	fprintf( stderr, "VC_GEE_corner: fatal error.\n");
	Seterr_and_terminate(CORNER_FAIL);
	}
tmp = VC_GEE_create_matrix( nr, nc, EPHEMERAL );
load = tmp->data;
for ( i = 0 ; i < nr ; i++ )
	{
	for ( j = 0 ; j < nc ; j++ )
		{
		*(load++) = MEL(mat, i, j);
		}
	}
free_if_ephemeral( mat );
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

for ( i = 0 ; i < rows_to_get ; i++ )
	{
	for ( j = 0 ; j < Source->ncols ; j++ )
		{
		*(ELREF(temp,i,j)) = *(ELREF(Source,VC_GEE_start,j));
		}
	VC_GEE_start++;
	}
/* DOES NOT CLEAN */
return temp;
}

static MATRIX *VC_GEE_extract_cols( x , VC_GEE_start , end )
MATRIX *x;
int VC_GEE_start, end;
{
MATRIX *tmp;
tmp = VC_GEE_transp(x);
tmp = VC_GEE_extract_rows( tmp, VC_GEE_start, end );
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
for ( i = 0 ; i < inmat->nrows ; i++ )
	{
	for ( j = 0 ; j < inmat->ncols ; j++ )
		{
		*(ELREF(outmat,i,j)) = *(ELREF(inmat,i,j));
		}
	}
/* DOES NOT CLEAN */
return outmat;
}
 
static int VC_GEE_split( matptr, discptr , matarrptr )
MATRIX *matptr, *discptr, *matarrptr[];
{   /* discriminator VC_GEE_vector assumed to be integer-valued dbls */
int i, iVC_GEE_start, k, VC_GEE_start, end;
if ( discptr->ncols != 1 )
	{
	fprintf(stderr,"VC_GEE_split: discriminator must be column vec.\n");
	fprintf(stderr,"VC_GEE_split: ncols = %d.\n", discptr->ncols);
	fprintf(stderr,"VC_GEE_split: fatal error.\n");
	Seterr_and_terminate(SPLIT_FAIL );
	}

k = 0;

iVC_GEE_start = (int)MEL( discptr , 0 , 0 );
VC_GEE_start = 0;
end = 0;
for ( i = 1 ; i <= discptr->nrows ; i++ )
	{
	if (i == discptr->nrows || MEL( discptr , i, 0 ) != iVC_GEE_start )
		{
		matarrptr[k] = VC_GEE_matcopy( VC_GEE_extract_rows( matptr , VC_GEE_start, end ) );
		make_permanent( matarrptr[k] );
		k++;
		VC_GEE_start = end+1;
		if (i < discptr->nrows) /* don't need iVC_GEE_start at end of loop */
			iVC_GEE_start = MEL( discptr, i, 0 );
		}
	if (VC_GEE_start < discptr->nrows ) end++ ;
	}
/* DOES NOT CLEAN */
return k;
}
 

static void VC_GEE_plug( VC_GEE_plugm, socket, row, col )
int row, col;
MATRIX *VC_GEE_plugm, *socket;  /* not a unix socket */
{
int pcol, prow;
double *sockload, *VC_GEE_plughead, *sockrow_VC_GEE_start;
int i,j;

pcol = VC_GEE_plugm->ncols;
prow = VC_GEE_plugm->nrows;

if ( pcol+col > socket->ncols || prow+row > socket->nrows )
	{
	fprintf( stderr,"M+-: VC_GEE_plug: socket too small. Dies.\n");
	Seterr_and_terminate(PLUG_FAIL);
	}

sockload = socket->data + col + row*(socket->ncols);
VC_GEE_plughead = VC_GEE_plugm->data;
sockrow_VC_GEE_start = sockload;

for ( i = 0 ; i < prow ; i++ )
	{
	sockload = sockrow_VC_GEE_start;
	for ( j = 0 ; j < pcol ; j++ )
		{
		*(sockload++) = *(VC_GEE_plughead++ );
		}
	sockrow_VC_GEE_start += socket->ncols;
	}
free_if_ephemeral( VC_GEE_plugm );
}
 
static MATRIX *VC_GEE_form_diag( vec )
MATRIX *vec;
{
MATRIX *tmp;
int i, ord;

ord = vec->nrows;
tmp = VC_GEE_create_matrix( ord,  ord, EPHEMERAL );
for ( i = 0 ; i < ord ; i++ )
	*(ELREF(tmp,i,i)) = MEL(vec,i,0);
free_if_ephemeral( vec );
return tmp;
}

static MATRIX *VC_GEE_band( in, wid )
MATRIX *in;
int wid;
{
MATRIX *tmp;
int i, j;
tmp = VC_GEE_matcopy( in );
for ( i = 0 ; i < in->nrows ; i++ )
	{
	for ( j = i+wid ; j < in->ncols ; j++ )
		{
		MEL( tmp, i, j ) = (double)0.;
		if (( i < in->ncols ) && ( j < in->nrows ))
			{
			MEL( tmp, j, i ) = (double)0.;
			}
		}
	}
free_if_ephemeral( in );
return tmp;
}

static MATRIX *VC_GEE_toeplitz( in )
MATRIX *in;
{
MATRIX *toep, *tin, *tmp;
int n, p, inrows, incols, i, j;

inrows = in->nrows;
incols = in->ncols;

if ( ( inrows > incols ) ? inrows % incols : incols % inrows )
	{
	fprintf(stderr,"M+-:VC_GEE_toeplitz: argument invalid. Dies.\n");
	Seterr_and_terminate(BAD_TOEPLITZ_ARG);
	}

if ( inrows > incols )
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

toep = VC_GEE_create_matrix( n*p, n*p, EPHEMERAL );

for ( i = 0 ; i < n ; i ++ )
	{
	tmp = VC_GEE_extract_rows( tin, i*p, (i*p)+p-1 );
	make_permanent(tmp);
	if ( i == 0 )
		{
		for ( j = 0 ; j < n ; j++ )
			{
			if ( inrows > incols )
				VC_GEE_plug( tmp, toep, j*p, j*p );
			else
				VC_GEE_plug( VC_GEE_transp(tmp), toep, j*p, j*p );
			}
		}
	else
		{
		for ( j = 0 ; j < n-i ; j++ )
			{
			VC_GEE_plug( VC_GEE_transp(tmp), toep, j*p, (j+i)*p );
			VC_GEE_plug( tmp, toep, (j+i)*p, j*p );
			}
		}
	VC_GEE_destroy_matrix(tmp);
	}
VC_GEE_destroy_matrix(tin);
return toep;
}


 
 

#define get_nelem( x ) (((x)->nrows) * ((x)->ncols))

static double VC_GEE_elsum( x )
MATRIX *x;
{
double t=0.;
double *loc;
int i, nelem;

nelem = get_nelem( x );
loc = x->data;
for ( i = 0 ; i < nelem ; i++ )
	t += *(loc++);
if ( is_ephemeral(x) ) VC_GEE_destroy_matrix(x);
return t;
}

static MATRIX *VC_GEE_matabs( x )
MATRIX *x;
{
double *load, *look;
MATRIX *tmp;
int nelem, i;

nelem = get_nelem( x );
tmp = VC_GEE_create_matrix( x->nrows, x->ncols , EPHEMERAL );
load = tmp->data;
look = x->data;
for ( i = 0 ; i < nelem ; i++ )
	*(load++) = fabs(*look++);
free_if_ephemeral(x);
return tmp ;
}

static double VC_GEE_matmax( x )
MATRIX *x;
{
double t;
double *loc;
int i, nelem;

nelem = get_nelem( x );
loc = x->data;
t = MEL(x,0,0);
for ( i = 0 ; i < nelem ; i++ )
	{
	if ( *(loc) > t ) t = *(loc);
	loc++;
	}
free_if_ephemeral( x );
return t;
}


static MATRIX *VC_GEE_matexp( x )
MATRIX *x;
{
double *load, *look;
MATRIX *tmp;
int nelem, i;

nelem = get_nelem( x );
tmp = VC_GEE_create_matrix( x->nrows, x->ncols , EPHEMERAL );
load = tmp->data;
look = x->data;
for ( i = 0 ; i < nelem ; i++ )
	*(load++) = exp(*look++);
free_if_ephemeral(x);
return tmp ;
}



static MATRIX *VC_GEE_matadd( mat1, mat2 )
MATRIX *mat1, *mat2;
{
MATRIX *result;
double *mat1base, *mat2base, *resbase;
int i, j;
if ( ( mat1->ncols != mat2->ncols ) || ( mat1->nrows != mat2->nrows ) )
	{
	fprintf(stderr,"VC_GEE_matadd: args (%dx%d) + (%dx%d) don't conform.\n",
		mat1->nrows, mat1->ncols, mat2->nrows,
		mat2->ncols );
	fprintf(stderr,"VC_GEE_matadd: fatal error.  exits. \n");
	Seterr_and_terminate(MATADD_NONCONFORMITY);
	}
result = VC_GEE_create_matrix( mat1->nrows , mat1->ncols, EPHEMERAL);
resbase = result->data;
mat1base = mat1->data;
mat2base = mat2->data;
for ( j = 0 ; j < result->nrows ; j++ )
	{
	for ( i = 0 ; i < result->ncols ; i++ )
		{
		*resbase = *mat1base + *mat2base ;
		resbase++ ; mat1base++ ; mat2base++ ;
		/* *(resbase++) = *(mat1base++) + *(mat2base++ ); */
		}
	}
if ( is_ephemeral( mat1 ) ) VC_GEE_destroy_matrix( mat1 );
if ( is_ephemeral( mat2 ) ) VC_GEE_destroy_matrix( mat2 );
return result;
}

static MATRIX *VC_GEE_matsub( mat1, mat2 )
MATRIX *mat1, *mat2;
{
MATRIX *result;
double *mat1base, *mat2base, *resbase;
int i, j;
if ( ( mat1->ncols != mat2->ncols ) || ( mat1->nrows != mat2->nrows ) )
	{
	fprintf(stderr,"VC_GEE_matsub: args (%dx%d) + (%dx%d) don't conform.\n",
		mat1->nrows, mat1->ncols, mat2->nrows,
		mat2->ncols );
	fprintf(stderr,"VC_GEE_matsub: fatal error.  exits. \n");
	Seterr_and_terminate(MATSUB_NONCONFORMITY);
	}
result = VC_GEE_create_matrix( mat1->nrows , mat1->ncols, EPHEMERAL);
resbase = result->data;
mat1base = mat1->data;
mat2base = mat2->data;
for ( j = 0 ; j < result->nrows ; j++ )
	{
	for ( i = 0 ; i < result->ncols ; i++ )
		{
		*resbase = *mat1base - *mat2base ;
		resbase++ ; mat1base++ ; mat2base++ ;
		/* *(resbase++) = *(mat1base++) - *(mat2base++ ); */
		}
	}
if ( is_ephemeral( mat1 ) ) VC_GEE_destroy_matrix( mat1 );
if ( is_ephemeral( mat2 ) ) VC_GEE_destroy_matrix( mat2 );
return result;
}

static MATRIX *VC_GEE_matmult( mat1, mat2 )
MATRIX *mat1, *mat2;
{
double *mat1base, *mat1loc, *mat2base, *mat2loc, *resbase;
MATRIX *result;
int i, rows, j;

if ( mat1->ncols != mat2->nrows )
	{
	fprintf(stderr,"VC_GEE_matmult: args (%dx%d) * (%dx%d) don't conform.\n",
		mat1->nrows, mat1->ncols, mat2->nrows,
		mat2->ncols );
	fprintf(stderr,"VC_GEE_matmult: fatal error.  exits. \n");
	Seterr_and_terminate(MATMULT_NONCONFORMITY);
	}

result = VC_GEE_create_matrix( mat1->nrows , mat2->ncols , EPHEMERAL);

resbase = result->data;
mat1base = mat1->data;
mat2base = mat2->data;

for ( j = 0 ; j < result->nrows ; j++ )
	{
	for ( i = 0 ; i < result->ncols ; i++ )
		{
		mat1loc = mat1base;
		mat2loc = mat2base;
		for ( rows = 0 ; rows < mat2->nrows ; rows++ )
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
if ( is_ephemeral( mat1 ) ) VC_GEE_destroy_matrix( mat1 );
if ( is_ephemeral( mat2 ) ) VC_GEE_destroy_matrix( mat2 );
return result;
}
 


static MATRIX *VC_GEE_px1_times_pxq( px1, pxq) /* mult elements of a colvec */
				/* across corresp row of mat */
MATRIX *px1, *pxq;
{
MATRIX *tmp;
double *load, colel;
int i, j;

if ( px1->ncols != 1 )
	{
	fprintf( stderr,"M+-: VC_GEE_px1_times_pxq: arg1 not a col-vec. Dies.\n");
	Seterr_and_terminate(PX1XPXQ_ARG1_BAD);
	}
if ( px1->nrows != pxq->nrows )
	{
	fprintf( stderr,"M+-: VC_GEE_px1_times_pxq: args not conforming.  Dies.\n");
	Seterr_and_terminate(PX1XPXQ_CONFORMITY);
	}
tmp = VC_GEE_matcopy( pxq );
load = tmp->data;
for ( i = 0 ; i < tmp->nrows ; i++ )
	{
	colel = MEL( px1, i, 0);
	for ( j = 0 ; j < tmp->ncols ; j++ )
		{
		*load *= colel ;
		load++ ;
		}
	}
free_if_ephemeral(px1);
free_if_ephemeral(pxq);
return tmp;
}

static MATRIX *VC_GEE_pxq_divby_px1( pxq, px1) /* divide elements of a colvec */
				/* into corresp row of mat */
MATRIX *px1, *pxq;
{
MATRIX *tmp;
double *load, colel;
int i, j;
if ( px1->ncols != 1 )
	{
	fprintf( stderr,"M+-: VC_GEE_pxq_divby_px1: arg2 not a col-vec. Dies.\n");
	Seterr_and_terminate(PXQDPX1_ARG1_BAD);
	}
if ( px1->nrows != pxq->nrows )
	{
	fprintf( stderr,"M+-: VC_GEE_pxq_divby_px1: args not conforming.  Dies.\n");
	Seterr_and_terminate(PXQDPX1_CONFORMITY);
	}

tmp = VC_GEE_matcopy( pxq );
load = tmp->data;
for (  i = 0 ; i < tmp->nrows ; i++ )
	{
	colel = MEL( px1, i, 0);
	for ( j = 0 ; j < tmp->ncols ; j++ )
		{
		*load = (*load) / colel ;
		load++ ;
		}
	}
free_if_ephemeral(px1);
free_if_ephemeral(pxq);
return tmp;
}
 
static MATRIX *VC_GEE_scalar_times_matrix( a , X )
double a;
MATRIX *X;
{
MATRIX *tmp;
double *tbase;
int i, nelem;
tmp = VC_GEE_matcopy(X);
nelem = get_nelem(tmp);
tbase = tmp->data;
for ( i = 0 ; i < nelem ; i++ ) {
	*tbase *= a ;
	tbase++ ;
}
free_if_ephemeral( X );
return tmp;
}

 
 


static void VC_GEE_matdump(mat)
MATRIX *mat;
{
double *curel;
int outtok = 0;
int nel;

nel = mat->nrows * mat->ncols;

for ( curel = mat->data ;  curel < mat->data + nel ; curel++ )
	{
	printf(  ((fabs(*curel)<.00001) && (fabs(*curel)>0.)) ? "%.4le%c" : "%.4lf%c" , *curel, 
			( outtok++%mat->ncols == mat->ncols-1 )
			? '\n' : ' ' );
	}
/* DOES NOT CLEAN */
}



 
 
static MATRIX *VC_GEE_luinv(X)
MATRIX *X;
{
  MATRIX *Y;
  doublereal det[2] , *work;
  integer job[1];
  static int VC_GEE_dgefaXXY_(), VC_GEE_dgediXXY_();

  doublereal *y;
  integer info[1], *ipvt;
  integer nrows, ncols;

  det[0] = 0.;
  det[1] = 0.;

  info[0] = 0;
  job[0] = 0;

  Y = VC_GEE_matcopy(X);  /* inversion in situ */

  nrows = X->nrows;
  ncols = X->ncols;

  ipvt = (integer *)malloc((unsigned)nrows*sizeof(integer));
  work = (doublereal *)malloc((unsigned)nrows*sizeof(doublereal));

  y = (doublereal *)Y->data;
  VC_GEE_dgefaXXY_(y,&nrows,&ncols,ipvt,info);

  job[0] = 11;
  VC_GEE_dgediXXY_(y,&nrows,&ncols,ipvt,det,work,job);

  free(ipvt);
  free(work);
  free_if_ephemeral(X);

  return Y;
}

static MATRIX *VC_GEE_covlag( inmat, lag, demean )
MATRIX *inmat;
int lag, demean;
{
MATRIX *xrows[MAX_COVLAG], *res, *temp;
int n, i, j, nv, q;
double nrec;

n = inmat->nrows;
nrec = (double)1./(double)n;
if ( n > MAX_COVLAG )
	{
	fprintf(stderr,"VC_GEE_covlag: arg has > MAX_COVLAG rows. Dies.\n");
	Seterr_and_terminate(EXCEED_MAX_COVLAG);
	}

nv = inmat->ncols;

res = VC_GEE_create_matrix( nv, lag*nv, EPHEMERAL );

for ( q = 0 ; q < n ; q++ )
	{
	xrows[q] = VC_GEE_extract_rows( inmat, q , q );
	make_permanent(xrows[q]);
	}


for ( i = 0 ; i < lag ; i++ )
	{
	temp = VC_GEE_create_matrix( nv, nv, EPHEMERAL );
	for ( j = i ; j < n ; j++ )
		{
		if ( (j-i) < n ) temp = VC_GEE_matadd( temp,
			VC_GEE_matmult(VC_GEE_transp(xrows[j]),xrows[j-i]));
		}
	VC_GEE_plug( VC_GEE_scalar_times_matrix( nrec, temp ), res, 0, i*nv );
	}

for ( q = 0 ; q < n ; q++ )
	{
	VC_GEE_destroy_matrix(xrows[q]);
	}
return res;
}




static MATRIX *VC_GEE_ident( ord )
int ord;
{
MATRIX *I;
int i;

I = VC_GEE_create_matrix( ord, ord, EPHEMERAL );
for ( i = 0 ; i < ord ; i++ )
	*(ELREF(I,i,i)) = (double)1.0;
return I;
}

static MATRIX *VC_GEE_col_1s( k )
int k;
{
MATRIX *tmp;
int i;
tmp = VC_GEE_create_matrix( k , 1 , EPHEMERAL );
for ( i = 0 ; i < k ; i++ )
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
	fprintf(stderr,"VC_GEE_nchanges:  must be column VC_GEE_vector; ncols = %d.\n",
				X->ncols);
	fprintf(stderr,"VC_GEE_nchanges: exiting.\n");
	exit(1);
	}

iVC_GEE_start = MEL( X , 0 , 0 );

for ( i = 1 ; i < X->nrows ; i++ )
	{
	if ( MEL ( X , i , 0 ) != iVC_GEE_start )
		{
		tmp++;
		iVC_GEE_start = MEL ( X , i , 0 );
		}
	}
return tmp;
}
 
static  MATRIX *VC_GEE_matanticlog( x )
 MATRIX *x;
 {
 double *load, *look;
 double exp();
 MATRIX *tmp;
 int nelem, i;
  
 nelem = get_nelem( x );
 tmp = VC_GEE_create_matrix( x->nrows, x->ncols , EPHEMERAL );
 load = tmp->data;
 look = x->data;
 for ( i = 0 ; i < nelem ; i++ )
         *(load++) = 1-exp(-exp(*look++));
 free_if_ephemeral(x);
 return tmp ;
 }


/* gee support @(#) clinluxxy.c 3.1 94/03/03 */



static /* Subroutine */ int VC_GEE_daxpyXXY_(n, da, dx, incx, dy, incy)
integer *n;
doublereal *da, *dx;
integer *incx;
doublereal *dy;
integer *incy;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i, m, ix, iy, mp1;


/*     constant times a VC_GEE_vector plus a VC_GEE_vector. */
/*     uses unrolled loops for increments equal to one. */
/*     jack dongarra, linpack, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*da == 0.) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	dy[iy] += *da * dx[ix];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 4;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i = 1; i <= i__1; ++i) {
	dy[i] += *da * dx[i];
/* L30: */
    }
    if (*n < 4) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= i__1; i += 4) {
	dy[i] += *da * dx[i];
	dy[i + 1] += *da * dx[i + 1];
	dy[i + 2] += *da * dx[i + 2];
	dy[i + 3] += *da * dx[i + 3];
/* L50: */
    }
    return 0;
} /* daxpy_ */


/* dgeco.f -- translated by f2c (version of 21 October 1993  13:46:10).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/


/* Table of constant values */

static integer c__1 = 1;


/* dgedi.f -- translated by f2c (version of 21 October 1993  13:46:10).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/


/* Table of constant values */


static /* Subroutine */ int VC_GEE_dgediXXY_(a, lda, n, ipvt, det, work, job)
doublereal *a;
integer *lda, *n, *ipvt;
doublereal *det, *work;
integer *job;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i, j, k, l;
    static doublereal t;
    static /* Subroutine */ int VC_GEE_dscalXXY_(), VC_GEE_dswapXXY_(), VC_GEE_daxpyXXY_();
    static integer kb, kp1, nm1;
    static doublereal ten;


/*     dgedi computes the determinant and inverse of a matrix */
/*     using the factors computed by dgeco or dgefa. */

/*     on entry */

/*        a       double precision(lda, n) */
/*                the output from dgeco or dgefa. */

/*        lda     integer */
/*                the leading dimension of the array  a . */

/*        n       integer */
/*                the order of the matrix  a . */

/*        ipvt    integer(n) */
/*                the pivot VC_GEE_vector from dgeco or dgefa. */

/*        work    double precision(n) */
/*                work VC_GEE_vector.  contents destroyed. */

/*        job     integer */
/*                = 11   both determinant and inverse. */
/*                = 01   inverse only. */
/*                = 10   determinant only. */

/*     on return */

/*        a       inverse of original matrix if requested. */
/*                otherwise unchanged. */

/*        det     double precision(2) */
/*                determinant of original matrix if requested. */
/*                otherwise not referenced. */
/*                determinant = det(1) * 10.0**det(2) */
/*                with  1.0 .le. dabs(det(1)) .lt. 10.0 */
/*                or  det(1) .eq. 0.0 . */

/*     error condition */

/*        a division by zero will occur if the input factor contains */
/*        a zero on the diagonal and the inverse is requested. */
/*        it will not occur if the subroutines are called correctly */
/*        and if dgeco has set rcond .gt. 0.0 or dgefa has set */
/*        info .eq. 0 . */

/*     linpack. this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     blas daxpy,dscal,dswap */
/*     fortran dabs,mod */

/*     internal variables */



/*     compute determinant */

    /* Parameter adjustments */
    --work;
    --det;
    --ipvt;
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    if (*job / 10 == 0) {
	goto L70;
    }
    det[1] = 1.;
    det[2] = 0.;
    ten = 10.;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	if (ipvt[i] != i) {
	    det[1] = -det[1];
	}
	det[1] = a[i + i * a_dim1] * det[1];
/*        ...exit */
	if (det[1] == 0.) {
	    goto L60;
	}
L10:
	if (abs(det[1]) >= 1.) {
	    goto L20;
	}
	det[1] = ten * det[1];
	det[2] += -1.;
	goto L10;
L20:
L30:
	if (abs(det[1]) < ten) {
	    goto L40;
	}
	det[1] /= ten;
	det[2] += 1.;
	goto L30;
L40:
/* L50: */
	;
    }
L60:
L70:

/*     compute inverse(u) */

    if (*job % 10 == 0) {
	goto L150;
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	a[k + k * a_dim1] = 1. / a[k + k * a_dim1];
	t = -a[k + k * a_dim1];
	i__2 = k - 1;
	VC_GEE_dscalXXY_(&i__2, &t, &a[k * a_dim1 + 1], &c__1);
	kp1 = k + 1;
	if (*n < kp1) {
	    goto L90;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    t = a[k + j * a_dim1];
	    a[k + j * a_dim1] = 0.;
	    VC_GEE_daxpyXXY_(&k, &t, &a[k * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &
		    c__1);
/* L80: */
	}
L90:
/* L100: */
	;
    }

/*        form inverse(u)*inverse(l) */

    nm1 = *n - 1;
    if (nm1 < 1) {
	goto L140;
    }
    i__1 = nm1;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n - kb;
	kp1 = k + 1;
	i__2 = *n;
	for (i = kp1; i <= i__2; ++i) {
	    work[i] = a[i + k * a_dim1];
	    a[i + k * a_dim1] = 0.;
/* L110: */
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    t = work[j];
	    VC_GEE_daxpyXXY_(n, &t, &a[j * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &
		    c__1);
/* L120: */
	}
	l = ipvt[k];
	if (l != k) {
	    VC_GEE_dswapXXY_(n, &a[k * a_dim1 + 1], &c__1, &a[l * a_dim1 + 1], &c__1);
	}
/* L130: */
    }
L140:
L150:
    return 0;
} /* dgedi_ */

/* dgefa.f -- translated by f2c (version of 21 October 1993  13:46:10).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/


/* Table of constant values */


static /* Subroutine */ int VC_GEE_dgefaXXY_(a, lda, n, ipvt, info)
doublereal *a;
integer *lda, *n, *ipvt, *info;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer j, k, l;
    static doublereal t;
    static /* Subroutine */ int VC_GEE_dscalXXY_(), VC_GEE_daxpyXXY_();
    static integer VC_GEE_idamaxXXY_();
    static integer kp1, nm1;

/*     dgefa factors a double precision matrix by gaussian elimination. */

/*     dgefa is usually called by dgeco, but it can be called */
/*     directly with a saving in time if  rcond  is not needed. */
/*     (time for dgeco) = (1 + 9/n)*(time for dgefa) . */

/*     on entry */

/*        a       double precision(lda, n) */
/*                the matrix to be factored. */

/*        lda     integer */
/*                the leading dimension of the array  a . */

/*        n       integer */
/*                the order of the matrix  a . */

/*     on return */

/*        a       an upper triangular matrix and the multipliers */
/*                which were used to obtain it. */
/*                the factorization can be written  a = l*u  where */
/*                l  is a product of permutation and unit lower */
/*                triangular matrices and  u  is upper triangular. */

/*        ipvt    integer(n) */
/*                an integer VC_GEE_vector of pivot indices. */

/*        info    integer */
/*                = 0  normal value. */
/*                = k  if  u(k,k) .eq. 0.0 .  this is not an error */
/*                     condition for this subroutine, but it does */
/*                     indicate that dgesl or dgedi will divide by zero */
/*                     if called.  use  rcond  in dgeco for a reliable */
/*                     indication of singularity. */

/*     linpack. this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     blas daxpy,dscal,idamax */

/*     internal variables */



/*     gaussian elimination with partial pivoting */

    /* Parameter adjustments */
    --ipvt;
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    *info = 0;
    nm1 = *n - 1;
    if (nm1 < 1) {
	goto L70;
    }
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	kp1 = k + 1;

/*        find l = pivot index */

	i__2 = *n - k + 1;
	l = VC_GEE_idamaxXXY_(&i__2, &a[k + k * a_dim1], &c__1) + k - 1;
	ipvt[k] = l;

/*        zero pivot implies this column already triangularized */

	if (a[l + k * a_dim1] == 0.) {
	    goto L40;
	}

/*           interchange if necessary */

	if (l == k) {
	    goto L10;
	}
	t = a[l + k * a_dim1];
	a[l + k * a_dim1] = a[k + k * a_dim1];
	a[k + k * a_dim1] = t;
L10:

/*           compute multipliers */

	t = -1. / a[k + k * a_dim1];
	i__2 = *n - k;
	VC_GEE_dscalXXY_(&i__2, &t, &a[k + 1 + k * a_dim1], &c__1);

/*           row elimination with column indexing */

	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    t = a[l + j * a_dim1];
	    if (l == k) {
		goto L20;
	    }
	    a[l + j * a_dim1] = a[k + j * a_dim1];
	    a[k + j * a_dim1] = t;
L20:
	    i__3 = *n - k;
	    VC_GEE_daxpyXXY_(&i__3, &t, &a[k + 1 + k * a_dim1], &c__1, &a[k + 1 + j * 
		    a_dim1], &c__1);
/* L30: */
	}
	goto L50;
L40:
	*info = k;
L50:
/* L60: */
	;
    }
L70:
    ipvt[*n] = *n;
    if (a[*n + *n * a_dim1] == 0.) {
	*info = *n;
    }
    return 0;
} /* dgefa_ */

/* dgidamax.f -- translated by f2c (version of 21 October 1993  13:46:10).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/


static integer VC_GEE_idamaxXXY_(n, dx, incx)
integer *n;
doublereal *dx;
integer *incx;
{
    /* System generated locals */
    integer ret_val, i__1;
    doublereal d__1;

    /* Local variables */
    static doublereal dmax_;
    static integer i, ix;


/*     finds the index of element having max. absolute value. */
/*     jack dongarra, linpack, 3/11/78. */


    /* Parameter adjustments */
    --dx;

    /* Function Body */
    ret_val = 0;
    if (*n < 1) {
	return ret_val;
    }
    ret_val = 1;
    if (*n == 1) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    ix = 1;
    dmax_ = abs(dx[1]);
    ix += *incx;
    i__1 = *n;
    for (i = 2; i <= i__1; ++i) {
	if ((d__1 = dx[ix], abs(d__1)) <= dmax_) {
	    goto L5;
	}
	ret_val = i;
	dmax_ = (d__1 = dx[ix], abs(d__1));
L5:
	ix += *incx;
/* L10: */
    }
    return ret_val;

/*        code for increment equal to 1 */

L20:
    dmax_ = abs(dx[1]);
    i__1 = *n;
    for (i = 2; i <= i__1; ++i) {
	if ((d__1 = dx[i], abs(d__1)) <= dmax_) {
	    goto L30;
	}
	ret_val = i;
	dmax_ = (d__1 = dx[i], abs(d__1));
L30:
	;
    }
    return ret_val;
} /* idamax_ */

/* dscal.f -- translated by f2c (version of 21 October 1993  13:46:10).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/


static /* Subroutine */ int VC_GEE_dscalXXY_(n, da, dx, incx)
integer *n;
doublereal *da, *dx;
integer *incx;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i, m, nincx, mp1;


/*     scales a VC_GEE_vector by a constant. */
/*     uses unrolled loops for increment equal to one. */
/*     jack dongarra, linpack, 3/11/78. */


    /* Parameter adjustments */
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    nincx = *n * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i = 1; i__2 < 0 ? i >= i__1 : i <= i__1; i += i__2) {
	dx[i] = *da * dx[i];
/* L10: */
    }
    return 0;

/*        code for increment equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__2 = m;
    for (i = 1; i <= i__2; ++i) {
	dx[i] = *da * dx[i];
/* L30: */
    }
    if (*n < 5) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__2 = *n;
    for (i = mp1; i <= i__2; i += 5) {
	dx[i] = *da * dx[i];
	dx[i + 1] = *da * dx[i + 1];
	dx[i + 2] = *da * dx[i + 2];
	dx[i + 3] = *da * dx[i + 3];
	dx[i + 4] = *da * dx[i + 4];
/* L50: */
    }
    return 0;
} /* dscal_ */

/* dswap.f -- translated by f2c (version of 21 October 1993  13:46:10).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/


static /* Subroutine */ int VC_GEE_dswapXXY_(n, dx, incx, dy, incy)
integer *n;
doublereal *dx;
integer *incx;
doublereal *dy;
integer *incy;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i, m;
    static doublereal dtemp;
    static integer ix, iy, mp1;


/*     interchanges two VC_GEE_vectors. */
/*     uses unrolled loops for increments equal one. */
/*     jack dongarra, linpack, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*       code for unequal increments or equal increments not equal */
/*         to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	dtemp = dx[ix];
	dx[ix] = dy[iy];
	dy[iy] = dtemp;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*       code for both increments equal to 1 */


/*       clean-up loop */

L20:
    m = *n % 3;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i = 1; i <= i__1; ++i) {
	dtemp = dx[i];
	dx[i] = dy[i];
	dy[i] = dtemp;
/* L30: */
    }
    if (*n < 3) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= i__1; i += 3) {
	dtemp = dx[i];
	dx[i] = dy[i];
	dy[i] = dtemp;
	dtemp = dx[i + 1];
	dx[i + 1] = dy[i + 1];
	dy[i + 1] = dtemp;
	dtemp = dx[i + 2];
	dx[i + 2] = dy[i + 2];
	dy[i + 2] = dtemp;
/* L50: */
    }
    return 0;
} /* dswap_ */

/* gee support @(#) VC_GEE_diag_as_vec.c 3.3 94/03/09 */


static MATRIX *VC_GEE_diag_as_vec(inmat)
MATRIX *inmat;
{
  int i;
  MATRIX *outmat;

  if(inmat->ncols!=inmat->nrows)
    {
      fprintf(stderr,"M+-: VC_GEE_diag_as_vec: arg is not a square matrix. Dies.\n");
      fprintf(stderr,"\nNumber of columns = %d",inmat->ncols);
      fprintf(stderr,"\nNumber of rows    = %d\n",inmat->nrows);
      Seterr_and_terminate(DIAG_AS_VEC_ARG_BAD);
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
  tmp= VC_GEE_matcopy(x);
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

#if defined(R_VERSION) && R_VERSION >= R_Version(0,99,0)
#include <R_ext/Mathlib.h>
#define dnorm1(x) dnorm4(x, 0.0, 1.0, 0)
#define pnorm1(x) pnorm5(x, 0.0, 1.0, 1, 0)
#else
#include <Mathlib.h>
#define dnorm1(x) dnorm(x, 0.0, 1.0)
#define pnorm1(x) pnorm(x, 0.0, 1.0)
#endif


/* following added by pj catalano to support probit link option in cgee.c */
   
static MATRIX *VC_GEE_matnpdf( x )
MATRIX *x;
{
double *load, *look;
double VC_GEE_npdf();
MATRIX *tmp;
int nelem, i;

nelem = get_nelem( x );
tmp = VC_GEE_create_matrix( x->nrows, x->ncols , EPHEMERAL );
load = tmp->data;
look = x->data;
for ( i = 0 ; i < nelem ; i++ )
	*(load++) = dnorm1(*look++);
free_if_ephemeral(x);
return tmp ;
}

static MATRIX *VC_GEE_matncdf( x )
MATRIX *x;
{
double *load, *look;
double VC_GEE_ncdf();
MATRIX *tmp;
int nelem, i;

nelem = get_nelem( x );
tmp = VC_GEE_create_matrix( x->nrows, x->ncols , EPHEMERAL );
load = tmp->data;
look = x->data;
for ( i = 0 ; i < nelem ; i++ )
	*(load++) = pnorm1(*look++);
free_if_ephemeral(x);
return tmp ;
}

/* end of pj catalano additions  */
