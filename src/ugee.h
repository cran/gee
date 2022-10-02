/* ugee.h -- constants, structs and macros for cgee */
/* previously gee support @(#) cgee.h */
/* 4.1 98/01/27 */



#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <setjmp.h>
#include <R_ext/Applic.h> /* BLAS */
extern void F77_NAME(dgefa)(double*,int*,int*,int*,int*);
extern void F77_NAME(dgedi)(double*,int*,int*,int*,double*,double*,int*);

#define PERMANENT 1
#define EPHEMERAL 0
#define MAX_COVLAG 30
#define NO_ERROR 0
#define UNKNOWN_FAILURE 1
#define NO_MEM_MATSTRUCT 2
#define NO_MEM_MATDATA 3
#define SPLIT_FAIL 4
#define MATMULT_NONCONFORMITY 5
#define MATADD_NONCONFORMITY 6
#define CCHOL_FAIL 7
#define CORNER_FAIL 8
#define EXCEED_MAX_COVLAG 9
#define BAD_TOEPLITZ_ARG 10
#define PLUG_FAIL 11
#define PX1XPXQ_ARG1_BAD 12
#define PX1XPXQ_CONFORMITY 13
#define PXQDPX1_ARG1_BAD 14
#define PXQDPX1_CONFORMITY 15
#define CCHOL_NOT_SQUARE 16
#define MATREAD_OPEN_FAIL 17
#define MATREAD_NOT_RECTANGLE 18

 

typedef struct matrix
		{
		int nrows, ncols;
		double *data;
		int permanence;
		} MATRIX;
/* element reference is handled principally by MEL */
#define ELREF( matp , s1, s2 ) ((matp)->data)+(s2)+((s1)*(matp->ncols))
#define MEL(X ,i, j) (*(ELREF( (X), (i), (j) ) ))
/* $ y = MEL(X,i,j) :\Rightarrow y \in {\cal R} \wedge y = x_{ij} $ */
#define get_nelem( x ) (((x)->nrows) * ((x)->ncols))
 

 
/* recommended allocation defs by B. Ripley, IX.1996 */

/* #ifdef DONT_USE_S_ALLOC
   #else
   extern void *S_alloc(void) ;
   #define malloc(n) S_alloc(n, 1)
   #define calloc S_alloc
   #define free Free
   #define cfree Free 
   #endif */

/* mem alloc macros supplied by Bill Dunlap, VII.1996 */

#define malloc(n) S_alloc(n, 1)
#define calloc S_alloc
#define free(p) {p;}
#define cfree(p) {p;} 

#define is_permanent( x ) (x)->permanence == PERMANENT
#define is_ephemeral( x ) (x)->permanence == EPHEMERAL
#define make_permanent( x ) (x)->permanence = PERMANENT;
#define make_ephemeral( x ) (x)->permanence = EPHEMERAL;

#define free_if_ephemeral( x ) if (is_ephemeral((x))) VC_GEE_destroy_matrix((x))
 
#define from_S( Sdblptr , Srowintptr , Scolintptr , Matptr ) \
Matptr = VC_GEE_create_matrix( (int)*Srowintptr, (int)*Scolintptr , EPHEMERAL ); \
{ \
int i, j, Scol, Srow; \
double *Sload; \
Scol = *Scolintptr; \
Srow = *Srowintptr; \
Sload = Sdblptr; \
for ( j = 0 ; j < Scol ; j++ ) \
	{ \
	for ( i = 0 ; i < Srow ; i++ ) \
		{ \
		MEL( Matptr , i , j ) = (double) * ( Sload ++ ); \
		} \
	} \
}
/* end define |from_S| */

#define to_S( Matptr, Sdblptr ) \
{ \
int i, j; \
double *Sload; \
Sload = Sdblptr; \
for ( j = 0 ; j < Matptr->ncols ; j++ ) \
	{ \
	for ( i = 0 ; i < Matptr->nrows ; i++ ) \
		{ \
		* ( Sload ++ ) = MEL( Matptr , i , j ); \
		} \
	} \
}
 


/* chanmatfuns.h */
/* gee support @(#) chanmatfuns.h 3.2 94/03/09 */


static MATRIX *VC_GEE_create_matrix(int, int, int),
    *VC_GEE_matcopy(MATRIX *),
    *VC_GEE_extract_rows(MATRIX *, int,  int),
    *VC_GEE_matadd(MATRIX *, MATRIX *),
    *VC_GEE_matsub(MATRIX *, MATRIX *),
    *VC_GEE_matmult(MATRIX *, MATRIX *),
    *VC_GEE_transp(MATRIX *),
    *VC_GEE_col_1s(int),
    *VC_GEE_matabs(MATRIX *),
    *VC_GEE_matexp(MATRIX *),
    *VC_GEE_px1_times_pxq(MATRIX *, MATRIX *),
    *VC_GEE_pxq_divby_px1(MATRIX *, MATRIX *),
    *VC_GEE_scalar_times_matrix(double, MATRIX *),
    *VC_GEE_ident(int),
    *VC_GEE_form_diag(MATRIX *),
    *VC_GEE_corner(MATRIX *, int, int),
    *VC_GEE_covlag(MATRIX *, int, int),
    *VC_GEE_toeplitz(MATRIX *),
    *VC_GEE_band(MATRIX *, int),
    *VC_GEE_extract_cols(MATRIX *, int, int),
     /* following two functions added by pj catalano  */
    *VC_GEE_matnpdf(MATRIX *),
    *VC_GEE_matncdf(MATRIX *),
    *VC_GEE_matanticlog(MATRIX *)
     ;

static double VC_GEE_matmax(MATRIX *), VC_GEE_elsum(MATRIX *);

static void VC_GEE_matdump(MATRIX *),
    VC_GEE_plug(MATRIX *, MATRIX *, int, int),
    VC_GEE_destroy_matrix(MATRIX *);

static MATRIX *VC_GEE_luinv(MATRIX *);


static int VC_GEE_split(MATRIX *, MATRIX *, MATRIX *[]),
    VC_GEE_nchanges(MATRIX *);

void Cgee(double *, double *, double *, double *, double *, int *, int *,
	  int *, int *, double *, double *, double *, double *, double *,
	  double *, int *, int *, int *, int *, int *, int *);


#ifndef MAX_NUM_CLUSTS
#define MAX_NUM_CLUSTS 5000
#endif

#define TOO_MANY_CLUSTS  100
#define UNKNOWN_LINK  101
#define UNKNOWN_VAR_MEAN_REL 102
#define M_DEP_SIZE_LT_M 103
#define MAXITER_EXCEEDED  104
#define LOGISTIC_DIVERGENCE  105

enum link_type { VC_GEE_identity , 
		 logarithm , 
		 logit , 
		 reciprocal,
		 probit,
		cloglog };

enum var_mean_rel_type { Gaussian , 
			 Poisson , 
			 Binomial , 
			 Gamma };

enum corstruct_type { independence, 
	    	      fixed , 
		      stat_M_dep , 
		      non_stat_M_dep , 
		      exchangeable , 
		      AR_M , 
		      unstructured };


