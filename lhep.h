//#ifndef _BN_ETX_H_
//#define _BN_ETX_H_
#include "bn_ext.h"
//#endif

//the public parameter for the lhe linearly homomorphic encryption scheme
typedef struct
{
	bn_t p;
	fmpz_t pf;
	bn_t q;
	fmpz_t qf;
	bn_t n;     
	int ni;
	fmpz_t nf;
	bn_t r;
	fmpz_t rf;

	bn_t r1;
    fmpz_t r1f;

	fq_ctx_t ctx;
    fq_ctx_t pctx;
	fq_poly_t modf;
	fq_poly_t modp;
	//g1_t g;
	//g2_t h;
	//gt_t gt;
	int t;
} lhe_par;

//generate the parameters for the encryption scheme
int lhep_new(lhe_par *par);
