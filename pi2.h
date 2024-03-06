#include "md2.h"

/*! \struct vc_k
 * m and d are the row and column of the F passing in and the vck->F going to generate.
 * F is the encrypted function
 * G and T is used to verify
 */
typedef struct
{
  int m;
  int d;
  fq_mat_t F;
  g1_t *G;
  g1_t *T;
} vc_k;


/*! \struct vc_p
 * C is the encrypted x
 * hk is the encrypted y
 */

typedef struct
{
    fq_poly_t a0;
    fq_poly_t b0;
} public_key;

typedef struct
{
  lhe_c **C;
  fq_poly_t s;
  gt_t pi;
  public_key pk;
} vc_p;

/**
 * reduce the coefficients of a polynomial
 * @param[out] m         the out polynomial
 * @param[in] par        some parameter in flint lib
 */
int msg_modp(fq_poly_t m, lhe_par *par);

/**
 * generate a random matrix F stored fq_t
 * @param[out]F        the random matrix
 * @param[in] m        the row of the random matrix
 * @param[in] d        the column of the random matrix
 * @param[in] par      some parameter in flint lib
 */
int vc_pgen(vc_p *vcp, vc_k * vck, fq_mat_t X, lhe_par *par);

/**
 * compute the encrypted output y
 * @param[out] nv		 the encrypted output y
 * @param[in] vck        the vck->F is needed here to compute
 * @param[in] vcp        the vcp->C is needed here to compute
 * @param[in] par        some parameter in flint lib   
 */
int vc_comp(lhe_c **nv, vc_k *vck, vc_p *vcp, lhe_par *par);

int polyToMatrixForm(vc_k *vck, vc_p *vcp, lhe_c **nv, lhe_par *par);

/**
 * verify the outcome. output 0 if the outcome is true, otherwise 2
 * @param[out] y         the decoded output y
 * @param[in] vck   	 the vck->G and vck->T is needed here to verify
 * @param[in] vcp   	 the vcp->C is needed here to verify
 * @param[in] nv		 the encrypted output y
 * @param[in] par   	 some parameter in flint lib  
 */
int vc_vrfy(vc_k *vck, vc_p *vcp, lhe_c **nv, lhe_par *par);

/**
 * dec the outcome.
 * @param[out] y         the decoded output y
 * @param[in] vck   	 the vck->G and vck->T is needed here to verify
 * @param[in] vcp   	 the vcp->C is needed here to verify
 * @param[in] nv		 the encrypted output y
 * @param[in] par   	 some parameter in flint lib  
 */
int vc_dec(fq_poly_t *y, vc_k *vck, vc_p *vcp, lhe_c **nv, lhe_par *par);
