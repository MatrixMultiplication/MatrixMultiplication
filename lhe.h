#include "lhep.h"


/*! \struct lhe_f
 * the function F is encrypted as two parts p0 and p1
 */
typedef struct
{
  fq_poly_t p0;
  fq_poly_t p1;
} lhe_c;


int fq_poly_rand(fq_poly_t pol, int deg, fq_ctx_t ctx, bn_t bound);

/**
 * generate a LHE key
 * @param[out]s          the secret key get from this function
 * @param[out]a0         the first term of the public key get from this function pk=(a0, b0)
 * @param[out]b0         the second term of the public key
 * @param[in] par        some parameter in flint lib
 */
int lhe_keygen(fq_poly_t s, fq_poly_t a0, fq_poly_t b0, lhe_par *par);

/**
 * encrype the information
 * @param[out]c          the out polynomial
 * @param[in] a0         the first term of the public key
 * @param[in] b0         the second term of the public key
 * @param[in] m          the information wanted to be encryped
 * @param[in] par        some parameter in flint lib
 */
int lhe_enc(lhe_c *c, fq_poly_t a0, fq_poly_t b0, fq_poly_t m,lhe_par *par);

/**
 * decode the information
 * @param[out]m1         the decoded result
 * @param[in] s          the secret key
 * @param[in] c          the encryped information
 * @param[in] par        some parameter in flint lib
 */
int lhe_dec(fq_poly_t m1, fq_poly_t s, lhe_c *c, lhe_par *par);
