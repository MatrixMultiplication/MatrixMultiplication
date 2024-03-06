#include "lhe.h"

//this function generate a random polynomial that has small coefficients
int fq_poly_rand(fq_poly_t pol, int deg, fq_ctx_t ctx, bn_t bound)
{

  fq_poly_init(pol,ctx);

  bn_t z1;
  bn_null(z1);

  fmpz_t z2;
  fmpz_init(z2);

  fq_t z3;
  fq_init(z3,ctx);

  for (int n=0;n<deg;n++)
  {
      bn_rand_mod(z1,bound);
      bn2fmpz(z2,z1);
      fq_set_fmpz(z3,z2,ctx);  
      fq_poly_set_coeff(pol,n,z3,ctx);
  }


  //release the memory
  bn_free(z1);
  fmpz_clear(z2);
  fq_clear(z3,ctx);
  return 0;

}

// the key generation algorithm of the linearly homomorphic encryption: LHE.KeyGen
int lhe_keygen(fq_poly_t s, fq_poly_t a0, fq_poly_t b0, lhe_par *par)
{
  //reduction modulo x^n+1
  fq_poly_t Q;
  fq_poly_init(Q,par->ctx);

  fq_poly_rand(s,par->ni,par->ctx,par->r);//Randomly select an element s from the discrete Gaussian distribution $\chi$

  fq_poly_t e0;
  //Randomly select an element a_0 from Rq
  fq_poly_rand(a0,par->ni,par->ctx,par->q);
  fq_poly_init(b0,par->ctx);
  //randomly select an element $e_0$ from the discrete Gaussian distribution $\chi$
  fq_poly_rand(e0,par->ni,par->ctx,par->r);

  fq_poly_t firstTerm, secondTerm;
  fq_t fqP;
  fq_poly_init(firstTerm,par->ctx);
  fq_poly_init(secondTerm,par->ctx);
  fq_init(fqP,par->ctx);

  //b_0=a_0s+pe_0
  fq_poly_mul(firstTerm,a0,s,par->ctx);//a_0s
  fq_poly_divrem(Q,firstTerm,firstTerm,par->modf,par->ctx);
  fq_set_fmpz(fqP,par->pf,par->ctx);
  fq_poly_scalar_mul_fq(secondTerm,e0,fqP,par->ctx);//pe_0
  fq_poly_divrem(Q,secondTerm,secondTerm,par->modf,par->ctx);
  fq_poly_add(b0,firstTerm,secondTerm,par->ctx);
  fq_poly_divrem(Q,b0,b0,par->modf,par->ctx);

  // release the memory
  fq_poly_clear(e0,par->ctx);
  fq_poly_clear(firstTerm,par->ctx);
  fq_poly_clear(secondTerm,par->ctx);
  fq_clear(fqP,par->ctx);
  fq_poly_clear(Q,par->ctx);

  return 0;
}

// the encryption algorithm of the linearly homomorphic encryption: LHE.Enc
int lhe_enc(lhe_c *c, fq_poly_t a0, fq_poly_t b0, fq_poly_t m,lhe_par *par)
{
  //reduction modulo x^n+1
  fq_poly_t Q;
  fq_poly_init(Q,par->ctx);

  fq_poly_t v, e1, e2;
  fq_poly_rand(v,par->ni,par->ctx,par->r);//randomly select an element v from the discrete Gaussian distribution $\chi$
  fq_poly_rand(e1,par->ni,par->ctx,par->r);//Randomly select an element e1 from the discrete Gaussian distribution $\chi$
  fq_poly_rand(e2,par->ni,par->ctx,par->r1);//Randomly select an element e1 from the discrete Gaussian distribution $\chi'$

  fq_poly_init(c->p0,par->ctx);
  fq_poly_init(c->p1,par->ctx);

  // get -1
  fq_t u1,u2;
  fq_init(u1,par->ctx);
  fq_zero(u1,par->ctx);
  fq_init(u2,par->ctx);
  fq_one(u2,par->ctx);
  fq_sub(u1,u1,u2,par->ctx);

  fq_poly_t firstTerm, secondTerm;
  fq_t fqP;
  fq_poly_init(firstTerm,par->ctx);
  fq_poly_init(secondTerm,par->ctx);
  fq_init(fqP,par->ctx);
  fq_set_fmpz(fqP,par->pf,par->ctx);

  //p1=-a_0v-pe'
  fq_poly_mul(firstTerm,a0,v,par->ctx);//a_0v
  fq_poly_divrem(Q,firstTerm,firstTerm,par->modf,par->ctx);
  fq_poly_scalar_mul_fq(secondTerm,e1,fqP,par->ctx);//pe1
  fq_poly_divrem(Q,secondTerm,secondTerm,par->modf,par->ctx);
  fq_poly_add(c->p1,firstTerm,secondTerm,par->ctx);
  fq_poly_scalar_mul_fq(c->p1,c->p1,u1,par->ctx);//-p1
  fq_poly_divrem(Q,c->p1,c->p1,par->modf,par->ctx);

  //p0=b_0v+pe''+m
  fq_poly_mul(firstTerm,b0,v,par->ctx);//b_0v
  fq_poly_divrem(Q,firstTerm,firstTerm,par->modf,par->ctx);
  fq_poly_scalar_mul_fq(secondTerm,e2,fqP,par->ctx);//pe2
  fq_poly_divrem(Q,secondTerm,secondTerm,par->modf,par->ctx);
  fq_poly_add(c->p0,firstTerm,secondTerm,par->ctx);
  fq_poly_divrem(Q,c->p0,c->p0,par->modf,par->ctx);
  fq_poly_add(c->p0,c->p0,m,par->ctx);//p0=p0+m
  fq_poly_divrem(Q,c->p0,c->p0,par->modf,par->ctx);

  // release the memory
  fq_poly_clear(v,par->ctx);
  fq_poly_clear(e1,par->ctx);
  fq_poly_clear(e2,par->ctx);
  fq_clear(u1,par->ctx);
  fq_clear(u2,par->ctx);
  fq_poly_clear(firstTerm,par->ctx);
  fq_poly_clear(secondTerm,par->ctx);
  fq_clear(fqP,par->ctx);
  fq_poly_clear(Q,par->ctx);

  return 0;   

}

// the decryption function of the linearly homomorphic encryption scheme: LEH.Dec
int lhe_dec(fq_poly_t m1, fq_poly_t s, lhe_c *c, lhe_par *par)
{
    fq_poly_init(m1,par->ctx);
    fq_poly_mul(m1,s,c->p1,par->ctx);
    fq_poly_add(m1,m1,c->p0,par->ctx);

    // reduction modulo x^n+1
    fq_poly_t Q;
    fq_poly_init(Q,par->ctx);
    fq_poly_divrem(Q,m1,m1,par->modf,par->ctx); 
   
    // release the memory
    fq_poly_clear(Q,par->ctx);
    
    return 0;   
}
