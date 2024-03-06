#include "pi2.h"


int msg_modp(fq_poly_t m, lhe_par *par)
{
   long len=fq_poly_length(m,par->ctx);

   fq_t zq;
   fq_init(zq,par->ctx);

   fmpz_t zf;
   fmpz_init(zf);

   for (int i=0;i<len;i++)
   {

     fq_poly_get_coeff(zq,m,i,par->ctx);
     fmpz_set_str(zf,fq_get_str_pretty(zq,par->ctx), 10);
     fmpz_mod(zf,zf,par->pf);
     fq_set_fmpz(zq,zf,par->ctx);
     fq_poly_set_coeff(m,i,zq,par->ctx);
   }

    fq_clear(zq,par->ctx);
    fmpz_clear(zf);

   return 0;
}

int vc_pgen(vc_p *vcp, vc_k * vck, fq_mat_t X, lhe_par *par)
{

   lhe_keygen(vcp->s,vcp->pk.a0,vcp->pk.b0, par);

   fmpz_t *zf;
   //fmpz_init(zf);
   fq_t zq;
   fq_init(zq,par->ctx);

   fq_poly_t zqp;
   fq_poly_init(zqp,par->ctx);
   vcp->C=malloc(sizeof(lhe_c*)*vck->d*par->t);//row(vck->d)Xcolumn(par->t)

   for (int j=0;j<vck->d;j++)//encryption the j-th row of x
   {      /*MatrixMatrix*/
     for(int i=0;i<par->t;i++)//encryption the i-th block in j-th row
     {
        int n = par->t*par->ni;
        for (slong k=0;k<par->ni;k++)
        {
           //bn2fmpz(zf,X[j*n+i*par->ni+k]);
           //zf=fq_mat_entry(X,j,i*par->ni+k)->coeffs;
           //fq_set_fmpz(zq,zf[0],par->ctx);
           fq_poly_set_coeff(zqp, k, fq_mat_entry(X,j,i*par->ni+k), par->ctx);
         }
        //printf("\nzqp\n");
        //fq_poly_print(zqp,par->ctx);
        //printf("\nendZqp\n");
        vcp->C[j*par->t+i]=malloc(sizeof(lhe_c));//(j,i)-th of matrix vcp->C
        lhe_enc(vcp->C[j*par->t+i],vcp->pk.a0,vcp->pk.b0,zqp,par);
     }
   }

   //fmpz_clear(zf);
   fq_clear(zq,par->ctx);
   fq_poly_clear(zqp,par->ctx);

   return 0;
}

// the vc_comp algorithm
int vc_comp(lhe_c **nv, vc_k *vck, vc_p *vcp, lhe_par *par)
{
    /*MatrixMatrix*/
    int block;
    int indexInBlock;
    fq_t zq;
    fq_init(zq,par->ctx);
    int elementsOfX=2*par->t*par->ni;

    fq_mat_t vcpC;
    fq_mat_init(vcpC,vck->d,elementsOfX,par->ctx);//vcpC is d*elementsOfX matrix

    fq_mat_t FC;
    fq_mat_init(FC,vck->m,elementsOfX,par->ctx);//vcpC is d*elementsOfX matrix
    for (int i=0;i<vck->d;i++)
    {
        //mpz_set_ui(tempVector[i],0);//set the value to 0
        for (int l=0;l<elementsOfX;l++)//i-th row of vcp->C
        {
            block = l / (2*par->ni);
            indexInBlock = l % (2*par->ni);

            if (indexInBlock < par->ni) fq_poly_get_coeff(zq,vcp->C[i*par->t+block]->p0,indexInBlock,par->ctx);
            else fq_poly_get_coeff(zq,vcp->C[i*par->t+block]->p1,indexInBlock-par->ni,par->ctx);

            fq_mat_entry_set(vcpC,i,l,zq,par->ctx);
        }
    }
    //FC=F*vcpC
    fq_mat_mul(FC,vck->F,vcpC,par->ctx);//FC(mx2tn)

    //fq_mat_print_pretty(FC,par->ctx);

    //FC into nv
    for (int i=0;i<vck->m;i++)
    {
        for (int l=0;l<elementsOfX;l++)//i-th row of vcp->C
        {
            block = l / (2*par->ni);
            indexInBlock = l % (2*par->ni);

            if (indexInBlock==0)
            {
                nv[i*par->t+block]=malloc(sizeof(lhe_c));
                fq_poly_init(nv[i*par->t+block]->p0,par->ctx);
                fq_poly_init(nv[i*par->t+block]->p1,par->ctx);
            }

            if (indexInBlock < par->ni) fq_poly_set_coeff(nv[i*par->t+block]->p0,indexInBlock, fq_mat_entry(FC,i,l), par->ctx);
            else fq_poly_set_coeff(nv[i*par->t+block]->p1,indexInBlock-par->ni, fq_mat_entry(FC,i,l), par->ctx);
        }
    }

    fq_mat_clear(vcpC,par->ctx);
    fq_mat_clear(FC,par->ctx);
  return 0;
}

fq_mat_t vcpC, nvMatrix;

int polyToMatrixForm(vc_k *vck, vc_p *vcp, lhe_c **nv, lhe_par *par)
{
    int elementsOfX=2*par->t*par->ni;

    fq_mat_init(vcpC,vck->d,elementsOfX,par->ctx);//vcpC is d*elementsOfX matrix
    fq_mat_init(nvMatrix,vck->m,elementsOfX,par->ctx);//nvMatrix is m*elementsOfX matrix

    int block;
    int indexInBlock;
    fq_t zq;
    fq_init(zq,par->ctx);
    //calculate left vector
    for (int i=0;i<vck->d;i++)
    {
        //mpz_set_ui(tempVector[i],0);//set the value to 0
        for (int l=0;l<elementsOfX;l++)//i-th row of vcp->C
        {
            block = l / (2*par->ni);
            indexInBlock = l % (2*par->ni);

            if (indexInBlock < par->ni) fq_poly_get_coeff(zq,vcp->C[i*par->t+block]->p0,indexInBlock,par->ctx);
            else fq_poly_get_coeff(zq,vcp->C[i*par->t+block]->p1,indexInBlock-par->ni,par->ctx);

            fq_mat_entry_set(vcpC,i,l,zq,par->ctx);
        }
    }

    //Next
    for (int i=0;i<vck->m;i++)
    {
        for (int l=0;l<elementsOfX;l++)//i-th row of nv
        {
            block = l / (2*par->ni);
            indexInBlock = l % (2*par->ni);

            if (indexInBlock < par->ni) fq_poly_get_coeff(zq,nv[i*par->t+block]->p0,indexInBlock,par->ctx);
            else fq_poly_get_coeff(zq,nv[i*par->t+block]->p1,indexInBlock-par->ni,par->ctx);

            fq_mat_entry_set(nvMatrix,i,l,zq,par->ctx);
        }
    }

}

int vc_vrfy(vc_k *vck, vc_p *vcp, lhe_c **nv, lhe_par *par)
{
    /*VerifyAlgorithm3-Monte Carlo: Freivald's matrix multiplication EXTENSION*/
    int flag=-1;
    int elementsOfX=2*par->t*par->ni;

    fq_mat_t x;
    flint_rand_t state;
    flint_randinit(state);
    fq_mat_init(x,elementsOfX,1,par->ctx);//x is elementsOfX*1 matrix
    fq_mat_randtest(x, state, par->ctx);
    flint_randclear(state);

    fq_mat_t tempVector, leftVector, rightVector;

    fq_mat_init(tempVector,vck->d,1,par->ctx);//tempVector is d*1 matrix
    //tempVector=(vcp->C*x)
    fq_mat_mul(tempVector,vcpC,x,par->ctx);
    fq_mat_clear(vcpC,par->ctx);

    fq_mat_init(leftVector,vck->m,1,par->ctx);//leftVector is m*1 matrix
    //leftVector=F*tempVector
    fq_mat_mul(leftVector,vck->F,tempVector,par->ctx);
    fq_mat_clear(tempVector,par->ctx);

    fq_mat_init(rightVector,vck->m,1,par->ctx);//rightVector is m*1 matrix
    //rightVector=nv*x
    fq_mat_mul(rightVector,nvMatrix,x,par->ctx);
    fq_mat_clear(nvMatrix,par->ctx);
    fq_mat_clear(x,par->ctx);

    flag = fq_mat_equal(leftVector,rightVector,par->ctx);
    fq_mat_clear(leftVector,par->ctx);
    fq_mat_clear(rightVector,par->ctx);

    if(flag==0) return 2;

    return 0;
}

int vc_dec(fq_poly_t *y, vc_k *vck, vc_p *vcp, lhe_c **nv, lhe_par *par)
{
	for (int i=0;i<vck->m;i++)
	{
		for (int j=0;j<par->t;j++)//(i,j)-th of matrix
		{
		     lhe_dec(y[i*par->t+j],vcp->s,nv[i*par->t+j],par);//(i,j)-th
             free(nv[i*par->t+j]);
		     msg_modp(y[i*par->t+j], par);//I add it!

//		     printf("\ndecryption\n");
//		     fq_poly_print(y[i*par->t+j],par->ctx);//print fq_poly_t
		}
  	 }
	return 0;
}