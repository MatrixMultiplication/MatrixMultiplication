#include "pi2.h"
#include <time.h>

void try(int mat_m, int mat_d, int t, lhe_par *par)//Note: all matrix are organized by row
{
	clock_t s1,s2,s3,s4,s5,s6,e1,e2,e3,e4,e5,e6;

    flint_rand_t state;
    flint_randinit(state);
	par->t = t;

        vc_k *vck=malloc(sizeof(*vck));
        vck->m=mat_m;
        vck->d=mat_d;
        fq_mat_init(vck->F,vck->m,vck->d,par->pctx);
        fq_mat_randtest(vck->F, state, par->pctx);
    flint_randclear(state);

	s1 = clock();

	e1 = clock();

    //fq_mat_clear(F,par->pctx);

        vc_p *vcp=malloc(sizeof(*vcp));

        //MatrixMatrix
        int n = par->t*par->ni;
        fq_mat_t X;
        fq_mat_init(X,mat_d,n,par->pctx);//X is dxn matrix
        fq_mat_randtest(X,state,par->pctx);

    fq_mat_t FX;
    fq_mat_init(FX,mat_m,n,par->pctx);

    s6 = clock();
    fq_mat_mul(FX,vck->F,X,par->pctx);
    e6 = clock();

//    printf("\nF\n");
//    fq_mat_print_pretty(vck->F,par->pctx);
//    printf("\nX\n");
//    fq_mat_print_pretty(X,par->pctx);
//    printf("\nFX\n");
//    fq_mat_print_pretty(FX,par->pctx);

    fq_mat_clear(FX,par->pctx);

	s2 = clock();
        vc_pgen(vcp,vck,X,par);
	e2 = clock();
    fq_mat_clear(X,par->pctx);

    fq_poly_t *y=malloc(sizeof(fq_poly_t)*(vck->m)*(par->t));

    lhe_c **nv=malloc(sizeof(lhe_c)*vck->m*par->t);
	s3 = clock();
        vc_comp(nv,vck,vcp,par);
	e3 = clock();

    polyToMatrixForm(vck, vcp, nv, par);

	s4 = clock();
	    int flag=vc_vrfy(vck,vcp,nv,par);
	e4 = clock();

    fq_mat_clear(vck->F,par->pctx);
    for (int j=0;j<vck->d;j++)
    {
        for(int i=0;i<par->t;i++)
        {
            free(vcp->C[j*par->t+i]=malloc(sizeof(lhe_c)));
        }
    }
    free(vcp->C);


	//fq_poly_t *y=malloc(sizeof(fq_poly_t)*(vck->m)*(par->t));
	s5 = clock();

	if (flag==0)
	{
		vc_dec(y,vck,vcp,nv,par);
	}else
    {
	    printf("\n Verify Error!\n");
    }
	e5 = clock();

    free(vcp);
//    for (int i=0;i<vck->m;i++)
//    {
//        for (int l=0;l<t;l++)
//        {
//            free(nv[i*par->t+l]);
//        }
//    }
    free(nv);
    free(vck);
    free(y);


        // output the running time of all algorithms
        double T_kg=(double) (e1-s1)/CLOCKS_PER_SEC;
        double T_pg=(double) (e2-s2)/CLOCKS_PER_SEC;
        double T_cm=(double) (e3-s3)/CLOCKS_PER_SEC;
        double T_vf=(double) (e4-s4)/CLOCKS_PER_SEC;
	    double T_dec=(double) (e5-s5)/CLOCKS_PER_SEC;
        double T_ntvCm=(double) (e6-s6)/CLOCKS_PER_SEC;
        //double T_us=T_pg+T_vf;
	printf("m=%d, d=%d, t=%d, T_kg=%fs, T_pg=%fs, T_vf=%fs, T_cm=%fs, T_dec=%fs, T_ntvCm=%fs, T_us=%fs, T_ntvCm/T_vf=%f\n", mat_m, mat_d, t, T_kg, T_pg, T_vf, T_cm, T_dec, T_ntvCm, T_vf+T_dec,T_ntvCm/T_vf);
}

 int main(int argc, char **argv)
{
        if (md2_init())
        {
                printf("Testing FAILED\n");
                printf("Problem initializing the library\n");
               return 1;
        }
        lhe_par *par=malloc(sizeof(*par));
        lhep_new(par);

        for (int t=1;t<=5;t=t+4)
          	try(t*par->ni,t*par->ni,t,par);
            //try(4,4,1,par);
        return 0;
}