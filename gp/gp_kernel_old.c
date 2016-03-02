#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>

#include "lbfgsb.h"
#include "clapack.h"
#include "cblas.h"
#include "gp.h"

//#define SAT_DEBUG

int log_lkhd(GP *gp, double *phi, double *fn, double *gr)
{
	int i, j, k;
	int offset;
	int one = 1;

	int ind;

	int err;
	char U = 'U';
	char N = 'N';
	char T = 'T';
	char onec = '1';

	double det = 0.0;

	double *deriv;
	double *re;

	double tmp;
	double tmptmp = 0.0;

	double anorm;
	double rcond;
	double *dwork;
	int *iwork;

	deriv = malloc(sizeof(double) * gp->ud);
	re = malloc(sizeof(double) * gp->ud);

	R_packed_U(gp->R, gp->phi, gp->ux, gp->ud, gp->dim);

/*
	for(i=0,offset=0;i<gp->ud;i++)
	{
		for(j=0;j<=i;j++)
		{
			printf("%e ", gp->R[offset+j]);
		}
		printf("\n");
		offset+=i+1;
	}
*/
#ifdef SAT_DEBUG
	dwork = malloc(sizeof(double) * gp->ud * 3);
	iwork = malloc(sizeof(double) * gp->ud);
	anorm = dlansp_(&onec, &U, &gp->ud, gp->R, dwork);
	printf("anorm %e\n", rcond);
#endif

	dpptrf_(&U, &gp->ud, gp->R, &err); // R^-1
	if(err!=0){printf("CHOLESKY ERROR: %d\n", err);return(err);}

#ifdef SAT_DEBUG
	dppcon_(&U, &gp->ud, gp->R, &anorm, &rcond, dwork, iwork, &err);
	printf("rcond %e\n", rcond);
	if(err!=0){printf("RCOND ERROR: %d\n", err);return(err);}

	free(dwork);
	free(iwork);
#endif

/*
	for(i=0,offset=0;i<gp->ud;i++)
	{
		for(j=0;j<=i;j++)
		{
			printf("%e ", gp->R[offset+j]);
		}
		printf("\n");
		offset+=i+1;
	}
*/

	for(i=0;i<gp->ud;i++)
	{
		gp->e[i] = gp->uy[i];
		gp->F[i] = 1.0;
	}
	//dpptrs_(&U, &gp->ud, &one, gp->R, gp->e, &gp->dd, &err); // R^-1 * y
	dtpsv_(&U, &T, &N, &gp->ud, gp->R, gp->e, &one); // R^-1 * y
	//if(err!=0){printf("error: %d\n", err);}
	//dpptrs_(&U, &gp->ud, &one, gp->R, gp->F, &gp->dd, &err); // R^-1 * 1
	dtpsv_(&U, &T, &N, &gp->ud, gp->R, gp->F, &one); // R^-1 * 1
	//if(err!=0){printf("error: %d\n", err);}
	*gp->FRF = ddot_(&gp->ud, gp->F, &one, gp->F, &one); //1^t * R^-1 * 1
	//printf("FRF %e\n", *gp->FRF);

	*gp->beta = (1.0/(*gp->FRF)) * ddot_(&gp->ud, gp->F, &one, gp->e, &one); //beta = (1^t * R^-1 * 1)^-1 * 1^t * R^-1 * y

	//printf("beta %e\n", *gp->beta);

	for(i=0;i<gp->ud;i++)
	{
		gp->e[i] = gp->uy[i] - *gp->beta;
	}
	for(i=0;i<gp->ud;i++)
	{
		//printf("R1 %e\n", gp->F[i]);
	}
	//dpptrs_(&U, &gp->ud, &one, gp->R, gp->e, &gp->dd, &err); // R^-1 * e
	dtpsv_(&U, &T, &N, &gp->ud, gp->R, gp->e, &one); // R^-1 * e
	//if(err!=0){printf("error: %d\n", err);}
	gp->sigma = (1.0/(double)gp->ud) * ddot_(&gp->ud, gp->e, &one, gp->e, &one); // sigma = 1/N * e^t * R^-1 * e

	for(i=0,ind=0;i<gp->ud;i++)
	{
		det += log(gp->R[ind+i]);
		ind += i+1;
	}
	det *= 2.0;

	*fn = det + gp->ud * gp->sigma;

	for(i=0;i<gp->dim;i++)
	{
		gr[i] = 0.0;
	}

	for(i=0;i<gp->ud;i++)
	{
		re[i] = gp->e[i];
	}
	//dpptrs_(&U, &gp->ud, &one, gp->R, re, &gp->ud, &err);
	dtpsv_(&U, &N, &N, &gp->ud, gp->R, re, &one); // R^-1 * e
	for(k=0;k<gp->ud;k++)
	{
		for(i=0;i<gp->ud;i++)
		{
			deriv[i] = -re[k] * (gp->uy[i] - *gp->beta);
		}
/*
		for(j=0;j<gp->ud;j++)
		{
			if(k==j){deriv[j] = 1.0;}else{deriv[j] = 0.0;}
		}
*/
		deriv[k] += 1.0;
		dpptrs_(&U, &gp->ud, &one, gp->R, deriv, &gp->dd, &err);
/*
		for(j=0;j<gp->ud;j++)
		{
			deriv[j] -= re[j] * re[k];
		}
*/
		//if(err){printf("err:%d\n", err);}
		//dtpsv_(&U, &N, &N, &gp->ud, gp->R, deriv, &one); // R^-1 * e
		//dtpsv_(&U, &T, &N, &gp->ud, gp->R, deriv, &one); // R^-1 * e
		for(j=0;j<gp->ud;j++)
		{
			for(i=0;i<gp->dim;i++)
			{
				//printf("bef %e\n", tmp);
				kernel_deriv(i, &gp->ux[j*gp->dim], &gp->ux[k*gp->dim], gp->phi, gp->dim, &tmp);
				//printf("aft %e\n", tmp);
				//printf("%d %d %e %e %e\n", k, j, gp->ux[j*gp->dim], gp->ux[k*gp->dim], tmp);
				gr[i] += deriv[j] * tmp;
				//printf("%d %d gr[%d] += %e * %e = %e\n", k, j, i, deriv[j], tmp, tmptmp);
			}
		}
	}
//	printf("%e\n", tmptmp);

/*
	for(i=0;i<gp->ud;i++)
	{
		for(j=0;j<gp->ud;j++)
		{
			deriv[i * gp->ud + j] = - gp->e[i] * gp->e[j];
		}
		deriv[i * gp->ud + i] += 1.0;
		dtpsv_(&U, &N, &N, &gp->ud, gp->R, &deriv[i*gp->ud], &one);
		dtpsv_(&U, &T, &N, &gp->ud, gp->R, &deriv[i*gp->ud], &one);
	}
	//dpptrs_(&U, &gp->ud, &gp-ud, gp->R, gp->deriv, &gp->ud, &err);
	//dpptrs_(&U, &gp->ud, &gp-ud, gp->R, gp->deriv, &gp->ud, &err);
	for(k=0;k<gp->dim;k++)
	{
		gr[k] = 0.0;
		for(j=0;j<gp->ud;j++)
		{
			for(i=0;i<gp->ud;i++)
			{
				gr[k] += deriv[i * gp->ud + j] * kernel_deriv(k, &gp->ux[i*gp->dim], &gp->ux[j*gp->dim], gp->phi, gp->dim);
			}
		}
	}
*/

	free(deriv);
	free(re);
}

void gp_mle(GP *gp)
{
	int i, iter;
	int zero = -1;
	int mem = 20;
	int *nbd;
	double *l, *u;

	double *dwork;
	int *iwork;

	double ftol = 1.e+7;
	double gtol = 1.e-7;

	char task[60];
	long tl = 60;
	char csave[60];
	int lsave[4];
	int isave[44];
	double dsave[29];


	double fn;
	double *gr;

	gr = malloc(sizeof(double) * gp->dim);

	nbd = malloc(sizeof(int) * gp->dim);
	l = malloc(sizeof(double) * gp->dim);
	u = malloc(sizeof(double) * gp->dim);

	//printf("%d\n", (2 * mem * gp->dim + 5 * gp->dim + 11 * mem * mem + 8 * mem));
	dwork = malloc(sizeof(double) * (2 * mem * gp->dim + 5 * gp->dim + 11 * mem * mem + 8 * mem));
	iwork = malloc(sizeof(int) * (3 * gp->dim));

	for(i=0;i<gp->dim;i++)
	{
		gr[i] = 0.0;
		nbd[i] = 2;
		l[i] = 1.e-2;
		u[i] = 1.e+5;
	}

	for(i=0;i<61;i++)
	{
		task[i] = '\x0';
		csave[i] = '\x0';
	}
	task[0] = 'S'; task[1] = 'T'; task[2] = 'A'; task[3] = 'R'; task[4] = 'T'; 

	iter = 0;
	while(!strncmp(task, "FG", 2)||!strncmp(task, "NEW_X", 5)||!strncmp(task, "START", 5))
	{
		setulb_(&gp->dim, &mem, gp->phi, l, u, nbd, &fn, gr, &ftol, &gtol, dwork, iwork, task, &zero, csave, lsave, isave, dsave, 60, 60);
		printf("%s\n", task);
		printf("iter %d phi: ", iter);
		for(i=0;i<gp->dim;i++)
		{
			printf("%e ", gp->phi[i]);
		}
		printf("\n");
		if(!strncmp(task, "FG", 2))
		{
			log_lkhd(gp, NULL, &fn, gr);
		}
		//printf("%e %e: %e %e %e\n", gp->phi[0], gp->phi[1], fn, gr[0], gr[1]);
		++iter;
	}

	free(iwork);
	free(dwork);
	free(u);
	free(l);
	free(nbd);
	free(gr);
}

void gp_surrogate(GP *gp)
{
	int i, j;
	int one = 1;
	double oned = 1.0;
	double zerod = 0.0;
	int err;
	char U = 'U';
	char N = 'N';
	char T = 'T';

	r_init_1d(gp->rr, gp->xx, gp->ux, gp->phi, gp->dd, gp->ud, gp->dim);

	//dpptrs_(&U, &gp->ud, &gp->dd, gp->R, gp->rr, &gp->dd, &err); // R^-1 * r 
	for(i=0;i<gp->dd;i++)
	{
		dtpsv_(&U, &T, &N, &gp->ud, gp->R, &gp->rr[i*gp->dd], &one); // R^-1 * r 
	}
	//if(err!=0){printf("error: %d\n", err);}

/*
	for(i=0;i<gp->dd;i++)
	{
		for(j=0;j<gp->ud;j++)
		{
			printf("%e ", gp->rr[i*gp->dd+j]);
		}
		printf("\n");
	}
*/

	dgemv_(&T, &gp->ud, &gp->dd, &oned, gp->rr, &gp->dd, gp->e, &one, &zerod, gp->mean, &one); //mean = r^t * R^-1 * e

	dgemv_(&T, &gp->ud, &gp->dd, &oned, gp->rr, &gp->dd, gp->F, &one, &zerod, gp->uu, &one); //u = 1^t * R^-1 * r

	for(i=0;i<gp->dd;i++)
	{
		gp->uu[i] -= 1.0;
		//gp->uu[i] /= sqrt(&gp->FRF);

		gp->mean[i] += *gp->beta; // mean += beta
		gp->mse[i] = 1.0 - ddot_(&gp->ud, &gp->rr[i*gp->dd], &one, &gp->rr[i*gp->dd], &one);
		//printf("1 - rRr %e\n", gp->mse[i]);
		gp->mse[i] += gp->uu[i] * (1.0/(*gp->FRF)) * gp->uu[i];
		//printf("1 - rRr - uFRFu %e\n", gp->mse[i]);
		//gp->mse[i] += ddot_(&gp->k, &gp->uu[i*gp->k], &one, &gp->uu[i*gp->k], &one);
		gp->mse[i] *= gp->sigma; // mse = sigma * (1 - r^t * R^-1 * r + u^t * F^t * R^-1 * F * u)
	}
}
