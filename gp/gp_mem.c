#include <stdlib.h>
#include <math.h>

#include "gp.h"

/*
void malloc_R(double *R, int n)
{
	R = malloc(sizeof(double) * ((n * (n+1)) / 2));
}

void malloc_FRF(double *FRF, int k)
{
	FRF = malloc(sizeof(double) * ((k * (k+1)) / 2));
}
*/

void gp_setdim(GP *gp, const int dim)
{
	int i;
	gp->dim = dim;
	gp->nd = malloc(sizeof(int) * dim);
	gp->phi = malloc(sizeof(double) * dim);
	gp->x = malloc(sizeof(double *) * dim);
	gp->k = 1;
	for(i=0;i<dim;i++)
	{
		gp->phi[i] = 10.0;
	}
}

void gp_setvec(GP *gp, const int size_buf)
{
	int i, j;
	gp->size_buf = size_buf;

	for(i=0;i<gp->dim;i++)
	{
		gp->x[i] = malloc(sizeof(double) * gp->nd[i]);
		for(j=0;j<gp->nd[i];j++)
		{
			gp->x[i][j] = (double)(j) / (double)(gp->nd[i] - 1.0);
		}
	}
	//malloc_R(gp->R, gp->dd);
	gp->R = malloc(sizeof(double) * ((size_buf * (size_buf + 1)) / 2));
	//malloc_FRF(gp->FRF, gp->k);
	gp->FRF = calloc(((gp->k * (gp->k+1)) / 2), sizeof(double));
	//gp->F = malloc(sizeof(double) * gp->k * gp->dd);
	gp->F = calloc(gp->k * gp->dd, sizeof(double));
	//gp->e = malloc(sizeof(double) * gp->dd);
	gp->e = calloc(size_buf, sizeof(double));
	gp->rr = calloc(size_buf, sizeof(double));
	// gp->uu is no longer used
	//gp->uu = calloc(gp->k * gp->dd, sizeof(double));
	//gp->uu = malloc(sizeof(double) * gp->k * gp->dd);

	gp->ux = malloc(sizeof(double) * gp->dim * size_buf);
	gp->uy = malloc(sizeof(double) * size_buf);

	//gp->xx = malloc(sizeof(double) * gp->dim * gp->dd);

	gp->mean = calloc(gp->dd, sizeof(double));
	gp->mse = calloc(gp->dd, sizeof(double));

	gp->beta = malloc(sizeof(double) * gp->k);
}

void gp_destroy(GP *gp)
{
	int i;
	for(i=0;i<gp->dim;i++)
	{
		free(gp->x[i]);
	}
	free(gp->x);
	free(gp->R);
	free(gp->FRF);
	free(gp->F);
	free(gp->e);
	free(gp->rr);
	//free(gp->uu);
	free(gp->ux);
	free(gp->uy);
	//free(gp->xx);
	free(gp->mean);
	free(gp->mse);
	free(gp->beta);
	free(gp->nd);
	free(gp->phi);
	//free(gp);
}
