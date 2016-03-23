#include <stdlib.h>
#include <stdio.h>

#include "sat.h"
#include "gp_kernel.h"
#include "gp_mem.h"
#include "gp_misc.h"

void sat_setdim(SAT *sat, const int dim)
{
	sat->dim = dim;

	sat->iter = 0;
	sat->cond = 0;
	sat->bestp = 0;
	sat->nextp = 0;
	sat->bestv = 1.e+7;
	
	sat->nd_whole = malloc(sizeof(int) * dim);

	//sat->x = malloc(sizeof(double *) * dim);

	sat->gp = malloc(sizeof(GP));

	gp_setdim(sat->gp, dim);
}

void sat_setvec(SAT *sat, const int size_buf)
{
	int i;

	for(i=0;i<sat->dim;i++)
	{
		//sat->x[i] = malloc(sizeof(double) * sat->nd_whole[i]);
		sat->gp->nd[i] = sat->nd_whole[i];
	}
	sat->gp->dd = sat->dd_whole;

	sat->id = malloc(sizeof(int) * sat->dd_init);
	sat->used = calloc(sat->dd_whole, sizeof(int));
	
	gp_setvec(sat->gp, size_buf);
}

/*
void sat_x_init(SAT * sat)
{
	gp_x_init(sat->gp, sat->x);
}
*/

void sat_destroy(SAT *sat)
{
	int i;

	gp_destroy(sat->gp);
	free(sat->gp);
/*
	for(i=0;i<sat->dim;i++)
	{
		free(sat->x[i]);
	}
*/

	//free(sat->x);
	free(sat->id);
	free(sat->used);
	free(sat->nd_whole);
}

int sat_save(const SAT *sat, const char *path)
{
	int i;
	FILE *fp;
	fp = fopen(path, "wb");
	if(fp == NULL){printf("fopen failed!!!\n");return(1);}

	fwrite(&sat->dim, sizeof(int), 1, fp);
	fwrite(&sat->dd_init, sizeof(int), 1, fp);
	fwrite(&sat->dd_whole, sizeof(int), 1, fp);
	fwrite(sat->nd_whole, sizeof(int), sat->dim, fp);

	fwrite(&sat->iter, sizeof(int), 1, fp);
	fwrite(&sat->cond, sizeof(int), 1, fp);
	fwrite(&sat->tol, sizeof(double), 1, fp);
	fwrite(&sat->nextp, sizeof(int), 1, fp);
	fwrite(&sat->bestp, sizeof(int), 1, fp);
	fwrite(&sat->bestv, sizeof(double), 1, fp);

	fwrite(sat->id, sizeof(int), sat->dd_init, fp);
	fwrite(sat->used, sizeof(int), sat->dd_whole, fp);


	fwrite(&sat->gp->dim, sizeof(int), 1, fp);
	fwrite(&sat->gp->dd, sizeof(int), 1, fp);
	fwrite(&sat->gp->ud, sizeof(int), 1, fp);
	fwrite(&sat->gp->k, sizeof(int), 1, fp);
	fwrite(sat->gp->nd, sizeof(int), sat->gp->dim, fp);
	fwrite(&sat->gp->size_buf, sizeof(int), 1, fp);

	fwrite(sat->gp->R, sizeof(double), (sat->gp->size_buf * (sat->gp->size_buf + 1) / 2), fp);
	fwrite(sat->gp->FRF, sizeof(double), (sat->gp->k * (sat->gp->k + 1) / 2), fp);
	fwrite(sat->gp->F, sizeof(double), sat->gp->size_buf * sat->gp->k, fp);
	fwrite(sat->gp->e, sizeof(double), sat->gp->size_buf, fp);
	fwrite(sat->gp->rr, sizeof(double), sat->gp->size_buf, fp);
	//fwrite(sat->gp->uu, sizeof(double), sat->gp->size_buf * sat->gp->k, fp);
	fwrite(sat->gp->ux, sizeof(double), sat->gp->size_buf * sat->gp->dim, fp);
	fwrite(sat->gp->uy, sizeof(double), sat->gp->size_buf, fp);
	fwrite(sat->gp->mean, sizeof(double), sat->gp->dd, fp);
	fwrite(sat->gp->mse, sizeof(double), sat->gp->dd, fp);
	fwrite(&sat->gp->sigma, sizeof(double), 1, fp);
	fwrite(sat->gp->beta, sizeof(double), sat->gp->k, fp);
	fwrite(sat->gp->phi, sizeof(double), sat->gp->dim, fp);
	for(i=0;i<sat->gp->dim;i++)
	{
		fwrite(sat->gp->x[i], sizeof(double), sat->gp->nd[i], fp);
	}
	fclose(fp);
	return(0);
}

int sat_load(SAT *sat, const char *path)
{
	int i;
	FILE *fp;

	fp = fopen(path, "rb");
	if(fp == NULL){printf("fopen failed!!!\n");return(1);}

	fread(&sat->dim, sizeof(int), 1, fp);
	fread(&sat->dd_init, sizeof(int), 1, fp);
	fread(&sat->dd_whole, sizeof(int), 1, fp);
	fread(sat->nd_whole, sizeof(int), sat->dim, fp);

	fread(&sat->iter, sizeof(int), 1, fp);
	fread(&sat->cond, sizeof(int), 1, fp);
	fread(&sat->tol, sizeof(double), 1, fp);
	fread(&sat->nextp, sizeof(int), 1, fp);
	fread(&sat->bestp, sizeof(int), 1, fp);
	fread(&sat->bestv, sizeof(double), 1, fp);

	fread(sat->id, sizeof(int), sat->dd_init, fp);
	fread(sat->used, sizeof(int), sat->dd_whole, fp);


	fread(&sat->gp->dim, sizeof(int), 1, fp);
	fread(&sat->gp->dd, sizeof(int), 1, fp);
	fread(&sat->gp->ud, sizeof(int), 1, fp);
	fread(&sat->gp->k, sizeof(int), 1, fp);
	fread(sat->gp->nd, sizeof(int), sat->gp->dim, fp);
	fread(&sat->gp->size_buf, sizeof(int), 1, fp);

	fread(sat->gp->R, sizeof(double), (sat->gp->size_buf * (sat->gp->size_buf + 1) / 2), fp);
	fread(sat->gp->FRF, sizeof(double), (sat->gp->k * (sat->gp->k + 1) / 2), fp);
	fread(sat->gp->F, sizeof(double), sat->gp->size_buf * sat->gp->k, fp);
	fread(sat->gp->e, sizeof(double), sat->gp->size_buf, fp);
	fread(sat->gp->rr, sizeof(double), sat->gp->size_buf, fp);
	//fread(sat->gp->uu, sizeof(double), sat->gp->size_buf * sat->gp->k, fp);
	fread(sat->gp->ux, sizeof(double), sat->gp->size_buf * sat->gp->dim, fp);
	fread(sat->gp->uy, sizeof(double), sat->gp->size_buf, fp);
	fread(sat->gp->mean, sizeof(double), sat->gp->dd, fp);
	fread(sat->gp->mse, sizeof(double), sat->gp->dd, fp);
	fread(&sat->gp->sigma, sizeof(double), 1, fp);
	fread(sat->gp->beta, sizeof(double), sat->gp->k, fp);
	fread(sat->gp->phi, sizeof(double), sat->gp->dim, fp);
	for(i=0;i<sat->gp->dim;i++)
	{
		fread(sat->gp->x[i], sizeof(double), sat->gp->nd[i], fp);
	}
	fclose(fp);
	return(0);
}
