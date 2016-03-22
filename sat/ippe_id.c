#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "sat.h"

void initial_grid(SAT *ip, int *nd_init)
{
	int j;

	int iter;

	int dim;
	int *nd;

	int coo,foo,tmp,tmptmp;
	int ind = 0;

	dim = ip->dim;
	nd = ip->nd_whole;

	for(iter=0;iter<ip->dd_init;iter++)
	{
		ip->id[iter] = 0;
		tmp = ip->dd_init;
		tmptmp = ip->dd_whole;
		foo = iter;
		for(j=dim-1;j>=0;j--)
		{
			tmp = tmp / nd_init[j];
			tmptmp = tmptmp / nd[j];
			coo = foo / tmp;
			foo = foo % tmp;
			ip->id[iter] += (int)(coo * (nd[j] - 1.0) / (nd_init[j] - 1.0)) * tmptmp;
			//ind += (int)(coo * (nd[j] - 1) / 3.0) * tmptmp;
		}
	}
	//return(ind);
}

void initial_user_2d(SAT *ip, char *path)
{
	int iter;
	FILE *fp = NULL;
	int num = ip->dd_init;
	char lhd[128];

	double *xx1, *xx2;
	int x1, x2;

	xx1 = malloc(sizeof(double) * num);
	xx2 = malloc(sizeof(double) * num);

	fp = fopen(path, "r");
	if(fp == NULL){printf("idopen fail!!!\n");}
	for(iter=0;iter<num;iter++)
	{
		fscanf(fp, "%lf %lf", &xx1[iter], &xx2[iter]);
	}
	fclose(fp);

	for(iter=0;iter<num;iter++)
	{
		x1 = (int)(xx1[iter] * ip->nd_whole[0]);
		x2 = (int)(xx2[iter] * ip->nd_whole[1]);
		ip->id[iter] = x2 * ip->nd_whole[0] + x1;
	}

	free(xx1);
	free(xx2);

	//return(x1 * ip->nd_whole[0] + x2);
}

void initial_user_3d(SAT *ip, char *path)
{
	int iter;
	FILE *fp;
	int num = ip->dd_init;
	char lhd[128];

	double *xx1, *xx2, *xx3;
	int x1, x2, x3;

	xx1 = malloc(sizeof(double) * num);
	xx2 = malloc(sizeof(double) * num);
	xx3 = malloc(sizeof(double) * num);

	fp = fopen(path, "r");
	if(fp == NULL){printf("idopen fail!!!\n");}
	for(iter=0;iter<num;iter++)
	{
		fscanf(fp, "%lf %lf %lf", &xx1[iter], &xx2[iter], &xx3[iter]);
	}
	fclose(fp);

	for(iter=0;iter<num;iter++)
	{
		x1 = (int)(xx1[iter] * ip->nd_whole[0]);
		x2 = (int)(xx2[iter] * ip->nd_whole[1]);
		x3 = (int)(xx3[iter] * ip->nd_whole[2]);
		ip->id[iter] = x3 * ip->nd_whole[1] * ip->nd_whole[0] + x2 * ip->nd_whole[0] + x1;
	}

	free(xx1);
	free(xx2);
	free(xx3);
}

void scan_2d(FILE *fp, double *xx, const int n)
{
	int iter;

	for(iter=0;iter<n;iter++)
	{
		fscanf(fp, "%lf %lf", &xx[iter], &xx[1*n+iter]);
	}
}

void scan_3d(FILE *fp, double *xx, const int n)
{
	int iter;

	for(iter=0;iter<n;iter++)
	{
		fscanf(fp, "%lf %lf %lf", &xx[iter], &xx[1*n+iter], &xx[2*n+iter]);
	}
}

void scan_4d(FILE *fp, double *xx, const int n)
{
	int iter;

	for(iter=0;iter<n;iter++)
	{
		fscanf(fp, "%lf %lf %lf %lf", &xx[iter], &xx[1*n+iter], &xx[2*n+iter], &xx[3*n+iter]);
	}
}

void scan_8d(FILE *fp, double *xx, const int n)
{
	int iter;

	for(iter=0;iter<n;iter++)
	{
		fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf", &xx[iter], &xx[1*n+iter], &xx[2*n+iter], &xx[3*n+iter], &xx[4*n+iter], &xx[5*n+iter], &xx[6*n+iter], &xx[7*n+iter]);
	}
}

void initial_user(SAT *ip, char *path)
{
	int i;
	int iter;
	int tmp;
	FILE *fp = NULL;
	int num = ip->dd_init;
	int dim = ip->dim;
	char lhd[128];

	double *xx;
	int *x;

	xx = malloc(sizeof(double) * num * dim);
	x = malloc(sizeof(int) * dim);

	printf("%s\n", path);

	fp = fopen(path, "r");
	if(fp == NULL){printf("idopen fail!!!\n");}
	switch(dim)
	{
		case 2:
			scan_2d(fp, xx, num);break;
		case 3:
			scan_3d(fp, xx, num);break;
		case 4:
			scan_4d(fp, xx, num);break;
		case 8:
			scan_8d(fp, xx, num);break;
	}
	fclose(fp);

	for(iter=0;iter<num;iter++)
	{
		ip->id[iter] = 0;
		tmp = 1;
		for(i=0;i<dim;i++)
		{
			x[i] = (int)(xx[i * num + iter] * ip->nd_whole[i]);
			ip->id[iter] += tmp * x[i];
			tmp *= ip->nd_whole[i];
			//ip->id[iter] = x3 * ip->nd_whole[1] * ip->nd_whole[0] + x2 * ip->nd_whole[0] + x1;
		}
	}

	free(x);
	free(xx);
}
