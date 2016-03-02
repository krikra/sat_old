#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "gp.h"
#include "gp_mem.h"
#include "gp_misc.h"
#include "gp_kernel.h"

int main(int argc, char **argv)
{
	int i, j, k;
	int n;
	int offset;

	int one = 1;
	double oned = 1.0;
	char U = 'U';
	char L = 'L';
	int err;

	double testx[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
	double testb[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	double x[6] = {0.0, 2.0, 4.0, 6.0, 8.0, 10.0};
	double y[6] = {};

	double r;
	double phi;
	double fn[100], gr[100]={};
	double tmp;


	GP *gp;

	for(i=0;i<6;i++)
	{
		y[i] = sin(x[i]);
	}

	gp = malloc(sizeof(GP));

	gp_setdim(gp, 1);
	gp->k = 1;
	gp->nd[0] = 100;
	gp->dd = 100;
	gp->ud = 6;
	*gp->phi = 1.0;
	gp_setvec(gp);

	//gp->xx = malloc(sizeof(double) * gp->dd);

	for(i=0;i<gp->dd;i++)
	{
		gp->xx[i] = ((double)(i) / (double)(gp->dd-1.0)) * 10.0;
		//printf("xx%e\n", gp->xx[i]);
	}

	//printf("pass\n");

	for(i=0;i<gp->ud;i++)
	{
		gp->ux[i] = x[i];
		gp->uy[i] = y[i];
		//printf("%e %e\n", gp->ux[i], gp->uy[i]);
	}

	//R_packed_U(gp->R, gp->phi, gp->ux, gp->ud, gp->dim);
	/*
	for(i=0,offset=0;i<gp->ud;i++)
	{
		for(j=0;j<=i;j++)
		{
			printf("%e ", kernel_deriv(0, &gp->ux[i], &gp->ux[j], gp->phi, gp->dim, &tmp));
		}
		printf("\n");
		offset += i+1;
	}
	*/

	for(i=0;i<100;i++)
	{
		r = (((double)i / 99.0) * 6.0) - 3.0;
		*gp->phi = pow(10.0, r);
		log_lkhd(gp, &phi, &fn[i], &gr[i]);
		printf("%e %e %e\n", r, fn[i], gr[i]);
	}

	//printf("lkhd: %e\n", lkhd);
	r_init_1d(gp->rr, gp->xx, gp->ux, &phi, gp->dd, 6, 1);
	for(i=0;i<gp->dd;i++)
	{
		for(j=0;j<gp->ud;j++)
		{
			//printf("%e ", gp->rr[i*gp->dd+j]);
		}
		//printf("\n");
	}

	//gp_surrogate(gp);

/*
	for(i=0;i<gp->dd;i++)
	{
		printf("%e %e\n", gp->xx[i], gp->mean[i]);
	}
*/

	return(0);
}
