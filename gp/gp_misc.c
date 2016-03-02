#include <stdio.h>
#include <math.h>

#include "gp.h"

double kernel(const double *x, const double *y, const double *phi, const int dim)
{
	int i;
	double tmp = 0.0;

	for(i=0;i<dim;i++)
	{
		tmp -= phi[i] * (x[i] - y[i]) * (x[i] - y[i]);
	}
	return(exp(tmp));
}

double kernel_deriv(const int cor, const double *x, const double *y, const double *phi, const int dim, double *deriv)
{
	//printf("%d %e %e %e %d\n", cor, *x, *y, *phi, dim);
	//printf("ker %e %e %e\n", (x[cor] - y[cor])*(x[cor] - y[cor]), kernel(x, y, phi, dim), -((x[cor] - y[cor]) * (x[cor] - y[cor])) * (kernel(x, y, phi, dim)));
	*deriv = -((x[cor] - y[cor]) * (x[cor] - y[cor])) * (kernel(x, y, phi, dim));
	return(-((x[cor] - y[cor]) * (x[cor] - y[cor])) * (kernel(x, y, phi, dim)));
}

/*
double kernel_deriv(const int cor, const double *x, const double *y, const double *phi, const int dim)
{
	int i;
	double tmp;

	for(i=0;i<dim;i++)
	{
		tmp -= phi[i] * (x[i] - y[i]) * (x[i] - y[i]);
	}
	return(-exp(tmp) * (x[cor] - y[cor]) * (x[cor] - y[cor]));
}
*/

void print_R_U(const double *R, const int ud)
{
	int i, j, offset;
	for(i=0,offset=0;i<ud;i++)
	{
		for(j=0;j<=i;j++)
		{
			printf("%e ", R[offset + j]);
		}
		printf("\n");
		offset += i+1;
	}
}

void R_packed_U(double *R, const double *phi, const double *ux, const int ud, const int dim)
{
	int i, j;
	int offset = 0;

	for(i=0;i<ud;i++)
	{
		for(j=0;j<=i;j++)
		{
			R[offset + j] = kernel(&ux[i*dim], &ux[j*dim], phi, dim);
		}
		offset += i+1;
	}
}

void R_packed_L(double *R, const double *phi, const double *ux, const int ud, const int dim)
{
	int i, j;
	int offset = 0;

	for(i=0;i<ud;i++)
	{
		for(j=i;j<ud;j++)
		{
			R[offset + j-i] = kernel(&ux[i*dim], &ux[j*dim], phi, dim);
		}
		offset += ud - i;
	}
}

void R_dense(double *R, const double *phi, const double *ux, const int ud, const int dim)
{
	int i, j;

	for(i=0;i<ud;i++)
	{
		for(j=0;j<ud;j++)
		{
			R[i * ud + j] = kernel(&ux[i*dim], &ux[j*dim], phi, dim);
		}
	}
}

/*
void sample_update(double *ux, double *uy, const double *pos, const double val)
{
	int i, j;
	for(i=0;i<dim;i++)
	{
		ux[iter*dim+i] = x[i];
	}
	uy[ud] = val;
}
*/

/*
void r_init_2d(double *rr, const double *xx, const double *ux, const double *phi, const int *nd, const int ud, const int dim)
{
	int i, j, h;
	int tmp;

	for(j=0;j<nd[1];j++)
	{
		for(i=0;i<nd[0];i++)
		{
			tmp = j * nd[0] + i;
			for(h=0;h<ud;h++)
			{
				rr[tmp*ud+h] = kernel(&xx[tmp*dim], &ux[h*dim], phi, dim);
			}
		}
	}
}
*/

/*
void r_init_1d(double *rr, const double *xx, const double *ux, const double *phi, const int dd, const int ud, const int dim)
{
	int i, h;
	int tmp;
	double x;

	for(i=0;i<dd;i++)
	{
		for(h=0;h<ud;h++)
		{
			rr[i*dd+h] = kernel(&xx[i*dim], &ux[h*dim], phi, dim);
		}
	}
}
*/

void r_init(double *rr, const double **xx, const double *ux, const double *phi, const int *n, const int ud, const int dim, const int size_buf, double *x, const int ii)
{
	int i, j;
	for(i=0;i<n[ii];i++)
	{
		x[ii] = xx[ii][i];
		if(ii > 0){r_init(rr, xx, ux, phi, n, ud, dim, size_buf, x, ii-1);}
		else
		{
			for(j=0;j<ud;j++)
			{
				rr[i*size_buf+j] = kernel(x, &ux[j*dim], phi, dim);
			}
		}
	}
}

/*
void gp_x_init(GP *gp, double **x)
{
	int i, j;
	int tmp, tmptmp, coo;

	for(i=0;i<gp->dd;i++)
	{
		tmp = i;
		tmptmp = gp->dd;
		for(j=gp->dim-1;j>=0;j--)
		{
			tmptmp /= gp->nd[j];
			coo = tmp / tmptmp;
			tmp = tmp % tmptmp;
			gp->xx[i*gp->dim+j] = x[j][coo];
		}
	}
}
*/
