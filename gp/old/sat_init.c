#include <math.h>

void malloc_R(double *R, int n)
{
	R = malloc(sizeof(double) * ((n * (n+1)) / 2));
}

void malloc_FRF(double *FRF, int k)
{
	FRF = malloc(sizeof(double) * ((k * (k+1)) / 2));
}

void gp_setvec(gp *gp)
{
	malloc_R(gp->R, gp->dd);
	malloc_FRF(gp->R, gp->k);
	gp->F = malloc(sizeof(double) * gp->k * gp->dd);
	gp->e = malloc(sizeof(double) * gp->dd);
	gp->rr = malloc(sizeof(double) * gp->dd * gp->dd);
	gp->uu = malloc(sizeof(double) * gp->k * gp->dd);

	gp->ux = malloc(sizeof(double) * gp->dim * gp->dd);
	gp->uy = malloc(sizeof(double) * gp->dd);

	gp->mean = malloc(sizeof(double) * gp->dd);
	gp->mse = malloc(sizeof(double) * gp->dd);

	gp->beta = malloc(sizeof(double) * gp->k);
	gp->phi = malloc(sizeof(double) * gp->dim);
}

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

double kernel_deriv(const double *x, const double *y, const double *phi, const int dim)
{
	int i;
	double tmp;

	for(i=0;i<dim;i++)
	{
		tmp -= phi[i] * (x[i] - y[i]) * (x[i] - y[i]);
	}
	return(-exp(tmp) * (x[i] - y[i]) * (x[i] - y[i]));
}

void R_packed_U(double *R, const double *phi, const double *ux, const int ud, const int dim)
{
	int i, j;
	int offset = 0;

	for(i=0;i<ud;i++)
	{
		for(j=0;j<=i;j++)
		{
			R[offset + j] = kernel(ux[i*ud], ux[j*ud], phi, dim);
		}
		offset += i;
	}
}

void R_packed_L(double *R, const double *phi, const double *ux, const int ud, const int dim)
{
	int i, j;
	int offset = 0;

	for(i=0;i<nd;i++)
	{
		for(j=i;j<nd;j++)
		{
			R[offset + j-i] = kernel(ux[i*ud], ux[j*ud], phi, dim);
		}
		offset += nd - i;
	}
}

void R_dense(double *R, const double *phi, const double *ux, const int ud, const int dim)
{
	int i, j;

	for(i=0;i<ud;i++)
	{
		for(j=0;j<ud;j++)
		{
			R[i * ud + j] = kernel(ux[i*ud], ux[j*ud], phi, dim);
		}
	}
}

void sample_update(double *ux, double *uy, const double *pos, const double val)
{
	for(i=0;i<dim;i++)
	{
		ux[iter*dim+i] = x[i];
	}
	uy[iter] = val;
}

void r_init_2d(double *rr, const double *ux, const double *phi, const int *nd, const int ud, const int dim)
{
	int i, j, h;
	int tmp;
	double x[2];

	for(j=0;j<nd[1];j++)
	{
		x[1] = j / (nd[1] - 1);
		for(i=0;i<nd[0];i++)
		{
			x[0] = i / (nd[0] - 1);
			tmp = j * nd[0] + i;
			for(h=0;h<ud;h++)
			{
				rr[tmp*ud+h] = kernel(x, &ux[h*ud], phi, dim);
			}
		}
	}
}

void r_init_1d(double *rr, const double *ux, const double *phi, const int *nd, const int ud, const int dim)
{
	int i, h;
	int tmp;
	double x;

	for(i=0;i<nd[0];i++)
	{
		x = i / (nd[0] - 1);
		for(h=0;h<ud;h++)
		{
			rr[i*ud+h] = kernel(&x, &ux[h*ud], phi, dim);
		}
	}
}

void set_e(double *e, const double *beta, const double *y, const int ud)
{
	int i;
	for(i=0;i<ud;i++)
	{
		e[i] = y[i] - beta;
	}
}
