#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "cblas.h"

void power(const double *R, const int ud, double *lambda)
{
	int one = 1;
	double oned = 1.0;
	double zerod = 0.0;
	char U = 'U';

	int i;

	double *p, *tmp, *x_old, *x_new;
	double xx, xrx, xrrx, res;
	double rho, rho_old;

	p = malloc(sizeof(double) * 2 * ud);
	x_old = p;
	x_new = (p + ud);

	for(i=0;i<ud;i++)
	{
		x_old[i] = 1.0;
	}
	rho = 0.0;

	for(i=0;i<ud;i++)
	{
		rho_old = rho;
		xx = ddot_(&ud, x_old, &one, x_old, &one);
		dspmv_(&U, &ud, &oned, R, x_old, &one, &zerod, x_new, &one);
		xrx = ddot_(&ud, x_new, &one, x_old, &one);
		rho = xrx / xx;
		res = fabs(rho - rho_old);
		//printf("%e %e\n", rho, rho_old);
		printf("%d RESIDUAL %e\n", i, res);
		if(res < 1.e-3)
		{
			break;
		}

		xrrx = 1.0 / ddot_(&ud, x_new, &one, x_new, &one);
		dscal_(&ud, &xrrx, x_new, &one);

		tmp = x_old;
		x_old = x_new;
		x_new = tmp;
	}
	*lambda = rho;
	free(p);
}

double nugget(const double *R, const int ud, double *nug)
{
	const double eps = exp(25);
	double lambda;

	power(R, ud, &lambda);
	*nug = lambda / (eps - 1.0);
	return(*nug);
}

