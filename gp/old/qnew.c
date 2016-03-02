void lbfgs_inv(const lbfgs *H, double *d, const double *g)
{
	int i, j, k;
	double *y, *s, *rho, *q;

	q = malloc(sizeof(double) * H->dim);

	for(i=0;i<n;i++)
	{
		q[i] = g[i];
	}

	for(j=H->k-1;j>=0;j--)
	{
		y = H->dy->read(j);
		s = H->dx->read(j);
		rho = H->rho->read(j);
		tmp = ddot_(&dim, s, &one, q, &one);
		for(i=0;i<dim;i++)
		{
			alpha[i] = tmp * rho[i];
			q[i] += tmp * rho[i] * y[i];
		}
	}

	for(j=0;j<H->k;j++)
	{
		//y = H->y->read(j);
		//s = H->s->read(j);
		//rho = H->rho->read(j);
		tmp = ddot_(&dim, y, &one, q, &one);
		for(i=0;i<dim;i++)
		{
			beta[i] = tmp * rho[i];
			d[i] = q[i] + (alpha[i] - tmp * rho[i]) * s[i];
		}
	}

	free(q);
}

void cauchy()
{
	int
}

void pqn_choose(double *x)
{
	int *a;
	double *d;
	double alpha;

	cauchy();
	direction(H, d, a);
	lnsearch(alpha);

	x += alpha * d;
}

int pqn_update(const double *x, const double y)
{
	if(r < rtol)
	{
		return(1);
	}
	else
	{
		H->dx->write(x);
		H->dy->write(y);
	}

	return(0);
}
