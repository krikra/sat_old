/*
 * do j=k-1, k-2, ... , k-m
 * 	a = rho * s^t * q
 * 	q = q - a * y
 * enddo
 * do j=k-m, k-m+1, ... , k-1
 * 	b = rho * y^t * z
 * 	z = z - s * (a - b)
 * enddo
 */

void write(lbfgs *A, const double *x)
{
	int goma;
	for(i=0;i<ncol;i++)
	{
		A->data[A->lu * nrow + i] = x[i];
	}
	A->head = A->lu;
	A->lu = A->lu+1%A->length;
}

double *read(lbfgs *A, const int i)
{
	return(A->data[(A->head+i)%A->nrow]);
}

/*
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
		y = H->y->read(j);
		s = H->s->read(j);
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
*/

/*
 *  H  -I I
 *  LI  G
 * -LI    G
 */

void ipm_choose(double *x)
{
	double *x, *l;
	double *d, *dx, *dl;
	double ax, al;
	double mu;
	double step;

	int one = 1;

	int moi, mii;

	for(i=0;i<moi;i++)
	{
		for(j=0;j<mii;j++)
		{
			dtpsv_(&U, &N, &N, &n, H, d, &one);
			lns(ax, al);
			x += ax * dx;
			l += al * dl;
			if(res < itol){break;}
		}
		mu *= step;
		if(mu < otol){break;}
	}

	for(i=0;i<n;i++)
	{
		x[i] = new[i];
	}
}

void lns(double *ax, double *al)
{
	
}

void ipm_choose()
{
	dtpsv_(&U, &N, &N, &n, H, d, &one);
	lns(ax, al);
	x += ax * dx;
	l += al * dl;
}

void ipm_update(const double *x, const double fn, const double *gr)
{
	if(res > itol)
	{
		H->dx->write(dx);
		H->dg->write(dg);
		ii++;
	}
	else
	{
		H->clean;
		ii = 0;
		mu *= step;
	}
	if(mu < otol){terminated == 1;}
}
