#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "sat.h"
#include "gp_kernel.h"

//#define SAT_DEBUG

double pnorm(double x)
{
	double pnorm;
	pnorm = 0.5 * (1.0 + erf(1.0 / sqrt(2.0) * x));
	return(pnorm);
}

double dnorm(double x)
{
	double dnorm;
	dnorm = 1.0 / sqrt(2.0 * M_PI) * exp(-0.5 * x * x);
	return(dnorm);
}

double ei(SAT *sat, int pos)
{
	double tmp, tmptmp, tmptmptmp;
	double ei;

	tmp = sat->bestv - sat->gp->mean[pos];
	tmptmp = sqrt(sat->gp->mse[pos]);
	tmptmptmp = tmp / tmptmp;
	ei = tmp * pnorm(tmptmptmp) / tmptmp + tmptmp * dnorm(tmptmptmp);

	return(ei);
}

void decomp(SAT *sat, int para, int *decomp)
{
	int i;
	int tmp, tmptmp;
	tmp = para;
	tmptmp = sat->dd_whole;
	for(i=sat->dim-1;i>=0;i--)
	{
		tmptmp /= sat->nd_whole[i];
		decomp[sat->dim - 1 - i] = tmp / tmptmp;
		tmp %= tmptmp;
	}
}
/*
void decomp(SAT *sat, int para, int *decomp)
{
	int i;
	int tmp, tmptmp;
	tmp = para;
	tmptmp = sat->dd_whole;
	for(i=sat->dim-1;i>=0;i--)
	{
		tmptmp /= sat->nd_whole[i];
		decomp[i] = tmp / tmptmp;
		tmp %= tmptmp;
	}
}
*/

int sat_choose(SAT *sat, int *para)
{
	int i, j;
	int where;
	int tmp;
	double t;
	double eimax = 0.0;
	double fn;
	double gr[2];

#ifdef SAT_DEBUG
	FILE *fp[3];
	char file[64];
	double *eiv;
#else
	double eiv;
#endif

	if(sat->cond){printf("estimation has been terminated.\n");return(1);}

	if(sat->iter < sat->dd_init)
	{
		*para = sat->id[sat->iter];
	}
	else
	{
#ifdef SAT_DEBUG
		eiv = malloc(sizeof(double) * sat->dd_whole);
#endif
		t = omp_get_wtime();
		gp_mle(sat->gp);
		t = omp_get_wtime() - t;
		printf("TIME_MLE %d %e\n", sat->iter, t);
		//log_lkhd(sat->gp, NULL, &fn, gr);
		//printf("%d PHI %e %e %e\n", sat->iter, sat->gp->phi[0], sat->gp->phi[1], sat->gp->phi[2]);
		t = omp_get_wtime();
		gp_surrogate(sat->gp);
		t = omp_get_wtime() - t;
		printf("TIME_SURROGATE %d %e\n", sat->iter, t);

		for(i=0;i<sat->dd_whole;i++)
		{
#ifdef SAT_DEBUG
			if(!sat->used[i])
			{
				eiv[i] = ei(sat, i);
			}else{eiv[i] = 0.0;}
			if(eimax < eiv[i]){eimax = eiv[i]; *para = i;/*printf("sat %e %d\n", eiv[i], i);*/}
#else
			if(!sat->used[i])
			{
				eiv = ei(sat, i);
			}else{eiv = 0.0;}
			if(eimax < eiv){eimax = eiv; *para = i;}
#endif
		}
#ifdef SAT_DEBUG
		sprintf(file, "./mean_%d.csv", sat->iter);
		fp[0] = fopen(file, "w");
		sprintf(file, "./mse_%d.csv", sat->iter);
		fp[1] = fopen(file, "w");
		sprintf(file, "./ei_%d.csv", sat->iter);
		fp[2] = fopen(file, "w");
		for(j=0;j<sat->dd_whole;j++)
		{
			fprintf(fp[0], "%d %d %e %e %e %e\n", j, sat->used[j], sat->gp->xx[j*sat->dim], sat->gp->xx[j*sat->dim+1], sat->gp->xx[j*sat->dim+2], sat->gp->mean[j]);
			fprintf(fp[1], "%d %d %e %e %e %e\n", j, sat->used[j], sat->gp->xx[j*sat->dim], sat->gp->xx[j*sat->dim+1], sat->gp->xx[j*sat->dim+2], sat->gp->mse[j]);
			fprintf(fp[2], "%d %d %e %e %e %e\n", j, sat->used[j], sat->gp->xx[j*sat->dim], sat->gp->xx[j*sat->dim+1], sat->gp->xx[j*sat->dim+2], eiv[j]);
		}
/*
		for(j=0;j<sat->nd_whole[1];j++)
		{
			for(i=0;i<sat->nd_whole[0];i++)
			{
				tmp = j * sat->nd_whole[0] + i;
				fprintf(fp[0], "%d %d %e %e %e\n", tmp, sat->used[tmp], sat->gp->xx[tmp*sat->dim], sat->gp->xx[tmp*sat->dim+1], sat->gp->mean[tmp]);
				fprintf(fp[1], "%d %d %e %e %e\n", tmp, sat->used[tmp], sat->gp->xx[tmp*sat->dim], sat->gp->xx[tmp*sat->dim+1], sat->gp->mse[tmp]);
				fprintf(fp[2], "%d %d %e %e %e\n", tmp, sat->used[tmp], sat->gp->xx[tmp*sat->dim], sat->gp->xx[tmp*sat->dim+1], eiv[tmp]);
			}
			fprintf(fp[0], "\n");
			fprintf(fp[1], "\n");
			fprintf(fp[2], "\n");
		}
*/
		fclose(fp[2]);
		fclose(fp[1]);
		fclose(fp[0]);
		free(eiv);
#endif
		printf("EI_MAX %e\n", eimax);
		if(eimax < sat->tol){sat->cond = 1; *para = sat->bestp; return(1);}
	}
	return(0);
}

int sat_update(SAT *sat, int para, double val)
{
	int i;
	int *p;

	if(sat->cond){printf("estimation has been terminated.\n");return(1);}

	p = malloc(sizeof(int) * sat->dim);

	decomp(sat, para, p);
	//printf("CHOSEN %e %e %e\n", sat->x[1][p[1]], sat->x[0][p[0]], val);
	for(i=0;i<sat->dim;i++)
	{
		//sat->gp->ux[sat->iter*sat->dim+i] = sat->x[i][p[i]];
		sat->gp->ux[sat->iter*sat->dim+i] = sat->gp->x[i][p[i]];
	}
	sat->gp->uy[sat->iter] = val;
	if(sat->bestv > val){sat->bestv = val;sat->bestp = para;}
	sat->used[para]++;

	sat->gp->ud = ++sat->iter;

	if(sat->iter >= sat->gp->size_buf){printf("BUFFER IS RAN OUT; ROUTINE IS TERMINATED\n");sat->cond = 2;return(2);}

	free(p);
	return(0);
}
