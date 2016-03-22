#include <stdlib.h>

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
