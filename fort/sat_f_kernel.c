#include "sat.h"
#include "sat_sys.h"
#include "ippe_id.h"

void *sat_f_setdim_(const int *dim)
{
	SAT *sat;
	sat = malloc(sizeof(SAT));
	sat_setdim(sat, dim);
	return((void *)sat);
}

void sat_f_setvec_(void **sat, const int *n, const int *num, const int *buff, char *id, int *id_num)
{
	int i;
	*(SAT **)sat->dd_whole = 1;
	for(i=0;i<*(SAT **)sat->dim;i++)
	{
		*(SAT **)sat->nd_whole[i] = n[i];
		*(SAT **)sat->dd_whole *= n[i];
	}
	*(SAT **)sat->tol = 0.01;
	*(SAT **)sat->dd_init = num;
	sat_setvec(*(SAT **)sat, *buff);
		
	initial_user(sat, id);
}

void sat_f_choose_(void **sat, int *para)
{
	sat_choose(*(SAT **)sat, para);
}

void sat_f_update_(void **sat, int *para, double *val)
{
	sat_update(*(SAT **)sat, *para, *val);
}

int sat_f_terminated(void **sat)
{
	return(sat->cond);
}

void sat_f_destroy_(void **sat)
{
	sat_destroy(*(SAT **)sat);
}
