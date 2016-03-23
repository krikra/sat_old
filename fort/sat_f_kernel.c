#include <stdlib.h>
#include "sat.h"
#include "sat_mem.h"
#include "sat_sys.h"
#include "ippe_id.h"

#define fx10

void *sat_f_setdim_(const int *dim)
{
	SAT *sat;
	sat = malloc(sizeof(SAT));
	sat_setdim(sat, *dim);
	return(sat);
}
#ifdef fx10
void sat_f_setvec_(void **s, const int *n, const int *num, const int *buff, char *id, int *id_num)
#else
void sat_f_setvec_(void **s, const int **n, const int *num, const int *buff, char *id, int *id_num)
#endif
{
	SAT *sat;
	sat = *(SAT **)s;

	int i;
	sat->dd_whole = 1;
	for(i=0;i<sat->dim;i++)
	{
#ifdef fx10
		sat->nd_whole[i] = n[i];
		sat->dd_whole *= n[i];
#else
		sat->nd_whole[i] = (*(n))[i];
		sat->dd_whole *= (*(n))[i];
#endif
	}
	printf("debug %d %d %d\n", n[0], n[1], n[2]);
	sat->tol = 0.01;
	sat->dd_init = *num;
	sat_setvec(sat, *buff);
		
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
	return((*(SAT **)sat)->cond);
}

void sat_f_destroy_(void **sat)
{
	sat_destroy(*(SAT **)sat);
}
