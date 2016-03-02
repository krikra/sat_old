#ifndef __SAT_MEM__
#define __SAT_MEM__

int sat_setdim(SAT *, const int);
int sat_setvec(SAT *, const int);
//int sat_x_init(SAT *);

int sat_destroy(SAT *);

#endif
