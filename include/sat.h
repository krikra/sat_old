#include "gp.h"

#ifndef __SAT__
#define __SAT__

struct _sat
{
	int dim;

	int *nd_whole;
	int dd_whole;
	int dd_init;

	int iter;
	int cond;
	double tol;

	int *id;
	int *used;

	double bestv;
	int bestp;
	int nextp;

	//double **x;

	GP *gp;

	void (* choose)(struct _sat *, int *);
	void (* update)(struct _sat *, int, double);
};

typedef struct _sat SAT;

#endif
