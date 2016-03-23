#ifndef __MY_GP__
#define __MY_GP__
struct _gp
{
	double *R; //matrix (packed) (N * N)
	double *FRF; //matrix (packed) (k * k)
	
	double *F; //matrix (k * N)
	double *e; //vector (N)
	double *rr; //vector (N)
	//double *uu; //matrix (k * n)

	double *ux; //matrix (dim * N)
	double *uy; //vector (N)

	double *mean;
	double *mse;

	double sigma;
	double *beta;
	double *phi;

	double **x;
	//double *xx; //matrix dim * n

	int size_buf;
	int dim;
	int *nd;
	int dd; //n
	int ud; //N
	int k; //k
};

typedef struct _gp GP;
#endif
