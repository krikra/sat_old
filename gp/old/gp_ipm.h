struct _hist_stack
{
	double *data;
	int nrow;
	int ncol;
	int lu;
	int head;
	double *(*read)(struct _hist_stack *, int);
	double *(*write)(struct _hist_stack *, double *);
};

typedef struct _hist_stack hist_stack;

struct _lbfgs
{
	int k;
	hist_stack *y; // g_k - g_k-1
	hist_stack *s; // x_k - x_k-1
	hist_stack *rho; // rho_k = y_k^t * g_k
};

typedef struct _lbfgs lbfgs;
