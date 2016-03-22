#ifndef __GP_MISC__
#define __GP_MISC__
double kernel(const double *, const double *, const double *, const int, double *);
double kernel_deriv(const int, const double *, const double *, const double *, const int, double *);
void R_packed_U(double *, const double *, const double *, const int, const int);
//void r_init_2d(double *, const double *, const double *, const double *, const int *, const int, const int);
//void r_init_1d(double *, const double *, const double *, const double *, const int, const int, const int);
void r_init(double *, const double **, const double *, const double *, const int *, const int, const int, const int, double *, const int);
//void gp_x_init(GP *, double ** /*, double *, int *, int*/);
#endif
