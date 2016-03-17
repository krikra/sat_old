#ifndef __CBLAS__
#define __CBLAS__
double ddot_(const int *, const double *, const int *, const double *, const int *);
void dspmv_(const char *, const int *, const double *, const double *, const double *, const int *, const double *, double *, const int *);
void dgemv_(const char *, const int *, const int *, const double *, const double *, const int *, const double *, const int *, const double *, double *, const int *);
void dgemm_(const char *, const char *, const int *, const int *, const int *, const double *, const double *, const int *, const double *, const int *, const double *, double *, const int *);
#endif
