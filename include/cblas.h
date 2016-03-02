#ifndef __CBLAS__
#define __CBLAS__
double ddot_(int *, double *, int *, double *, int *);
void dspmv_(char *, char *, char *, int *, double *, double *, int *);
void dgemv_(char *, int *, int *, double *, double *, int *, double *, int *, double *, double *, int *);
void dgemm_(char *, char *, int *, int *, int *, double *, double *, int *, double *, int *, double *, double *, int *);
#endif
