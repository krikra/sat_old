#ifndef __CLAPACK__
#define __CLAPACK__
double dlansp_(char *, char *, int *, double *, double *);
double dppcon_(char *, int *, double *, double *, double *, double *, int *, int *);
void dpptrf_(const char *, const int *, double *, int *);
void dpptrs_(const char *, const int *, const int *, const double *, double *, const int *, int *);
void dtpsv_(const char *,const  char *,const  char *, const int *, const double *, double *, const int *);
#endif
