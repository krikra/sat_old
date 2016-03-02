#ifndef __CLAPACK__
#define __CLAPACK__
double dlansp_(char *, char *, int *, double *, double *);
double dppcon_(char *, int *, double *, double *, double *, double *, int *, int *);
void dpptrf_(char *, int *, double *, int *);
void dpptrs_(char *, int *, int *, double *, double *, int *, int *);
#endif
