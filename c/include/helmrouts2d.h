
#include <complex.h>

// these are header files for some routines that needed to be called from C

// forms a multipole expansion from charges only
void h2dformmpc_(int *nd, double complex *zk, double *rscale,  double *source,
                 int *ns, double complex *charge, double *center,
                 int *nterms, double complex *mpole);

// performs a multipole-to-multipole translation
void h2dmpmp_(int *nd, double complex *zk, double *rscale1, double *center1,
              double complex *hexp1, int *nterms1, double *rscale2,
              double *center2,  double complex *hexp2,  int *nterms2);

