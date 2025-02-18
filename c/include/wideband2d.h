

#include <complex.h>

// performs a diagonal high-frequency MP-to-MP translation
void h2dmpmphf_(int *nd, double _Complex *zk, double *rscale0, double *center0,
                  double _Complex *mpole0, int *nterms0,
                  double *rscale1, double *center1,
                  double _Complex *mpole1, int *nterms1);

// performs a diagonal high-frequency MP-to-Local translation
void h2dmplochf_(int *nd, double _Complex *zk, double *rscale0, double *center0,
                  double _Complex *mpole0, int *nterms0,
                  double *rscale1, double *center1,
                  double _Complex *mpole1, int *nterms1);

// performs a diagonal high-frequency Local-to-Local translation
void h2dloclochf_(int *nd, double _Complex *zk, double *rscale0, double *center0,
                  double _Complex *mpole0, int *nterms0,
                  double *rscale1, double *center1,
                  double _Complex *mpole1, int *nterms1);

